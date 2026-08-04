// Submodules: OCCT BOPDS data structures
pub mod common_block;
pub mod face_info;
pub mod pave;
pub mod iterator;
pub use iterator::BOPDS_Iterator;
pub mod topods_builder;

// OCCT BOPDS_DS 1:1 translation.
// `BOPDS_DS.hxx` ?Data Structure for Boolean Operations.
//
// Maps:
//   TopoDS_Shape          ?Shape (rcad_kernel::topo_shape::Shape)
//   TopAbs_ShapeEnum      ?ShapeType
//   BOPDS_ShapeInfo       ?ShapeInfo (defined below)
//   BOPDS_IndexRange      ?IndexRange
//   BOPDS_CommonBlock     ?CommonBlock
//   BOPDS_PaveBlock       ?PaveBlock via SharedPB
//   BOPDS_FaceInfo        ?FaceInfo
//   BOPDS_InterfVV     ?InterferenceVV
//   NCollection_DynamicArray ?Vec
//   NCollection_DataMap   ?HashMap
//   NCollection_Map       ?HashSet
//   NCollection_List      ?Vec

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

/// Empty vertex data for placeholder vertices in push_edge/push_wire.
fn empty_vertex_data() -> rcad_kernel::topods::TVertexData {
    rcad_kernel::topods::TVertexData {
        my_shapes: Vec::new(), flags: 0,
        point: glam::DVec3::ZERO, tolerance: 0.0, points: Vec::new(),
    }
}
use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve3, Surface3};
use rcad_kernel::topods::{self, Orientation, ShapeType, TShape, TVertexData};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::base::bnd_lib::surface_bounding_box;
use rcad_kernel::curve_bounding_box;
use rcad_kernel::curve_bounding_box_range;
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::CurveEval;
use rcad_kernel::topology;
use rcad_kernel::{is_negative_infinite_value, is_positive_infinite_value};
use crate::bop::ds::face_info::FaceInfo;
use crate::bop::ds::pave::{Pave, PaveBlock, SharedPB};
use crate::bop::ds::common_block::CommonBlock;

/// Identifies which input shape a sub-shape came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShapeOrigin {
    ShapeA,
    ShapeB,
}

/// OCCT TopAbs_State: classification result against a shape.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Classification {
    In,
    Out,
    On,
}

/// =BOPDS_PassKey =sorted (index1, index2) pair key.
/// OCCT BOPDS_PassKey.hxx =wraps two integers with index1 <= index2.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PassKey {
    pub i1: usize,
    pub i2: usize,
}

/// =BOPTools_ConnexityBlock =connected component with IsRegular flag.
/// OCCT BOPTools_ConnexityBlock.hxx
#[derive(Debug, Clone)]
pub struct ConnexityBlock {
    pub shapes: Vec<usize>,
    pub is_regular: bool,
    pub loops: Vec<Vec<usize>>,
}

impl ConnexityBlock {
    pub fn new() -> Self {
        ConnexityBlock {
            shapes: Vec::new(),
            is_regular: false,
            loops: Vec::new(),
        }
    }
    pub fn is_regular(&self) -> bool {
        self.is_regular
    }
    pub fn set_regular(&mut self, r: bool) {
        self.is_regular = r;
    }
    pub fn shapes(&self) -> &[usize] {
        &self.shapes
    }
    pub fn add_shape(&mut self, s: usize) {
        self.shapes.push(s);
    }
    pub fn loops(&self) -> &[Vec<usize>] {
        &self.loops
    }
    pub fn change_shapes(&mut self) -> &mut Vec<usize> {
        &mut self.shapes
    }
    pub fn change_loops(&mut self) -> &mut Vec<Vec<usize>> {
        &mut self.loops
    }
}

impl PassKey {
    pub fn new(a: usize, b: usize) -> Self {
        if a <= b {
            PassKey { i1: a, i2: b }
        } else {
            PassKey { i1: b, i2: a }
        }
    }
}

/// =lightweight pair iterator over shape indices.
/// OCCT BOPDS_Iterator =produces sorted (i,j) pairs, optionally cross-group (A ).
pub struct PairIterator {
    i: usize,
    j: usize,
    a_end: usize,
    b_end: usize,
    done: bool,
    cross: bool, // true = cross-group (A ), false = all pairs (0..n)
}

impl PairIterator {
    /// OCCT: BOPDS_Iterator =iterate all pairs over [0, count).
    pub fn new(count: usize) -> Self {
        PairIterator {
            i: 0,
            j: 1,
            a_end: count,
            b_end: count,
            done: count < 2,
            cross: false,
        }
    }

    /// OCCT: BOPDS_Iterator::Prepare =iterate cross-group pairs A[end_a]  ?B[end_b..].
    /// For rcad: A = [0, a_end), B = [a_end, b_end).
    /// This matches the PaveFiller's A  cross-shape pair iteration pattern.
    pub fn prepare_ab(a_end: usize, b_end: usize) -> Self {
        let has_pairs = a_end > 0 && b_end > a_end;
        PairIterator {
            i: 0,
            j: a_end,
            a_end,
            b_end,
            done: !has_pairs,
            cross: true,
        }
    }

    pub fn more(&self) -> bool {
        !self.done
    }
    pub fn value(&self) -> PassKey {
        PassKey {
            i1: self.i,
            i2: self.j,
        }
    }

    pub fn next(&mut self) {
        if self.cross {
            self.j += 1;
            if self.j >= self.b_end {
                self.i += 1;
                self.j = self.a_end;
            }
            if self.i >= self.a_end {
                self.done = true;
            }
        } else {
            self.j += 1;
            if self.j >= self.b_end {
                self.i += 1;
                self.j = self.i + 1;
            }
            if self.i >= self.b_end - 1 || self.i >= self.a_end {
                self.done = true;
            }
        }
    }
}

/// =BOPDS_ShapeSD =same-domain shape mappings.
/// Wraps SharedTopologyInfo data with OCCT-style IsSubShape/HasSource queries.
/// OCCT BOPDS_ShapeSD.hxx, BOPDS_ShapeSD.cxx
#[derive(Debug, Clone)]
pub struct ShapeSD {
    sd_vertices: std::collections::HashSet<(usize, usize)>,
    sd_edges: std::collections::HashSet<(usize, usize)>,
    sd_faces: std::collections::HashSet<(usize, usize)>,
}

impl ShapeSD {
    pub fn new(a_count: usize, shared: &SharedTopologyInfo) -> Self {
        let mut sv: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
        let mut se: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
        let mut sf: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
        for &(a, b) in &shared.shared_vertices {
            // Store both (a,b) and (b,a) for bidirectional lookup.
            sv.insert((a, b));
            sv.insert((b, a));
        }
        for &(a, b) in &shared.shared_edges {
            se.insert((a, b));
            se.insert((b, a));
        }
        for &(a, b) in &shared.shared_faces {
            sf.insert((a, b));
            sf.insert((b, a));
        }
        ShapeSD {
            sd_vertices: sv,
            sd_edges: se,
            sd_faces: sf,
        }
    }

    /// OCCT: HasSource(sd, src) =true if sd has a same-domain counterpart src.
    pub fn has_source_vertex(&self, v: usize) -> bool {
        self.sd_vertices.contains(&(v, usize::MAX))
    }
    pub fn has_source_edge(&self, e: usize) -> bool {
        self.sd_edges.contains(&(e, usize::MAX))
    }
    pub fn has_source_face(&self, f: usize) -> bool {
        self.sd_faces.contains(&(f, usize::MAX))
    }

    /// OCCT: IsSubShape(shape) =true if shape participates in any SD mapping.
    pub fn is_sub_vertex(&self, v: usize) -> bool {
        self.sd_vertices.iter().any(|(a, _)| *a == v)
    }
    pub fn is_sub_edge(&self, e: usize) -> bool {
        self.sd_edges.iter().any(|(a, _)| *a == e)
    }
    pub fn is_sd_face(&self, fi: usize) -> bool {
        self.sd_faces.iter().any(|(a, _)| *a == fi)
    }

    pub fn has_sd_vertex(&self, a: usize, b: usize) -> bool {
        self.sd_vertices.contains(&(a, b))
    }
    pub fn has_sd_edge(&self, a: usize, b: usize) -> bool {
        self.sd_edges.contains(&(a, b))
    }
    pub fn has_sd_face(&self, a: usize, b: usize) -> bool {
        self.sd_faces.contains(&(a, b))
    }

    /// OCCT: ShapesSD iterator =(source, same_domain) pairs.
    pub fn sd_vertices_iter(&self) -> impl Iterator<Item = &(usize, usize)> {
        self.sd_vertices.iter()
    }

    /// AddShapeSD =register a dynamic same-domain vertex pair.
    pub fn add_sd_vertex(&mut self, a: usize, b: usize) {
        self.sd_vertices.insert((a, b));
        self.sd_vertices.insert((b, a));
    }

    /// HasShapeSD(n, nSD) =find the SD partner for a vertex.
    pub fn find_sd_partner(&self, v: usize) -> Option<usize> {
        self.sd_vertices
            .iter()
            .filter(|(a, _)| *a == v)
            .map(|(_, b)| *b)
            .min()
    }
}

/// Information about shared topology between the two input shapes.
///
/// This is used by the glue path to skip interference detection for
/// sub-shapes that are already coincident between the two inputs.
#[derive(Debug, Clone, Default)]
pub struct SharedTopologyInfo {
    /// Pairs of vertex indices (v_a, v_b) that are coincident.
    /// v_a is from ShapeA (index < a_vertex_count), v_b from ShapeB.
    pub shared_vertices: Vec<(usize, usize)>,
    /// Pairs of edge indices (e_a, e_b) that share the same geometry.
    /// e_a is from ShapeA (index < a_edge_count), e_b from ShapeB.
    pub shared_edges: Vec<(usize, usize)>,
    /// Pairs of face indices (f_a, f_b) that have shared topology.
    /// This includes both fully-overlapping faces and faces with partial overlap.
    pub shared_faces: Vec<(usize, usize)>,
    /// Face pairs with full boundary overlap (can be skipped entirely).
    pub fully_glued_faces: Vec<(usize, usize)>,
    /// Face pairs with partial edge sharing.
    pub partially_glued_faces: Vec<(usize, usize)>,
}

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
        range1: [f64; 2],
        range2: [f64; 2],
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

/// type-specific interference records replacing the flat Vec<Interference>.
/// OCCT BOPDS_DS stores interferences per-type in separate IndexedDataMaps,
/// which provide O(log n) lookup by shape index and natural pair dedup.
/// These are used by the new TypedInterferences container.
#[derive(Debug, Clone)]
pub struct InterferenceVV {
    pub v1: usize,
    pub v2: usize,
    pub merged_vertex: usize,
}

#[derive(Debug, Clone)]
pub struct InterferenceVE {
    pub vertex: usize,
    pub edge: usize,
    pub param: f64,
    /// OCCT IndexNew: new vertex index after tolerance-based fusing (UpdateVertex).
    pub index_new: usize,
}

#[derive(Debug, Clone)]
pub struct InterferenceEE {
    pub e1: usize,
    pub e2: usize,
    pub point: DVec3,
    pub param1: f64,
    pub param2: f64,
    pub new_vertex: usize,
    /// OCCT BOPDS_InterfEE::myCommonPart::Range1
    pub range1: [f64; 2],
    /// OCCT BOPDS_InterfEE::myCommonPart::Range2
    pub range2: [f64; 2],
}

#[derive(Debug, Clone)]
pub struct InterferenceVF {
    pub vertex: usize,
    pub face: usize,
    /// BOPDS_InterfVF::myU / myV (Interf.hxx L360-362).
    ///   UV coordinates of the vertex projection on the face surface.
    ///   Stored as (f64, f64) ?use u/v getters for clarity.
    pub u: f64,
    pub v: f64,
    /// BOPDS_Interf::myIndexNew (Interf.hxx L203).
    ///   Set via SetIndexNew() when UpdateVertex produces a new vertex
    ///   (SD resolution during VF processing).
    pub index_new: Option<usize>,
}

#[derive(Debug, Clone)]
pub struct InterferenceEF {
    pub edge: usize,
    pub face: usize,
    pub point: DVec3,
    pub edge_param: f64,
    pub new_vertex: usize,
}

/// BOPDS_Point (BOPDS_Point.hxx L29-72).
/// Stores intersection point data before it is promoted to a DS vertex.
#[derive(Debug, Clone)]
pub struct FFPoint {
    /// 3D intersection point (=BOPDS_Point::myPnt)
    pub pnt: DVec3,
    /// UV on face 1 (=BOPDS_Point::myPnt2D1)
    pub uv1: DVec2,
    /// UV on face 2 (=BOPDS_Point::myPnt2D2)
    pub uv2: DVec2,
    /// Index of the associated DS vertex, or usize::MAX if not yet assigned (=BOPDS_Point::myIndex, -1)
    pub vertex_index: usize,
}

impl FFPoint {
    pub fn new(pnt: DVec3, uv1: DVec2, uv2: DVec2) -> Self {
        Self {
            pnt,
            uv1,
            uv2,
            vertex_index: usize::MAX,
        }
    }
}

/// FF entry: keyed by (Fmin,Fmax) pair with all curves and touch points merged.
/// BOPDS_InterfFF (Interf.hxx L445-495).
#[derive(Debug, Clone)]
pub struct InterferenceFF {
    pub f1: usize,
    pub f2: usize,
    pub curves: Vec<usize>,
    /// BOPDS_Point array stored inline (not as DS vertex indices).
    /// OCCT equivalent: BOPDS_InterfFF::myPoints
    pub points: Vec<FFPoint>,
    /// BOPDS_InterfFF::myTangentFaces (Interf.hxx L490-492).
    ///   True when the two faces are tangent at the intersection curve(s).
    ///   Used by FillSameDomainFaces to decide whether faces can be merged.
    pub tangent_faces: bool,
}

/// BOPDS_InterfVZ (Interf.hxx L497-510).
///   Interference between a Vertex and a Solid (vertex is inside/on the solid).
#[derive(Debug, Clone)]
pub struct InterferenceVZ {
    pub vertex: usize,
    pub solid: usize,
}

/// BOPDS_InterfEZ (Interf.hxx L512-525).
///   Interference between an Edge and a Solid.
#[derive(Debug, Clone)]
pub struct InterferenceEZ {
    pub edge: usize,
    pub solid: usize,
}

/// BOPDS_InterfFZ (Interf.hxx L527-540).
///   Interference between a Face and a Solid.
#[derive(Debug, Clone)]
pub struct InterferenceFZ {
    pub face: usize,
    pub solid: usize,
}

/// BOPDS_InterfZZ (Interf.hxx L542-555).
///   Interference between two Solids.
#[derive(Debug, Clone)]
pub struct InterferenceZZ {
    pub s1: usize,
    pub s2: usize,
}

/// An intersection curve from F-F intersection, bounded by vertices.
/// =BOPDS_Curve (hxx:31-119).
#[derive(Debug, Clone)]
pub struct IntersectionCurve {
    pub curve: Curve3,
    /// Sampled points from numerical marching (non-empty for marched curves).
    /// When non-empty this takes priority over `curve` for face splitting.
    pub polyline: Vec<DVec3>,
    pub start_vertex: usize,
    pub end_vertex: usize,
    pub t_range: [f64; 2],
    /// PCurve (2D parametric curve) of this intersection on surface A (populated in Task 3+).
    pub pcurve_on_a: Option<Curve2d>,
    /// PCurve (2D parametric curve) of this intersection on surface B (populated in Task 3+).
    pub pcurve_on_b: Option<Curve2d>,
    /// tolerance of this section edge (CorrectToleranceOfSE).
    pub geom_tol: f64,
    /// =BOPDS_Curve::myPaveBlocks (hxx:115).
    /// Sub-segments of this intersection curve, created by splitting at paves.
    pub pave_blocks: Vec<crate::bop::ds::pave::SharedPB>,
    /// =BOPDS_Curve / IntTools_Curve extra fields.
    pub curve_extra: CurveExtra,
}

/// =IntTools_Curve (tangential_tol) + BOPDS_Curve
/// (techno_vertices, my_box) fields.
#[derive(Debug, Clone)]
pub struct CurveExtra {
    pub tangential_tol: f64,
    pub techno_vertices: Vec<usize>,
    pub my_box: Option<(glam::DVec3, glam::DVec3)>,
}

impl Default for CurveExtra {
    fn default() -> Self {
        CurveExtra {
            tangential_tol: 0.0,
            techno_vertices: Vec::new(),
            my_box: None,
        }
    }
}

impl IntersectionCurve {
    /// =BOPDS_Curve::ChangePaveBlock1 (lxx:96-100).
    pub fn change_pave_block1(pb_indices: &[crate::bop::ds::pave::SharedPB]) -> Option<usize> {
        pb_indices.first().map(|_| 0)
    }
}

/// Type of near-tangency between faces (used by glue detection).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NearTangentType {
    /// Planes that are nearly parallel.
    PlaneParallel,
    /// Cylinder tangent to plane.
    CylinderPlane,
    /// Sphere tangent to plane.
    SpherePlane,
    /// Two cylinders tangent along a generator.
    CylinderCylinder,
    /// Cone tangent to plane.
    ConePlane,
    /// General surface tangency.
    General,
}

/// =BOPDS_ShapeInfo =per-shape metadata in the flat DS index.
/// Corresponds to one entry in BOPDS_DS::myLines.

// ===
// BOPDS_IndexRange ?index range for an argument's shapes
// ===
#[derive(Debug, Clone, Copy)]
pub struct IndexRange {
    pub first: usize,
    pub last: usize,
}
impl IndexRange {
    pub fn new(f: usize, l: usize) -> Self { IndexRange { first: f, last: l } }
    pub fn contains(&self, i: usize) -> bool { i >= self.first && i <= self.last }
}

// ===
// BOPDS_ShapeInfo ?type, bounding box, sub-shapes, reference, flag
// ===
#[derive(Debug, Clone)]
pub struct ShapeInfo {
    pub shape: Shape,
    pub shape_type: ShapeType,
    pub bbox: BndBox, // OCCT: Bnd_Box myBox
    pub sub_shapes: Vec<usize>,
    pub reference: i64,
    pub flag: i64,
}
impl ShapeInfo {
    pub fn shape_type(&self) -> ShapeType { self.shape_type }
    pub fn shape(&self) -> &Shape { &self.shape }
    pub fn has_brep(&self) -> bool {
        matches!(self.shape_type, ShapeType::Vertex | ShapeType::Edge
            | ShapeType::Wire | ShapeType::Face | ShapeType::Shell)
    }
    pub fn is_interfering(&self) -> bool { self.has_brep() || self.shape_type == ShapeType::Solid }
    pub fn has_reference(&self) -> bool { self.reference >= 0 }
    pub fn reference(&self) -> i64 { self.reference }
    pub fn set_reference(&mut self, r: i64) { self.reference = r; }
    pub fn has_flag(&self) -> bool { self.flag >= 0 }
    pub fn flag(&self) -> i64 { self.flag }
    pub fn set_flag(&mut self, f: i64) { self.flag = f; }
    pub fn has_sub_shape(&self, i: usize) -> bool { self.sub_shapes.contains(&i) }
    pub fn sub_shapes(&self) -> &[usize] { &self.sub_shapes }
}

// ===
// BOPDS_DS ?Data Structure
// ===
#[derive(Debug)]
pub struct DS {
    // BOPDS_DS.hxx fields ?1:1 mapping
    pub arguments: Vec<Shape>,
    pub nb_source_shapes: usize,
    pub ranges: Vec<IndexRange>,
    pub shapes: Vec<ShapeInfo>,
    // (ptr_id, location) ?flat index  (TopoDS_Shape ?int map)
    pub map_shape_index: HashMap<(u64, u32), usize>,
    pub pave_blocks_pool: Vec<Vec<SharedPB>>,
    pub map_pb_cb: HashMap<u64, usize>,
    pub face_info_pool: Vec<FaceInfo>,
    pub shapes_sd: HashMap<usize, usize>,
    pub map_ve: HashMap<usize, Vec<usize>>,
    pub interf_tb: HashSet<(usize, usize)>,
    pub interf_vv: Vec<InterferenceVV>,  pub interf_ve: Vec<InterferenceVE>,
    pub interf_vf: Vec<InterferenceVF>,  pub interf_ee: Vec<InterferenceEE>,
    pub interf_ef: Vec<InterferenceEF>,  pub interf_ff: Vec<InterferenceFF>,
    pub interf_vz: Vec<InterferenceVZ>,  pub interf_ez: Vec<InterferenceEZ>,
    pub interf_fz: Vec<InterferenceFZ>,  pub interf_zz: Vec<InterferenceZZ>,
    pub interfered: HashSet<usize>,
    pub intersection_curves: Vec<crate::bop::int_tools::face_face::IntersectionCurve>,
    // CommonBlock storage (OCCT: myMapPBCB, but stored as Vec for index-based access)
    pub common_blocks: Vec<CommonBlock>,
}

impl DS {
    // ================================================================    // Construction / initialisation
    // ================================================================
    pub fn new() -> Self {
        DS {
            arguments: Vec::new(), nb_source_shapes: 0, ranges: Vec::new(),
            shapes: Vec::new(), map_shape_index: HashMap::new(),
            pave_blocks_pool: Vec::new(), map_pb_cb: HashMap::new(),
            face_info_pool: Vec::new(), shapes_sd: HashMap::new(), map_ve: HashMap::new(),
            interf_tb: HashSet::new(),
            interf_vv: Vec::new(), interf_ve: Vec::new(), interf_vf: Vec::new(),
            interf_ee: Vec::new(), interf_ef: Vec::new(), interf_ff: Vec::new(),
            interf_vz: Vec::new(), interf_ez: Vec::new(), interf_fz: Vec::new(),
            interf_zz: Vec::new(), interfered: HashSet::new(),
            common_blocks: Vec::new(),
            intersection_curves: Vec::new(),
        }
    }

    pub fn clear(&mut self) {
        self.nb_source_shapes = 0;
        self.arguments.clear();
        self.ranges.clear();
        self.shapes.clear();
        self.map_shape_index.clear();
        self.pave_blocks_pool.clear();
        self.face_info_pool.clear();
        self.shapes_sd.clear();
        self.map_ve.clear();
        self.map_pb_cb.clear();
        self.interf_tb.clear();
        self.interf_vv.clear();
        self.interf_ve.clear();
        self.interf_vf.clear();
        self.interf_ee.clear();
        self.interf_ef.clear();
        self.interf_ff.clear();
        self.interf_vz.clear();
        self.interf_ez.clear();
        self.interf_fz.clear();
        self.interf_zz.clear();
        self.interfered.clear();
        self.common_blocks.clear();
    }

    // ================================================================    // Arguments
    // ================================================================
    pub fn set_arguments(&mut self, a: Vec<Shape>) { self.arguments = a; }
    pub fn arguments(&self) -> &[Shape] { &self.arguments }

    // ================================================================    // Init
    // ================================================================
    /// BOPDS_DS::Init ?builds shape index, ranges, and bounding boxes.
    // OCCT BOPDS_DS.cxx L285-324
    pub fn init(&mut self, fuzz: f64) {
        if self.arguments.is_empty() { return; }
        let args = self.arguments.clone();
        let mut i1 = 0usize;
        for s in &args {
            if self.map_shape_index.contains_key(&(s.ptr_id(), s.location)) { continue; }
            let idx = self.append_shape(s.clone());
            self.init_shape(idx, s);
            let i2 = self.nb_shapes() - 1;
            self.ranges.push(IndexRange::new(i1, i2));
            i1 = i2 + 1;
        }
        if std::env::var("RCAD_EE_DEBUG").is_ok() {
            for i in 0..self.nb_shapes() {
                if self.shapes[i].shape_type == ShapeType::Face {
                    let surf = self.shapes[i].shape.as_face().and_then(|f| f.surface.clone());
                    let surf_desc = match &surf {
                        Some(rcad_kernel::geom::Surface3::Plane(p)) => format!("Plane n=({:.1},{:.1},{:.1})", p.normal.x, p.normal.y, p.normal.z),
                        Some(other) => format!("{:?}", std::mem::discriminant(other)),
                        None => "None".to_string(),
                    };
                    let wires: Vec<String> = self.shapes[i].sub_shapes.iter().map(|&w| {
                        let edges: Vec<String> = self.shapes.get(w).map(|ws| ws.sub_shapes.iter().map(|&e| e.to_string()).collect()).unwrap_or_default();
                        format!("w{}:{}", w, edges.join(","))
                    }).collect();
                    eprintln!("[DS-FACE] shape {} surf={} wires=[{}]", i, surf_desc, wires.join(" "));
                }
            }
        }
        self.nb_source_shapes = self.nb_shapes();
        // OCCT L312: max(theFuzz, Precision::Confusion()) * 0.5
        let tol = fuzz.max(1e-7) * 0.5;
        // OCCT L313-316: prepare
        self.prepare_vertices(tol);
        let an_edge_count = self.prepare_edges(tol);
        let a_face_count  = self.prepare_faces(tol);
        self.prepare_solids();
        // OCCT L319: buildVertexEdgeMap
        self.build_vertex_edge_map();
        // OCCT L322-323: prepare pools
        self.pave_blocks_pool.reserve(an_edge_count);
        self.face_info_pool.reserve(a_face_count);
    }

    /// OCCT BOPDS_DS::InitShape — add sub-shapes recursively.
    fn init_shape(&mut self, idx: usize, s: &Shape) {
        self.shapes[idx].shape_type = s.shape_type();
        // OCCT: no dedup — closed edges need duplicate vertex entries.
        let children = sub_shapes_of(s);
        for child in children {
            let pk = (child.ptr_id(), child.location);
            let ci = match self.map_shape_index.get(&pk) {
                Some(&e) => e,
                None => {
                    let ci = self.append_shape(child.clone());
                    self.init_shape(ci, &child);
                    ci
                }
            };
            self.shapes[idx].sub_shapes.push(ci);
        }
    }

    // ================================================================    // Queries ?shape count, range, rank
    // ================================================================
    pub fn nb_shapes(&self) -> usize { self.shapes.len() }
    pub fn nb_source_shapes(&self) -> usize { self.nb_source_shapes }
    pub fn nb_ranges(&self) -> usize { self.ranges.len() }
    pub fn range(&self, i: usize) -> &IndexRange { &self.ranges[i] }

    /// BOPDS_DS::Rank ?returns which argument (0-based) the shape belongs to.
    pub fn rank(&self, i: usize) -> isize {
        for ri in 0..self.nb_ranges() {
            if self.range(ri).contains(i) { return ri as isize; }
        }
        -1
    }

    pub fn is_new_shape(&self, i: usize) -> bool { i >= self.nb_source_shapes }

     /// Append with pre-built ShapeInfo.
    pub fn append(&mut self, si: ShapeInfo) -> usize {
        let pk = (si.shape.ptr_id(), si.shape.location);
        self.shapes.push(si);
        let idx = self.shapes.len() - 1;
        self.map_shape_index.insert(pk, idx);
        idx
    }

    /// Append shape, create default ShapeInfo.
    pub fn append_shape(&mut self, s: Shape) -> usize {
        let pk = (s.ptr_id(), s.location);
        let st = s.shape_type();
        self.shapes.push(ShapeInfo {
            shape: s, shape_type: st,
            bbox: BndBox::new(),
            sub_shapes: Vec::new(), reference: -1, flag: -1,
        });
        let idx = self.shapes.len() - 1;
        self.map_shape_index.insert(pk, idx);
        idx
    }

    pub fn shape_info(&self, i: usize) -> &ShapeInfo { &self.shapes[i] }
    pub fn change_shape_info(&mut self, i: usize) -> &mut ShapeInfo { &mut self.shapes[i] }
    pub fn shape(&self, i: usize) -> &Shape { &self.shapes[i].shape }
    pub fn index(&self, s: &Shape) -> isize {
        match self.map_shape_index.get(&(s.ptr_id(), s.location)) {
            Some(&i) => i as isize,
            None => -1,
        }
    }

    /// OCCT mutates TShapes in place (e.g. BRep_Builder::UpdateVertex), so the
    /// shape-index map keyed by TShape pointer never goes stale. rcad uses
    /// Arc::make_mut (clone-on-write) when the shape's TShape Arc is shared
    /// (referenced by multiple shapes): the clone changes the TShape pointer,
    /// so the map key would no longer match the shape's current pointer. This
    /// re-maps the current pointer to the same index, restoring the OCCT
    /// invariant that Index(shape) finds a shape by its current TShape.
    pub fn remap_shape_idx(&mut self, idx: usize) {
        if idx < self.shapes.len() {
            let (pk, loc) = {
                let si = &self.shapes[idx];
                (si.shape.ptr_id(), si.shape.location)
            };
            self.map_shape_index.insert((pk, loc), idx);
        }
    }

    // Pave blocks pool
    pub fn pave_blocks_pool(&self) -> &[Vec<SharedPB>] { &self.pave_blocks_pool }
    pub fn change_pave_blocks_pool(&mut self) -> &mut Vec<Vec<SharedPB>> { &mut self.pave_blocks_pool }
    pub fn has_pave_blocks(&self, i: usize) -> bool { self.shapes[i].has_reference() }
    /// OCCT: myPaveBlocksMap(theIndex) — pave blocks of edge by shape index.
    pub fn pave_blocks(&self, i: usize) -> &[SharedPB] {
        if self.has_pave_blocks(i) {
            &self.pave_blocks_pool[self.shapes[i].reference as usize]
        } else {
            &[]
        }
    }
    /// OCCT BOPDS_DS::ChangePaveBlocks — returns mutable ref to existing pave blocks.
    /// OCCT assumes InitPaveBlocks was called first. rcad: panics if not initialized.
    pub fn change_pave_blocks(&mut self, i: usize) -> &mut Vec<SharedPB> {
        let idx = self.shapes[i].reference;
        assert!(idx >= 0, "change_pave_blocks({}): not initialized, call init_pave_blocks first", i);
        &mut self.pave_blocks_pool[idx as usize]
    }

    /// Get vertex parameters on an edge (OCCT: BRep_Tool::Parameter).
    fn edge_vertex_params(&self, edge_idx: usize, v1_ds_idx: usize, v2_ds_idx: usize) -> (f64, f64) {
        let ed = match self.shapes[edge_idx].shape.as_edge() {
            Some(e) => e,
            None => return (0.0, 0.0),
        };
        let v1_brep_idx = self.shapes[v1_ds_idx].shape.index;
        let v2_brep_idx = self.shapes[v2_ds_idx].shape.index;
        // OCCT: BRep_Tool::Parameter reads stored param (set during edge construction).
        (ed.vertex_params.get(&v1_brep_idx).copied().unwrap_or(0.0),
         ed.vertex_params.get(&v2_brep_idx).copied().unwrap_or(0.0))
    }

    // ================================================================    // Common block map
    // ================================================================
    pub fn is_common_block(&self, pb: &SharedPB) -> bool {
        let ptr = std::sync::Arc::as_ptr(&pb.0) as u64;
        self.map_pb_cb.contains_key(&ptr)
    }
    pub fn common_block(&self, pb: &SharedPB) -> Option<usize> {
        let ptr = std::sync::Arc::as_ptr(&pb.0) as u64;
        self.map_pb_cb.get(&ptr).copied()
    }
    pub fn set_common_block(&mut self, pb: &SharedPB, cb: usize) {
        let ptr = std::sync::Arc::as_ptr(&pb.0) as u64;
        self.map_pb_cb.insert(ptr, cb);
    }
    /// OCCT BOPDS_DS::RealPaveBlock 鈥?if PB is in a CommonBlock, return CB's first PB.
    pub fn real_pave_block(&self, pb: &SharedPB) -> SharedPB {
        if let Some(cb_idx) = self.common_block(pb) {
            if let Some(cb) = self.common_blocks.get(cb_idx) {
                if let Some(fpb) = cb.pave_block1() {
                    return fpb;
                }
            }
        }
        pb.clone()
    }
    pub fn is_common_block_on_edge(&self, pb: &SharedPB) -> bool { self.common_block(pb).is_some() }

    /// Pool index of a SharedPB (OCCT: PaveBlock handle → pool position).
    /// Returns `(pool_index, position_within_pool)`.
    pub fn pb_pool_index(&self, pb: &SharedPB) -> Option<(usize, usize)> {
        let ptr = std::sync::Arc::as_ptr(&pb.0);
        for (pi, pool) in self.pave_blocks_pool.iter().enumerate() {
            for (li, spb) in pool.iter().enumerate() {
                if std::sync::Arc::as_ptr(&spb.0) == ptr {
                    return Some((pi, li));
                }
            }
        }
        None
    }

    /// OCCT BOPDS_DS::AddCommonBlock.
    /// Creates a new CommonBlock containing `the_pbs` and associates all PBs with it.
    /// The specific PaveBlock handles are stored so `PaveBlock1()` returns the
    /// exact block (an edge split into several blocks cannot be identified by a
    /// pool index alone).
    pub fn add_common_block(&mut self, the_pbs: &[SharedPB]) -> usize {
        let mut a_cb = CommonBlock::new();
        for pb in the_pbs {
            a_cb.add_pave_block(pb.clone(), 0); // face_idx = 0 placeholder
            pb.0.write().unwrap().common_block_idx = Some(self.common_blocks.len());
        }
        let cb_idx = self.common_blocks.len();
        self.common_blocks.push(a_cb);
        for pb in the_pbs {
            let ptr = std::sync::Arc::as_ptr(&pb.0) as u64;
            self.map_pb_cb.insert(ptr, cb_idx);
        }
        cb_idx
    }

    // ================================================================    // Face info pool
    // ================================================================
    pub fn face_info_pool(&self) -> &[FaceInfo] { &self.face_info_pool }
    pub fn change_face_info_pool(&mut self) -> &mut Vec<FaceInfo> { &mut self.face_info_pool }
    pub fn has_face_info(&self, i: usize) -> bool { self.shapes[i].has_reference() }
    pub fn face_info(&self, i: usize) -> &FaceInfo {
        if self.has_face_info(i) {
            &self.face_info_pool[self.shapes[i].reference as usize]
        } else {
            use std::sync::LazyLock;
            static E: LazyLock<FaceInfo> = LazyLock::new(FaceInfo::default);
            &E
        }
    }
    pub fn change_face_info(&mut self, i: usize) -> &mut FaceInfo {
        if !self.has_face_info(i) {
            let pi = self.face_info_pool.len();
            self.face_info_pool.push(FaceInfo::default());
            self.shapes[i].reference = pi as i64;
            // OCCT BOPDS_DS::InitFaceInfo (BOPDS_DS.cxx L738-747): InitFaceInfoIn +
            // UpdateFaceInfoOn. Note: only the direct VERTEX sub-shapes go into
            // VerticesIn here (OCCT InitFaceInfoIn L751-769); the boundary-edge
            // PBs populate VerticesOn/PaveBlocksOn via FaceInfoOn.
            self.face_info_pool[pi].set_index(i);
            let sub_shapes = self.shapes[i].sub_shapes.clone();
            for sub in sub_shapes {
                if sub >= self.nb_shapes() {
                    continue;
                }
                if self.shapes[sub].shape_type == ShapeType::Vertex {
                    let sd = self.get_same_domain_index(sub as isize);
                    if sd >= 0 {
                        self.face_info_pool[pi].vertices_in.insert(sd as usize);
                    }
                }
            }
            self.update_face_info_on(i);
        }
        &mut self.face_info_pool[self.shapes[i].reference as usize]
    }
    /// OCCT BOPDS_DS::UpdateFaceInfoIn (BOPDS_DS.cxx L773-791) + FaceInfoIn
    /// (BOPDS_DS.cxx L837-890): rebuild the face's PaveBlocksIn/VerticesIn.
    ///   step 1: pure internal (direct VERTEX) sub-shapes -> VerticesIn (SD)
    ///   step 2: VF interferences containing the face -> VerticesIn (SD vertex)
    ///   step 3: EF interferences containing the face -> VerticesIn (SD new
    ///           vertex) or, without a new vertex, the common-block PaveBlock1
    ///           -> PaveBlocksIn
    pub fn update_face_info_in(&mut self, the_index: usize) {
        if self.shapes[the_index].reference < 0 {
            return;
        }
        let mut pb_marks: Vec<usize> = Vec::new();
        let mut vertex_marks: Vec<usize> = Vec::new();
        // FaceInfoIn step 1: pure internal (free) vertices on the face.  OCCT
        // uses TopoDS_Iterator on the face BRep (direct children = wires + free
        // vertices), so the edge-endpoint vertices are NOT included.  The DS
        // face sub-shapes are flattened to edges + vertices (prepare_faces), so
        // a vertex that is a sub-shape of one of the face's edges is a boundary
        // vertex, not a free one.
        let sub_shapes = self.shapes[the_index].sub_shapes.clone();
        let mut edge_vertex_set: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &si in &sub_shapes {
            if si >= self.nb_shapes() {
                continue;
            }
            if self.shapes[si].shape_type == ShapeType::Edge {
                for &vi in &self.shapes[si].sub_shapes {
                    if vi < self.nb_shapes() {
                        edge_vertex_set.insert(vi);
                    }
                }
            }
        }
        for &si in &sub_shapes {
            if si >= self.nb_shapes() {
                continue;
            }
            if self.shapes[si].shape_type == ShapeType::Vertex && !edge_vertex_set.contains(&si) {
                let sd = self.get_same_domain_index(si as isize);
                if sd >= 0 {
                    vertex_marks.push(sd as usize);
                }
            }
        }
        // FaceInfoIn step 2: Vertex-Face interferences.
        for vf in &self.interf_vf {
            if vf.face == the_index {
                let sd = self.get_same_domain_index(vf.vertex as isize);
                if sd >= 0 {
                    vertex_marks.push(sd as usize);
                }
            }
        }
        // FaceInfoIn step 3: Edge-Face interferences.
        for ef in &self.interf_ef {
            if ef.face == the_index {
                if ef.new_vertex != usize::MAX {
                    let sd = self.get_same_domain_index(ef.new_vertex as isize);
                    if sd >= 0 {
                        vertex_marks.push(sd as usize);
                    }
                } else {
                    for pb in self.edge_pave_blocks(ef.edge) {
                        if let Some(cb_idx) = self.common_block(pb) {
                            if let Some(cb) = self.common_blocks.get(cb_idx) {
                                if cb.faces().contains(&the_index) {
                                    if let Some(pb1) = cb.pave_block1() {
                                        if let Some((pi, _)) = self.pb_pool_index(&pb1) {
                                            pb_marks.push(pi);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        let pfi = self.change_face_info(the_index);
        pfi.pave_blocks_in.clear();
        pfi.vertices_in.clear();
        for idx in pb_marks {
            pfi.pave_blocks_in.insert(idx);
        }
        for v in vertex_marks {
            pfi.vertices_in.insert(v);
        }
    }
    pub fn update_face_info_on(&mut self, the_index: usize) {
        // OCCT BOPDS_DS::UpdateFaceInfoOn (BOPDS_DS.cxx L792-807) + FaceInfoOn (L811-833):
        //   for each boundary edge, add its PBs' endpoint vertices to VerticesOn and the
        //   PBs to PaveBlocksOn; for each boundary vertex, add its same-domain index.
        if self.shapes[the_index].reference < 0 {
            return;
        }
        let sub_shapes = self.shapes[the_index].sub_shapes.clone();
        let mut pb_marks: Vec<usize> = Vec::new();
        let mut vertex_marks: Vec<usize> = Vec::new();
        for &si in &sub_shapes {
            if si >= self.nb_shapes() {
                continue;
            }
            match self.shapes[si].shape_type {
                ShapeType::Edge => {
                    let pbs = self.edge_pave_blocks(si).to_vec();
                    for pb in &pbs {
                        // OCCT FaceInfoOn: VerticesOn gets the ORIGINAL PB's endpoints;
                        // PaveBlocksOn gets the REAL pave block (resolved through the
                        // common block), so a coincident edge's PB is found there.
                        let (n_v1, n_v2) = { let r = pb.0.read().unwrap(); r.indices() };
                        vertex_marks.push(n_v1);
                        vertex_marks.push(n_v2);
                        let pbr = self.real_pave_block(pb);
                        let pptr = std::sync::Arc::as_ptr(&pbr.0) as u64;
                        if let Some(pidx) = self.pave_blocks_pool.iter().position(|pool| {
                            pool.iter().any(|spb| std::sync::Arc::as_ptr(&spb.0) as u64 == pptr)
                        }) {
                            pb_marks.push(pidx);
                        }
                    }
                }
                ShapeType::Vertex => {
                    let sd = self.get_same_domain_index(si as isize);
                    if sd >= 0 {
                        vertex_marks.push(sd as usize);
                    }
                }
                _ => {}
            }
        }
        let pfi = self.change_face_info(the_index);
        pfi.pave_blocks_on.clear();
        pfi.vertices_on.clear();
        for &idx in &pb_marks {
            pfi.pave_blocks_on.insert(idx);
        }
        for &v in &vertex_marks {
            pfi.vertices_on.insert(v);
        }
    }

    // ================================================================    // Same-domain shapes
    // ================================================================
    pub fn shapes_sd(&mut self) -> &mut HashMap<usize, usize> { &mut self.shapes_sd }
    pub fn add_shape_sd(&mut self, i: usize, sd: usize) {
        if i != sd { self.shapes_sd.insert(i, sd); }
    }
    pub fn has_shape_sd(&self, i: usize, sd: &mut usize) -> bool {
        let mut p = self.shapes_sd.get(&i);
        let mut f = false;
        while let Some(&n) = p { *sd = n; f = true; p = self.shapes_sd.get(&n); }
        f
    }
    /// OCCT BOPDS_DS::GetSameDomainIndex (BOPDS_DS.cxx L1244-1253): follows the
    /// same-domain chain unconditionally (myShapesSD.Seek until null), regardless
    /// of whether the SD index is smaller or larger than the query index. The
    /// SD vertex is usually a newly appended vertex (higher index), so resolving
    /// only downward (as the previous guard did) failed to unify vertices.
    pub fn get_same_domain_index(&self, i: isize) -> isize {
        let mut r = i;
        let mut guard = 0;
        loop {
            match self.shapes_sd.get(&(r as usize)) {
                Some(&n) => {
                    r = n as isize;
                    guard += 1;
                    if guard > 1000 { break; }
                }
                None => break,
            }
        }
        r
    }

    // ================================================================    // Interferences ?typed accessors
    // ================================================================
    pub fn interf_vv(&mut self) -> &mut Vec<InterferenceVV> { &mut self.interf_vv }
    pub fn interf_ve(&mut self) -> &mut Vec<InterferenceVE> { &mut self.interf_ve }
    pub fn interf_vf(&mut self) -> &mut Vec<InterferenceVF> { &mut self.interf_vf }
    pub fn interf_ee(&mut self) -> &mut Vec<InterferenceEE> { &mut self.interf_ee }
    pub fn interf_ef(&mut self) -> &mut Vec<InterferenceEF> { &mut self.interf_ef }
    pub fn interf_ff(&mut self) -> &mut Vec<InterferenceFF> { &mut self.interf_ff }
    pub fn interf_vz(&mut self) -> &mut Vec<InterferenceVZ> { &mut self.interf_vz }
    pub fn interf_ez(&mut self) -> &mut Vec<InterferenceEZ> { &mut self.interf_ez }
    pub fn interf_fz(&mut self) -> &mut Vec<InterferenceFZ> { &mut self.interf_fz }
    pub fn interf_zz(&mut self) -> &mut Vec<InterferenceZZ> { &mut self.interf_zz }

    pub fn nb_interf_types() -> usize { 10 }

    /// BOPDS_DS::AddInterf ?register an interference pair.
    pub fn add_interf(&mut self, i1: usize, i2: usize) -> bool {
        let k = if i1 < i2 { (i1, i2) } else { (i2, i1) };
        if self.interf_tb.insert(k) {
            self.interfered.insert(i1);
            self.interfered.insert(i2);
            true
        } else {
            false
        }
    }

    /// BOPDS_DS::HasInterf (single shape) ?true if shape has any interference.
    pub fn has_interf_single(&self, i: usize) -> bool { self.interfered.contains(&i) }

    /// BOPDS_DS::HasInterf (pair) ?true if the two shapes interfere.
    pub fn has_interf(&self, i1: usize, i2: usize) -> bool {
        let k = if i1 < i2 { (i1, i2) } else { (i2, i1) };
        self.interf_tb.contains(&k)
    }

    pub fn has_interf_shape_sub_shapes(&self, i1: usize, i2: usize, any: bool) -> bool {
        let s = &self.shapes[i2].sub_shapes;
        if s.is_empty() { return false; }
        if any { s.iter().any(|&ss| self.has_interf(i1, ss)) }
        else { s.iter().all(|&ss| self.has_interf(i1, ss)) }
    }

    pub fn has_interf_sub_shapes(&self, i1: usize, i2: usize) -> bool {
        self.shapes[i1].sub_shapes.iter().any(|&ss| self.has_interf_shape_sub_shapes(ss, i2, true))
    }

    pub fn interferences(&self) -> &HashSet<(usize, usize)> { &self.interf_tb }

    // ================================================================    // Dump
    // ================================================================
    pub fn dump(&self) -> String {
        let mut s = String::new();
        s.push_str(" *** DS ***\n");
        s.push_str(&format!(" Ranges: {}\n", self.nb_ranges()));
        for i in 0..self.nb_ranges() {
            let r = self.range(i);
            s.push_str(&format!("  range[{}]: [{},{}]\n", i, r.first, r.last));
        }
        s.push_str(&format!(" Shapes: {}\n", self.nb_shapes()));
        for i in 0..self.nb_shapes() {
            let si = self.shape_info(i);
            s.push_str(&format!("  {}: type={:?} ref={} flag={}\n",
                i, si.shape_type, si.reference, si.flag));
            if i == self.nb_source_shapes() - 1 { s.push_str(" ****** adds\n"); }
        }
        s.push_str(" ******\n");
        s
    }

    // ================================================================    // Sub-shape / topology queries
    // ================================================================
    pub fn is_sub_shape(&self, c: usize, p: usize) -> bool {
        self.shapes[p].sub_shapes.iter().any(|&s| s == c)
    }

    /// BOPDS_DS::Paves ?collect sorted paves for an edge.
    pub fn paves(&self, e: usize, lp: &mut Vec<Pave>) {
        let pbs = self.pave_blocks(e);
        if pbs.is_empty() { return; }
        let mut r: Vec<Pave> = Vec::new();
        for pb in pbs {
            let x = pb.0.read().unwrap();
            for pv in [&x.pave1, &x.pave2] {
                if !r.iter().any(|p: &Pave| p.vertex_idx == pv.vertex_idx && p.param == pv.param) {
                    r.push(*pv);
                }
            }
        }
        r.sort_by(|a, b| a.param.partial_cmp(&b.param).unwrap_or(std::cmp::Ordering::Equal));
        lp.extend(r);
    }

    /// Vertex count in source shapes.
    pub fn vertex_count(&self) -> usize {
        self.shapes.iter().filter(|s| s.shape_type == ShapeType::Vertex).count()
    }
    /// Edge count in source shapes.
    pub fn edge_count(&self) -> usize {
        self.shapes.iter().filter(|s| s.shape_type == ShapeType::Edge).count()
    }
    /// Face count in source shapes.
    pub fn face_count(&self) -> usize {
        self.shapes.iter().filter(|s| s.shape_type == ShapeType::Face).count()
    }
    /// Vertex count from shape A (first operand).
    pub fn a_vertex_count(&self) -> usize {
        self.shapes[..self.nb_source_shapes].iter().filter(|s| s.shape_type == ShapeType::Vertex).count()
    }
    /// Edge count from shape A.
    pub fn a_edge_count(&self) -> usize {
        self.shapes[..self.nb_source_shapes].iter().filter(|s| s.shape_type == ShapeType::Edge).count()
    }
    /// Face count from shape A.
    pub fn a_face_count(&self) -> usize {
        self.shapes[..self.nb_source_shapes].iter().filter(|s| s.shape_type == ShapeType::Face).count()
    }

    // ================================================================    // Update* methods
    // ================================================================
    pub fn update_pave_blocks_with_sd_vertices(&mut self) {
        if std::env::var("RCAD_EE_DEBUG").is_ok() {
            eprintln!("[EE-DBG] update_sd: shapes_sd={:?}", self.shapes_sd);
        }
        for list in self.pave_blocks_pool.clone() {
            for pb in &list { self.update_pave_block_with_sd_vertices(pb); }
        }
    }
    pub fn update_pave_block_with_sd_vertices(&self, pb: &SharedPB) {
        let mut w = pb.0.write().unwrap();
        w.pave1.vertex_idx = self.get_same_domain_index(w.pave1.vertex_idx as isize) as usize;
        w.pave2.vertex_idx = self.get_same_domain_index(w.pave2.vertex_idx as isize) as usize;
    }
    pub fn update_common_block_with_sd_vertices(&self, _cb: &CommonBlock) {
        // OCCT BOPDS_DS::UpdateCommonBlockWithSDVertices
    }

    /// OCCT BOPDS_DS::InitPaveBlocks — creates initial pave block for an edge.
    /// Must be called before ChangePaveBlocks (OCCT: Init → Change sequence).
    pub fn init_pave_blocks(&mut self, edge_idx: usize) {
        if self.has_pave_blocks(edge_idx) { return; }
        let spb = if self.shapes[edge_idx].shape_type == ShapeType::Edge
            && self.shapes[edge_idx].sub_shapes.len() >= 2
        {
            let mut n_v1 = self.shapes[edge_idx].sub_shapes[0];
            let mut n_v2 = self.shapes[edge_idx].sub_shapes[1];
            // OCCT BOPDS_DS::InitPaveBlocks (BOPDS_DS.cxx L437-485) appends each
            // vertex pave in edge order, then BOPDS_PaveBlock::Update (L291)
            // std::sorts the paves by parameter and creates one PB per consecutive
            // pair. For a closed edge (single shared vertex, e.g. a full circle)
            // the two paves are at range[0] and range[1]. For an open edge whose
            // vertex order is reversed relative to its curve (e.g. the sphere
            // seam, whose vertices are (north, south) but whose curve runs
            // south -> north), the geometric vertex parameters are sorted to keep
            // the PB range increasing, mirroring the std::sort.
            let (mut p1, mut p2): (f64, f64) = if n_v1 == n_v2 {
                self.shapes[edge_idx].shape.as_edge()
                    .map(|ed| (ed.range[0], ed.range[1]))
                    .unwrap_or((0.0, 0.0))
            } else {
                self.edge_vertex_params(edge_idx, n_v1, n_v2)
            };
            if p1 > p2 {
                std::mem::swap(&mut p1, &mut p2);
                std::mem::swap(&mut n_v1, &mut n_v2);
            }
            let pb = PaveBlock::new(edge_idx,
                Pave { vertex_idx: n_v1, param: p1 },
                Pave { vertex_idx: n_v2, param: p2 },
            );
            SharedPB::new(pb)
        } else {
            let p0 = Pave { vertex_idx: 0, param: 0.0 };
            SharedPB::new(PaveBlock::new(edge_idx, p0, p0))
        };
        self.pave_blocks_pool.push(vec![spb]);
        self.shapes[edge_idx].reference = (self.pave_blocks_pool.len() - 1) as i64;
    }

    pub fn init_pave_blocks_for_vertex(&mut self, v: usize) {
        let e: Vec<usize> = self.map_ve.get(&v).cloned().unwrap_or_default();
        for &ei in &e { self.init_pave_blocks(ei); }
    }

    pub fn release_pave_blocks(&mut self) {
        for i in 0..self.pave_blocks_pool.len() {
            if self.pave_blocks_pool[i].len() != 1 { continue; }
            let pb = &self.pave_blocks_pool[i][0];
            if self.is_common_block(pb) { continue; }
            let (v1, v2) = {
                let r = pb.0.read().unwrap();
                (r.pave1.vertex_idx, r.pave2.vertex_idx)
            };
            if !self.is_new_shape(v1) && !self.is_new_shape(v2) {
                let oe = pb.0.read().unwrap().original_edge;
                if oe < self.nb_shapes() { self.shapes[oe].reference = -1; }
                let ptr = std::sync::Arc::as_ptr(&pb.0) as u64;
                for e in &mut self.pave_blocks_pool {
                    e.retain(|spb| std::sync::Arc::as_ptr(&spb.0) as u64 != ptr);
                }
            }
        }
    }

    // OCCT BOPDS_DS::IsValidShrunkData (BOPDS_DS.cxx L1547-1585).
    // Compares the distances from the bounds of the shrunk range to the vertices
    // with the tolerance values of vertices.
    pub fn is_valid_shrunk_data(&self, pb: &PaveBlock) -> bool {
        if !pb.has_shrunk_data() { return false; }
        let (ts1, ts2, _) = pb.shrunk_data();
        let a_parameter = [ts1, ts2];
        let (v1i, v2i) = pb.indices();
        let a_vertex_index = [v1i, v2i];
        if v1i >= self.nb_shapes() || v2i >= self.nb_shapes() { return false; }

        // OCCT L1565-1567: BRepAdaptor_Curve on edge + epsilon = edge tolerance * 0.01
        let n_e = pb.original_edge;
        let curve = match self.edge_curve(n_e) {
            Some(c) => c,
            None => return false,
        };
        let an_epsilon = self.edge_tolerance(n_e) * 0.01;

        for i in 0..2 {
            // OCCT L1571-1572: aTol = BRep_Tool::Tolerance(aVertex) + Precision::Confusion()
            let a_vertex_point = self.vertex_point_by_idx(a_vertex_index[i]);
            let a_tol = self.vertex_tolerance_by_idx(a_vertex_index[i]) + rcad_kernel::CONFUSION;
            // OCCT L1577: gp_Pnt aPointOnEdge = aCurveAdaptor.Value(aParameter[i]);
            let a_point_on_edge = curve.point_at(a_parameter[i]);
            // OCCT L1579: if (aTol - aVertexPoint.Distance(aPointOnEdge) > anEpsilon) return false;
            let a_dist = a_vertex_point.distance(a_point_on_edge);
            if a_tol - a_dist > an_epsilon {
                return false;
            }
        }
        true
    }

    // ================================================================    // BuildBndBoxSolid ?compute solid bounding box from sub-shapes
    // ================================================================
    pub fn build_bnd_box_solid(&mut self, idx: usize, the_box: &mut (DVec3, DVec3, f64), _ci: bool) {
        let subs: Vec<usize> = self.shapes[idx].sub_shapes.clone();
        let mut faces: Vec<usize> = Vec::new();
        for &shi in &subs {
            if shi < self.nb_shapes() && self.shapes[shi].shape_type == ShapeType::Shell {
                faces.extend(self.shapes[shi].sub_shapes.clone());
            }
        }
        for &fi in &faces {
            if fi < self.nb_shapes() && self.shapes[fi].shape_type == ShapeType::Face {
                if let Some(b) = self.build_bnd_box(fi) {
                    if the_box.0.x.is_infinite() {
                        the_box.0 = b.0; the_box.1 = b.1; the_box.2 = b.2;
                    } else {
                        the_box.0 = the_box.0.min(b.0);
                        the_box.1 = the_box.1.max(b.1);
                        the_box.2 = the_box.2.max(b.2);
                    }
                }
                if self.shapes[fi].bbox.is_void() {
                    // open face ?solid is unbounded
                    the_box.0 = DVec3::splat(f64::NEG_INFINITY);
                    the_box.1 = DVec3::splat(f64::INFINITY);
                    return;
                }
            }
        }
    }

    // ================================================================    // Helpers ?DS internal
    // ================================================================
    /// Prepare vertex shape info: compute bounding boxes from geometry.
    // OCCT BOPDS_DS.cxx L1589-1610
    fn prepare_vertices(&mut self, tol: f64) -> usize {
        let mut a_vertex_count = 0;
        for a_vertex_index in 0..self.nb_source_shapes {
            if self.shapes[a_vertex_index].shape_type != ShapeType::Vertex { continue; }
            a_vertex_count += 1;
            // OCCT L1603-1606: SetGap(Tolerance + tol) + Add(point) → degenerate box [pt, pt]
            let vt = self.vertex_tolerance(&self.shapes[a_vertex_index].shape);
            if let Some(pt) = self.vertex_point_on_shape(&self.shapes[a_vertex_index].shape) {
                self.shapes[a_vertex_index].bbox = BndBox::from_point(pt);
                self.shapes[a_vertex_index].bbox.set_gap(vt + tol);
            }
        }
        a_vertex_count
    }

    // OCCT BOPDS_DS.cxx L1614-1692
    // OCCT BOPDS_DS.cxx L1614-1692
    fn prepare_edges(&mut self, tol: f64) -> usize {
        let mut an_edge_count = 0;
        for an_edge_index in 0..self.nb_source_shapes {
            // OCCT L1620: BOPDS_ShapeInfo& anEdgeInfo = ChangeShapeInfo(anEdgeIndex);
            // rcad: self.shapes[i] access — Rust borrow checker prevents
            // holding &mut ShapeInfo and &ShapeInfo simultaneously.
            if self.shapes[an_edge_index].shape_type != ShapeType::Edge { continue; }
            an_edge_count += 1;

            // OCCT L1627: const TopoDS_Edge& anEdge = TopoDS::Edge(anEdgeInfo.Shape());
            let shape = self.shapes[an_edge_index].shape.clone();
            // OCCT L1628: const double anEdgeTolerance = BRep_Tool::Tolerance(anEdge);
            let an_edge_tolerance = shape.as_edge().map_or(0.0, |ed| ed.tolerance);
            let is_degenerated = shape.as_edge().map_or(false, |ed| ed.degenerated);

            // OCCT L1630: if (!BRep_Tool::Degenerated(anEdge))
            if !is_degenerated {
                // OCCT L1636: NCollection_List<int>& aVertexIndices = anEdgeInfo.ChangeSubShapes();
                let a_forward_edge_curve_opt =
                    shape.as_edge().and_then(|ed| ed.curve.as_ref().map(|c| c.clone()));
                let a_range_opt = shape.as_edge().map(|ed| ed.range);
                if let (Some(ref a_forward_edge_curve), Some(a_range)) =
                    (a_forward_edge_curve_opt.as_ref(), a_range_opt)
                {
                    let (a_curve_start, a_curve_end) = (a_range[0], a_range[1]);
                    let mut new_vtx_indices: Vec<usize> = Vec::new();
                    if is_negative_infinite_value(a_curve_start) {
                        let a_point = a_forward_edge_curve.point_at(a_curve_start);
                        let a_vertex = Shape::new(
                            Arc::new(TShape::Vertex(TVertexData {
                                my_shapes: Vec::new(), flags: 0,
                                point: a_point, tolerance: an_edge_tolerance,
                                points: Vec::new(),
                            })),
                            0, Orientation::Forward,
                        );
                        let vi = self.append_shape(a_vertex);
                        self.shapes[vi].flag = 1;
                        new_vtx_indices.push(vi);
                    }
                    if is_positive_infinite_value(a_curve_end) {
                        let a_point = a_forward_edge_curve.point_at(a_curve_end);
                        let a_vertex = Shape::new(
                            Arc::new(TShape::Vertex(TVertexData {
                                my_shapes: Vec::new(), flags: 0,
                                point: a_point, tolerance: an_edge_tolerance,
                                points: Vec::new(),
                            })),
                            0, Orientation::Forward,
                        );
                        let vi = self.append_shape(a_vertex);
                        self.shapes[vi].flag = 1;
                        new_vtx_indices.push(vi);
                    }
                    self.shapes[an_edge_index].sub_shapes.extend(new_vtx_indices);
                }
            } else {
                self.shapes[an_edge_index].flag = an_edge_index as i64;
            }

            // OCCT L1677-1679: BRepBndLib::Add(anEdge, anEdgeBoundBox) adds curve box.
            // BRepBndLib::Add → BndLib_Add3dCurve::Add(BC, Tol, B)
            //   → GeomBndLib_Curve(C).Add(U1, U2, Tol, B)
            //   → per-type dispatch (each Box applies Enlarge(Tol) internally):
            //       Circle   → GeomBndLib_Circle::Box(gp_Circ, U1, U2, Tol)
            //       Ellipse  → GeomBndLib_Ellipse::Box(gp_Elips, U1, U2, Tol)
            //       Line     → GeomBndLib_Line::Box(gp_Lin, U1, U2, Tol)
            //       BSpline  → GeomBndLib_BSplineCurve::Box(U1, U2, Tol)
            //       Bezier   → GeomBndLib_BezierCurve::Box(U1, U2, Tol)
            //       Other    → GeomBndLib_OtherCurve::Box(U1, U2, Tol) (33pt sampling)
            let (mut an_edge_bx_min, mut an_edge_bx_max) =
                if let (Some(curve), Some(range)) = (
                    shape.as_edge().and_then(|ed| ed.curve.as_ref()),
                    shape.as_edge().map(|ed| ed.range),
                ) {
                    if let Some([cmn, cmx]) = curve_bounding_box_range(
                        curve, range[0], range[1], an_edge_tolerance,
                    ) {
                        (cmn, cmx)
                    } else {
                        (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY))
                    }
                } else {
                    (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY))
                };

            // OCCT L1681-1686: anEdgeBoundBox.Add(aVertexInfo.Box())
            // OCCT Add(Box) adds the full box INCLUDING the vertex tolerance gap.
            let sub_shapes = self.shapes[an_edge_index].sub_shapes.clone();
            for &vi in &sub_shapes {
                if !self.shapes[vi].bbox.is_void() {
                    if let Some((vxmn, vymn, vzmn, vxmx, vymx, vzmx)) = self.shapes[vi].bbox.get() {
                        if an_edge_bx_min.x.is_finite() {
                            an_edge_bx_min = an_edge_bx_min.min(DVec3::new(vxmn, vymn, vzmn));
                            an_edge_bx_max = an_edge_bx_max.max(DVec3::new(vxmx, vymx, vzmx));
                        } else {
                            an_edge_bx_min = DVec3::new(vxmn, vymn, vzmn);
                            an_edge_bx_max = DVec3::new(vxmx, vymx, vzmx);
                        }
                    }
                }
            }

            // OCCT L1688: SetGap(existing_gap + tol).
            // NOTE: no outer Enlarge here — the curve box already applied
            // Enlarge(anEdgeTolerance) inside GeomBndLib_*::Box, matching OCCT.
            if an_edge_bx_min.x.is_finite() {
                self.shapes[an_edge_index].bbox = BndBox::from_corners(
                    an_edge_bx_min.x, an_edge_bx_min.y, an_edge_bx_min.z,
                    an_edge_bx_max.x, an_edge_bx_max.y, an_edge_bx_max.z);
                let cur_gap = self.shapes[an_edge_index].bbox.get_gap();
                self.shapes[an_edge_index].bbox.set_gap(cur_gap + tol);
            }
        }
        an_edge_count
    }

    // OCCT BOPDS_DS.cxx L1696-1779
    fn prepare_faces(&mut self, tol: f64) -> usize {
        let mut a_face_count = 0;
        for a_face_index in 0..self.nb_source_shapes {
            if self.shapes[a_face_index].shape_type != ShapeType::Face { continue; }
            a_face_count += 1;

            let mut a_new_sub_shape_indices = HashSet::new();
            let shape = self.shapes[a_face_index].shape.clone();

            // OCCT L1715-1717: BRepBndLib::Add(aFace, aFaceBoundBox)
            let face_tolerance = shape.as_face().map_or(0.0, |fd| fd.tolerance);
            let (mut mn, mut mx) =
                if let Some(surface) = shape.as_face().and_then(|fd| fd.surface.as_ref()) {
                    let mut verts: Vec<topology::Vertex> = Vec::new();
                    for &wi in &self.shapes[a_face_index].sub_shapes {
                        if wi >= self.nb_shapes() { continue; }
                        for &ei in &self.shapes[wi].sub_shapes {
                            if ei >= self.nb_shapes() { continue; }
                            if self.shapes[ei].shape_type == ShapeType::Edge {
                                for &vi in &self.shapes[ei].sub_shapes {
                                    if vi < self.nb_shapes() {
                                        if let Some(pt) = self.vertex_point_on_shape(&self.shapes[vi].shape) {
                                            verts.push(topology::Vertex { point: pt });
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if let Some([cmn, cmx]) = surface_bounding_box(surface, &verts) {
                        (cmn - DVec3::splat(face_tolerance),
                         cmx + DVec3::splat(face_tolerance))
                    } else {
                        (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY))
                    }
                } else {
                    (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY))
                };

            // OCCT L1725-1752: iterate wires → edges
            let wire_indices = self.shapes[a_face_index].sub_shapes.clone();
            for &wi in &wire_indices {
                if wi >= self.nb_shapes() { continue; }
                let edge_indices = self.shapes[wi].sub_shapes.clone();
                for &ei in &edge_indices {
                    if ei >= self.nb_shapes() { continue; }
                    if self.shapes[ei].shape_type != ShapeType::Edge { continue; }

                    // OCCT L1729-1732: Add edge bounding box to face box
                    if !self.shapes[ei].bbox.is_void() {
                        mn = mn.min(self.shapes[ei].bbox.raw_min());
                        mx = mx.max(self.shapes[ei].bbox.raw_max());
                    }

                    // OCCT L1742-1746: Mark degenerated edges
                    if self.shapes[ei].shape.as_edge().map_or(false, |ed| ed.degenerated) {
                        self.shapes[ei].flag = a_face_index as i64;
                    }

                    // OCCT L1735 + L1748-1752: Add edge + vertices to map
                    a_new_sub_shape_indices.insert(ei);
                    for &vi in &self.shapes[ei].sub_shapes {
                        if vi < self.nb_shapes() {
                            a_new_sub_shape_indices.insert(vi);
                        }
                    }
                }
            }

            // OCCT L1756-1764: Add standalone face vertices (TopoDS_Iterator)
            if let Some(fd) = shape.as_face() {
                for iv in &fd.internal_vertices {
                    let pk = (iv.ptr_id(), iv.location);
                    if let Some(&i) = self.map_shape_index.get(&pk) {
                        a_new_sub_shape_indices.insert(i);
                    }
                }
            }

            // OCCT L1767-1773: Replace wire indices with edge+vertex indices
            self.shapes[a_face_index].sub_shapes = a_new_sub_shape_indices.into_iter().collect();

            if mn.x.is_finite() {
                self.shapes[a_face_index].bbox = BndBox::from_corners(
                    mn.x, mn.y, mn.z, mx.x, mx.y, mx.z);
                // OCCT L1775: SetGap(existing_gap + tol)
                let cur_gap = self.shapes[a_face_index].bbox.get_gap();
                self.shapes[a_face_index].bbox.set_gap(cur_gap + tol);
            }
        }
        a_face_count
    }

    // OCCT BOPDS_DS.cxx L1783-1852
    fn prepare_solids(&mut self) -> usize {
        // OCCT L1789-1792: check mode — only single argument
        if self.arguments.len() != 1 { return 0; }
        let mut a_solid_count = 0;
        for a_solid_index in 0..self.nb_source_shapes {
            if self.shapes[a_solid_index].shape_type != ShapeType::Solid { continue; }
            a_solid_count += 1;

            // OCCT L1807-1808: BuildBndBoxSolid
            let mut a_solid_bound_box =
                (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY), 0.0);
            self.build_bnd_box_solid(a_solid_index, &mut a_solid_bound_box, false);
            if a_solid_bound_box.0.x.is_finite() {
                self.shapes[a_solid_index].bbox = BndBox::from_corners(
                    a_solid_bound_box.0.x, a_solid_bound_box.0.y, a_solid_bound_box.0.z,
                    a_solid_bound_box.1.x, a_solid_bound_box.1.y, a_solid_bound_box.1.z);
            }

            // OCCT L1803-1804: map of sub-shape indices
            let mut a_new_sub_shape_indices = HashSet::new();

            // OCCT L1814-1839: iterate shells → faces → edges
            for &a_shell_index in &self.shapes[a_solid_index].sub_shapes.clone() {
                if a_shell_index >= self.nb_shapes() { continue; }
                if self.shapes[a_shell_index].shape_type != ShapeType::Shell { continue; }
                for &a_face_index in &self.shapes[a_shell_index].sub_shapes {
                    if a_face_index >= self.nb_shapes() { continue; }
                    if self.shapes[a_face_index].shape_type != ShapeType::Face { continue; }
                    a_new_sub_shape_indices.insert(a_face_index);
                    for &an_edge_index in &self.shapes[a_face_index].sub_shapes {
                        a_new_sub_shape_indices.insert(an_edge_index);
                    }
                }
            }

            // OCCT L1841-1848: replace shell indices with face+edge indices
            self.shapes[a_solid_index].sub_shapes = a_new_sub_shape_indices.into_iter().collect();
        }
        a_solid_count
    }

    fn build_vertex_edge_map(&mut self) {
        for ei in 0..self.nb_source_shapes {
            if self.shapes[ei].shape_type != ShapeType::Edge { continue; }
            for &vi in &self.shapes[ei].sub_shapes {
                if vi >= self.nb_shapes() { continue; }
                let e = self.map_ve.entry(vi).or_default();
                if !e.contains(&ei) { e.push(ei); }
            }
        }
    }

    fn build_bnd_box(&mut self, i: usize) -> Option<(DVec3, DVec3, f64)> {
        if !self.shapes[i].bbox.is_void() {
            let g = self.shapes[i].bbox.get_gap();
            return Some((self.shapes[i].bbox.raw_min(), self.shapes[i].bbox.raw_max(), g));
        }
        match self.shapes[i].shape_type {
            ShapeType::Vertex => {
                let shape = self.shapes[i].shape.clone();
                let p = self.vertex_point_on_shape(&shape);
                let t = self.vertex_tolerance(&shape);
                if let Some(pt) = p {
                    let tol = t.max(1e-10);
                    let b = (pt - DVec3::splat(tol), pt + DVec3::splat(tol), tol);
                    self.shapes[i].bbox = BndBox::from_corners(
                        b.0.x, b.0.y, b.0.z, b.1.x, b.1.y, b.1.z);
                    self.shapes[i].bbox.set_gap(b.2);
                    Some(b)
                } else { None }
            }
            _ => {
                let mut mn = DVec3::splat(f64::INFINITY);
                let mut mx = DVec3::splat(f64::NEG_INFINITY);
                let mut gap = 0.0f64;
                for &c in &self.shapes[i].sub_shapes.clone() {
                    if c < self.nb_shapes() {
                        if let Some(b) = self.build_bnd_box(c) {
                            mn = mn.min(b.0); mx = mx.max(b.1); gap = gap.max(b.2);
                        }
                    }
                }
                if mn.x.is_finite() {
                    self.shapes[i].bbox = BndBox::from_corners(
                        mn.x, mn.y, mn.z, mx.x, mx.y, mx.z);
                    self.shapes[i].bbox.set_gap(gap);
                    Some((mn, mx, gap))
                } else { None }
            }
        }
    }

    pub fn vertex_tolerance(&self, s: &Shape) -> f64 {
        s.as_vertex().map_or(0.0, |vd| vd.tolerance)
    }
    pub fn vertex_point_on_shape(&self, s: &Shape) -> Option<DVec3> {
        s.as_vertex().map(|vd| vd.point)
    }

    // ================================================================    // BRep_Tool-style query helpers
    // ================================================================
    /// Edge curve by shape index.
    pub fn edge_curve(&self, i: usize) -> Option<&rcad_kernel::geom::Curve3> {
        self.shapes.get(i).and_then(|si| {
            if si.shape_type != ShapeType::Edge { return None; }
            si.shape.as_edge().and_then(|e| e.curve.as_ref())
        })
    }

    /// Face tolerance by shape index.
    pub fn face_tolerance(&self, i: usize) -> f64 {
        self.shapes.get(i).and_then(|si| {
            si.shape.as_face().map(|f| f.tolerance)
        }).unwrap_or(0.0)
    }

    /// Edge tolerance by shape index.
    pub fn edge_tolerance(&self, i: usize) -> f64 {
        self.shapes.get(i).and_then(|si| {
            si.shape.as_edge().map(|e| e.tolerance)
        }).unwrap_or(0.0)
    }

    /// Edge parameter range by shape index.
    pub fn edge_range(&self, i: usize) -> [f64; 2] {
        self.shapes.get(i).and_then(|si| {
            si.shape.as_edge().map(|e| e.range)
        }).unwrap_or([0.0, 0.0])
    }

    /// Face surface by shape index.
    pub fn face_surface(&self, i: usize) -> Option<rcad_kernel::geom::Surface3> {
        self.shapes.get(i).and_then(|si| {
            if si.shape_type != ShapeType::Face { return None; }
            si.shape.as_face().and_then(|f| f.surface.clone())
        })
    }

    /// True if the edge is degenerate (its 3D curve collapses to a point).
    /// OCCT BRep_Tool::Degenerated — NOT the same as a closed edge (a circle's
    /// first/last vertices are the same shape but the edge is not degenerate).
    pub fn is_edge_degenerated(&self, i: usize) -> bool {
        self.shapes.get(i).and_then(|si| {
            if si.shape_type != ShapeType::Edge { return None; }
            si.shape.as_edge().map(|e| e.degenerated)
        }).unwrap_or(false)
    }

    /// True if the vertex has internal flag set.
    pub fn vertex_is_internal(&self, i: usize) -> bool {
        self.shapes.get(i).and_then(|si| {
            if si.shape_type != ShapeType::Vertex { return None; }
            si.shape.as_vertex().map(|v| v.flags != 0)
        }).unwrap_or(false)
    }

    /// Natural restriction (infinite face bounds).
    pub fn face_natural_restriction(&self, i: usize) -> bool {
        // Returns true if the face is bounded (has explicit wire bounds)
        // rather than being an infinite face (no bounds = natural restriction false).
        self.shapes.get(i).and_then(|si| {
            if si.shape_type != ShapeType::Face { return None; }
            Some(!si.sub_shapes.is_empty())
        }).unwrap_or(false)
    }

    /// Push a vertex into the DS (BRep_Builder equivalent).
    pub fn push_vertex(&mut self, pt: glam::DVec3, tol: f64) -> usize {
        let vd = rcad_kernel::topods::TVertexData {
            my_shapes: Vec::new(), flags: 0,
            point: pt, tolerance: tol, points: Vec::new(),
        };
        let s = Shape::new(Arc::new(TShape::Vertex(vd)), 0, topods::Orientation::Forward);
        self.append_shape(s)
    }

    /// Push an edge into the DS (BRep_Builder equivalent).
    pub fn push_edge(&mut self, curve: rcad_kernel::geom::Curve3, range: [f64; 2],
        first: usize, last: usize) -> usize {
        let empty_vertex = Shape::new(
            Arc::new(TShape::Vertex(empty_vertex_data())), 0, Orientation::Forward,
        );
        let v_first = self.shapes.get(first)
            .map(|s| Shape::new(s.shape.data.clone(), 0, Orientation::Forward))
            .unwrap_or_else(|| empty_vertex.clone());
        let v_last = self.shapes.get(last)
            .map(|s| Shape::new(s.shape.data.clone(), 0, Orientation::Forward))
            .unwrap_or(empty_vertex);
        let ed = rcad_kernel::topods::TEdgeData {
            curve: Some(curve.clone()), range,
            first: v_first, last: v_last,
            tolerance: 0.0, same_parameter: true, same_range: true,
            degenerated: false, pcurves: HashMap::new(),
            representations: Vec::new(), vertex_params: HashMap::new(),
            my_shapes: Vec::new(), flags: 0,
        };
        let s = Shape::new(Arc::new(TShape::Edge(ed)), 0, topods::Orientation::Forward);
        self.append_shape(s)
    }

    /// Push a wire into the DS (BRep_Builder equivalent).
    pub fn push_wire(&mut self, edges: Vec<(usize, topods::Orientation)>) -> usize {
        let edge_shapes: Vec<Shape> = edges.iter().map(|&(ei, orient)| {
            self.shapes.get(ei).map(|si| {
                Shape::new(si.shape.data.clone(), 0, orient)
            }).unwrap_or_else(|| {
                Shape::new(Arc::new(TShape::Vertex(empty_vertex_data())), 0, orient)
            })
        }).collect();
        let wd = rcad_kernel::topods::TWireData { edges: edge_shapes, my_shapes: Vec::new(), flags: 0 };
        let s = Shape::new(Arc::new(TShape::Wire(wd)), 0, topods::Orientation::Forward);
        self.append_shape(s)
    }

    /// Push a face into the DS (BRep_Builder equivalent).
    pub fn push_face(&mut self, surface: rcad_kernel::geom::Surface3,
        outer_wire: usize, inner_wires: Vec<usize>, natural_restriction: bool) -> usize {
        let empty_vertex = Shape::new(
            Arc::new(TShape::Vertex(empty_vertex_data())), 0, Orientation::Forward,
        );
        let ow = self.shapes.get(outer_wire)
            .map(|si| Shape::new(si.shape.data.clone(), 0, Orientation::Forward))
            .unwrap_or_else(|| empty_vertex.clone());
        let iw: Vec<Shape> = inner_wires.iter().filter_map(|&wi| {
            self.shapes.get(wi).map(|si| Shape::new(si.shape.data.clone(), 0, Orientation::Forward))
        }).collect();
        let fd = rcad_kernel::topods::TFaceData {
            surface: Some(surface.clone()),
            outer_wire: ow, inner_wires: iw,
            tolerance: 0.0, natural_restriction,
            sample_point: None, uv_domain: None,
            internal_vertices: Vec::new(),
            surface_location: 0,
            my_shapes: Vec::new(), flags: 0,
        };
        let s = Shape::new(Arc::new(TShape::Face(fd)), 0, Orientation::Forward);
        self.append_shape(s)
    }

    /// Push a wire into the DS (BRep_Builder equivalent).
    pub fn push_wire_edges(&mut self, edges: Vec<(usize, Orientation)>) -> usize {
        let edge_shapes: Vec<Shape> = edges.iter().map(|&(ei, orient)| {
            self.shapes.get(ei)
                .map(|si| Shape::new(si.shape.data.clone(), 0, orient))
                .unwrap_or_else(|| Shape::new(
                    Arc::new(TShape::Vertex(empty_vertex_data())), 0, orient))
        }).collect();
        let wd = rcad_kernel::topods::TWireData {
            edges: edge_shapes, my_shapes: Vec::new(), flags: 0,
        };
        let s = Shape::new(Arc::new(TShape::Wire(wd)), 0, Orientation::Forward);
        self.append_shape(s)
    }

    /// Source face index for an image face.
    pub fn source_face_idx(&self, i: usize) -> usize {
        if i < self.nb_source_shapes { i } else { 0 }
    }

    /// Vertex tolerance by shape index.
    pub fn vertex_tolerance_by_idx(&self, i: usize) -> f64 {
        self.shapes.get(i).and_then(|si| {
            if si.shape_type != ShapeType::Vertex { return None; }
            si.shape.as_vertex().map(|v| v.tolerance)
        }).unwrap_or(0.0)
    }

    /// Vertex point (alias for vertex_point_by_idx).
    pub fn vertex_point(&self, i: usize) -> glam::DVec3 { self.vertex_point_by_idx(i) }

    /// Vertex point by shape index (returns DVec3::ZERO if not a vertex).
    pub fn vertex_point_by_idx(&self, i: usize) -> glam::DVec3 {
        self.shapes.get(i).and_then(|si| {
            if si.shape_type != ShapeType::Vertex { return None; }
            si.shape.as_vertex().map(|v| v.point)
        }).unwrap_or(glam::DVec3::ZERO)
    }

    /// Face count (number of shapes with type Face).
    pub fn face_count_by_type(&self) -> usize {
        self.shapes.iter().filter(|s| s.shape_type == ShapeType::Face).count()
    }

    /// Shape index of the fi-th face in the shapes array.
    pub fn face_shape_idx(&self, fi: usize) -> usize {
        let mut n = 0;
        for i in 0..self.nb_shapes() {
            if self.shapes[i].shape_type == ShapeType::Face {
                if n == fi { return i; }
                n += 1;
            }
        }
        0
    }

    /// Edge shape index of the ei-th edge.
    pub fn edge_shape_idx(&self, ei: usize) -> usize {
        let mut n = 0;
        for i in 0..self.nb_shapes() {
            if self.shapes[i].shape_type == ShapeType::Edge {
                if n == ei { return i; }
                n += 1;
            }
        }
        0
    }

    /// Vertex shape index of the vi-th vertex.
    pub fn vertex_shape_idx(&self, vi: usize) -> usize {
        let mut n = 0;
        for i in 0..self.nb_shapes() {
            if self.shapes[i].shape_type == ShapeType::Vertex {
                if n == vi { return i; }
                n += 1;
            }
        }
        0
    }

    /// DS index of a vertex Shape, by pointer identity (OCCT TopoDS::IsSame —
    /// the vertex is shared iff the two edges reference the same TShape).
    fn vertex_ds_index(&self, v: &Shape) -> usize {
        let pid = v.ptr_id();
        for (i, si) in self.shapes.iter().enumerate() {
            if si.shape_type == ShapeType::Vertex && si.shape.ptr_id() == pid {
                return i;
            }
        }
        usize::MAX
    }

    /// Edge start vertex DS index, by edge shape index (same convention as
    /// edge_curve / edge_tolerance / edge_pave_blocks).
    pub fn edge_start_vertex_ds(&self, ei: usize) -> usize {
        self.shapes.get(ei).and_then(|si| {
            if si.shape_type != ShapeType::Edge { return None; }
            si.shape.as_edge().map(|ed| self.vertex_ds_index(&ed.first))
        }).unwrap_or(usize::MAX)
    }

    /// Edge end vertex DS index, by edge shape index.
    pub fn edge_end_vertex_ds(&self, ei: usize) -> usize {
        self.shapes.get(ei).and_then(|si| {
            if si.shape_type != ShapeType::Edge { return None; }
            si.shape.as_edge().map(|ed| self.vertex_ds_index(&ed.last))
        }).unwrap_or(usize::MAX)
    }

    /// Same as pave_blocks — OCCT: myPaveBlocksMap(theIndex)
    pub fn edge_pave_blocks(&self, ei: usize) -> &[SharedPB] {
        self.pave_blocks(ei)
    }

    /// True if the edge is a geometric (non-degenerate) edge.
    pub fn edge_is_geometric(&self, i: usize) -> bool {
        !self.is_edge_degenerated(i)
    }

    /// Boundary edges of a face.
    pub fn face_boundary_edges(&self, fi: usize) -> Vec<usize> {
        // `fi` is a flat shape index (like face_surface).  The face's sub_shapes
        // are its wires; each wire stores its ordered edges as TShape::Wire::edges
        // (Shape refs mapped via map_shape_index).
        if fi >= self.nb_shapes() || self.shapes[fi].shape_type != ShapeType::Face {
            return vec![];
        }
        let mut out = Vec::new();
        for &ss in &self.shapes[fi].sub_shapes {
            if ss >= self.nb_shapes() {
                continue;
            }
            match self.shapes[ss].shape_type {
                ShapeType::Edge => out.push(ss),
                ShapeType::Wire => {
                    if let TShape::Wire(w) = &*self.shapes[ss].shape.data {
                        for eshape in &w.edges {
                            if let Some(&e) = self
                                .map_shape_index
                                .get(&(eshape.ptr_id(), eshape.location))
                            {
                                if e < self.nb_shapes()
                                    && self.shapes[e].shape_type == ShapeType::Edge
                                {
                                    out.push(e);
                                }
                            }
                        }
                    }
                }
                _ => {}
            }
        }
        out
    }

    /// UV boundary of a face.
    pub fn face_uv_boundary(&self, fi: usize) -> [f64; 4] {
        // OCCT BRepAdaptor_Surface(aF) — the face surface's natural UV bounds
        // (Geom_Surface::FirstUParameter/LastUParameter/FirstVParameter/
        // LastVParameter). A Geom_Plane is unbounded, a Geom_Sphere is
        // [0,2*PI]x[-PI/2,PI/2], a Geom_Cylinder is [0,2*PI]x(-INF,+INF).
        use rcad_kernel::geom::SurfaceEval;
        if let Some(surf) = self.face_surface(fi) {
            let d = surf.default_domain();
            [d[0], d[1], d[2], d[3]]
        } else {
            [0.0, 1.0, 0.0, 1.0]
        }
    }

    /// OCCT BRep_Tool::UVBounds — the face's actual UV bounds computed by
    /// sampling the boundary edges' pcurves. rcad faces build pcurves
    /// incrementally (MakePCurves runs after VF/EF/FF), so the boundary edges'
    /// 3D curves are projected onto the face surface instead (each sample point
    /// is on the surface, so projection recovers the UV parameter).
    pub fn face_actual_uv_bounds(&self, fi: usize) -> [f64; 4] {
        let Some(surf) = self.face_surface(fi) else {
            return [f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY];
        };
        let face_data = match &*self.shapes[fi].shape.data {
            TShape::Face(fd) => fd,
            _ => return [f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY],
        };
        // OCCT IntTools_Context::UVBounds L1029-1040: for a natural-restriction
        // face BRepAdaptor_Surface returns the surface natural domain, NOT the
        // boundary-sampled rect.  Sampling the boundary of a closed face (e.g. a
        // sphere's seam) yields a degenerate V range and breaks the FF domain
        // classification.
        //
        // rcad note: the rcad primitives mark the lateral cylinder face as
        // natural-restriction even though it is trimmed by its caps, so the
        // override is restricted to genuinely closed surfaces (sphere/torus),
        // whose boundary never restricts the UV domain.
        if face_data.natural_restriction
            && matches!(surf, Surface3::Sphere(_) | Surface3::Torus(_))
        {
            use rcad_kernel::geom::SurfaceEval;
            return surf.default_domain();
        }
        let mut umin = f64::INFINITY;
        let mut umax = f64::NEG_INFINITY;
        let mut vmin = f64::INFINITY;
        let mut vmax = f64::NEG_INFINITY;
        let wire_shapes: Vec<Shape> = std::iter::once(face_data.outer_wire.clone())
            .chain(face_data.inner_wires.iter().cloned())
            .collect();
        let mut any = false;
        for ws in &wire_shapes {
            let Some(&wi) = self.map_shape_index.get(&(ws.ptr_id(), ws.location)) else { continue };
            if wi >= self.nb_shapes() || self.shapes[wi].shape_type != ShapeType::Wire {
                continue;
            }
            let wire_edge_shapes = match &*self.shapes[wi].shape.data {
                TShape::Wire(w) => w.edges.clone(),
                _ => Vec::new(),
            };
            for eshape in &wire_edge_shapes {
                let Some(&ei) = self.map_shape_index.get(&(eshape.ptr_id(), eshape.location)) else { continue };
                if ei >= self.nb_shapes() {
                    continue;
                }
                let edge_data = match &*self.shapes[ei].shape.data {
                    TShape::Edge(ed) => ed,
                    _ => continue,
                };
                let Some(c3d) = edge_data.curve.clone() else { continue };
                let t0 = edge_data.range[0];
                let t1 = edge_data.range[1];
                const NS: usize = 8;
                for k in 0..=NS {
                    let t = t0 + (t1 - t0) * (k as f64 / NS as f64);
                    let p = c3d.point_at(t);
                    let (uv, _) = crate::bop::closest_point_on_surface(&surf, p);
                    umin = umin.min(uv.x);
                    umax = umax.max(uv.x);
                    vmin = vmin.min(uv.y);
                    vmax = vmax.max(uv.y);
                    any = true;
                }
            }
        }
        if !any {
            return [f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY];
        }
        [umin, umax, vmin, vmax]
    }

    /// Inner boundary of a face (wire inner loops).
    pub fn face_inner_boundary(&self, _fi: usize) -> Vec<Vec<usize>> {
        Vec::new()
    }

    /// Boundary vertices of a face (sub-shapes that are vertices).
    pub fn face_boundary_verts(&self, fi: usize) -> Vec<usize> {
        let si = self.face_shape_idx(fi);
        if si < self.nb_shapes() {
            self.shapes[si].sub_shapes.iter().filter(|&&ss| {
                ss < self.nb_shapes() && self.shapes[ss].shape_type == ShapeType::Vertex
            }).copied().collect()
        } else {
            Vec::new()
        }
    }
}

// ===
// Free function: extract sub-shapes of a Shape (TopExp_Explorer equivalent)
// ===
fn sub_shapes_of(s: &Shape) -> Vec<Shape> {
    // Preserve original BRep index so edge_vertex_params can look up vertex_params.
    let cp = |sr: &Shape| Shape::from_parts(sr.data.clone(), sr.index, sr.location, sr.orientation);
    match &*s.data {
        TShape::Vertex(_) => vec![],
        TShape::Edge(ed) => vec![cp(&ed.first), cp(&ed.last)],
        TShape::Wire(wd) => wd.edges.iter().map(cp).collect(),
        TShape::Face(fd) => {
            let mut v = vec![cp(&fd.outer_wire)];
            v.extend(fd.inner_wires.iter().map(cp));
            v
        }
        TShape::Shell(sd) => sd.faces.iter().map(cp).collect(),
        TShape::Solid(sd) => sd.shells.iter().map(cp).collect(),
        TShape::CompSolid(cd) => cd.iter().map(cp).collect(),
        TShape::Compound(cd) => cd.iter().map(cp).collect(),
    }
}

/// ShapeSource (topalgo) for the boolean DS — lets the BRepClass/BRepClass3d
/// classifiers read shape data without depending on the DS type (OCCT:
/// BRepClass uses BRepAdaptor, never BOPDS).
impl crate::topalgo::shape_source::ShapeSource for DS {
    fn nb_shapes(&self) -> usize {
        self.shapes.len()
    }
    fn shape_at(&self, i: usize) -> Shape {
        self.shapes[i].shape.clone()
    }
    fn shape_type(&self, i: usize) -> ShapeType {
        self.shapes[i].shape_type
    }
    fn sub_shapes(&self, i: usize) -> &[usize] {
        &self.shapes[i].sub_shapes
    }
    fn map_shape_index(&self, ptr_id: u64, location: u32) -> Option<usize> {
        self.map_shape_index.get(&(ptr_id, location)).copied()
    }
    fn map_ve(&self, vertex: usize) -> Option<&Vec<usize>> {
        self.map_ve.get(&vertex)
    }
    fn face_surface(&self, i: usize) -> Option<Surface3> {
        DS::face_surface(self, i)
    }
    fn vertex_tolerance(&self, i: usize) -> f64 {
        DS::vertex_tolerance_by_idx(self, i)
    }
    fn is_edge_degenerated(&self, i: usize) -> bool {
        DS::is_edge_degenerated(self, i)
    }
}

impl Default for DS { fn default() -> Self { Self::new() } }







