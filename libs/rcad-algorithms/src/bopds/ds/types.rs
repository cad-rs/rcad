use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, Line2d, Line3, Plane, Surface3, any_perpendicular};
use rcad_kernel::topods;

use crate::bopds::common_block::CommonBlock;
use crate::bopds::face_info::FaceInfo;
use crate::bopds::pave::{Pave, PaveBlock, NO_EDGE};
use crate::tolerance::*;
use std::collections::HashMap;

/// Identifies which input shape a sub-shape came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShapeOrigin {
    ShapeA,
    ShapeB,
}

/// 鉁?OCCT-aligned: BOPDS_PassKey 鈥?sorted (index1, index2) pair key.
/// OCCT BOPDS_PassKey.hxx 鈥?wraps two integers with index1 <= index2.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct PassKey {
    pub i1: usize,
    pub i2: usize,
}

/// 鉁?BOPTools_ConnexityBlock 鈥?connected component with IsRegular flag.
/// OCCT BOPTools_ConnexityBlock.hxx
#[derive(Debug, Clone)]
pub struct ConnexityBlock {
    pub shapes: Vec<usize>,
    pub is_regular: bool,
    pub loops: Vec<Vec<usize>>,
}

impl ConnexityBlock {
    pub fn new() -> Self { ConnexityBlock { shapes: Vec::new(), is_regular: false, loops: Vec::new() } }
    pub fn is_regular(&self) -> bool { self.is_regular }
    pub fn set_regular(&mut self, r: bool) { self.is_regular = r; }
    pub fn shapes(&self) -> &[usize] { &self.shapes }
    pub fn add_shape(&mut self, s: usize) { self.shapes.push(s); }
    pub fn loops(&self) -> &[Vec<usize>] { &self.loops }
    pub fn change_shapes(&mut self) -> &mut Vec<usize> { &mut self.shapes }
    pub fn change_loops(&mut self) -> &mut Vec<Vec<usize>> { &mut self.loops }
}

impl PassKey {
    pub fn new(a: usize, b: usize) -> Self {
        if a <= b { PassKey { i1: a, i2: b } } else { PassKey { i1: b, i2: a } }
    }
}

/// 鉁?OCCT-aligned: lightweight pair iterator over shape indices.
/// OCCT BOPDS_Iterator 鈥?produces sorted (i,j) pairs, optionally cross-group (A脳B).
pub struct PairIterator {
    i: usize, j: usize, a_end: usize, b_end: usize, done: bool,
    cross: bool,  // true = cross-group (A脳B), false = all pairs (0..n)
}

impl PairIterator {
    /// OCCT: BOPDS_Iterator 鈥?iterate all pairs over [0, count).
    pub fn new(count: usize) -> Self {
        PairIterator { i: 0, j: 1, a_end: count, b_end: count, done: count < 2, cross: false }
    }

    /// OCCT: BOPDS_Iterator::Prepare 鈥?iterate cross-group pairs A[end_a] 脳 B[end_b..].
    /// For rcad: A = [0, a_end), B = [a_end, b_end).
    /// This matches the PaveFiller's A脳B cross-shape pair iteration pattern.
    pub fn prepare_ab(a_end: usize, b_end: usize) -> Self {
        let has_pairs = a_end > 0 && b_end > a_end;
        PairIterator { i: 0, j: a_end, a_end, b_end, done: !has_pairs, cross: true }
    }

    pub fn more(&self) -> bool { !self.done }
    pub fn value(&self) -> PassKey { PassKey { i1: self.i, i2: self.j } }

    pub fn next(&mut self) {
        if self.cross {
            self.j += 1;
            if self.j >= self.b_end { self.i += 1; self.j = self.a_end; }
            if self.i >= self.a_end { self.done = true; }
        } else {
            self.j += 1;
            if self.j >= self.b_end { self.i += 1; self.j = self.i + 1; }
            if self.i >= self.b_end - 1 || self.i >= self.a_end { self.done = true; }
        }
    }
}

/// 鉁?OCCT-aligned: BOPDS_ShapeSD 鈥?same-domain shape mappings.
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
            sv.insert((a, b)); sv.insert((b, a));
        }
        for &(a, b) in &shared.shared_edges {
            se.insert((a, b)); se.insert((b, a));
        }
        for &(a, b) in &shared.shared_faces {
            sf.insert((a, b)); sf.insert((b, a));
        }
        ShapeSD { sd_vertices: sv, sd_edges: se, sd_faces: sf }
    }

    /// OCCT: HasSource(sd, src) 鈥?true if sd has a same-domain counterpart src.
    pub fn has_source_vertex(&self, v: usize) -> bool { self.sd_vertices.contains(&(v, usize::MAX)) }
    pub fn has_source_edge(&self, e: usize) -> bool { self.sd_edges.contains(&(e, usize::MAX)) }
    pub fn has_source_face(&self, f: usize) -> bool { self.sd_faces.contains(&(f, usize::MAX)) }

    /// OCCT: IsSubShape(shape) 鈥?true if shape participates in any SD mapping.
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

    /// OCCT: ShapesSD iterator 鈥?(source, same_domain) pairs.
    pub fn sd_vertices_iter(&self) -> impl Iterator<Item = &(usize, usize)> {
        self.sd_vertices.iter()
    }

    /// OCCT-aligned: AddShapeSD 鈥?register a dynamic same-domain vertex pair.
    pub fn add_sd_vertex(&mut self, a: usize, b: usize) {
        self.sd_vertices.insert((a, b));
        self.sd_vertices.insert((b, a));
    }

    /// OCCT-aligned: HasShapeSD(n, nSD) 鈥?find the SD partner for a vertex.
    pub fn find_sd_partner(&self, v: usize) -> Option<usize> {
        self.sd_vertices.iter()
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

/// A vertex in the DS pool.
#[derive(Debug, Clone)]
pub struct DSVertex {
    pub point: DVec3,
    /// None for vertices created at intersections.
    pub origin: Option<ShapeOrigin>,
    /// Model tolerance at this vertex (`vertex_tolerance` from source BRep when loaded;
    /// [`TOLERANCE_ABS`](crate::tolerance::TOLERANCE_ABS) for vertices added by the DS).
    pub geom_tol: f64,
    /// OCCT-aligned: TopAbs_INTERNAL orientation marker.
    ///   True when this vertex is INTERNAL to its source solid volume
    ///   (not on the boundary).  Used by FillInternalShapes.
    pub is_internal: bool,
    /// OCCT-aligned: TopLoc_Location index into DS.locations[]; 0 = identity.
    ///   Populated when loading from topods::BRep with non-identity Location.
    ///   Used by emit_wire_face_topods to create ShapeRefs with correct Location.
    pub location: u32,
}

/// 鉁?OCCT-aligned: edge's pcurve on one face (BRep_CurveRepresentation equivalent).
/// Mirrors OCCT's BRep_TEdge per-face pcurve storage with PCurve/PCurve2.
#[derive(Debug, Clone)]
pub struct DSCurveRepOnFace {
    pub face_idx: usize,
    pub pcurve: Curve2d,
    pub pcurve2: Option<Curve2d>,
    pub pcurve_range: [f64; 2],
    pub start_param: f64,
    pub end_param: f64,
}

/// An edge in the DS pool with curve reference.
#[derive(Debug, Clone)]
pub struct DSEdge {
    /// Index into DS.vertices.
    pub start_vertex: usize,
    pub end_vertex: usize,
    pub curve: Curve3,
    /// Parametric range `[t_start, t_end]` on the curve.
    pub t_range: [f64; 2],
    pub origin: ShapeOrigin,
    /// Model tolerance on the source edge (`edge_tolerance` from BRep when loaded).
    pub geom_tol: f64,
    /// Paves inserted on this edge by intersection passes (unsorted until build_split_edges).
    pub paves: Vec<Pave>,
    /// After `build_split_edges`, the edge is represented by these sub-segments.
    pub pave_blocks: Vec<PaveBlock>,
    /// 鉁?OCCT-aligned: per-face pcurve representations (BRep_CurveRepresentation).
    ///   Populated by DS::build_face_reps() after edges and faces are loaded.
    pub face_reps: Vec<DSCurveRepOnFace>,
    /// 鉁?OCCT-aligned: TopAbs_INTERNAL orientation marker.
    ///   True when this edge is INTERNAL to its source solid volume.
    pub is_internal: bool,
    /// 鉁?OCCT-aligned: BRep_Tool::Parameter(aV, aE) 鈥?per-vertex parameter on this edge's curve.
    ///   Populated during DS loading; for intersection/section edges, set to t_range bounds.
    pub     vertex_params: std::collections::HashMap<usize, f64>,
}

impl DSEdge {
    /// 鉁?OCCT-aligned: BRep_Tool::Parameter(aV, this edge) equivalent.
    pub fn vertex_param(&self, v: usize) -> Option<f64> {
        self.vertex_params.get(&v).copied()
    }
}

/// A wire in the DS pool — first-class entity matching OCCT TopoDS_Wire.
#[derive(Debug, Clone)]
pub struct DSWire {
    /// Edge indices in traversal order (into DS.edges).
    pub edges: Vec<usize>,
}

/// A shell in the DS pool 鈥?first-class entity matching OCCT TopoDS_Shell.
#[derive(Debug, Clone)]
pub struct DSShell {
    /// Face indices forming this shell (into DS.faces).
    pub faces: Vec<usize>,
}

/// A face in the DS pool with surface reference.
#[derive(Debug, Clone)]
pub struct DSFace {
    pub surface: Surface3,
    /// Boundary vertex indices (ordered, into DS.vertices) 鈥?outer wire.
    pub boundary_verts: Vec<usize>,
    /// Boundary edge indices (into DS.edges) 鈥?outer wire.
    pub boundary_edges: Vec<usize>,
    /// FORWARD/REVERSED orientation for each boundary edge, matching OCCT's
    /// edge orientation in the face wire. Same length as boundary_edges.
    pub boundary_edge_forwards: Vec<bool>,
    /// Inner wire edges (TopExp_Explorer iterates outer wire first, then inner wires).
    /// Each entry is one inner wire: Vec<(edge_idx, forward_in_wire)>.
    pub inner_boundary_edges: Vec<Vec<(usize, bool)>>,
    /// Outer wire index into DS.wires (OCCT TopAbs_WIRE reference).
    pub outer_wire_idx: Option<usize>,
    /// Inner wire indices into DS.wires.
    pub inner_wire_idxs: Vec<usize>,
    pub normal: DVec3,
    pub origin: ShapeOrigin,
    pub face_info: FaceInfo,
    /// Original face index within the source BRep's flattened face list.
    pub source_face_idx: usize,
    /// Model tolerance on the source face (`face_tolerance` from BRep when loaded).
    pub geom_tol: f64,
    /// UV-space boundary polygon on this face's surface (populated in Task 3+).
    pub uv_boundary: Option<Vec<DVec2>>,
    /// 鉁?OCCT-aligned: natural_restriction 鈥?true when the face surface has
    ///   natural boundaries (full untrimmed sphere, cylinder, cone, etc.).
    ///   BRep_Tool::NaturalRestriction in OCCT, used by BuilderFace::PerformAreas
    ///   to decide whether an empty wire produces the whole surface face.
    pub natural_restriction: bool,
    /// 鉁?OCCT-aligned: source shell index within the source BRep.
    ///   OCCT: BOPDS_ShapeInfo tracks TopAbs_SHELL hierarchy; each source face
    ///   knows which shell it belongs to via its ShapeInfo parent pointer.
    ///   rcad: shell index (0-based, counting all shells across all solids in
    ///   the source BRep) assigned during load_brep.  Used by
    ///   fill_images_containers_shells to group result faces by source shell
    ///   boundary (OCCT FillImagesContainer preserves source shell structure).
    pub source_shell_idx: Option<usize>,
    /// 鉁?OCCT-aligned: source solid index within the source BRep.
    ///   OCCT: each TopAbs_SOLID in the DS has its own ShapeInfo entry.
    ///   rcad: solid index (0-based within the source BRep's flat solids Vec)
    ///   assigned during load_brep.  Used by fill_images_compounds to
    ///   reconstruct the compound hierarchy in the result BRep.
    pub source_solid_idx: Option<usize>,
    /// 鉁?OCCT-aligned: source compsolid index (0-based). OCCT BOPDS_DS
    ///   tracks TopAbs_COMPSOLID in ShapeInfo hierarchy.  rcad: assigned
    ///   during load_brep when solid belongs to a CompSolid (else None).
    ///   Used by fill_images_containers_compsolid to preserve compsolid
    ///   boundaries in the result (FillImagesContainer, Builder_1.cxx L221-276).
    pub source_compsolid_idx: Option<usize>,
}

/// Record of an intersection between two sub-shapes.
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

/// OCCT-aligned: type-specific interference records replacing the flat Vec<Interference>.
///   OCCT BOPDS_DS stores interferences per-type in separate IndexedDataMaps,
///   which provide O(log n) lookup by shape index and natural pair dedup.
///   These are used by the new TypedInterferences container.
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
}

#[derive(Debug, Clone)]
pub struct InterferenceEE {
    pub e1: usize,
    pub e2: usize,
    pub point: DVec3,
    pub param1: f64,
    pub param2: f64,
    pub new_vertex: usize,
}

#[derive(Debug, Clone)]
pub struct InterferenceVF {
    pub vertex: usize,
    pub face: usize,
}

#[derive(Debug, Clone)]
pub struct InterferenceEF {
    pub edge: usize,
    pub face: usize,
    pub point: DVec3,
    pub edge_param: f64,
    pub new_vertex: usize,
}

/// FF entry: keyed by (Fmin,Fmax) pair with all curves and touch points merged.
#[derive(Debug, Clone)]
pub struct InterferenceFF {
    pub f1: usize,
    pub f2: usize,
    pub curves: Vec<usize>,
    pub points: Vec<usize>,
}

/// An intersection curve from F-F intersection, bounded by vertices.
/// 鉁?OCCT-aligned: BOPDS_Curve (hxx:31-119).
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
    /// OCCT-aligned: tolerance of this section edge (CorrectToleranceOfSE).
    pub geom_tol: f64,
    /// 鉁?OCCT-aligned: BOPDS_Curve::myPaveBlocks (hxx:115).
    ///   Sub-segments of this intersection curve, created by splitting at paves.
    pub pave_blocks: Vec<PaveBlock>,
    /// 鉁?OCCT-aligned: BOPDS_Curve / IntTools_Curve extra fields.
    pub curve_extra: CurveExtra,
}

/// 鉁?OCCT-aligned: IntTools_Curve (tangential_tol) + BOPDS_Curve
///   (techno_vertices, my_box) fields.
#[derive(Debug, Clone)]
pub struct CurveExtra {
    pub tangential_tol: f64,
    pub techno_vertices: Vec<usize>,
    pub my_box: Option<(glam::DVec3, glam::DVec3)>,
}

impl Default for CurveExtra {
    fn default() -> Self {
        CurveExtra { tangential_tol: 0.0, techno_vertices: Vec::new(), my_box: None }
    }
}

impl IntersectionCurve {
    /// 鉁?OCCT-aligned: BOPDS_Curve::InitPaveBlock1 (lxx:85-92).
    /// OCCT only pushes an empty PB to the list. PB vertices are set by
    /// PutPavesOnCurve (ext_paves) -> Update(false) (sub-PBs from ext_paves).
    pub fn init_pave_block1(&mut self) {
        if self.pave_blocks.is_empty() {
            self.pave_blocks.push(PaveBlock::new_curve_block());
        }
    }

    /// 鉁?OCCT-aligned: BOPDS_Curve::ChangePaveBlock1 (lxx:96-100).
    pub fn change_pave_block1(&mut self) -> Option<&mut PaveBlock> {
        self.pave_blocks.first_mut()
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

/// 鉁?OCCT-aligned: BOPDS_ShapeInfo 鈥?per-shape metadata in the flat DS index.
///   Corresponds to one entry in BOPDS_DS::myLines.
#[derive(Debug, Clone)]
pub struct ShapeInfo {
    pub shape_type: rcad_kernel::topods::ShapeType,
    pub sub_shapes: Vec<usize>,
    pub flag: i64,
    pub reference: i64,
    pub has_brep: bool,
    /// Bounding box min corner (Bnd_Box equivalent). None = not computed.
    pub box_min: Option<DVec3>,
    /// Bounding box max corner. None = not computed.
    pub box_max: Option<DVec3>,
    /// Box gap expansion (Bnd_Box::SetGap equivalent).
    pub box_gap: f64,
    /// True for shapes created during intersection (not from source BRep).
    pub is_new: bool,
    /// Operand rank: 0 = ShapeA, 1 = ShapeB.
    pub rank: usize,
    /// Original source index within the type-specific array.
    pub source_idx: usize,
}

impl ShapeInfo {
    pub fn new(shape_type: rcad_kernel::topods::ShapeType) -> Self {
        ShapeInfo {
            shape_type, sub_shapes: Vec::new(), flag: -1, reference: -1, has_brep: true,
            box_min: None, box_max: None, box_gap: 0.0,
            is_new: true, rank: 0, source_idx: 0,
        }
    }
    pub fn has_flag(&self) -> bool { self.flag >= 0 }
    pub fn has_reference(&self) -> bool { self.reference >= 0 }
    pub fn has_sub_shape(&self, idx: usize) -> bool { self.sub_shapes.contains(&idx) }
    pub fn has_brep(&self) -> bool { self.has_brep }

    /// OCCT-aligned: Bnd_Box::IsOut check — true if the two boxes do not overlap.
    ///   Uses box_gap expansion (Bnd_Box.SetGap equivalent).
    pub fn box_is_out(&self, other: &Self) -> bool {
        let (Some(mn1), Some(mx1)) = (self.box_min, self.box_max) else { return false };
        let (Some(mn2), Some(mx2)) = (other.box_min, other.box_max) else { return false };
        let g = self.box_gap + other.box_gap;
        mx1.x + g < mn2.x - g || mx1.y + g < mn2.y - g || mx1.z + g < mn2.z - g
            || mx2.x + g < mn1.x - g || mx2.y + g < mn1.y - g || mx2.z + g < mn1.z - g
    }
}

/// Central data structure (OCCT: BOPDS_DS).
#[derive(Debug)]
pub struct DS {
    pub vertices: Vec<DSVertex>,
    pub edges: Vec<DSEdge>,
    pub wires: Vec<DSWire>,
    pub shells: Vec<DSShell>,
    pub faces: Vec<DSFace>,
    /// OCCT-aligned: type-specific interference vecs (BOPDS_DS myInterfVV/VE/VF/EE/EF/FF).
    /// Replaces the generic Vec<Interference> enum — each variant has its own typed Vec.
    pub interf_vv: Vec<InterferenceVV>,
    pub interf_ve: Vec<InterferenceVE>,
    pub interf_vf: Vec<InterferenceVF>,
    pub interf_ee: Vec<InterferenceEE>,
    pub interf_ef: Vec<InterferenceEF>,
    pub interf_ff: Vec<InterferenceFF>,

    pub intersection_curves: Vec<IntersectionCurve>,
    /// Mapping: intersection curve index -> DSEdge indices created by make_section_edges_from_curve_pbs.
    /// Populated during PaveFiller::make_section_edges_from_curve_pbs.
    /// Used by ds_to_brep to skip ICs already converted to DSEdges (Step 2, A2).
    pub section_edge_refs: Vec<Vec<usize>>,
    /// Fuzzy tolerance used during interference detection.
    ///
    /// Vertices/edges within this distance are considered coincident.
    /// When set to a value larger than `TOLERANCE_ABS`, approximate
    /// near-miss intersections (analogous to OCCT `BOPAlgo_Options::SetFuzzyValue`).
    pub fuzzy_tol: f64,
    /// Number of vertices loaded from shape A (first shape). Shape A DS vertex indices are 0..a_vertex_count.
    pub a_vertex_count: usize,
    /// Number of edges loaded from shape A. Shape A DS edge indices are 0..a_edge_count.
    pub a_edge_count: usize,
    /// Number of faces loaded from shape A. Shape A DS face indices are 0..a_face_count.
    pub a_face_count: usize,
    /// Shared topology information for glue path optimization.
    pub shared_topology: SharedTopologyInfo,
    /// 鉁?OCCT-aligned: BOPDS_ShapeSD 鈥?same-domain shape mapping (built from shared_topology).
    pub shape_sd: ShapeSD,
    /// Pre-computed overlap polygons for same-domain (coplanar) face pairs.
    /// Each entry is (face_a_index, face_b_index, overlap_boundary_in_3d).
    /// Populated during PaveFiller's coplanar analysis, consumed by Builder.
    pub same_domain_overlaps: Vec<(usize, usize, Vec<DVec3>)>,

    /// Common blocks grouping geometrically coincident PaveBlocks
    /// (OCCT: BOPDS_CommonBlock). Populated by the PaveFiller
    /// (`ForceInterfEE`) and consumed by the Builder (`FillSameDomainFaces`).
    pub common_blocks: Vec<CommonBlock>,

    /// Edge image mapping (OCCT: BOPAlgo_Builder::myImages).
    /// Indexed by original edge index, each entry lists sub-edge indices
    /// created by `build_edge_images()`.
    pub my_images: Vec<Vec<usize>>,
    /// Edge origin mapping (OCCT: BOPAlgo_Builder::myOrigins).
    /// Indexed by sub-edge index, value is the original edge index.
    pub my_origins: Vec<usize>,

    /// OCCT FillImagesContainers(WIRE): pre-built edge lists for wires whose
    /// edges were split by the PaveFiller.  Each entry corresponds to one
    /// original wire (flat index across all solids/shells of the source BRep).
    /// None = wire unchanged (no image needed).
    pub wire_images: Vec<Option<Vec<(usize, bool)>>>,

    /// OCCT FillImagesContainers(SHELL): placeholder for shell-level images.
    /// Populated by checking if any wire in the shell has split edges.
    pub shell_images: Vec<bool>,

    /// OCCT FillImagesSolids: placeholder for solid-level images.
    pub solid_images: Vec<bool>,

    /// OCCT-aligned: TopLoc_Location storage. Index 0 = identity (implicit), 1+ stored here.
    /// Populated by load_brep when loading from topods::BRep with non-identity Location.
    pub locations: Vec<glam::DAffine3>,

    /// Global PaveBlock array (OCCT: BOPDS_DS::myPaveBlocks).
    /// Indices in FaceInfo::pave_blocks_on / pave_blocks_in refer to this array.
    pub pave_blocks: Vec<PaveBlock>,
    /// 鉁?OCCT-aligned: myIncreasedSS 鈥?vertices whose tolerance was increased
    ///   during intersection processing.  Read by RepeatIntersection to determine
    ///   which vertices need VV/VE/VF re-checks.
    pub increased_ss: std::collections::HashSet<usize>,
    /// 鉁?OCCT-aligned: BOPDS_DS::myLines 鈥?flat array of BOPDS_ShapeInfo for all shapes.
    ///   Each entry records shape_type, sub_shapes, flag, reference, has_brep.
    ///   First nb_source_shapes entries are original source shapes; entries beyond
    ///   are shapes created during intersection.
    pub shape_info: Vec<ShapeInfo>,
    /// 鉁?OCCT-aligned: myNbSourceShapes 鈥?count of original source shapes.
    pub nb_source_shapes: usize,
}

