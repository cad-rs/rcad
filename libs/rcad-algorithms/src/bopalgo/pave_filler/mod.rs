use std::collections::HashMap;
use std::collections::HashSet;
use std::sync::{Arc, RwLock};

use glam::{DVec2, DVec3};
use rcad_kernel::PCurve;
use rcad_kernel::geom::*;
use rcad_kernel::topods;

use crate::bopalgo::{Alert, GlueEnum, Report};
use crate::bopds::common_block::CommonBlock;
use crate::bopds::ds::face_aabb;
use crate::bopds::ds::{
    DS, DSCurveRepOnFace, DSEdge, DSVertex, Interference, InterferenceEE, InterferenceEF,
    InterferenceFF, InterferenceVE, InterferenceVF, InterferenceVV, IntersectionCurve, ShapeOrigin,
};
use crate::bopds::pave::*;
use crate::boptools::bvh::BoxTree;
use crate::boptools::bvh::Bvh;
use crate::inttools;
use crate::inttools::context::Context as IntToolsContext;
use crate::inttools::fclass2d::{FClass2d, State};
use crate::tolerance::*;
use indexmap::IndexSet;
use rcad_kernel::closest_point_on_curve;
pub mod helpers;
use self::helpers::*;

// =IntPatch_Intersection surface category (L1264-1294).
// GeomGeom = ts1==ts2==1 (both analytic, ImpImpIntersection)
// ParamParam = ts1==ts2==0 (both parametric, PrmPrmIntersection)

// Re-export NearTangentType from bopds::ds for use in this module's public types
pub use crate::bopds::ds::NearTangentType;

/// Minimum total face count before BVH acceleration is used.
/// Below this threshold, brute-force O(n ? is faster due to BVH build overhead.
const BVH_THRESHOLD: usize = 20;

pub(crate) mod analytics;
mod config;
mod ff_intersect;
///  BOPAlgo_PaveFiller =six intersection passes
/// (PaveFiller.hxx L106-107, PaveFiller.cxx L234-355).
mod glue;
mod interf;
mod intersection;
mod make_blocks;
pub(crate) mod marching;
pub(crate) mod p_walking;
mod paves;
pub(crate) mod polyhedron;
pub(crate) mod polyhedron_bvh;
mod posttreat;
pub(crate) mod prm_prm_intersection;
mod tolerances;

/// BOPAlgo_SectionAttribute =controls approximation and
/// pcurve computation for section edges (BOPAlgo_SectionAttribute.hxx).
#[derive(Debug, Clone)]
pub(crate) struct SectionAttribute {
    pub approximation: bool,
    pub pcurve_on_s1: bool,
    pub pcurve_on_s2: bool,
}

impl Default for SectionAttribute {
    fn default() -> Self {
        Self {
            approximation: true,
            pcurve_on_s1: true,
            pcurve_on_s2: true,
        }
    }
}

/// PaveFiller::EdgeRangeDistance =stores minimal distance
/// between an edge range and a face that don't geometrically intersect.
/// Used by PostTreatFF to re-check E-F pairs after tolerance updates.
#[derive(Debug, Clone)]
pub(crate) struct EdgeRangeDistance {
    pub first: f64,
    pub last: f64,
    pub distance: f64,
}

/// OCCT BOPAlgo_PaveFiller::UpdateExistingPaveBlocks (PaveFiller_6.cxx L3311-3529).
/// Replaces old pave blocks (from aPBf / its CommonBlock) with new ones from
/// PostTreatFF (aLPB). Handles CommonBlock splitting and re-projects edges to
/// faces from thePBFacesMap.
fn update_existing_pave_blocks(
    ds: &mut DS,
    context: &mut IntToolsContext,
    a_pbf: usize,
    a_lpb: &[usize],
    the_pb_faces_map: &HashMap<usize, Vec<usize>>,
    fuzzy_value: f64,
) {
    if a_lpb.is_empty() {
        return;
    }

    // 1. Determine the set of old PBs (from aPBf's CB, or just aPBf)
    let b_cb: bool;
    // Preserve (pb_idx, face_idx) pairs from the CommonBlock.
    let old_pbs: Vec<(usize, usize)>;
    {
        let spb = &ds.pave_blocks[a_pbf];
        let pb = spb.0.read().unwrap();
        b_cb = pb.common_block_idx.is_some();
        if let Some(cb_idx) = pb.common_block_idx {
            if let Some(cb) = ds.common_blocks.get(cb_idx) {
                old_pbs = cb.pave_blocks().to_vec();
            } else {
                old_pbs = vec![(a_pbf, 0)];
            }
        } else {
            old_pbs = vec![(a_pbf, 0)];
        }
    }

    // 2. Remove old PBs from edge PB lists
    for &(old_pb_idx, _) in &old_pbs {
        let orig_e = {
            let spb = &ds.pave_blocks[old_pb_idx];
            spb.0.read().unwrap().original_edge
        };
        if orig_e < ds.edge_count() {
            ds.edges[orig_e].pave_blocks.retain(|spb| {
                let ptr = Arc::as_ptr(&spb.0);
                let old_ptr = Arc::as_ptr(&ds.pave_blocks[old_pb_idx].0);
                ptr != old_ptr
            });
        }
    }

    // 3. If CB: create new CommonBlocks per replacement PB
    //    If not CB: just append new PBs to the edge
    let mut new_pb_list: Vec<usize> = Vec::new();

    if b_cb {
        let orig_faces: Vec<usize> = {
            let spb = &ds.pave_blocks[a_pbf];
            let pb = spb.0.read().unwrap();
            if let Some(cb_idx) = pb.common_block_idx {
                if let Some(cb) = ds.common_blocks.get(cb_idx) {
                    cb.faces().to_vec()
                } else {
                    Vec::new()
                }
            } else {
                Vec::new()
            }
        };

        for &rp_idx in a_lpb {
            let mut a_cb = crate::bopds::common_block::CommonBlock::new();
            let rp_pave1;
            let rp_pave2;
            {
                let rp_pb = ds.pave_blocks[rp_idx].0.read().unwrap();
                rp_pave1 = rp_pb.pave1.clone();
                rp_pave2 = rp_pb.pave2.clone();
            }

            for &(old_pb_idx, old_face_idx) in &old_pbs {
                let orig_e = {
                    let spb = &ds.pave_blocks[old_pb_idx];
                    spb.0.read().unwrap().original_edge
                };
                let mut pb2n =
                    crate::bopds::pave::PaveBlock::new(orig_e, rp_pave1.clone(), rp_pave2.clone());
                pb2n.new_edge = {
                    let rp_pb = ds.pave_blocks[rp_idx].0.read().unwrap();
                    rp_pb.new_edge
                };
                // Register in DS
                let new_g_idx = ds.pave_blocks.len();
                ds.pave_blocks.push(crate::bopds::pave::SharedPB::new(pb2n));
                a_cb.add_pave_block(new_g_idx, old_face_idx);
                // Append to edge's PB list
                if orig_e < ds.edge_count() {
                    ds.edges[orig_e]
                        .pave_blocks
                        .push(ds.pave_blocks[new_g_idx].clone());
                }
            }
            a_cb.set_faces(orig_faces.clone());
            let cb_idx = ds.common_blocks.len();
            // Capture first PB before a_cb is moved into common_blocks
            let first_pb = a_cb.pave_blocks().first().copied().map(|(pbi, _)| pbi);
            // Set common_block_idx on all PBs in this CB
            for &(pbi, _) in a_cb.pave_blocks() {
                if pbi < ds.pave_blocks.len() {
                    ds.pave_blocks[pbi].0.write().unwrap().common_block_idx = Some(cb_idx);
                }
            }
            ds.common_blocks.push(a_cb);
            if let Some(fp) = first_pb {
                new_pb_list.push(fp);
            }
        }
    } else {
        let orig_e = {
            let spb = &ds.pave_blocks[a_pbf];
            spb.0.read().unwrap().original_edge
        };
        if orig_e < ds.edge_count() {
            for &rp_idx in a_lpb {
                ds.edges[orig_e]
                    .pave_blocks
                    .push(ds.pave_blocks[rp_idx].clone());
                new_pb_list.push(rp_idx);
            }
        }
    }

    // 4. Project replacement PBs to faces in thePBFacesMap
    // OCCT L3481-3528: IntTools_EdgeFace check for each replacement PB
    if let Some(pb_faces) = the_pb_faces_map.get(&a_pbf) {
        let a_tol_f = crate::tolerance::CONFUSION;
        for &fi in pb_faces {
            if fi >= ds.face_count() {
                continue;
            }
            let face_surface = ds.faces[fi].surface.clone();
            for &new_pb_idx in &new_pb_list {
                if new_pb_idx >= ds.pave_blocks.len() {
                    continue;
                }
                // Check if already registered in this face's ON/IN sets (OCCT L3498)
                {
                    let face_info = &ds.faces[fi].face_info;
                    if face_info.pave_blocks_on.contains(&new_pb_idx)
                        || face_info.pave_blocks_in.contains(&new_pb_idx)
                    {
                        continue;
                    }
                }
                // OCCT L3503-3511: IntTools_EdgeFace check
                let (ei, t1, t2) = {
                    let spb = &ds.pave_blocks[new_pb_idx];
                    let pb = spb.0.read().unwrap();
                    let e_idx = pb.new_edge.unwrap_or(pb.original_edge);
                    let (tt1, tt2) = pb.range();
                    (e_idx, tt1, tt2)
                };
                if ei >= ds.edge_count() {
                    continue;
                }
                let a_tol_e = ds.edge_tolerance(ei).max(CONFUSION);
                // OCCT L3514: bCoincide = (aCPrts.Length() == 1 && aCPrts(1).Type() == TopAbs_EDGE)
                let b_coincide = {
                    crate::inttools::edge_face::is_coincident_edge_face(
                        &ds.edges[ei].curve,
                        [t1, t2],
                        &face_surface,
                        a_tol_f,
                        a_tol_e,
                        context,
                        ds,
                        fi,
                    )
                };
                if b_coincide {
                    // OCCT L3517-3524: create/find CommonBlock + AddFace
                    let spb = &ds.pave_blocks[new_pb_idx];
                    let pb_ptr = Arc::as_ptr(&spb.0);
                    let mut found_cb = false;
                    for cb in &mut ds.common_blocks {
                        if cb.pave_blocks().iter().any(|&(pbi, _)| {
                            pbi < ds.pave_blocks.len()
                                && Arc::as_ptr(&ds.pave_blocks[pbi].0) == pb_ptr
                        }) {
                            cb.add_face(fi);
                            found_cb = true;
                            break;
                        }
                    }
                    if !found_cb {
                        // OCCT L3520-3522: create new CommonBlock
                        let mut new_cb = crate::bopds::common_block::CommonBlock::new();
                        new_cb.add_pave_block(new_pb_idx, fi);
                        new_cb.add_face(fi);
                        let cb_idx = ds.common_blocks.len();
                        ds.pave_blocks[new_pb_idx]
                            .0
                            .write()
                            .unwrap()
                            .common_block_idx = Some(cb_idx);
                        ds.common_blocks.push(new_cb);
                    }
                    // OCCT L3525: aFI.ChangePaveBlocksIn().Add(aPB)
                    ds.face_info_mut(fi).pave_blocks_in.insert(new_pb_idx);
                }
            }
        }
    }
}

pub struct PaveFiller<'a> {
    pub ds: &'a mut DS,
    /// =myIterator (BOPAlgo_PaveFiller.hxx) — BOPDS_Iterator for pair enumeration.
    /// Uses 'static lifetime via unsafe transmute (see config.rs for safety).
    pub(crate) my_iterator: crate::bopds::ds::BOPDS_Iterator<'static>,
    bvh_a: Option<&'a Bvh>,
    bvh_b: Option<&'a Bvh>,
    /// DS-based face BVH for FF pair detection. Uses DS face indices directly,
    /// matching OCCT's BOPTools_BoxTree which operates on source shape indices.
    pub(crate) glue: GlueEnum,
    glue_tolerance: f64,
    /// convenience  ?true when glue is active (not GlueOff).
    /// =BOPAlgo_Options::SetFuzzyValue
    fuzzy_tolerance: f64,
    /// =BOPAlgo_Algo::myRunParallel
    run_parallel: bool,
    /// =BOPAlgo_PaveFiller::myNonDestructive
    non_destructive: bool,
    /// =BOPAlgo_Algo::myUseOBB
    use_obb: bool,
    /// =IntTools_Context (PaveFiller::Init L203)
    context: IntToolsContext,
    /// =myArguments =original input shapes (BOPAlgo_PaveFiller.hxx L639).
    /// rcad: carries the original BRep operands for OCCT-API compatibility.
    my_arguments: Vec<rcad_kernel::topods::BRep>,
    /// =mySectionAttribute (BOPAlgo_SectionAttribute.hxx)
    section_attribute: SectionAttribute,
    /// =myIsPrimary (BOPAlgo_PaveFiller.cxx L62)
    is_primary: bool,
    /// =myAvoidBuildPCurve (BOPAlgo_PaveFiller.cxx L63)
    avoid_build_pcurve: bool,
    /// =myFPBDone =fence map tracking processed (face, pave_block) pairs.
    /// Map: face_idx =set of pave_block indices already processed in PostTreatFF.
    fpbdone: std::collections::HashMap<usize, std::collections::HashSet<usize>>,
    /// =myVertsToAvoidExtension =vertices that should NOT have
    /// their tolerance extended further (near EE/EF intersection points).
    verts_to_avoid_extension: std::collections::HashSet<usize>,
    /// =myIncreasedSS (PaveFiller.hxx L651) — sub-shapes with increased tolerance.
    /// OCCT: NCollection_Map<int> on PaveFiller, NOT on DS.
    my_increased_ss: std::collections::HashSet<usize>,
    /// =myDistances =minimal edge-face distances for non-intersecting
    /// pairs.  Map: (edge_idx, face_idx) =Vec<EdgeRangeDistance>.
    distances: std::collections::HashMap<(usize, usize), Vec<EdgeRangeDistance>>,
    ///  myReport  ?collects alerts during PaveFiller execution.
    my_report: Report,
    /// Pipeline stage dump context (rcad PF stages). Created from env vars;
    /// disabled when RCAD_DUMP_PIPELINE is not set.
    pub dump_ctx: crate::pipeline_dump::DumpCtx,
    /// Stage to stop after (for stage-by-stage testing). When set, perform()
    /// returns early after completing the named stage.
    pub stop_after: Option<String>,
}

impl<'a> PaveFiller<'a> {
    /// IsGlue  ?true when glue mode is active (not GlueOff).
    pub fn use_glue(&self) -> bool {
        self.glue != GlueEnum::GlueOff
    }

    /// GetGlue  ?return current glue mode.
    pub fn glue_mode(&self) -> GlueEnum {
        self.glue
    }
}

/// GetFullShapeMap (PaveFiller_6.cxx L2941-2958).
/// Builds a set of all sub-shape indices belonging to face `fi`:
/// the face itself, its boundary edges, and their endpoint vertices.
pub(crate) fn build_face_shape_map(ds: &DS, fi: usize) -> std::collections::HashSet<usize> {
    let mut aMI = std::collections::HashSet::new();
    aMI.insert(fi);
    if fi < ds.face_count() {
        for &ei in ds.face_boundary_edges(fi) {
            aMI.insert(ei);
            if ei < ds.edge_count() {
                let v_start = ds.edge_start_vertex_ds(ei);
                let v_end = ds.edge_end_vertex_ds(ei);
                aMI.insert(v_start);
                aMI.insert(v_end);
                // OCCT: resolve SD vertices — if a boundary vertex has an SD
                // partner, include it so IsSubShape checks match.
                if let Some(n_sd) = ds.has_shape_sd(v_start) {
                    aMI.insert(n_sd);
                }
                if let Some(n_sd) = ds.has_shape_sd(v_end) {
                    aMI.insert(n_sd);
                }
            }
        }
    }
    aMI
}

/// =Propagate IC vertices to all faces sharing boundary edges
/// (OCCT BOPDS_FaceInfo::AppendBlock equivalent).
/// OCCT BOPAlgo_PaveFiller propagates pave block vertices to all faces
/// referencing the split edge. rcad's add_curve only adds vertices to the
/// two FF-interference faces (f1, f2), but the vertex may lie on boundary
/// edges of other faces (e.g. side-face tangent-line IC endpoints on the
/// top face's boundary edge).
fn propagate_ic_vertices_to_shared_faces(
    ds: &mut DS,
    ic_vertices: &[usize],
    skip_faces: &[usize; 2],
) {
    let vtol = TOLERANCE_ABS * 1000.0; // 1e-4 geometric tolerance for on-edge check
    let vtol_sq = vtol * vtol;
    for fi in 0..ds.face_count() {
        if fi == skip_faces[0] || fi == skip_faces[1] {
            continue;
        }
        for &vi in ic_vertices {
            if ds.face_info(fi).vertices_in.contains(&vi) {
                continue;
            }
            let vp = ds.vertex_point(vi);
            for &ei in ds.face_boundary_edges(fi) {
                let Some(edge) = ds.edges.get(ei) else {
                    continue;
                };
                let a = ds.vertex_point(edge.start_vertex);
                let b = ds.vertex_point(edge.end_vertex);
                let ab = b - a;
                let ab_len2 = ab.length_squared();
                if ab_len2 < TOLERANCE_LEN_SQ_DIV_SAFE {
                    continue;
                }
                let ap = vp - a;
                let t = ap.dot(ab) / ab_len2;
                if t > -0.01 && t < 1.01 {
                    let proj = a + ab * t.clamp(0.0, 1.0);
                    if (vp - proj).length_squared() < vtol_sq {
                        ds.face_info_mut(fi).vertices_in.insert(vi);
                        break;
                    }
                }
            }
        }
    }
}

impl<'a> PaveFiller<'a> {
    ///  Prepare (PaveFiller_7.cxx L850-929).
    /// Build 2D pcurves for edges on planar faces.
    /// OCCT: iterate all V/E, E/E, E/F pairs to find planar faces, collect edge-face pairs,
    /// compute pcurves in parallel, update edges.  rcad: DS::build_face_reps already computes
    /// pcurves for all face types; this step ensures planar-face pcurves exist for any
    /// edges that build_face_reps may have missed (non-boundary or intersection-relevant).
    // OCCT PaveFiller_7.cxx L850-932
    fn prepare(&mut self) {
        // OCCT L850: void BOPAlgo_PaveFiller::Prepare
        // OCCT L852-856: in non-destructive mode, do not modify original edges
        if self.non_destructive {
            return;
        }
        // OCCT L857-879: Iterate V/E/F vs F via Iterator to find planar faces
        use crate::bopds::ds::BOPDS_Iterator;
        use rcad_kernel::topods::ShapeType;
        // OCCT L858: bool bIsBasedOnPlane
        let mut a_mf: std::collections::HashSet<usize> = std::collections::HashSet::new();
        // OCCT L857: TopAbs_ShapeEnum aType[] = {TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE}
        let a_type = [ShapeType::Vertex, ShapeType::Edge, ShapeType::Face];
        // OCCT L863: aNb = 3
        let a_nb = 3;
        // OCCT L864: Message_ProgressScope aPSOuter (rcad: sequential, no progress)
        for i in 0..a_nb {
            // OCCT L867: myIterator->Initialize(aType[i], aType[2])
            let mut it = BOPDS_Iterator::new(self.ds);
            it.prepare();
            let pairs = it.pairs(a_type[i], ShapeType::Face);
            // OCCT L868-878: iterate pairs
            for &(_n1, n_f) in pairs.iter() {
                if n_f < self.ds.face_count()
     // OCCT L873: IsBasedOnPlane(aF) — checks Geom_Plane with TrimmedSurface unwrap
     // rcad: Surface3 has no TrimmedSurface wrapper; locate_surface applies Location
     && matches!(self.ds.locate_surface(n_f), Surface3::Plane(_))
                {
                    // OCCT L876: aMF.Add(aF)
                    a_mf.insert(n_f);
                }
            }
        }
        // OCCT L881: aNbF = aMF.Extent()
        let a_nb_f = a_mf.len();
        // OCCT L882-884: if (!aNbF) { return; }
        if a_nb_f == 0 {
            return;
        }

        // OCCT L888-901: collect edge-face pairs from planar faces' boundary topology
        // OCCT L888: BOPAlgo_VectorOfBPC aVBPC (rcad: Vec<(usize, usize)> pairs)
        use crate::bopds::ds::DSCurveRepOnFace;
        let mut a_vbpc: Vec<(usize, usize)> = Vec::new();
        // OCCT L890: for (i = 1; i <= aNbF; ++i)
        for i in 0..a_nb_f {
            // OCCT L891: const TopoDS_Face& aF = *(TopoDS_Face*)&aMF(i)
            // Architecture diff: HashSet iteration order is arbitrary (not IndexedMap order)
            let fi = *a_mf.iter().nth(i).unwrap();
            let f = &self.ds.faces[fi];
            // OCCT L893: aExp.Init(aF, aType[1]) — explore EDGE sub-shapes
            // rcad: boundary_edges + inner_boundary_edges (same semantics)
            for &ei in &f.boundary_edges {
                if ei < self.ds.edge_count() && self.ds.edge_on_face(ei, fi).is_none() {
                    // OCCT L898: aBPC.SetEdge(aE); aBPC.SetFace(aF)
                    a_vbpc.push((ei, fi));
                }
            }
            for w in &f.inner_boundary_edges {
                for &(ei, _) in w {
                    if ei < self.ds.edge_count() && self.ds.edge_on_face(ei, fi).is_none() {
                        a_vbpc.push((ei, fi));
                    }
                }
            }
        }
        // OCCT L903-910: Build pcurves (BOPTools_Parallel)
        // OCCT L904-908: prepare BPC list  (rcad: sequential, compute inline)
        // OCCT L909: BOPTools_Parallel::Perform (rcad: sequential loop)
        let mut pcurve_results: Vec<(usize, usize, rcad_kernel::geom::Curve2d, f64)> = Vec::new();
        for &(ei, fi) in &a_vbpc {
            let surf = self.ds.faces[fi].surface.clone();
            if let Some(edge) = self.ds.edges.get(ei) {
                if let Some((pcurve, span)) = DS::compute_edge_pcurve(&edge.curve, &surf, None) {
                    pcurve_results.push((ei, fi, pcurve, span));
                }
            }
        }
        // OCCT L916-931: Update edges with pcurves
        // OCCT L917: BRep_Builder aBB (rcad: direct DS update)
        for &(ei, fi, ref pcurve, span) in &pcurve_results {
            // OCCT L926: if (aBPC.IsToUpdate())
            if let Some(edge) = self.ds.edges.get_mut(ei) {
                if edge.face_reps.iter().any(|r| r.face_idx == fi) {
                    continue;
                }
                // OCCT L928: double aTolE = BRep_Tool::Tolerance(aBPC.GetEdge())
                #[allow(unused_variables)]
                let a_tol_e = edge.geom_tol;
                // OCCT L929: aBB.UpdateEdge(aBPC.GetEdge(), aBPC.GetCurve2d(), aBPC.GetFace(), aTolE)
                edge.face_reps.push(DSCurveRepOnFace {
                    face_idx: fi,
                    pcurve: pcurve.clone(),
                    pcurve2: None,
                    pcurve_range: [0.0, span],
                    start_param: 0.0,
                    end_param: span,
                });
                // OCCT: UpdateEdge also updates edge tolerance (here aTolE == current, no change)
            }
        }
    }

    //  ?OCCT L248 Prepare: build pcurves on planar faces
    // RepeatInt->ForceEE->ForceEF->FF->UpdBlk->RefFI->MkSEdges->MkBlks->
    // ChkSI->RefFO->RmvME->MkPCurves->ProcDE

    /// BOPAlgo_PaveFiller::Init (PaveFiller.cxx L176-213).
    /// Populates the DS from operand BReps, matching what
    /// BOPDS_DS::Init + BOPDS_DS::SetArguments does in OCCT.
    /// Architecture diff: rcad borrows DS from caller (not new + SetArguments).
    /// Guard: allows pre-populated DS (rcad pattern), matching OCCT Init fresh.
    fn init(&mut self, a: &topods::BRep, b: &topods::BRep, fuzzy_tol: f64) {
        // OCCT L178-182: check arguments non-empty
        if !self.ds.faces.is_empty() {
            return;
        }
        // OCCT L196: Clear() — reset report and state
        self.clear();
        // OCCT L199-201: myDS = new BOPDS_DS; myDS->SetArguments(myArguments); myDS->Init(myFuzzyValue)
        let tol = fuzzy_tol.max(TOLERANCE_ABS);
        self.ds.fuzzy_tol = tol;
        // OCCT BOPDS_DS::Init: prepareVertices → prepareEdges → prepareFaces → prepareSolids
        self.prepare_vertices(a, b);
        self.prepare_edges(a, b);
        self.prepare_faces(a, b);
        self.prepare_solids(a, b);
        // Post-processing (shared by all prepare steps)
        self.ds.compute_uv_boundaries();
        self.ds.build_face_reps();
        self.ds.nb_source_shapes = self.ds.shape_info.len();
        self.ds.build_map_ve();
        // OCCT L204: myContext = new IntTools_Context — rcad: resize context
        if self.ds.face_count() != self.context.num_faces {
            self.context.resize(self.ds.face_count());
        }
        // OCCT L213: SetNonDestructive()
        self.set_non_destructive_auto();
    }

    /// OCCT BOPAlgo_PaveFiller::Clear (PaveFiller.cxx L185-192).
    /// Resets PaveFiller state while keeping the DS intact.
    pub fn clear(&mut self) {
        self.my_report.clear();
        self.fpbdone.clear();
        self.verts_to_avoid_extension.clear();
        self.my_increased_ss.clear();
        self.distances.clear();
    }

    /// OCCT BOPDS_DS::prepareVertices — traverse both operands and load all shapes into DS.
    /// Architecture diff: rcad loads all shape types in a single pass. Captures per-operand A counts
    /// between the two operand loads (matching original init() interleaving).
    fn prepare_vertices(&mut self, a: &topods::BRep, b: &topods::BRep) {
        use crate::bopds::ds::topods_builder::load_vertices_from_brep;
        load_vertices_from_brep(&mut self.ds, a, ShapeOrigin::ShapeA);
        // Capture A counts before loading B (original init() interleaving behavior)
        self.ds.a_vertex_count = self.ds.vertex_count();
        self.ds.a_edge_count = self.ds.edge_count();
        self.ds.a_face_count = self.ds.face_count();
        load_vertices_from_brep(&mut self.ds, b, ShapeOrigin::ShapeB);
    }

    /// OCCT BOPDS_DS::prepareEdges — structural no-op (loading done in prepareVertices).
    fn prepare_edges(&mut self, _a: &topods::BRep, _b: &topods::BRep) {
        use crate::bopds::ds::topods_builder::load_edges_from_brep;
        load_edges_from_brep(&mut self.ds, _a, ShapeOrigin::ShapeA);
        load_edges_from_brep(&mut self.ds, _b, ShapeOrigin::ShapeB);
    }

    /// OCCT BOPDS_DS::prepareFaces — structural no-op (loading done in prepareVertices).
    fn prepare_faces(&mut self, _a: &topods::BRep, _b: &topods::BRep) {
        use crate::bopds::ds::topods_builder::load_faces_from_brep;
        load_faces_from_brep(&mut self.ds, _a, ShapeOrigin::ShapeA);
        load_faces_from_brep(&mut self.ds, _b, ShapeOrigin::ShapeB);
    }

    /// OCCT BOPDS_DS::prepareSolids — structural no-op (loading done in prepareVertices).
    fn prepare_solids(&mut self, _a: &topods::BRep, _b: &topods::BRep) {
        use crate::bopds::ds::topods_builder::load_solids_from_brep;
        load_solids_from_brep(&mut self.ds, _a, ShapeOrigin::ShapeA);
        load_solids_from_brep(&mut self.ds, _b, ShapeOrigin::ShapeB);
    }

    /// Check whether to stop after the given stage name.
    /// Returns true if perform() should return early.
    fn check_stop(&self, stage: &str) -> bool {
        if self.stop_after.as_deref() == Some(stage) {
            return true;
        }
        if let Ok(s) = std::env::var("RCAD_STOP_AFTER") {
            return s == stage;
        }
        false
    }

    // OCCT BOPAlgo_PaveFiller.cxx L235-372: PerformInternal
    pub fn perform(&mut self, a: &topods::BRep, b: &topods::BRep) {
        // Init (PaveFiller.cxx L176-213).
        self.init(a, b, self.fuzzy_tolerance);
        self.dump_ctx.snapshot("after_Init", self.ds, None);
        if self.check_stop("after_Init") {
            return;
        }

        // OCCT L251-366: pipeline body (Prepare → ProcessDE)
        self.perform_body();
    }

    /// Run the pipeline body (Prepare through ProcessDE) on an already-initialized DS.
    /// This is the equivalent of PerformInternal after Init, and is used by
    /// PostTreatFF for the nested PaveFiller (OCCT PaveFiller_6.cxx L1392: aPF.Perform()).
    pub(crate) fn perform_body(&mut self) {
        // OCCT L251: Prepare =build pcurves on planar faces.
        self.prepare();
        if self.my_report.has_errors() {
            return;
        }
        self.dump_ctx.snapshot("after_Prepare", self.ds, None);
        if self.check_stop("after_Prepare") {
            return;
        }

        // OCCT L265: myIterator->Intersect(/*...*/) → BOPDS_Iterator with BVH
        // rcad: prepare iterator with all cross-operand pairs (stored in my_lists)
        use crate::bopds::ds::BOPDS_Iterator;
        use rcad_kernel::topods::ShapeType;
        self.my_iterator.prepare();

        self.perform_vv();
        if self.my_report.has_errors() {
            return;
        }
        self.dump_ctx.snapshot("after_PerformVV", self.ds, None);
        if self.check_stop("after_PerformVV") {
            return;
        }

        // OCCT: BOPDS_Iterator::Initialize(VERTEX, EDGE)
        self.perform_ve();
        if self.my_report.has_errors() {
            return;
        }
        self.dump_ctx.snapshot("after_PerformVE", self.ds, None);
        if self.check_stop("after_PerformVE") {
            return;
        }
        // OCCT: UpdatePaveBlocksWithSDVertices (after PerformVE, after dump)
        self.ds.update_pave_blocks_with_sd_vertices();

        // OCCT: BOPDS_Iterator::Initialize(EDGE, EDGE)
        self.perform_ee();
        if self.my_report.has_errors() {
            return;
        }
        self.dump_ctx.snapshot("after_PerformEE", self.ds, None);
        if self.check_stop("after_PerformEE") {
            return;
        }
        // OCCT: UpdatePaveBlocksWithSDVertices (after PerformEE, after dump)
        self.ds.update_pave_blocks_with_sd_vertices();

        // OCCT: BOPDS_Iterator::Initialize(VERTEX, FACE)
        self.perform_vf();
        if self.my_report.has_errors() {
            return;
        }
        self.dump_ctx.snapshot("after_PerformVF", self.ds, None);
        if self.check_stop("after_PerformVF") {
            return;
        }
        // OCCT: UpdatePaveBlocksWithSDVertices (after PerformVF, after dump)
        self.ds.update_pave_blocks_with_sd_vertices();

        // OCCT: BOPDS_Iterator::Initialize(EDGE, FACE) with BVH.
        let ef_pairs = {
            let mut ef_iterator = BOPDS_Iterator::new(self.ds);
            ef_iterator.prepare();
            ef_iterator.pairs(ShapeType::Edge, ShapeType::Face).to_vec()
        };
        self.perform_ef(&ef_pairs);
        if self.my_report.has_errors() {
            return;
        }
        // OCCT L295: UpdatePaveBlocksWithSDVertices
        self.ds.update_pave_blocks_with_sd_vertices();
        // OCCT L296: UpdateInterfsWithSDVertices
        self.update_interfs_with_sd_vertices();
        self.dump_ctx.snapshot("after_PerformEF", self.ds, None);
        if self.check_stop("after_PerformEF") {
            return;
        }

        // OCCT L300: RepeatIntersection
        self.repeat_intersection();
        if self.my_report.has_errors() {
            return;
        }
        self.dump_ctx
            .snapshot("after_RepeatIntersection", self.ds, None);
        if self.check_stop("after_RepeatIntersection") {
            return;
        }

        // OCCT L308: ForceInterfEE
        self.force_interf_ee();
        if self.my_report.has_errors() {
            return;
        }
        self.dump_ctx.snapshot("after_ForceInterfEE", self.ds, None);
        if self.check_stop("after_ForceInterfEE") {
            return;
        }

        // OCCT L316: ForceInterfEF
        self.force_interf_ef();
        if self.my_report.has_errors() {
            return;
        }
        self.dump_ctx.snapshot("after_ForceInterfEF", self.ds, None);
        if self.check_stop("after_ForceInterfEF") {
            return;
        }

        // OCCT L324: PerformFF
        self.perform_ff();
        if self.my_report.has_errors() {
            return;
        }
        self.dump_ctx.snapshot("after_PerformFF", self.ds, None);
        if self.check_stop("after_PerformFF") {
            return;
        }

        // OCCT L331: UpdateBlocksWithSharedVertices
        self.update_blocks_with_shared_vertices();

        // OCCT L333: RefineFaceInfoIn before MakeSplitEdges
        for fi in 0..self.ds.face_count() {
            self.ds.refine_face_info_in(fi);
        }

        // OCCT L335: MakeSplitEdges
        self.make_split_edges();
        if self.my_report.has_errors() {
            return;
        }
        self.dump_ctx
            .snapshot("after_MakeSplitEdges", self.ds, None);
        if self.check_stop("after_MakeSplitEdges") {
            return;
        }

        // OCCT L342: UpdatePaveBlocksWithSDVertices
        self.ds.update_pave_blocks_with_sd_vertices();

        // OCCT L344: MakeBlocks
        self.make_blocks();
        if self.my_report.has_errors() {
            return;
        }
        self.dump_ctx.snapshot("after_MakeBlocks", self.ds, None);
        if self.check_stop("after_MakeBlocks") {
            return;
        }

        // OCCT L351: CheckSelfInterference
        let _si_warnings = self.check_self_interference();

        // OCCT L353: UpdateInterfsWithSDVertices
        self.update_interfs_with_sd_vertices();

        // OCCT L354: ReleasePaveBlocks
        self.ds.release_pave_blocks();

        // OCCT L355: RefineFaceInfoOn =after ReleasePaveBlocks, remove
        // zero-length On pave blocks (BOPDS_DS::RefineFaceInfoOn).
        for fi in 0..self.ds.face_count() {
            self.ds.refine_face_info_on(fi);
        }

        // OCCT L357: RemoveMicroEdges =after MakeBlocks, before MakePCurves
        self.remove_micro_edges();

        // OCCT L359: MakePCurves =after RemoveMicroEdges
        self.make_pcurves();
        if self.my_report.has_errors() {
            return;
        }
        self.dump_ctx.snapshot("after_MakePCurves", self.ds, None);
        if self.check_stop("after_MakePCurves") {
            return;
        }

        // OCCT L366: ProcessDE =after MakePCurves
        self.process_de();
        if self.my_report.has_errors() {
            return;
        }
        self.dump_ctx.snapshot("after_ProcessDE", self.ds, None);
        if self.check_stop("after_ProcessDE") {
            return;
        }
    }

    /// BOPAlgo_PaveFiller::AddIntersectionFailedWarning (PaveFiller_2.cxx).
    /// Adds a warning that intersection between two shapes failed.
    pub(crate) fn add_intersection_failed_warning(&self, s1_idx: usize, s2_idx: usize) {
        // OCCT: creates BOPAlgo_AlertIntersectionFailed alert with shape pair info.
        // rcad: non-fatal warning logged to my_report. Shape indices logged for debugging.
        if std::env::var("RCAD_DEBUG_PF").is_ok() {
            eprintln!(
                "[PF] Intersection failed between shapes {} and {}",
                s1_idx, s2_idx
            );
        }
    }

    /// Helper: true when edge n_e is an intersection-created shape (not a source edge).
    /// OCCT equivalent: myDS->IsNewShape(nE) for an edge shape index.
    fn is_new_edge(&self, n_e: usize) -> bool {
        let si = if n_e < self.ds.edge_shape_idx.len() {
            self.ds.edge_shape_idx[n_e]
        } else {
            self.ds.vertex_count() + n_e
        };
        if si < self.ds.shape_info.len() {
            self.ds.shape_info[si].is_new
        } else {
            // no ShapeInfo entry -> created by push_edge during intersection
            true
        }
    }

    // OCCT BOPAlgo_PaveFiller_10.cxx L63-101
    pub(crate) fn update_edge_tolerance(&mut self, n_e: usize, a_tol_new: f64) {
        if n_e >= self.ds.edge_count() {
            return;
        }
        // OCCT L68-85: avoid modifying input shapes in safe (non-destructive) mode
        if self.non_destructive {
            // OCCT L71-74: if edge is not a new shape, return
            if !self.is_new_edge(n_e) {
                return;
            }
            // OCCT L76-84: if any vertex is old and has no SD, return
            let sv = self.ds.edge_start_vertex_ds(n_e);
            let ev = self.ds.edge_end_vertex_ds(n_e);
            for &n_v in &[sv, ev] {
                if !self.ds.is_new_vertex(n_v) && self.ds.has_shape_sd(n_v).is_none() {
                    return;
                }
            }
        }
        // OCCT L87-89: update edge tolerance (rcad: no TopoDS bounding box)
        let a_tol_e = self.ds.edge_tolerance(n_e);
        if a_tol_new > a_tol_e {
            self.ds.edge_data_mut(n_e).tolerance = a_tol_new;
        }
        // OCCT L94-100: update vertex tolerances
        let sv = self.ds.edge_start_vertex_ds(n_e);
        let ev = self.ds.edge_end_vertex_ds(n_e);
        self.update_vertex(sv, a_tol_new);
        self.update_vertex(ev, a_tol_new);
    }

    // OCCT BOPAlgo_PaveFiller_10.cxx L105-162
    pub(crate) fn update_vertex(&mut self, n_v: usize, a_tol_new: f64) -> usize {
        if n_v >= self.ds.vertex_count() {
            return n_v;
        }
        // OCCT L111: nVNew = nV
        let mut n_v_new = n_v;
        // OCCT L112: check is_new, has_sd, or non-destructive not in force
        let sd_opt = self.ds.has_shape_sd(n_v);
        // OCCT L112: IsNewShape || HasShapeSD || !NonDestructive
        if self.ds.is_new_vertex(n_v) || sd_opt.is_some() || !self.non_destructive {
            // OCCT L115: if HasShapeSD(nV, nVNew), nVNew becomes the SD partner
            if let Some(n_sd) = sd_opt {
                n_v_new = n_sd;
            }
            // OCCT L116: get current tolerance of (possibly SD-partner) vertex
            let a_tol_v = self.ds.vertex_tolerance(n_v_new);
            // OCCT L117-125: increase tolerance if needed
            if a_tol_v < a_tol_new {
                self.ds.vertex_data_mut(n_v_new).tolerance = a_tol_new;
                // OCCT L120-123: update bounding box in shape_info
                //   BOPDS_ShapeInfo& aSIV = myDS->ChangeShapeInfo(nVNew);
                //   Bnd_Box& aBoxV = aSIV.ChangeBox();
                //   BRepBndLib::Add(aVSD, aBoxV);
                //   aBoxV.SetGap(aBoxV.GetGap() + Precision::Confusion());
                let si_idx = self.ds.vertex_shape_idx.get(n_v_new).copied();
                if let Some(si) = si_idx {
                    if si < self.ds.shape_info.len() {
                        let pt = self.ds.vertex_point(n_v_new);
                        let new_gap = a_tol_new + self.fuzzy_tolerance * 0.5;
                        self.ds.shape_info[si].box_min = Some(pt - DVec3::splat(a_tol_new));
                        self.ds.shape_info[si].box_max = Some(pt + DVec3::splat(a_tol_new));
                        self.ds.shape_info[si].box_gap = new_gap + crate::tolerance::CONFUSION;
                    }
                }
                // OCCT L124: myIncreasedSS.Add(nV) — adds the original vertex
                self.my_increased_ss.insert(n_v);
            }
            return n_v_new;
        }
        // OCCT L129-159: nV is old vertex — create new SD vertex (non-destructive mode)
        // rcad: non_destructive is always false, so this path is never reached.
        // The branch above handles all cases (is_new || has_sd || !non_destructive).
        return n_v;
    }

    /// BOPAlgo_PaveFiller::UpdateCommonBlocksWithSDVertices (PaveFiller_10.cxx L173-221).
    pub(crate) fn update_common_blocks_with_sd_vertices(&mut self) {
        if !self.non_destructive {
            self.ds.update_pave_blocks_with_sd_vertices();
            return;
        }
        // Collect CB indices first to avoid borrow conflicts
        let mut cb_indices: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for ei in 0..self.ds.edge_count() {
            for spb in &self.ds.edges[ei].pave_blocks {
                let pb = spb.0.read().unwrap();
                if let Some(cb_idx) = pb.common_block_idx {
                    cb_indices.insert(cb_idx);
                }
            }
        }
        let a_tol_v = crate::tolerance::TOLERANCE_ABS;
        for &cb_idx in &cb_indices {
            // Find the first PB associated with this CB to get its vertices
            let vertices: Vec<usize> = {
                let mut verts = Vec::new();
                for ei in 0..self.ds.edge_count() {
                    for spb in &self.ds.edges[ei].pave_blocks {
                        let pb = spb.0.read().unwrap();
                        if pb.common_block_idx == Some(cb_idx) {
                            let (nv1, nv2) = pb.indices();
                            verts.push(nv1);
                            verts.push(nv2);
                            break;
                        }
                    }
                    if !verts.is_empty() {
                        break;
                    }
                }
                verts
            };
            for &nv in &vertices {
                self.update_vertex(nv, a_tol_v);
            }
        }
        self.ds.update_pave_blocks_with_sd_vertices();
    }

    /// BOPAlgo_PaveFiller::UpdateVerticesOfCB (PaveFiller_3.cxx L959-993).
    /// Updates vertices of CommonBlocks with the CommonBlock's tolerance.
    pub(crate) fn update_vertices_of_cb(&mut self) {
        // Collect (vertex_idx, tolerance) pairs first to avoid borrow conflicts
        let mut updates: Vec<(usize, f64)> = Vec::new();
        let mut a_mpb_fence: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for ei in 0..self.ds.edge_count() {
            for spb in &self.ds.edges[ei].pave_blocks {
                let pb = spb.0.read().unwrap();
                if let Some(cb_idx) = pb.common_block_idx {
                    if cb_idx < self.ds.common_blocks.len() {
                        let a_tol_cb = self.ds.common_blocks[cb_idx].tolerance();
                        if a_tol_cb > 0.0 {
                            let fence_key = pb.new_edge.unwrap_or(pb.original_edge);
                            if a_mpb_fence.insert(fence_key) {
                                updates.push((pb.pave1.vertex_idx, a_tol_cb));
                                updates.push((pb.pave2.vertex_idx, a_tol_cb));
                            }
                        }
                    }
                }
            }
        }
        for (vi, tol) in updates {
            self.update_vertex(vi, tol);
        }
    }
    /// Removes all PaveBlocks belonging to the given edge indices from:
    /// 1. the PaveBlocks Pool (edge_pave_blocks)
    /// 2. section curve PB lists
    /// 3. FaceInfo ON/IN/SC sets
    pub(crate) fn remove_pave_blocks(&mut self, the_edges: &std::collections::HashSet<usize>) {
        if the_edges.is_empty() {
            return;
        }
        // 1. From edge pave blocks
        for &ei in the_edges {
            if ei < self.ds.edge_count() {
                self.ds.edge_pave_blocks_mut(ei).clear();
            }
        }
        // 2. From section curves
        for ic in &mut self.ds.intersection_curves {
            ic.pave_blocks.retain(|spb| {
                let e = spb
                    .0
                    .read()
                    .unwrap()
                    .new_edge
                    .unwrap_or(spb.0.read().unwrap().original_edge);
                !the_edges.contains(&e)
            });
        }
        // 3. From FaceInfo
        for fi in 0..self.ds.face_count() {
            let fi_copy = fi; // avoid borrow conflict
            // Collect PB indices to remove
            let to_remove: std::collections::HashSet<usize> = {
                let face_info = &self.ds.faces[fi_copy].face_info;
                face_info
                    .pave_blocks_on
                    .iter()
                    .chain(face_info.pave_blocks_in.iter())
                    .chain(face_info.pave_blocks_sc.iter())
                    .filter(|&&pb_idx| {
                        if pb_idx >= self.ds.pave_blocks.len() {
                            return true;
                        }
                        let pb = &self.ds.pave_blocks[pb_idx];
                        let e =
                            pb.0.read()
                                .unwrap()
                                .new_edge
                                .unwrap_or(pb.0.read().unwrap().original_edge);
                        the_edges.contains(&e)
                    })
                    .copied()
                    .collect()
            };
            let fi_mut = fi;
            for &pb_idx in &to_remove {
                self.ds.face_info_mut(fi_mut).pave_blocks_on.remove(&pb_idx);
                self.ds.face_info_mut(fi_mut).pave_blocks_in.remove(&pb_idx);
                self.ds.face_info_mut(fi_mut).pave_blocks_sc.remove(&pb_idx);
            }
        }
    }

    /// BOPAlgo_PaveFiller::RemoveMicroSectionEdges (PaveFiller_6.cxx L4341-4417).
    /// Identifies micro section edges (too short / no valid shrunk data) and
    /// removes them from the section edge map, adding them to theMicroPB set
    /// for vertex unification in PostTreatFF.
    pub(crate) fn remove_micro_section_edges(
        &mut self,
        a_mscpb: &mut std::collections::HashMap<usize, (usize, usize)>,
        a_micro_pb: &mut Vec<crate::bopds::pave::PaveBlock>,
    ) {
        if a_mscpb.is_empty() {
            return;
        }
        let mut a_sepb_map: std::collections::HashMap<usize, (usize, usize)> =
            std::collections::HashMap::new();
        let keys: Vec<usize> = a_mscpb.keys().copied().collect();
        for &edge_or_vertex in &keys {
            let Some(&cpb) = a_mscpb.get(&edge_or_vertex) else {
                continue;
            };
            if edge_or_vertex < self.ds.edge_count() {
                // It's an edge — check if it's a micro edge
                let ei = edge_or_vertex;
                let is_micro = {
                    let (sv, ev) = {
                        let e = &self.ds.edges[ei];
                        (e.start_vertex, e.end_vertex)
                    };
                    if sv < self.ds.vertex_count() && ev < self.ds.vertex_count() {
                        let v1 = self.ds.vertex_point(sv);
                        let v2 = self.ds.vertex_point(ev);
                        v1.distance(v2) < TOLERANCE_ABS * 10.0
                    } else {
                        false
                    }
                };
                if !is_micro {
                    a_sepb_map.insert(edge_or_vertex, cpb);
                } else {
                    // Micro edge: add PB to theMicroPB for PostTreatFF
                    let pb = &self.ds.pave_blocks[edge_or_vertex];
                    a_micro_pb.push(pb.0.read().unwrap().clone());
                    // Remove from section curves
                    if cpb.0 < self.ds.intersection_curves.len() {
                        let ci = cpb.0;
                        self.ds.intersection_curves[ci].pave_blocks.retain(|spb| {
                            let e = spb
                                .0
                                .read()
                                .unwrap()
                                .new_edge
                                .unwrap_or(spb.0.read().unwrap().original_edge);
                            e != ei
                        });
                    }
                }
            } else {
                // Not an edge — pass through
                a_sepb_map.insert(edge_or_vertex, cpb);
            }
        }
        *a_mscpb = a_sepb_map;
    }

    // OCCT BOPAlgo_PaveFiller::UpdatePaveBlocks (PaveFiller_6.cxx L3712-3844)
    pub(crate) fn update_pave_blocks(&mut self, a_dm_new_sd: &HashMap<usize, usize>) {
        if a_dm_new_sd.is_empty() {
            return;
        }
        let mut a_micro_edges: HashSet<usize> = HashSet::new();
        let mut a_mpb: HashSet<usize> = HashSet::new();

        // Collect all PBs: section curves + pool
        let mut an_all_pbs: Vec<SharedPB> = Vec::new();

        // OCCT L3728-3746: section curve PBs
        let a_nb_ff = self.ds.interf_ff.len();
        for i in 0..a_nb_ff {
            let ff = &self.ds.interf_ff[i];
            for &ci in &ff.curves {
                if ci >= self.ds.intersection_curves.len() {
                    continue;
                }
                let ic = &self.ds.intersection_curves[ci];
                for spb in &ic.pave_blocks {
                    an_all_pbs.push(spb.clone());
                }
            }
        }

        // OCCT L3748-3760: pool PBs (all edge PBs)
        for ei in 0..self.ds.edge_count() {
            for spb in &self.ds.edges[ei].pave_blocks {
                an_all_pbs.push(spb.clone());
            }
        }

        // OCCT L3762-3837: process all PBs
        for spb in &an_all_pbs {
            let a_pb = spb.clone();
            // OCCT L3767: handle<CommonBlock>& aCB = myDS->CommonBlock(aPB);
            let orig_cb_idx = a_pb.0.read().unwrap().common_block_idx;
            let b_cb = orig_cb_idx.is_some();

            // OCCT L3769-3772: if (bCB) { aPB = aCB->PaveBlock1(); }
            let a_pb_primary = if let Some(cb_idx) = orig_cb_idx {
                if cb_idx < self.ds.common_blocks.len() {
                    let cb = &self.ds.common_blocks[cb_idx];
                    let pb1_local = cb.pave_block1();
                    let edge_idx = cb.edge();
                    if let (Some(pli), Some(ei)) = (pb1_local, edge_idx) {
                        if ei < self.ds.edge_count() && pli < self.ds.edges[ei].pave_blocks.len() {
                            self.ds.edges[ei].pave_blocks[pli].clone()
                        } else {
                            a_pb.clone()
                        }
                    } else {
                        a_pb.clone()
                    }
                } else {
                    a_pb.clone()
                }
            } else {
                a_pb.clone()
            };

            // OCCT L3774: if (!aMPB.Add(aPB)) { continue; }
            let ptr = Arc::as_ptr(&a_pb_primary.0) as usize;
            if !a_mpb.insert(ptr) {
                continue;
            }

            // OCCT L3776-3778: aPB->Indices(nV[0], nV[1]); aPB->Range(aT[0], aT[1]);
            let (mut n_v, a_t) = {
                let pb_ref = a_pb_primary.0.read().unwrap();
                (
                    [pb_ref.pave1.vertex_idx, pb_ref.pave2.vertex_idx],
                    [pb_ref.pave1.param, pb_ref.pave2.param],
                )
            };

            // OCCT L3780: bool wasRegularEdge = (nV[0] != nV[1]);
            let was_regular_edge = n_v[0] != n_v[1];
            let mut b_rebuild = false;

            // OCCT L3782-3801: replace vertices via aDMNewSD
            for j in 0..2 {
                if let Some(&new_v) = a_dm_new_sd.get(&n_v[j]) {
                    n_v[j] = new_v;
                    b_rebuild = true;
                    let a_pave = Pave {
                        vertex_idx: new_v,
                        param: a_t[j],
                    };
                    if j == 0 {
                        a_pb_primary.0.write().unwrap().pave1 = a_pave;
                    } else {
                        a_pb_primary.0.write().unwrap().pave2 = a_pave;
                    }
                }
            }

            // OCCT L3804: if (bRebuild) { ... }
            if !b_rebuild {
                continue;
            }

            // OCCT L3806-3812: int nE = aPB->Edge(); if (nE < 0) nE = aPB->OriginalEdge();
            let n_e = {
                let pb_ref = a_pb_primary.0.read().unwrap();
                if let Some(ne) = pb_ref.new_edge {
                    ne
                } else {
                    pb_ref.original_edge
                }
            };
            if n_e >= self.ds.edge_count() {
                continue;
            }

            // OCCT L3813: bool isDegEdge = myDS->ShapeInfo(nE).HasFlag();
            let is_degen_edge = self.ds.edge_has_flag(n_e);

            // OCCT L3814-3824: micro edge check
            if was_regular_edge && !is_degen_edge && n_v[0] == n_v[1] {
                // OCCT L3818: FillShrunkData(aPB);
                // OCCT L3819: if (!aPB->HasShrunkData())
                let has_shrunk = a_pb_primary.0.read().unwrap().has_shrunk_data();
                if !has_shrunk {
                    a_micro_edges.insert(n_e);
                    continue;
                }
            }

            // OCCT L3826: nSp = SplitEdge(nE, nV[0], aT[0], nV[1], aT[1]);
            let n_sp = self.split_edge(n_e, n_v[0], a_t[0], n_v[1], a_t[1]);

            // OCCT L3827-3834: if (bCB) aCB->SetEdge(nSp); else aPB->SetEdge(nSp);
            if b_cb {
                if let Some(cb_idx) = orig_cb_idx {
                    if cb_idx < self.ds.common_blocks.len() {
                        self.ds.common_blocks[cb_idx].set_edge(n_sp);
                    }
                }
            } else {
                a_pb_primary.0.write().unwrap().new_edge = Some(n_sp);
            }
        } // for spb in an_all_pbs

        // OCCT L3840-3843: if (aMicroEdges.Extent()) RemovePaveBlocks(aMicroEdges);
        if !a_micro_edges.is_empty() {
            self.remove_pave_blocks(&a_micro_edges);
        }
    }

    // OCCT BOPAlgo_PaveFiller::UpdateFaceInfo (PaveFiller_6.cxx L1705-1978)
    #[allow(non_snake_case)]
    pub(crate) fn update_face_info_post(
        &mut self,
        theDME: &HashMap<usize, Vec<usize>>,
        theDMV: &HashMap<usize, usize>,
        thePBFacesMap: &HashMap<usize, Vec<usize>>,
    ) {
        // OCCT L1715: anEdgeLPB — map from edge index to list of PBs
        let mut an_edge_lpb: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut a_mf: HashSet<usize> = HashSet::new();

        // OCCT L1720: process all FF interferences
        let a_nb_ff = self.ds.interf_ff.len();
        for i in 0..a_nb_ff {
            let (n_f1, n_f2) = {
                let ff = &self.ds.interf_ff[i];
                (ff.f1, ff.f2)
            };
            let a_fi1 = n_f1;
            let a_fi2 = n_f2;

            // OCCT L1729-1778: 1.1 Section edges
            let curve_idxs: Vec<usize> = self.ds.interf_ff[i].curves.clone();
            for &ci in &curve_idxs {
                if ci >= self.ds.intersection_curves.len() {
                    continue;
                }
                // Collect curve's PBs (need to handle removal during iteration)
                let curve_pbs: Vec<usize> = {
                    let ic = &self.ds.intersection_curves[ci];
                    (0..ic.pave_blocks.len()).collect()
                };
                let mut to_keep: Vec<usize> = Vec::new();

                for &pb_local in &curve_pbs {
                    let pb_idx = pb_local; // global PB index

                    // OCCT L1744: if (theDME.IsBound(aPB))
                    if let Some(replacements) = theDME.get(&pb_idx) {
                        // OCCT L1747: UpdateExistingPaveBlocks(aPB, aLPB, thePBFacesMap);
                        update_existing_pave_blocks(
                            &mut self.ds,
                            &mut self.context,
                            pb_idx,
                            replacements,
                            thePBFacesMap,
                            self.fuzzy_tolerance,
                        );

                        // OCCT L1749-1758: add replacement PBs to anEdgeLPB
                        for &rp in replacements {
                            let n_e = {
                                let pb_r = self.ds.pave_blocks[rp].0.read().unwrap();
                                pb_r.new_edge.unwrap_or(pb_r.original_edge)
                            };
                            an_edge_lpb.entry(n_e).or_default().push(rp);
                        }

                        // OCCT L1761: aLPBC.Remove(aItPB) — don't add to keep list
                        continue;
                    }

                    // OCCT L1765-1766: normal section PB → add to both faces' pave_blocks_sc
                    self.ds.face_info_mut(a_fi1).pave_blocks_sc.insert(pb_idx);
                    self.ds.face_info_mut(a_fi2).pave_blocks_sc.insert(pb_idx);

                    // OCCT L1768-1774: add to anEdgeLPB
                    let n_e = {
                        let pb_r = self.ds.pave_blocks[pb_idx].0.read().unwrap();
                        pb_r.new_edge.unwrap_or(pb_r.original_edge)
                    };
                    an_edge_lpb.entry(n_e).or_default().push(pb_idx);
                    to_keep.push(pb_idx);
                }
                // OCCT L1761: actual removal from curve PB list
                if to_keep.len() < curve_pbs.len() {
                    let ic = &mut self.ds.intersection_curves[ci];
                    // Build a set of Arc pointers for PBs to keep
                    let keep_ptrs: HashSet<*const RwLock<PaveBlock>> = to_keep
                        .iter()
                        .filter_map(|&i| ic.pave_blocks.get(i))
                        .map(|spb| Arc::as_ptr(&spb.0))
                        .collect();
                    ic.pave_blocks
                        .retain(|spb| keep_ptrs.contains(&Arc::as_ptr(&spb.0)));
                }
            } // for each curve

            // OCCT L1781-1793: 1.2 Section vertices (point contacts)
            let (f1, f2, points) = {
                let ff = &self.ds.interf_ff[i];
                (ff.f1, ff.f2, ff.points.clone())
            };
            for ffp in &points {
                if ffp.vertex_index < self.ds.vertex_count() {
                    let n_v = ffp.vertex_index;
                    if !self.ds.faces[f1].face_info.vertices_in.contains(&n_v) {
                        self.ds.faces[f1].face_info.vertices_in.insert(n_v);
                    }
                    if !self.ds.faces[f2].face_info.vertices_in.contains(&n_v) {
                        self.ds.faces[f2].face_info.vertices_in.insert(n_v);
                    }
                }
            }

            // OCCT L1795-1796: track faces
            a_mf.insert(a_fi1);
            a_mf.insert(a_fi2);
        } // for each FF

        // OCCT L1799-1889: unify PBs on the same edge (anEdgeLPB) via CommonBlocks
        let mut b_new_cb = false;
        {
            for (_n_e, pb_list) in an_edge_lpb.iter() {
                if pb_list.len() <= 1 {
                    continue;
                }

                b_new_cb = true;

                let mut a_cb_idx: Option<usize> = None;
                let mut a_m_faces: HashSet<usize> = HashSet::new();
                let mut a_mpave_blocks: IndexSet<usize> = IndexSet::new();

                for &pb_idx in pb_list {
                    a_mpave_blocks.insert(pb_idx);

                    // OCCT L1828-1850: if PB has a CommonBlock, collect its PBs and faces
                    let spb = &self.ds.pave_blocks[pb_idx];
                    if let Some(cb_idx) = spb.0.read().unwrap().common_block_idx {
                        if cb_idx < self.ds.common_blocks.len() {
                            let cb = &self.ds.common_blocks[cb_idx];
                            // Collect all PBs from this CB
                            for &(local_pb, _fi) in cb.pave_blocks() {
                                a_mpave_blocks.insert(local_pb);
                            }
                            // Collect faces
                            for &f in cb.faces() {
                                a_m_faces.insert(f);
                            }
                            if a_cb_idx.is_none() {
                                a_cb_idx = Some(cb_idx);
                            }
                        }
                    }
                }

                if let Some(cb_idx) = a_cb_idx {
                    // OCCT L1868-1888: extend existing CommonBlock
                    let all_pbs: Vec<(usize, usize)> = a_mpave_blocks
                        .iter()
                        .map(|&pb| (pb, 0)) // face_idx placeholder
                        .collect();
                    self.ds.common_blocks[cb_idx].set_pave_blocks(all_pbs);
                    let faces: Vec<usize> = a_m_faces.iter().copied().collect();
                    self.ds.common_blocks[cb_idx].set_faces(faces);
                    // Set CommonBlock on each PB
                    // OCCT L1875: myDS->SetCommonBlock(aPB, aCB);
                    // In rcad, CB stores local indices; set common_block_idx on each
                } else {
                    // OCCT L1857-1864: create new CommonBlock
                    let mut a_cb = CommonBlock::new();
                    let all_pbs: Vec<(usize, usize)> = pb_list.iter().map(|&pb| (pb, 0)).collect();
                    a_cb.set_pave_blocks(all_pbs);
                    let cb_idx = self.ds.common_blocks.len();
                    // Set common_block_idx on each PB
                    for &pb_idx in pb_list {
                        if pb_idx < self.ds.pave_blocks.len() {
                            self.ds.pave_blocks[pb_idx]
                                .0
                                .write()
                                .unwrap()
                                .common_block_idx = Some(cb_idx);
                        }
                    }
                    self.ds.common_blocks.push(a_cb);
                }
            }
        }

        // OCCT L1892-1897: early return if no changes needed
        let b_verts = !theDMV.is_empty();
        let b_edges = !theDME.is_empty() || b_new_cb;
        if !b_verts && !b_edges {
            return;
        }

        // OCCT L1906-1977: update face info for each face in aMF
        for &n_f1 in &a_mf {
            let a_fi = n_f1;

            // OCCT L1914-1934: 2.1 Update vertex ON/IN sets
            if b_verts {
                for (&n_v1, &n_v2) in theDMV.iter() {
                    if self.ds.face_info(a_fi).vertices_on.contains(&n_v1) {
                        self.ds.face_info_mut(a_fi).vertices_on.remove(&n_v1);
                        self.ds.face_info_mut(a_fi).vertices_on.insert(n_v2);
                    }
                    if self.ds.face_info(a_fi).vertices_in.contains(&n_v1) {
                        self.ds.face_info_mut(a_fi).vertices_in.remove(&n_v1);
                        self.ds.face_info_mut(a_fi).vertices_in.insert(n_v2);
                    }
                }
            }

            // OCCT L1938-1975: 2.2 Update PaveBlock ON/IN/SC sets
            if b_edges {
                let mut a_mpb_fence: HashSet<usize> = HashSet::new();

                // Copy the three PB sets before clearing (OCCT: aMPBCopy = *pMPB[i])
                let on_copy: Vec<usize> = self
                    .ds
                    .face_info(a_fi)
                    .pave_blocks_on
                    .iter()
                    .copied()
                    .collect();
                let in_copy: Vec<usize> = self
                    .ds
                    .face_info(a_fi)
                    .pave_blocks_in
                    .iter()
                    .copied()
                    .collect();
                let sc_copy: Vec<usize> = self
                    .ds
                    .face_info(a_fi)
                    .pave_blocks_sc
                    .iter()
                    .copied()
                    .collect();

                // Clear all three sets (OCCT: pMPB[i]->Clear())
                self.ds.face_info_mut(a_fi).pave_blocks_on.clear();
                self.ds.face_info_mut(a_fi).pave_blocks_in.clear();
                self.ds.face_info_mut(a_fi).pave_blocks_sc.clear();

                // Rebuild all three sets using aMPBFence for dedup (OCCT L1940-1974)
                let mut new_on: Vec<usize> = Vec::new();
                let mut new_in: Vec<usize> = Vec::new();
                let mut new_sc: Vec<usize> = Vec::new();

                for &pb_idx in &on_copy {
                    if let Some(replacements) = theDME.get(&pb_idx) {
                        for &rp in replacements {
                            // OCCT L1959: RealPaveBlock(aPB1) — in rcad the index IS the real block
                            if a_mpb_fence.insert(rp) {
                                new_on.push(rp);
                            }
                        }
                    } else {
                        if a_mpb_fence.insert(pb_idx) {
                            new_on.push(pb_idx);
                        }
                    }
                }
                for &pb_idx in &in_copy {
                    if let Some(replacements) = theDME.get(&pb_idx) {
                        for &rp in replacements {
                            if a_mpb_fence.insert(rp) {
                                new_in.push(rp);
                            }
                        }
                    } else {
                        if a_mpb_fence.insert(pb_idx) {
                            new_in.push(pb_idx);
                        }
                    }
                }
                for &pb_idx in &sc_copy {
                    if let Some(replacements) = theDME.get(&pb_idx) {
                        for &rp in replacements {
                            if a_mpb_fence.insert(rp) {
                                new_sc.push(rp);
                            }
                        }
                    } else {
                        if a_mpb_fence.insert(pb_idx) {
                            new_sc.push(pb_idx);
                        }
                    }
                }

                // Insert rebuilt sets
                for &pb in &new_on {
                    self.ds.face_info_mut(a_fi).pave_blocks_on.insert(pb);
                }
                for &pb in &new_in {
                    self.ds.face_info_mut(a_fi).pave_blocks_in.insert(pb);
                }
                for &pb in &new_sc {
                    self.ds.face_info_mut(a_fi).pave_blocks_sc.insert(pb);
                }
            }
        }
    }

    // ===== BVH-based pair enumeration (OCCT BOPDS_Iterator) =====

    // = = = =Pass 1: Vertex-Vertex = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    // = = = =Pass 2: Vertex-Edge = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    // = = = =Pass 3: Edge-Edge = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    // = = = =Pass 4: Vertex-Face = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    // = = = =Pass 5: Edge-Face = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    // = = = =Pass 6: Face-Face = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    // OCCT BOPAlgo_PaveFiller.cxx L331: UpdateBlocksWithSharedVertices
    fn update_blocks_with_shared_vertices(&mut self) {
        if !self.non_destructive {
            return;
        }
        let has_ff = !self.ds.interf_ff.is_empty();
        if !has_ff {
            return;
        }

        // Collect face pairs with shared (old) vertices
        let ff_entries: Vec<(usize, usize, Vec<usize>)> = self
            .ds
            .interf_ff
            .iter()
            .filter_map(|ff| {
                if ff.curves.is_empty() {
                    return None;
                }
                let fi1 = ff.f1;
                let fi2 = ff.f2;
                let on1 = &self.ds.faces[fi1].face_info.vertices_on;
                let in1 = &self.ds.face_info(fi1).vertices_in;
                let on2 = &self.ds.faces[fi2].face_info.vertices_on;
                let in2 = &self.ds.face_info(fi2).vertices_in;

                let shared: Vec<usize> = on1
                    .iter()
                    .chain(in1.iter())
                    .filter(|&&vi| {
                        if self.ds.is_new_vertex(vi) {
                            return false;
                        }
                        on2.contains(&vi) || in2.contains(&vi)
                    })
                    .copied()
                    .collect();

                if shared.is_empty() {
                    return None;
                }
                Some((fi1, fi2, ff.curves.clone()))
            })
            .collect();
        for (f1, f2, curves) in &ff_entries {
            // rcad: not needed =FaceInfo data is already populated.

            for &ci in curves {
                if ci >= self.ds.intersection_curves.len() {
                    continue;
                }
                let ic_curve = self.ds.intersection_curves[ci].curve.clone();
                let ic_geom_tol = self.ds.intersection_curves[ci].geom_tol;
                let _a_tol_r3d = ic_geom_tol;
                let f1 = *f1;
                let f2 = *f2;

                // Collect shared vertices for this face pair
                let on1 = &self.ds.faces[f1].face_info.vertices_on;
                let in1 = &self.ds.face_info(f1).vertices_in;
                let on2 = &self.ds.faces[f2].face_info.vertices_on;
                let in2 = &self.ds.face_info(f2).vertices_in;

                let shared: Vec<usize> = on1
                    .iter()
                    .chain(in1.iter())
                    .filter(|&&vi| {
                        if self.ds.is_new_vertex(vi) {
                            return false;
                        }
                        if !on2.contains(&vi) && !in2.contains(&vi) {
                            return false;
                        }
                        // rcad: check shape_sd
                        if self.ds.shape_sd.is_sub_vertex(vi) {
                            return false;
                        }
                        true
                    })
                    .copied()
                    .collect();

                for &n_v in &shared {
                    if self.estimate_pave_on_curve(ci, n_v).is_none() {
                        continue;
                    }
                    let v_tol = self.ds.vertex_tolerance(n_v);
                    // UpdateVertex: increase tolerance if the projection distance is larger
                    let t_result =
                        self.project_vertex_on_curve(n_v, &self.ds.intersection_curves[ci]);
                    if let Some(t) = t_result {
                        let pt_on_curve = ic_curve.point_at(t);
                        let dist = self.ds.vertex_point(n_v).distance(pt_on_curve);
                        if dist > v_tol {
                            self.ds.vertex_data_mut(n_v).tolerance = dist;
                            self.my_increased_ss.insert(n_v);
                        }
                    }
                    // InitPaveBlocksForVertex: collect edge indices + params, then apply
                    let mut new_paves: Vec<(usize, f64)> = Vec::new();
                    for (ei, edge) in self.ds.edges.iter().enumerate() {
                        if edge.start_vertex == n_v {
                            let has = edge.paves.iter().any(|p| p.vertex_idx == n_v);
                            if !has {
                                new_paves.push((ei, edge.t_range[0]));
                            }
                        } else if edge.end_vertex == n_v {
                            let has = edge.paves.iter().any(|p| p.vertex_idx == n_v);
                            if !has {
                                new_paves.push((ei, edge.t_range[1]));
                            }
                        }
                    }
                    for (ei, param) in new_paves {
                        self.add_pave_to_edge(
                            ei,
                            Pave {
                                vertex_idx: n_v,
                                param,
                            },
                        );
                    }
                }
            }
        }
        self.ds.update_pave_blocks_with_sd_vertices();
    }

    // OCCT BOPAlgo_PaveFiller.cxx L296: UpdateInterfsWithSDVertices
    fn update_interfs_with_sd_vertices(&mut self) {
        // Build vertex =SD vertex lookup (OCCT HasShapeSD equivalent)
        let sd_for: std::collections::HashMap<usize, usize> = self
            .ds
            .shape_sd
            .sd_vertices_iter()
            .filter_map(|&(a, b)| {
                // Stored symmetrically; only process (a,b) where a < b
                // to avoid double-insert.  Both directions work since
                // all SD pairs have symmetric entries.
                if a < b { Some((a, b)) } else { None }
            })
            .collect();
        for inf in &mut self.ds.interf_ee {
            if let Some(&sd) = sd_for.get(&inf.new_vertex) {
                inf.new_vertex = sd;
            }
        }
        for inf in &mut self.ds.interf_ef {
            if let Some(&sd) = sd_for.get(&inf.new_vertex) {
                inf.new_vertex = sd;
            }
        }
        for inf in &mut self.ds.interf_vv {
            if let Some(&sd) = sd_for.get(&inf.v1) {
                inf.merged_vertex = sd;
            } else if let Some(&sd) = sd_for.get(&inf.v2) {
                inf.merged_vertex = sd;
            }
        }
    }
    // Dead code removed: make_section_edges was always no-op.
    // Section edge creation is now handled per-curve inside make_blocks (OCCT form alignment).

    // OCCT BOPAlgo_PaveFiller_6.cxx L4341-4417 - RemoveMicroEdges
    // Removes edges whose PaveBlocks have no valid shrunk range.
    // For each edge with PBs: if any PB has zero-length range, no ShrunkData,
    // or shrunk range < edge tolerance, the edge's PBs are cleared.
    fn remove_micro_edges(&mut self) {
        let a_nb_e = self.ds.edge_count();
        for i in 0..a_nb_e {
            // OCCT L4343: skip removed/degenerated edges
            if self.ds.edge_has_flag(i) || self.ds.is_edge_degenerated(i) {
                continue;
            }
            // OCCT L4345: aLPB = myDS->PaveBlocks(i)
            let a_lpb: Vec<crate::bopds::pave::SharedPB> = self.ds.edge_pave_blocks(i).to_vec();
            if a_lpb.is_empty() {
                continue;
            }
            let mut b_to_remove = false;
            for spb in &a_lpb {
                let pb = spb.0.read().unwrap();
                // OCCT L4351: aPB->Range(aT1, aT2)
                let (a_t1, a_t2) = pb.range();
                // OCCT L4352: if (Abs(aT2 - aT1) < Precision::PConfusion())
                if (a_t2 - a_t1).abs() < rcad_kernel::tolerance::CONFUSION {
                    b_to_remove = true;
                    break;
                }
                // OCCT L4357: if (!aPB->HasShrunkData())
                if !pb.has_shrunk_data() {
                    b_to_remove = true;
                    break;
                }
                // OCCT L4362: aPB->ShrunkData(aTS1, aTS2, aIsSplittable)
                let (a_ts1, a_ts2, _a_is_splittable) = pb.shrunk_data();
                // OCCT L4365: aTolE = Max(Tolerance(i, EDGE), Precision::Confusion())
                let a_tol_e = self.ds.edge_tolerance(i).max(crate::tolerance::CONFUSION);
                // OCCT L4367: if (Abs(aTS2 - aTS1) < aTolE)
                if (a_ts2 - a_ts1).abs() < a_tol_e {
                    b_to_remove = true;
                    break;
                }
            }
            // OCCT L4374: if bToRemove -> clear all PBs for this edge
            if b_to_remove {
                self.ds.edge_pave_blocks_mut(i).clear();
            }
        }
    }

    // Missing WIP methods
    pub(crate) fn faces_of(&self, origin: ShapeOrigin) -> Vec<usize> {
        self.ds
            .faces
            .iter()
            .enumerate()
            .filter(|(_, f)| f.origin == origin)
            .map(|(i, _)| i)
            .collect()
    }
    pub(crate) fn verts_of(&self, origin: ShapeOrigin) -> Vec<usize> {
        self.ds
            .vertices
            .iter()
            .enumerate()
            .filter(|(_, v)| v.origin == Some(origin))
            .map(|(i, _)| i)
            .collect()
    }
    pub(crate) fn edges_of(&self, origin: ShapeOrigin) -> Vec<usize> {
        self.ds
            .edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.origin == origin)
            .map(|(i, _)| i)
            .collect()
    }
    // OCCT BOPAlgo_PaveFiller_7.cxx L371-549 MakeSplitEdges
    pub(crate) fn make_split_edges(&mut self) {
        // L391: UpdateCommonBlocksWithSDVertices
        self.ds.update_common_blocks_with_sd_vertices();
        //
        let n_edges = self.ds.edge_count();
        // L377-380: return if no edges (aNbPBP == 0 equivalent)
        if n_edges == 0 {
            return;
        }
        //
        // L386: aMCB — dedup set for CommonBlocks (store CB indices)
        let mut a_mcb: HashSet<usize> = HashSet::new();
        //
        // ----- Phase 1 (L396-498): collect split edge tasks -----
        // ----- Phase 1 (L396-498): collect split edge tasks -----
        //
        struct SplitEdgeTask {
            edge_idx: usize,
            n_v1: usize,
            n_v2: usize,
            a_t1: f64,
            a_t2: f64,
            // For non-CB: the PB to call SetEdge on; for CB: PaveBlock1
            pb_shared: crate::bopds::pave::SharedPB,
            cb_idx: Option<usize>,
        }
        //
        struct DeferredCBUpdate {
            cb_idx: usize,
            edge_idx: usize,
        }
        //
        let mut tasks: Vec<SplitEdgeTask> = Vec::new();
        let mut deferred_cb_updates: Vec<DeferredCBUpdate> = Vec::new();
        //
        for ei in 0..n_edges {
            // L398-401: UserBreak check omitted in rcad (no progress range)
            //
            // aLPB = aPBP(i) — the list of PBs for this edge (L402)
            let n_e = ei;
            //
            // L408-410: aSIE.HasFlag() → skip degenerated edges
            if self.ds.is_edge_degenerated(n_e) {
                continue;
            }
            //
            // L402: aLPB = self.ds.edges[n_e].pave_blocks
            let pb_count = self.ds.edges[n_e].pave_blocks.len();
            //
            // Clone PB references for iteration (avoids borrow during mutable Phase 1.5)
            let edge_pbs: Vec<crate::bopds::pave::SharedPB> =
                self.ds.edges[n_e].pave_blocks.clone();
            //
            for spb in &edge_pbs {
                // L407-423: extract PB data with short-lived guard (avoids borrow conflicts)
                let pb_orig_edge: usize;
                let pb_cb_idx: Option<usize>;
                let b_cb: bool;
                let n_v1: usize;
                let n_v2: usize;
                {
                    let pb = spb.0.read().unwrap();
                    pb_orig_edge = pb.original_edge;
                    pb_cb_idx = pb.common_block_idx;
                    b_cb = pb_cb_idx.is_some();
                    let (nv1, nv2) = pb.indices();
                    n_v1 = nv1;
                    n_v2 = nv2;
                } // pb guard dropped
                //
                // L409-414: skip if original edge is degenerated
                if self.ds.is_edge_degenerated(pb_orig_edge) {
                    continue;
                }
                //
                // L416-421: Check CommonBlock, dedup CBs
                if b_cb && !a_mcb.insert(pb_cb_idx.unwrap()) {
                    // CB already processed → skip (OCCT: aMCB.Add returns false)
                    continue;
                }
                //
                // L425-468: Check if split is necessary
                let mut b_to_split = true;
                let b_v1 = self.ds.is_new_vertex(n_v1);
                let b_v2 = self.ds.is_new_vertex(n_v2);
                //
                if !b_v1 && !b_v2 {
                    // L430: no new vertices — may avoid splitting
                    if !self.non_destructive || !b_cb {
                        if b_cb {
                            // L436-455: CB — find the PB whose original edge has 1 PB
                            let cb_idx = pb_cb_idx.unwrap();
                            let cb_pbs: Vec<(usize, usize)> = {
                                let cb = &self.ds.common_blocks[cb_idx];
                                cb.pave_blocks().to_vec()
                            };
                            //
                            let mut a_found_it = false;
                            let mut a_found_n_e = n_e;
                            for &(cb_pb_idx, _) in &cb_pbs {
                                if cb_pb_idx < self.ds.pave_blocks.len() {
                                    let oe = self.ds.pave_blocks[cb_pb_idx]
                                        .0
                                        .read()
                                        .unwrap()
                                        .original_edge;
                                    if oe < self.ds.edge_count()
                                        && self.ds.edges[oe].pave_blocks.len() == 1
                                    {
                                        a_found_it = true;
                                        a_found_n_e = oe;
                                        break;
                                    }
                                }
                            }
                            //
                            if a_found_it {
                                // L447-456: edge has only 1 PB → no split needed
                                b_to_split = false;
                                // L449: aCB->SetRealPaveBlock(it.Value())
                                // Reorder CB PBs so the matched PB is first.
                                if let Some(pos) = cb_pbs.iter().position(|&(cb_pb_idx, _)| {
                                    cb_pb_idx < self.ds.pave_blocks.len()
                                        && self.ds.pave_blocks[cb_pb_idx]
                                            .0
                                            .read()
                                            .unwrap()
                                            .original_edge
                                            == a_found_n_e
                                }) {
                                    if pos != 0 {
                                        let mut reordered = cb_pbs;
                                        reordered.swap(0, pos);
                                        self.ds.common_blocks[cb_idx].set_pave_blocks(reordered);
                                    }
                                }
                                // L450-454: SetEdge + ComputeTolerance + UpdateEdgeTol
                                deferred_cb_updates.push(DeferredCBUpdate {
                                    cb_idx,
                                    edge_idx: a_found_n_e,
                                });
                            }
                        } else if pb_count == 1 {
                            // L457-461: single PB, no new vertices → keep original edge
                            b_to_split = false;
                            // L460: aPB->SetEdge(nE) — short-lived guard for write
                            let orig = {
                                let p = spb.0.read().unwrap();
                                p.original_edge
                            };
                            spb.0.write().unwrap().new_edge = Some(orig);
                        }
                        // L462-465: if !bToSplit → skip (continue)
                        if !b_to_split {
                            continue;
                        }
                    }
                }
                //
                // L470-496: This PB needs splitting
                let (task_edge_idx, task_n_v1, task_n_v2, task_t1, task_t2, task_pb) = {
                    if b_cb {
                        // L471-476: use CB's PaveBlock1
                        let cb_idx = pb_cb_idx.unwrap();
                        let cb = &self.ds.common_blocks[cb_idx];
                        if let Some(pb1_idx) = cb.pave_block1() {
                            if pb1_idx < self.ds.pave_blocks.len() {
                                let pb1 = &self.ds.pave_blocks[pb1_idx];
                                let pb1g = pb1.0.read().unwrap();
                                let w_edge = pb1g.original_edge;
                                let (w_v1, w_v2) = pb1g.indices();
                                let (w_t1, w_t2) = pb1g.range();
                                drop(pb1g);
                                (w_edge, w_v1, w_v2, w_t1, w_t2, pb1.clone())
                            } else {
                                continue;
                            }
                        } else {
                            continue;
                        }
                    } else {
                        // L477: use current PB's data — short-lived guard
                        let (t1, t2) = {
                            let p = spb.0.read().unwrap();
                            p.range()
                        };
                        (pb_orig_edge, n_v1, n_v2, t1, t2, spb.clone())
                    }
                };
                //
                // L488-496: create SplitEdge task
                tasks.push(SplitEdgeTask {
                    edge_idx: task_edge_idx,
                    n_v1: task_n_v1,
                    n_v2: task_n_v2,
                    a_t1: task_t1,
                    a_t2: task_t2,
                    pb_shared: task_pb,
                    cb_idx: pb_cb_idx,
                });
            }
        }
        // OCCT: process intersection curve PBs (pool continuation after edges)
        // ----- Phase 1.5: Apply deferred CB tolerance updates -----
        // (OCCT L449-454: CB edge with 1 PB — no split, just tolerance)
        for update in &deferred_cb_updates {
            // L450: aCB->SetEdge(nE)
            self.ds.common_blocks[update.cb_idx].set_edge(update.edge_idx);
            // L453: ComputeToleranceOfCB
            let a_tol = crate::bopds::tools::compute_tolerance_of_cb(self.ds, update.cb_idx);
            // L454: UpdateEdgeTolerance
            self.update_edge_tolerance(update.edge_idx, a_tol);
        }
        //
        // ----- Phase 2 (L500-548): create split edges and register in DS -----
        // OCCT L500-508: BOPTools_Parallel::Perform (sequential in rcad)
        //
        for task in &tasks {
            // L521-523: get split edge and box from BOPAlgo_SplitEdge result
            let orig = &self.ds.edges[task.edge_idx];
            //
            // Build vertex_params map
            let mut vertex_params = std::collections::HashMap::new();
            vertex_params.insert(task.n_v1, task.a_t1);
            vertex_params.insert(task.n_v2, task.a_t2);
            //
            // L529-537: create ShapeInfo for the new edge and Append to DS
            let n_sp = self.ds.push_edge(
                DSEdge { shape_idx: 0,
                    start_vertex: task.n_v1,
                    end_vertex: task.n_v2,
                    curve: orig.curve.clone(),
                    t_range: [task.a_t1, task.a_t2],
                    origin: orig.origin,
                    geom_tol: orig.geom_tol,
                    paves: Vec::new(),
                    pave_blocks: Vec::new(),
                    face_reps: Vec::new(),
                    is_internal: orig.is_internal,
                    vertex_params,
                    face_tolerances: Vec::new(),
                    is_geometric: orig.is_geometric,
                    location: orig.location,
                },
                None,
            );
            //
            // L539-547: register new edge on CB or PB
            if let Some(cb_idx) = task.cb_idx {
                // OCCT L541: UpdateEdgeTolerance(nSp, aBSE.Tolerance())
                // OCCT's BOPAlgo_SplitEdge computes tolerance from actual
                // vertex-to-curve distance at the split parameters.
                let (v1_p, v2_p, orig_curve, orig_tol) = {
                    let o = &self.ds.edges[task.edge_idx];
                    (
                        self.ds.vertex_point(task.n_v1),
                        self.ds.vertex_point(task.n_v2),
                        o.curve.clone(),
                        o.geom_tol,
                    )
                };
                let (_t1, cp1) = crate::extrema::closest_point_on_curve(&orig_curve, v1_p);
                let (_t2, cp2) = crate::extrema::closest_point_on_curve(&orig_curve, v2_p);
                let v1_dist = (v1_p - cp1).length();
                let v2_dist = (v2_p - cp2).length();
                let a_tol_e = orig_tol.max(crate::tolerance::CONFUSION);
                let a_tol_v1 = self.ds.vertex_tolerance(task.n_v1);
                let a_tol_v2 = self.ds.vertex_tolerance(task.n_v2);
                let se_tol = a_tol_e.max(a_tol_v1 + v1_dist).max(a_tol_v2 + v2_dist);
                self.update_edge_tolerance(n_sp, se_tol);
                // OCCT L542: aCBk->SetEdge(nSp)
                self.ds.common_blocks[cb_idx].set_edge(n_sp);
            } else {
                // L546: aPBk->SetEdge(nSp)
                task.pb_shared.0.write().unwrap().new_edge = Some(n_sp);
            }
        }
    }

    // OCCT BOPAlgo_PaveFiller_11.cxx L1-126 CheckSelfInterference
    pub(crate) fn check_self_interference(&self) -> Vec<crate::bopalgo::Alert> {
        if self.my_arguments.len() <= 1 {
            return Vec::new();
        }

        let mut a_alerts: Vec<crate::bopalgo::Alert> = Vec::new();

        for a_rank in 0..2 {
            let mut a_mcsi: std::collections::HashMap<usize, indexmap::IndexSet<usize>> =
                std::collections::HashMap::new();
            let mut a_cb_fence: std::collections::HashSet<usize> = std::collections::HashSet::new();

            // Process EDGES from this operand
            for ei in 0..self.ds.edge_count() {
                let e_origin = self.ds.edge_origin(ei);
                let e_rank = match e_origin {
                    ShapeOrigin::ShapeA => 0,
                    ShapeOrigin::ShapeB => 1,
                };
                if e_rank != a_rank {
                    continue;
                }
                if self.ds.edge_pave_blocks(ei).is_empty() {
                    continue;
                }
                if self.ds.edge_has_flag(ei) {
                    continue;
                }

                // Sub-shape vertices with SD resolution
                let sv = self.ds.edge_start_vertex_ds(ei);
                let ev = self.ds.edge_end_vertex_ds(ei);
                let mut a_sub_s: std::collections::HashSet<usize> =
                    std::collections::HashSet::new();
                for n_v in [sv, ev] {
                    let n_v = self.ds.has_shape_sd(n_v).unwrap_or(n_v);
                    a_sub_s.insert(n_v);
                }

                let a_lpb = self.ds.edge_pave_blocks(ei);
                let b_analyze_v = a_lpb.len() > 1;

                for spb in a_lpb {
                    let pb = spb.0.read().unwrap();

                    if b_analyze_v {
                        let (nv1, nv2) = pb.indices();
                        for &n_v in &[nv1, nv2] {
                            let n_v = self.ds.has_shape_sd(n_v).unwrap_or(n_v);
                            let v_in_range = match self.ds.vertex_origin(n_v) {
                                Some(ShapeOrigin::ShapeA) => a_rank == 0,
                                Some(ShapeOrigin::ShapeB) => a_rank == 1,
                                None => false,
                            };
                            if !v_in_range && !a_sub_s.contains(&n_v) {
                                a_mcsi.entry(n_v).or_default().insert(ei);
                            }
                        }
                    }

                    if let Some(cb_idx) = pb.common_block_idx {
                        if a_cb_fence.insert(cb_idx) {
                            if let Some(cb) = self.ds.common_blocks.get(cb_idx) {
                                let mut a_le: Vec<usize> = Vec::new();
                                for &(pb_gi, _) in cb.pave_blocks() {
                                    if pb_gi < self.ds.pave_blocks.len() {
                                        let n_e_or = self.ds.pave_blocks[pb_gi]
                                            .0
                                            .read()
                                            .unwrap()
                                            .original_edge;
                                        let eo = self.ds.edge_origin(n_e_or);
                                        let eor_rank = match eo {
                                            ShapeOrigin::ShapeA => 0,
                                            ShapeOrigin::ShapeB => 1,
                                        };
                                        if eor_rank == a_rank {
                                            a_le.push(n_e_or);
                                        }
                                    }
                                }
                                if a_le.len() > 1 {
                                    a_alerts.push(crate::bopalgo::Alert::AcquiredSelfIntersection(
                                        a_le,
                                    ));
                                }
                            }
                        }
                    }
                }
            }

            // Process FACES from this operand
            for fi in 0..self.ds.face_count() {
                let f_origin = self.ds.face_origin(fi);
                let f_rank = match f_origin {
                    ShapeOrigin::ShapeA => 0,
                    ShapeOrigin::ShapeB => 1,
                };
                if f_rank != a_rank {
                    continue;
                }

                let a_fi = self.ds.face_info(fi);

                for vertices_set in [&a_fi.vertices_in, &a_fi.vertices_sc] {
                    for &n_v in vertices_set {
                        a_mcsi.entry(n_v).or_default().insert(fi);
                    }
                }
                for pb_set in [&a_fi.pave_blocks_in, &a_fi.pave_blocks_sc] {
                    for &pb_gi in pb_set {
                        if pb_gi < self.ds.pave_blocks.len() {
                            let n_e = {
                                let pb = self.ds.pave_blocks[pb_gi].0.read().unwrap();
                                pb.new_edge.unwrap_or(pb.original_edge)
                            };
                            a_mcsi.entry(n_e).or_default().insert(fi);
                        }
                    }
                }
            }

            // Analyze connections
            for (_sub_shape, shapes) in &a_mcsi {
                if shapes.len() > 1 {
                    a_alerts.push(crate::bopalgo::Alert::AcquiredSelfIntersection(
                        shapes.iter().copied().collect(),
                    ));
                }
            }
        }

        a_alerts
    }
}
