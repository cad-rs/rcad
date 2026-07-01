use std::collections::HashSet;

use glam::{DVec2, DVec3};
use rcad_kernel::geom::*;

use crate::bopds::ds::{
    DS, DSEdge, DSRepOnFace, Interference, IntersectionCurve, ShapeOrigin,
};
use crate::bopds::pave::*;
use crate::bvh::Bvh;
use crate::inttools;
use crate::inttools::context::Context as IntToolsContext;
use crate::inttools::fclass2d::{FClass2d, State};
use crate::tolerance::*;
use rcad_kernel::closest_point_on_curve;
pub mod helpers;
use self::helpers::*;

/// 锟?OCCT-aligned: IntPatch_Intersection surface category (L1264-1294).
///   GeomGeom  = ts1==ts2==1 锟?ImpImpIntersection (analytic-analytic)
///   GeomParam = ts1!=ts2     锟?ImpPrmIntersection (analytic-parametric)
///   ParamParam = ts1==ts2==0 锟?PrmPrmIntersection (parametric-parametric)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SurfaceCategory { GeomGeom, ParamParam }

fn classify_surface_type(surf: &Surface3) -> SurfaceCategory {
    match surf {
        Surface3::Plane(_) | Surface3::Cylinder(_) | Surface3::Sphere(_)
        | Surface3::Cone(_) | Surface3::Torus(_) => SurfaceCategory::GeomGeom,
        _ => SurfaceCategory::ParamParam,
    }
}

// Re-export NearTangentType from bopds::ds for use in this module's public types
pub use crate::bopds::ds::NearTangentType;

/// Minimum total face count before BVH acceleration is used.
/// Below this threshold, brute-force O(n虏) is faster due to BVH build overhead.
const BVH_THRESHOLD: usize = 20;

/// 鉁?OCCT-aligned: BOPAlgo_PaveFiller 鈥?six intersection passes
///   (PaveFiller.hxx L106-107, PaveFiller.cxx L234-355).
mod glue;
mod intersection;
pub(crate) mod analytics;
pub(crate) mod analytic_plane;
pub(crate) mod analytic_sphere;
pub(crate) mod analytic_cylinder;
pub(crate) mod analytic_cone;
pub(crate) mod analytic_torus;
pub(crate) mod marching;
pub(crate) mod prm_prm_intersection;
pub(crate) mod p_walking;
mod make_blocks;
mod config;
mod tolerances;
        mod interf;
mod paves;
mod ff_intersect;

/// OCCT-aligned: BOPAlgo_SectionAttribute 鈥?controls approximation and
/// pcurve computation for section edges (BOPAlgo_SectionAttribute.hxx).
#[derive(Debug, Clone)]
pub(crate) struct SectionAttribute {
    pub approximation: bool,
    pub pcurve_on_s1: bool,
    pub pcurve_on_s2: bool,
}

impl Default for SectionAttribute {
    fn default() -> Self {
        Self { approximation: true, pcurve_on_s1: true, pcurve_on_s2: true }
    }
}

/// OCCT-aligned: PaveFiller::EdgeRangeDistance 鈥?stores minimal distance
/// between an edge range and a face that don't geometrically intersect.
/// Used by PostTreatFF to re-check E-F pairs after tolerance updates.
#[derive(Debug, Clone)]
pub(crate) struct EdgeRangeDistance {
    pub first: f64,
    pub last: f64,
    pub distance: f64,
}

pub struct PaveFiller<'a> {
    pub ds: &'a mut DS,
    /// Optional BRep for direct output (dual-write mode). When set, PaveFiller
    /// populates the BRep on completion, eliminating the need for ds_to_brep.
    pub brep: Option<&'a mut rcad_kernel::topods::BRep>,
    /// Output: face_refs by ds_face_idx (populated by export_to_brep).
    pub face_refs: Vec<rcad_kernel::topods::ShapeRef>,
    /// Output: ic_edge_map: ci -> BRep edge ShapeRef (populated by export_to_brep).
    pub ic_edge_map: Vec<Option<rcad_kernel::topods::ShapeRef>>,
    bvh_a: Option<&'a Bvh>,
    bvh_b: Option<&'a Bvh>,
    use_glue: bool,
    glue_tolerance: f64,
    /// 锟?OCCT-aligned: BOPAlgo_Options::SetFuzzyValue
    fuzzy_tolerance: f64,
    /// 锟?OCCT-aligned: PaveFiller_6.cxx L393-479 seam edge shift tolerance
    seam_shift_tol: f64,
    /// 锟?OCCT-aligned: BOPAlgo_Algo::myRunParallel
    run_parallel: bool,
    /// 锟?OCCT-aligned: BOPAlgo_PaveFiller::myNonDestructive
    non_destructive: bool,
    /// 锟?OCCT-aligned: BOPAlgo_Algo::myUseOBB
    use_obb: bool,
    /// 锟?OCCT-aligned: IntTools_Context (PaveFiller::Init L203)
    context: IntToolsContext,
    /// 锟?OCCT-aligned: myArguments 鈥?original input shapes (BOPAlgo_PaveFiller.hxx L639).
    ///   rcad: carries the original BRep operands for OCCT-API compatibility.
    my_arguments: Vec<rcad_kernel::BRep>,
    /// 锟?OCCT-aligned: mySectionAttribute (BOPAlgo_SectionAttribute.hxx)
    section_attribute: SectionAttribute,
    /// 锟?OCCT-aligned: myIsPrimary (BOPAlgo_PaveFiller.cxx L62)
    is_primary: bool,
    /// 锟?OCCT-aligned: myAvoidBuildPCurve (BOPAlgo_PaveFiller.cxx L63)
    avoid_build_pcurve: bool,
    /// 锟?OCCT-aligned: myFPBDone 鈥?fence map tracking processed (face, pave_block) pairs.
    ///   Map: face_idx 鈫?set of pave_block indices already processed in PostTreatFF.
    fpbdone: std::collections::HashMap<usize, std::collections::HashSet<usize>>,
    /// 锟?OCCT-aligned: myVertsToAvoidExtension 鈥?vertices that should NOT have
    ///   their tolerance extended further (near EE/EF intersection points).
    verts_to_avoid_extension: std::collections::HashSet<usize>,
    /// 锟?OCCT-aligned: myDistances 鈥?minimal edge-face distances for non-intersecting
    ///   pairs.  Map: (edge_idx, face_idx) 鈫?Vec<EdgeRangeDistance>.
    distances: std::collections::HashMap<(usize, usize), Vec<EdgeRangeDistance>>,
}

/// 锟?OCCT-aligned:Propagate IC vertices to all faces sharing boundary edges
///    (OCCT BOPDS_FaceInfo::AppendBlock equivalent).
///    OCCT BOPAlgo_PaveFiller propagates pave block vertices to all faces
///    referencing the split edge. rcad's add_curve only adds vertices to the
///    two FF-interference faces (f1, f2), but the vertex may lie on boundary
///    edges of other faces (e.g. side-face tangent-line IC endpoints on the
///    top face's boundary edge).
fn propagate_ic_vertices_to_shared_faces(
    ds: &mut DS,
    ic_vertices: &[usize],
    skip_faces: &[usize; 2],
) {
    let vtol = TOLERANCE_ABS * 1000.0; // 1e-4 geometric tolerance for on-edge check
    let vtol_sq = vtol * vtol;
    for fi in 0..ds.faces.len() {
        if fi == skip_faces[0] || fi == skip_faces[1] {
            continue;
        }
        for &vi in ic_vertices {
            if ds.faces[fi].face_info.vertices_in.contains(&vi) {
                continue;
            }
            let vp = ds.vertices[vi].point;
            for &ei in &ds.faces[fi].boundary_edges {
                let Some(edge) = ds.edges.get(ei) else { continue };
                let a = ds.vertices[edge.start_vertex].point;
                let b = ds.vertices[edge.end_vertex].point;
                let ab = b - a;
                let ab_len2 = ab.length_squared();
                if ab_len2 < 1e-30 {
                    continue;
                }
                let ap = vp - a;
                let t = ap.dot(ab) / ab_len2;
                if t > -0.01 && t < 1.01 {
                    let proj = a + ab * t.clamp(0.0, 1.0);
                    if (vp - proj).length_squared() < vtol_sq {
                        ds.faces[fi].face_info.vertices_in.insert(vi);
                        break;
                    }
                }
            }
        }
    }
}

impl<'a> PaveFiller<'a> {

    /// 鉁?OCCT-aligned: Prepare (PaveFiller_7.cxx L850-929).
    ///   Build 2D pcurves for edges on planar faces.
    ///   OCCT: iterate all V/E, E/E, E/F pairs to find planar faces, collect edge-face pairs,
    ///   compute pcurves in parallel, update edges.  rcad: DS::build_face_reps already computes
    ///   pcurves for all face types; this step ensures planar-face pcurves exist for any
    ///   edges that build_face_reps may have missed (non-boundary or intersection-relevant).
    fn prepare(&mut self) {
        let mut planar_faces: Vec<usize> = Vec::new();
        for (fi, f) in self.ds.faces.iter().enumerate() {
            if matches!(f.surface, Surface3::Plane(_)) {
                planar_faces.push(fi);
            }
        }
        if planar_faces.is_empty() { return; }

        let surf: Vec<Surface3> = planar_faces.iter().map(|&fi| self.ds.faces[fi].surface.clone()).collect();

        // Collect all edge indices from each planar face's boundary
        let mut face_edges: Vec<Vec<usize>> = Vec::with_capacity(planar_faces.len());
        for (pos, &fi) in planar_faces.iter().enumerate() {
            let f = &self.ds.faces[fi];
            let mut eids: Vec<usize> = f.boundary_edges.clone();
            for w in &f.inner_boundary_edges {
                eids.extend(w.iter().map(|&(ei, _)| ei));
            }
            face_edges.push(eids);
        }

        // Compute pcurves for edge-face pairs that don't already have one
        for (pos, &fi) in planar_faces.iter().enumerate() {
            for &ei in &face_edges[pos] {
                if self.ds.edge_on_face(ei, fi).is_some() { continue; }
                let Some(edge) = self.ds.edges.get_mut(ei) else { continue; };
                if let Some((pcurve, span)) = DS::compute_edge_pcurve(&edge.curve, &surf[pos]) {
                    edge.face_reps.push(DSRepOnFace {
                        face_idx: fi,
                        pcurve,
                        pcurve2: None,
                        pcurve_range: [0.0, span],
                        start_param: 0.0,
                        end_param: span,
                    });
                }
            }
        }
    }


    // 鉁?OCCT L248 Prepare: build pcurves on planar faces
    // OCCT L234-355 PerformInternal: Init->Prepare->VV->VE->EE->VF->EF->
    //   RepeatInt->ForceEE->ForceEF->FF->UpdBlk->RefFI->MkSEdges->MkBlks->
    //   ChkSI->RefFO->RmvME->MkPCurves->ProcDE
    pub fn perform(&mut self) {
        // 鉁?OCCT L248: Prepare 鈥?build pcurves on planar faces.
        self.prepare();

        // Detect shared topology before interference passes when glue is enabled
        if self.use_glue {
            self.ds.detect_shared_topology(self.glue_tolerance);
        }

        // Skip redundant interference passes when glue is enabled and shared topology is detected
        let skip_ve = self.should_skip_ve_pass();
        let skip_ee = self.should_skip_ee_pass();
        let skip_vf = self.should_skip_vf_pass();
        let skip_ef = self.should_skip_ef_pass();
        let skip_ff = self.should_skip_ff_pass();

        self.perform_vv();

        // OCCT L145: BOPDS_Iterator 鈥?single combined BVH per shape type.
        let bvh_all_verts = self.build_ds_bvh_combined(false);
        let bvh_all_edges = self.build_ds_bvh_combined(true);
        let bvh_all_faces = self.build_ds_bvh_face_all();

        if !skip_ve {
            // OCCT: BOPDS_Iterator::Initialize(VERTEX, EDGE) 鈥?single traversal.
            self.perform_ve_bvh(&bvh_all_verts, &bvh_all_edges);
        }
        // 锟?OCCT-aligned: UpdatePaveBlocksWithSDVertices (PerformInternal L266)
        self.ds.update_pave_blocks_with_sd_vertices();

        let ee_survivors: Vec<usize> = if !skip_ee {
            // OCCT L145-267: PerformEE — build intersections
            let ee_modified = self.perform_ee_bvh(&bvh_all_edges, &bvh_all_edges);
            // OCCT L558-565: PerformCommonBlocks + PerformNewVertices
            let survivors = self.treat_new_vertices();
            // OCCT L571-585: SplitPaveBlocks for remaining modified edges
            //   (edges with new vertices already handled by treat_new_vertices)
            if !ee_modified.is_empty() {
                let sd_vertices: std::collections::HashSet<usize> =
                    self.ds.shape_sd.sd_vertices_iter()
                        .map(|&(a, _)| a)
                        .collect();
                // Edges to split = modified - edges whose new vertex is in SD
                let remaining: std::collections::HashSet<usize> = ee_modified.iter()
                    .filter(|&&ei| {
                        let e = &self.ds.edges[ei];
                        !sd_vertices.contains(&e.start_vertex) &&
                        !sd_vertices.contains(&e.end_vertex)
                    })
                    .copied()
                    .collect();
                if !remaining.is_empty() {
                    self.split_pave_blocks(&remaining, false);
                }
            }
            // OCCT L273: UpdatePaveBlocksWithSDVertices
            self.ds.update_pave_blocks_with_sd_vertices();
            survivors
        } else { vec![] };

        if !skip_vf {
            // OCCT: BOPDS_Iterator::Initialize(VERTEX, FACE) 鈥?single traversal.
            self.perform_vf_bvh(&bvh_all_verts, &bvh_all_faces);
        }
        // 锟?OCCT-aligned: UpdatePaveBlocksWithSDVertices (PerformInternal L280)
        self.ds.update_pave_blocks_with_sd_vertices();

        if !skip_ef {
            // OCCT-aligned: PerformEF uses BOPDS_Iterator for pair enumeration,
            //   not a pre-built BVH.  rcad perform_ef iterates all A脳B pairs.
            self.perform_ef();
            // 锟?OCCT-aligned: TreatNewVertices 锟?merge new vertices created by EF intersection.
            //    OCCT PaveFiller_5.cxx L570: PerformNewVertices(aMVCPB, ..., false)
            let ef_survivors = self.treat_new_vertices();

            // 锟?OCCT-aligned: RepeatIntersection (PaveFiller.cxx L296-299, L359-420).
            //    After EF, before FF, re-run VV/VE/VF for vertices with increased tolerance.
            //    OCCT reads from myIncreasedSS (populated by TreatNewVertices).
            //    rcad: ds.increased_ss is populated by treat_new_vertices above.
            self.ds.update_pave_blocks_with_sd_vertices();
            self.update_interfs_with_sd_vertices();
            self.repeat_intersection();
        }

        // 锟?OCCT-aligned: ForceInterfEE (PaveFiller_3.cxx L978-1276)
        //    OCCT L302: ForceInterfEE 锟?after RepeatIntersection, force intersection
        //    of edge pairs sharing a vertex with increased tolerance, detecting
        //    collinear/coincident edges (common block).
        //    锟?rcad: simplified, only checks line-line edge pairs sharing a pave vertex.
        if !skip_ee {
            self.force_interf_ee();
        }

        // 锟?OCCT-aligned: ForceInterfEF (PaveFiller_5.cxx L764-1099+)
        //    OCCT L309: ForceInterfEF 锟?after ForceInterfEE, force intersection of
        //    edges whose both endpoints are on a face with increased tolerance.
        //    锟?rcad: simplified, only checks edge-face pairs where both endpoints are on the face.
        if !skip_ef {
            self.force_interf_ef();
        }

        if !skip_ff {
            self.perform_ff();

            // 鉁?OCCT-aligned: InitPaveBlock1 for each IC (BOPDS_Curve::InitPaveBlock1
            //   is called during PerformFF in OCCT).  rcad ICs are created with empty
            //   pave_blocks; this gives them the default PB needed by
            //   make_section_edges (inside make_blocks).
            for ci in 0..self.ds.intersection_curves.len() {
                self.ds.intersection_curves[ci].init_pave_block1();
            }

            // 锟?OCCT-aligned: MakeSDVerticesFF (PaveFiller_6.cxx L1113)
            //    After FF, create shared SD vertices for same-domain (coplanar) face
            //    overlap boundaries so that overlap polygon vertices are shared between
            //    both faces and registered in face_info.vertices_in.
            self.make_sd_vertices_ff();
        }

        // 鉁?OCCT L318: UpdateBlocksWithSharedVertices
        self.update_blocks_with_shared_vertices();

        // 鉁?OCCT L320: RefineFaceInfoIn
        for fi in 0..self.ds.faces.len() {
            self.ds.refine_face_info_in(fi);
        }

        // 鉁?OCCT L322: MakeSplitEdges
        self.make_split_edges();

        // 鉁?OCCT L328: UpdatePaveBlocksWithSDVertices
        self.ds.update_pave_blocks_with_sd_vertices();

        // 鉁?OCCT L330: MakeBlocks
        self.make_blocks();

        // 鉁?OCCT L336: CheckSelfInterference (BOPAlgo_PaveFiller_11.cxx L28-221).
        //    OCCT uses AddWarning 鈥?non-fatal, the operation continues.
        let _ = self.check_self_interference();

        // 鉁?OCCT-aligned: UpdateInterfsWithSDVertices (PerformInternal L338)
        self.update_interfs_with_sd_vertices();

        // 鉁?OCCT-aligned: ReleasePaveBlocks (PerformInternal L339) 鈥?OCCT frees
        //   UNUSED PBs.  rcad: PBs remain in pool (PaveBlocksSc indices stay valid).
        //   Clearing the pool here would invalidate PaveBlocksSc indices.

        // 鉁?OCCT-aligned: RefineFaceInfoOn 鈥?after ReleasePaveBlocks, remove
        //    zero-length On pave blocks (PerformInternal L340, BOPDS_DS::RefineFaceInfoOn).
        for fi in 0..self.ds.faces.len() {
            self.ds.refine_face_info_on(fi);
        }

        // 锟?OCCT-aligned: RemoveMicroEdges 锟?after MakeBlocks, before MakePCurves
        //    (PerformInternal L342, PaveFiller_6.cxx L4229-4270).
        self.remove_micro_edges();

        // 锟?OCCT-aligned: MakePCurves 锟?after RemoveMicroEdges (PerformInternal L344)
        self.make_pcurves();

        // 鉁?OCCT-aligned: ProcessDE 鈥?after MakePCurves (PerformInternal L350)
        self.process_de();

        // Export to BRep if direct output is enabled (A3 dual-write).
        if let Some(ref mut brep) = self.brep {
            let (face_refs, ic_edge_map) = crate::ds_to_brep::export_to_brep(&self.ds, brep);
            self.face_refs = face_refs;
            self.ic_edge_map = ic_edge_map;
        }
    }

    // ===== BVH-based pair enumeration (OCCT BOPDS_Iterator) =====













    // 閳光偓閳光偓閳光偓 Pass 1: Vertex-Vertex 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓


    // 閳光偓閳光偓閳光偓 Pass 2: Vertex-Edge 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓



    // 閳光偓閳光偓閳光偓 Pass 3: Edge-Edge 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓




    // 閳光偓閳光偓閳光偓 Pass 4: Vertex-Face 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓





// 閳光偓閳光偓閳光偓 Pass 5: Edge-Face 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓






    // 閳光偓閳光偓閳光偓 Pass 6: Face-Face 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓


    fn update_blocks_with_shared_vertices(&mut self) {
        // OCCT L3948-3951: non-destructive guard
        if !self.non_destructive { return; }

        // OCCT L3953-3960: check FF interferences
        let has_ff = self.ds.interferences.iter().any(|inf|
            matches!(inf, Interference::FaceFace { curves, .. } if !curves.is_empty())
        );
        if !has_ff { return; }

        // Collect face pairs with shared (old) vertices
        // OCCT L3967-4049: iterate each FF
        let ff_entries: Vec<(usize, usize, Vec<usize>)> = self.ds.interferences.iter()
            .filter_map(|inf| {
                if let Interference::FaceFace { f1, f2, curves, .. } = inf {
                    if curves.is_empty() { return None; }

                    // OCCT L3996-4017: collect shared old vertices
                    let fi1 = *f1;
                    let fi2 = *f2;
                    let on1 = &self.ds.faces[fi1].face_info.vertices_on;
                    let in1 = &self.ds.faces[fi1].face_info.vertices_in;
                    let on2 = &self.ds.faces[fi2].face_info.vertices_on;
                    let in2 = &self.ds.faces[fi2].face_info.vertices_in;

                    let shared: Vec<usize> = on1.iter()
                        .chain(in1.iter())
                        .filter(|&&vi| {
                            // OCCT L4008-4011: skip new vertices (only old shapes)
                            if self.ds.is_new_vertex(vi) { return false; }
                            // OCCT L4012: must be ON or IN in the other face too
                            on2.contains(&vi) || in2.contains(&vi)
                        })
                        .copied()
                        .collect();

                    if shared.is_empty() { return None; }
                    Some((fi1, fi2, curves.clone()))
                } else { None }
            })
            .collect();

        // OCCT L4020-4048: for each FF entry, try shared vertices on each curve
        for (f1, f2, curves) in &ff_entries {
            // OCCT L3980-3987: UpdateFaceInfoOn equivalent
            //   rcad: not needed 锟?FaceInfo data is already populated.

            for &ci in curves {
                if ci >= self.ds.intersection_curves.len() { continue; }

                // OCCT L4023: aTolR3D = max(curve.Tolerance(), curve.TangentialTolerance())
                let ic = &self.ds.intersection_curves[ci];
                let _a_tol_r3d = ic.geom_tol;
                let f1 = *f1;
                let f2 = *f2;

                // Collect shared vertices for this face pair
                let on1 = &self.ds.faces[f1].face_info.vertices_on;
                let in1 = &self.ds.faces[f1].face_info.vertices_in;
                let on2 = &self.ds.faces[f2].face_info.vertices_on;
                let in2 = &self.ds.faces[f2].face_info.vertices_in;

                let shared: Vec<usize> = on1.iter()
                    .chain(in1.iter())
                    .filter(|&&vi| {
                        if self.ds.is_new_vertex(vi) { return false; }
                        // OCCT L4012: present in the other face's On/In sets
                        if !on2.contains(&vi) && !in2.contains(&vi) { return false; }
                        // OCCT L4030-4034: skip if already has SD mapping
                        //   rcad: check shape_sd
                        if self.ds.shape_sd.is_sub_vertex(vi) { return false; }
                        true
                    })
                    .copied()
                    .collect();

                for &n_v in &shared {
                    // OCCT L4036: EstimatePaveOnCurve
                    if self.estimate_pave_on_curve(ci, n_v).is_none() { continue; }

                    // OCCT L4042-4046: UpdateVertex + InitPaveBlocksForVertex
                    let v_tol = self.ds.vertices[n_v].geom_tol;
                    // UpdateVertex: increase tolerance if the projection distance is larger
                    if let Some(t) = self.project_vertex_on_curve(n_v, ic) {
                        let pt_on_curve = ic.curve.point_at(t);
                        let dist = self.ds.vertices[n_v].point.distance(pt_on_curve);
                        if dist > v_tol {
                            self.ds.vertices[n_v].geom_tol = dist;
                            self.ds.increased_ss.insert(n_v);
                        }
                    }
                    // InitPaveBlocksForVertex: collect edge indices + params, then apply
                    let mut new_paves: Vec<(usize, f64)> = Vec::new();
                    for (ei, edge) in self.ds.edges.iter().enumerate() {
                        if edge.start_vertex == n_v {
                            let has = edge.paves.iter().any(|p| p.vertex_idx == n_v);
                            if !has { new_paves.push((ei, edge.t_range[0])); }
                        } else if edge.end_vertex == n_v {
                            let has = edge.paves.iter().any(|p| p.vertex_idx == n_v);
                            if !has { new_paves.push((ei, edge.t_range[1])); }
                        }
                    }
                    for (ei, param) in new_paves {
                        self.ds.edges[ei].paves.push(Pave { vertex_idx: n_v, param });
                    }
                }
            }
        }
        // OCCT L4051: UpdateCommonBlocksWithSDVertices()
        self.ds.update_pave_blocks_with_sd_vertices();
    }

    fn update_interfs_with_sd_vertices(&mut self) {
        // Build vertex 锟?SD vertex lookup (OCCT HasShapeSD equivalent)
        let sd_for: std::collections::HashMap<usize, usize> = self.ds.shape_sd
            .sd_vertices_iter()
            .filter_map(|&(a, b)| {
                // Stored symmetrically; only process (a,b) where a < b
                // to avoid double-insert.  Both directions work since
                // all SD pairs have symmetric entries.
                if a < b { Some((a, b)) } else { None }
            })
            .collect();
        for inf in &mut self.ds.interferences {
            match inf {
                Interference::EdgeEdge { new_vertex, .. }
                | Interference::EdgeFace { new_vertex, .. } => {
                    if let Some(&sd) = sd_for.get(new_vertex) {
                        *new_vertex = sd;
                    }
                }
                Interference::VertexVertex { v1, v2, merged_vertex } => {
                    // OCCT: find SD partner for either v1 or v2
                    if let Some(&sd) = sd_for.get(v1) {
                        *merged_vertex = sd;
                    } else if let Some(&sd) = sd_for.get(v2) {
                        *merged_vertex = sd;
                    }
                }
                _ => {}
            }
        }
    }
    fn make_section_edges(&mut self) {
        // Section edges are now created per-curve inside make_blocks (OCCT form alignment).
        // This function is kept for reference but is no-op.
        return;
        // Collect section edge data per curve to avoid borrow conflicts
        struct SECurve { curve_idx: usize, sv: usize, ev: usize, curve: Curve3, geom_tol: f64, t_range: [f64; 2], pbs: Vec<PaveBlock> }
        let mut se_data: Vec<SECurve> = Vec::new();

        // 鉁?OCCT-aligned: build position鈫抳ertex map for IC endpoint remapping (ShapesSD equivalent).
        //   OCCT's PaveFiller records same-domain vertices via AddShapeSD.
        //   rcad: find the minimum-index vertex at each distinct position from ALL DS vertices.
        let bv_positions: Vec<(DVec3, usize)> = self.ds.vertices.iter().enumerate()
            .map(|(vi, v)| (v.point, vi)).collect();
        let remap_ds_v = |v: usize| -> usize {
            if v >= self.ds.vertices.len() { return v; }
            let p = self.ds.vertices[v].point;
            let tol = TOLERANCE_ABS * 1000.0;
            bv_positions.iter()
                .filter(|(bp, _)| (bp - p).length_squared() <= tol * tol)
                .map(|&(_, bv)| bv)
                .min()
                .unwrap_or(v)
        };
        // OCCT L854-875 (M6): Build BVH tree of existing PBs for IsExistingPaveBlock lookup.
        // rcad: key by (face1, face2, nV1, nV2) 鈥?same geometric edge from same face pair reuses.
        // Different face pairs produce separate edges even with same vertices (sphere pole case).
        let mut existing_edge_map: std::collections::HashMap<(usize, usize, usize, usize), usize> = std::collections::HashMap::new();

        for ci in 0..self.ds.intersection_curves.len() {
            let ic = &self.ds.intersection_curves[ci];

            // Find the two faces for IsValidBlockForFaces check (OCCT L906-918)
            let face_ids = find_face_idxs_for_curve(&self.ds, ci);
            let ff_tol = if face_ids[0] != usize::MAX && face_ids[1] != usize::MAX {
                self.ff_tol(face_ids[0], face_ids[1])
            } else { ic.geom_tol };
            // Pre-extract surface references for borrow-free comparison
            let surf0 = if face_ids[0] != usize::MAX { Some(self.ds.faces[face_ids[0]].surface.clone()) } else { None };
            let surf1 = if face_ids[1] != usize::MAX { Some(self.ds.faces[face_ids[1]].surface.clone()) } else { None };

            let mut sub_with_edge: Vec<PaveBlock> = Vec::new();

            for pbi in 0..ic.pave_blocks.len() {
                let pb = &ic.pave_blocks[pbi];

                // Clone all data before mutable access
                let mut pb_clone = pb.clone();
                let sub_pbs = if pb_clone.is_to_update() {
                    pb_clone.update(true) // flag=true: include boundary paves, matching OCCT Update() usage
                } else {
                    // OCCT-aligned: curves without ext_paves produce a single section edge
                    // spanning the entire IC range (OCCT uses Curve.StartVertex/EndVertex).
                    vec![PaveBlock::new(
                        crate::bopds::pave::NO_EDGE,
                        Pave { vertex_idx: ic.start_vertex, param: ic.t_range[0] },
                        Pave { vertex_idx: ic.end_vertex, param: ic.t_range[1] },
                    )]
                };
            for mut sub_pb in sub_pbs {
                let (nV1_raw, nV2_raw) = sub_pb.indices();
                // 鉁?OCCT-aligned: remap IC endpoint vertices to canonical boundary vertices
                //   (ShapesSD equivalent).  OCCT records SD during PaveFiller vertex creation;
                //   rcad does it here so section edges connect boundary vertices, not orphan IC vertices.
                let nV1 = remap_ds_v(nV1_raw);
                let nV2 = remap_ds_v(nV2_raw);
                if nV1 != nV1_raw || nV2 != nV2_raw {
                    sub_pb.pave1.vertex_idx = nV1;
                    sub_pb.pave2.vertex_idx = nV2;
                }
                let (aT1, aT2) = sub_pb.range();
                if (aT2 - aT1).abs() < crate::tolerance::TOLERANCE_ABS {
                    if std::env::var("RCAD_DEBUG_PB").is_ok() { eprintln!("[PB_FAIL] ci={} RANGE_TOO_SMALL", ci); }
                    continue;
                }
                // OCCT L906-918: IsValidBlockForFaces 锟?check midpoint of sub-PB against both faces
                if surf0.is_some() && surf1.is_some() {
                    let s0 = surf0.as_ref().unwrap();
                    let s1 = surf1.as_ref().unwrap();
                    let mid_t = (aT1 + aT2) * 0.5;
                    let mid_pt = ic.curve.point_at(mid_t);
                    let check_tol = ff_tol.max(TOLERANCE_ABS);
                    let mut b_flag = true;
                    for (i, &fi) in [face_ids[0], face_ids[1]].iter().enumerate() {
                        if fi == usize::MAX { continue; }
                        let pcurve = if i == 0 { ic.pcurve_on_a.as_ref() } else { ic.pcurve_on_b.as_ref() };
                        if let Some(pc) = pcurve {
                            let uv = pc.point_at(mid_t);
                            let in_on = self.context.is_point_in_on_face(self.ds, fi, uv);
                            if std::env::var("RCAD_DEBUG_ISVALID").is_ok() && (fi == 3 || fi == 0) {
                                eprintln!("[IV] ci={} fi={} uv=({:.6},{:.6}) in_on={}", ci, fi, uv.x, uv.y, in_on);
                            }
                            if !in_on {
                                // 3D fallback: check distance from midpoint to face surface.
                                let surf = if i == 0 { surf0.as_ref().unwrap() } else { surf1.as_ref().unwrap() };
                                let (_, proj_pt) = crate::extrema::closest_point_on_surface(surf, mid_pt);
                                let dist_3d = proj_pt.distance(mid_pt);
                                if dist_3d > check_tol {
                                    b_flag = false; break;
                                }
                            }
                        } else {
                            let surf = if i == 0 { surf0.as_ref().unwrap() } else { surf1.as_ref().unwrap() };
                            let (_, proj) = crate::extrema::closest_point_on_surface(surf, mid_pt);
                            if proj.distance(mid_pt) > check_tol {
                                b_flag = false; break;
                            }
                        }
                    }
                    if !b_flag { continue; }
                } // end IsValidBlockForFaces
                // OCCT L936-947: FindValidRange check 锟?skip micro-edges where vertex tolerance
                //   spheres cover the entire parameter range.
                if nV1 < self.ds.vertices.len() && nV2 < self.ds.vertices.len() {
                    let v1_pt = self.ds.vertices[nV1].point;
                    let v2_pt = self.ds.vertices[nV2].point;
                    let v1_tol = ff_tol.max(self.ds.vertices[nV1].geom_tol);
                    let v2_tol = ff_tol.max(self.ds.vertices[nV2].geom_tol);
                    if find_valid_range(&ic.curve, aT1, aT2, v1_pt, v1_tol, v2_pt, v2_tol).is_none() {
                        if std::env::var("RCAD_DEBUG_PB").is_ok() {
                            eprintln!("[PB] ci={} BLOCKED FindValidRange nV=({},{}) v1_tol={:.12} v2_tol={:.12} v1_pt=({:.4},{:.4},{:.4}) v2_pt=({:.4},{:.4},{:.4})",
                                ci, nV1, nV2, v1_tol, v2_tol,
                                v1_pt.x, v1_pt.y, v1_pt.z, v2_pt.x, v2_pt.y, v2_pt.z);
                        }
                        continue;
                    }
                }
                // OCCT L920-963 (M7e): IsExistingPaveBlock 鈥?check if this sub-PB already has
                //   a DSEdge from another curve (via BVH tree / existing_edge_map).
                // key includes face pair so different FF pairs create separate edges.
                let (v1, v2) = if nV1 < nV2 { (nV1, nV2) } else { (nV2, nV1) };
                let f1 = face_ids[0].min(face_ids[1]);
                let f2 = face_ids[0].max(face_ids[1]);
                let edge_key = (f1, f2, v1, v2);
                if let Some(&existing_ei) = existing_edge_map.get(&edge_key) {
                    // OCCT L924-928: UpdateEdgeTolerance + UpdateSavedTolerance for reused edge
                    sub_pb.new_edge = Some(existing_ei);
                    sub_with_edge.push(sub_pb);
                    if std::env::var("RCAD_DEBUG_PB").is_ok() && face_ids[0] == 0 { eprintln!("[PB_PASS] ci={} REUSE edge={}", ci, existing_ei); }
                    continue;
                }
                // Create new DSEdge for this sub-PB
                let new_ei = self.ds.edges.len();
                // OCCT-aligned: propagate pcurves from IC to section DSEdge face_reps.
                let mut sec_face_reps = Vec::new();
                if let Some(ref pca) = ic.pcurve_on_a {
                    sec_face_reps.push(DSRepOnFace {
                        face_idx: face_ids[0],
                        pcurve: pca.clone(),
                        pcurve2: None,
                        pcurve_range: [aT1, aT2],
                        start_param: aT1, end_param: aT2,
                    });
                }
                if let Some(ref pcb) = ic.pcurve_on_b {
                    sec_face_reps.push(DSRepOnFace {
                        face_idx: face_ids[1],
                        pcurve: pcb.clone(),
                        pcurve2: None,
                        pcurve_range: [aT1, aT2],
                        start_param: aT1, end_param: aT2,
                    });
                }
                self.ds.edges.push(DSEdge {
                    start_vertex: nV1, end_vertex: nV2,
                    curve: ic.curve.clone(),
                    t_range: [aT1, aT2],
                    origin: ShapeOrigin::ShapeA,
                    geom_tol: ic.geom_tol,
                    paves: Vec::new(),
                    pave_blocks: vec![sub_pb.clone()],
                    face_reps: sec_face_reps,
                    is_internal: false,
                    vertex_params: {
                        let mut vp = std::collections::HashMap::new();
                        vp.insert(nV1, aT1);
                        vp.insert(nV2, aT2);
                        vp
                    },
                });
                // Set new_edge in the PB stored inside the edge AND in sub_with_edge
                if let Some(epb) = self.ds.edges.last_mut().and_then(|e| e.pave_blocks.first_mut()) {
                    epb.new_edge = Some(new_ei);
                }
                sub_pb.new_edge = Some(new_ei);
                self.ds.section_edge_refs[ci].push(new_ei);
                existing_edge_map.insert(edge_key, new_ei);
                sub_with_edge.push(sub_pb);
            }
            } // end for pbi in 0..ic.pave_blocks.len()

            if !sub_with_edge.is_empty() {
                se_data.push(SECurve {
                    curve_idx: ci,
                    sv: ic.start_vertex, ev: ic.end_vertex,
                    curve: ic.curve.clone(), geom_tol: ic.geom_tol,
                    t_range: ic.t_range,
                    pbs: sub_with_edge,
                });
            }
        }

        // Register section edge PBs into global pool and pave_blocks_sc
        // OCCT-aligned: each section edge belongs only to the TWO faces of its FF pair.
        for se in &se_data {
            // Find the two faces referencing this curve
            let face_ids = find_face_idxs_for_curve(&self.ds, se.curve_idx);
            for pb in &se.pbs {
                if pb.new_edge.is_some() {
                    let g_pb_idx = self.ds.allocate_pave_block(pb.clone());
                    for &fi in &face_ids {
                        if fi != usize::MAX {
                            self.ds.faces[fi].face_info.pave_blocks_sc.insert(g_pb_idx);
                        }
                    }
                }
            }
        }
    }

    fn remove_micro_edges(&mut self) {
        let mut micro_edges: Vec<usize> = Vec::new();

        // OCCT L4398-4433: iterate edges in PaveBlocks pool
        for ei in 0..self.ds.edges.len() {
            // OCCT L4401-4403: skip edges with <2 PaveBlocks (no splits)
            if self.ds.edges[ei].pave_blocks.len() < 2 {
                continue;
            }

            // OCCT L4407-4410: skip degenerated edges (HasFlag)
            if self.ds.is_edge_degenerated(ei) {
                continue;
            }

            // OCCT L4412-4432: iterate PaveBlocks, find nV1==nV2
            for pb in &self.ds.edges[ei].pave_blocks {
                let nv1 = pb.pave1.vertex_idx;
                let nv2 = pb.pave2.vertex_idx;
                if nv1 == nv2 {
                    // OCCT L4425-4426: FillShrunkData + HasShrunkData check
                    // 锟?rcad: shrunk data not available 锟?skip valid-range check
                    micro_edges.push(ei);
                    break;
                }
            }
        }

        // OCCT L4434: RemovePaveBlocks(aMicroEdges)
        // rcad: clear pave_blocks, restore single-span representation
        for &ei in &micro_edges {
            let sv = self.ds.edges[ei].start_vertex;
            let ev = self.ds.edges[ei].end_vertex;
            let t0 = self.ds.edges[ei].t_range[0];
            let t1 = self.ds.edges[ei].t_range[1];
            self.ds.edges[ei].pave_blocks = vec![PaveBlock::new(
                ei,
                Pave { vertex_idx: sv, param: t0 },
                Pave { vertex_idx: ev, param: t1 },
            )];
        }
    }

}

// Second empty impl block (kept for OCCT alignment)
impl<'a> PaveFiller<'a> {
}

#[cfg(test)]
mod tests;

