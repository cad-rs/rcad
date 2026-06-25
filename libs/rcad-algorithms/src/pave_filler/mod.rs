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

/// �?OCCT-aligned: IntPatch_Intersection surface category (L1264-1294).
///   GeomGeom  = ts1==ts2==1 �?ImpImpIntersection (analytic-analytic)
///   GeomParam = ts1!=ts2     �?ImpPrmIntersection (analytic-parametric)
///   ParamParam = ts1==ts2==0 �?PrmPrmIntersection (parametric-parametric)
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
/// Below this threshold, brute-force O(n�? is faster due to BVH build overhead.
const BVH_THRESHOLD: usize = 20;

/// �?OCCT-aligned: BOPAlgo_PaveFiller �?six intersection passes
///   (PaveFiller.hxx L106-107, PaveFiller.cxx L234-355).
mod glue;
mod intersection;
pub struct PaveFiller<'a> {
    pub ds: &'a mut DS,
    bvh_a: Option<&'a Bvh>,
    bvh_b: Option<&'a Bvh>,
    use_glue: bool,
    glue_tolerance: f64,
    /// �?OCCT-aligned: BOPAlgo_Options::SetFuzzyValue
    fuzzy_tolerance: f64,
    /// �?OCCT-aligned: PaveFiller_6.cxx L393-479 seam edge shift tolerance
    seam_shift_tol: f64,
    /// �?OCCT-aligned: BOPAlgo_Algo::myRunParallel
    run_parallel: bool,
    /// �?OCCT-aligned: BOPAlgo_PaveFiller::myNonDestructive
    non_destructive: bool,
    /// �?OCCT-aligned: BOPAlgo_Algo::myUseOBB
    use_obb: bool,
    /// �?OCCT-aligned: IntTools_Context (PaveFiller::Init L203)
    context: IntToolsContext,
}

/// �?OCCT-aligned:Propagate IC vertices to all faces sharing boundary edges
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
    pub fn new(ds: &'a mut DS) -> Self {
        let n_faces = ds.faces.len();
        let context = IntToolsContext::new(n_faces, TOLERANCE_ABS * 100.0);
        Self {
            ds,
            bvh_a: None,
            bvh_b: None,
            use_glue: false,
            glue_tolerance: TOLERANCE_ABS,
            fuzzy_tolerance: 0.0,
            seam_shift_tol: 0.0,
            // �?OCCT-aligned: RunParallel (default false)
            run_parallel: false,
            // �?OCCT-aligned: NonDestructive (default false)
            non_destructive: false,
            // �?OCCT-aligned: UseOBB (default false)
            use_obb: false,
            // �?OCCT-aligned: IntTools_Context with FClass2d cache
            // OCCT PaveFiller.cxx L203: myContext = new IntTools_Context
            context,
        }
    }

    pub fn with_bvh(ds: &'a mut DS, bvh_a: &'a Bvh, bvh_b: &'a Bvh) -> Self {
        let total_faces = ds.faces.len();
        let use_bvh = total_faces >= BVH_THRESHOLD;
        let context = IntToolsContext::new(total_faces, TOLERANCE_ABS * 100.0);
        Self {
            ds,
            bvh_a: if use_bvh { Some(bvh_a) } else { None },
            bvh_b: if use_bvh { Some(bvh_b) } else { None },
            use_glue: false,
            glue_tolerance: TOLERANCE_ABS,
            fuzzy_tolerance: 0.0,
            seam_shift_tol: 0.0,
            // �?OCCT-aligned: RunParallel (default false)
            run_parallel: false,
            // �?OCCT-aligned: NonDestructive (default false)
            non_destructive: false,
            // �?OCCT-aligned: UseOBB (default false)
            use_obb: false,
            // �?OCCT-aligned: IntTools_Context with FClass2d cache
            context,
        }
    }

    pub fn configure_glue(&mut self, enable: bool, tolerance: f64) {
        self.use_glue = enable;
        self.glue_tolerance = tolerance.max(TOLERANCE_ABS);
    }

    pub fn configure_glue_adaptive(&mut self, enable: bool, base_tolerance: f64, adaptive: bool) -> f64 {
        if !enable {
            self.use_glue = false;
            return TOLERANCE_ABS;
        }

        self.use_glue = true;

        if !adaptive {
            self.glue_tolerance = base_tolerance.max(TOLERANCE_ABS);
            return self.glue_tolerance;
        }

        // Compute adaptive tolerance based on geometry
        let adaptive_tol = self.compute_adaptive_glue_tolerance(base_tolerance);
        self.glue_tolerance = adaptive_tol;
        adaptive_tol
    }

    pub fn configure_fuzzy(&mut self, fuzzy: f64) {
        self.fuzzy_tolerance = fuzzy.max(0.0);
    }

    pub fn set_run_parallel(&mut self, parallel: bool) {
        self.run_parallel = parallel;
    }

    pub fn set_non_destructive(&mut self, nd: bool) {
        self.non_destructive = nd;
    }

    pub fn set_non_destructive_auto(&mut self) {
        // OCCT: checks if any argument has a locked sub-shape.
        // rcad does not support locked shapes.
        self.non_destructive = false;
    }

    pub fn set_use_obb(&mut self, use_obb: bool) {
        self.use_obb = use_obb;
    }

    pub fn effective_tolerance(&self, base: f64) -> f64 {
        base.max(self.fuzzy_tolerance)
    }

    fn compute_adaptive_glue_tolerance(&self, base_tolerance: f64) -> f64 {
        let mut min_feature_size = f64::INFINITY;
        let mut min_edge_length = f64::INFINITY;
        let mut min_face_area = f64::INFINITY;

        // Analyze edge lengths
        for edge in &self.ds.edges {
            let p1 = edge.curve.point_at(edge.t_range[0]);
            let p2 = edge.curve.point_at(edge.t_range[1]);
            let length = (p2 - p1).length();
            if length > TOLERANCE_LINEAR_ULTRA_STRICT {
                min_edge_length = min_edge_length.min(length);
            }
        }

        // Analyze face areas (approximate from bounding box)
        for face in &self.ds.faces {
            let pts = self.ds.face_boundary_points(
                self.ds.faces.iter().position(|f| std::ptr::eq(f, face)).unwrap_or(0)
            );
            if pts.len() >= 3 {
                // Compute bounding box diagonal as area proxy
                let mut min_pt = pts[0];
                let mut max_pt = pts[0];
                for p in &pts[1..] {
                    min_pt = min_pt.min(*p);
                    max_pt = max_pt.max(*p);
                }
                let diag = (max_pt - min_pt).length();
                if diag > TOLERANCE_LINEAR_ULTRA_STRICT {
                    min_face_area = min_face_area.min(diag * diag);
                }
            }
        }

        // Use minimum feature size to bound tolerance
        if min_edge_length.is_finite() {
            min_feature_size = min_feature_size.min(min_edge_length);
        }
        if min_face_area.is_finite() {
            min_feature_size = min_feature_size.min(min_face_area.sqrt());
        }

        // Compute adaptive tolerance
        let adaptive_tol = if min_feature_size.is_finite() && min_feature_size > 0.0 {
            // Use a fraction of minimum feature size, but at least base tolerance
            let feature_based = min_feature_size * 0.01;
            base_tolerance.max(feature_based).min(min_feature_size * 0.1)
        } else {
            base_tolerance
        };

        adaptive_tol.max(TOLERANCE_ABS)
    }

    // (Extreme geometry detection removed �?rcad invention.
    //  OCCT Prepare builds pcurves on planar faces; rcad DS::build_face_reps subsumes this.)

    #[inline]
    fn tol(&self) -> f64 {
        self.ds.fuzzy_tol
    }

    #[inline]
    fn vv_pair_tol(&self, vi: usize, vj: usize) -> f64 {
        self.ds.vertices[vi].geom_tol
            + self.ds.vertices[vj].geom_tol
            + self.tol()
    }

    #[inline]
    fn ve_tol(&self, vi: usize, ei: usize) -> f64 {
        self.ds.vertices[vi].geom_tol
            + self.ds.edges[ei].geom_tol
            + self.tol()
    }

    #[inline]
    fn ee_tol(&self, e1: usize, e2: usize) -> f64 {
        self.ds.edges[e1].geom_tol
            + self.ds.edges[e2].geom_tol
            + self.tol()
    }

    #[inline]
    fn vf_tol(&self, vi: usize, fi: usize) -> f64 {
        self.ds.vertices[vi].geom_tol
            + self.ds.faces[fi].geom_tol
            + self.tol()
    }

    #[inline]
    fn ef_tol(&self, ei: usize, fi: usize) -> f64 {
        self.ds.edges[ei].geom_tol
            + self.ds.faces[fi].geom_tol
            + self.tol()
    }

    #[inline]
    fn ff_tol(&self, f1: usize, f2: usize) -> f64 {
        self.tol()
            .max(self.ds.faces[f1].geom_tol)
            .max(self.ds.faces[f2].geom_tol)
            .max(self.seam_shift_tol)
    }

    fn find_face_face_curve_indices(&self, f1: usize, f2: usize) -> Option<Vec<usize>> {
        for inf in &self.ds.interferences {
            if let Interference::FaceFace { f1: a, f2: b, curves, .. } = inf {
                if *a == f1 && *b == f2 {
                    return Some(curves.clone());
                }
            }
        }
        None
    }

    fn sampled_face_boundary_points(&self, face_idx: usize, samples_per_edge: usize) -> Vec<DVec3> {
        let mut pts = Vec::new();
        for &ei in &self.ds.faces[face_idx].boundary_edges {
            if let Some(edge) = self.ds.edges.get(ei) {
                let [t0, t1] = edge.t_range;
                let n = samples_per_edge.max(1);
                for k in 0..=n {
                    let t = t0 + (t1 - t0) * k as f64 / n as f64;
                    let p = edge.curve.point_at(t);
                    if p.is_finite() {
                        pts.push(p);
                    }
                }
            }
        }
        if pts.is_empty() {
            self.ds.face_boundary_points(face_idx)
        } else {
            pts
        }
    }

    fn closest_point_on_boundary_samples(&self, point: DVec3, samples: &[DVec3]) -> DVec3 {
        samples
            .iter()
            .copied()
            .min_by(|a, b| {
                let da = (*a - point).length_squared();
                let db = (*b - point).length_squared();
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(point)
    }

    fn snap_polyline_endpoints_to_face_boundaries(
        &self,
        chain: &mut Vec<DVec3>,
        f1: usize,
        f2: usize,
    ) {
        if chain.len() < 2 {
            return;
        }

        let boundary_a = self.sampled_face_boundary_points(f1, 12);
        let boundary_b = self.sampled_face_boundary_points(f2, 12);
        if boundary_a.is_empty() || boundary_b.is_empty() {
            return;
        }

        let snap_start_a = self.closest_point_on_boundary_samples(chain[0], &boundary_a);
        let snap_start_b = self.closest_point_on_boundary_samples(chain[0], &boundary_b);
        let snap_end_a = self.closest_point_on_boundary_samples(chain[chain.len() - 1], &boundary_a);
        let snap_end_b = self.closest_point_on_boundary_samples(chain[chain.len() - 1], &boundary_b);

        let choose_better = |orig: DVec3, p1: DVec3, p2: DVec3| {
            let d1 = (p1 - orig).length_squared();
            let d2 = (p2 - orig).length_squared();
            if d1 <= d2 { p1 } else { p2 }
        };

        let start = choose_better(chain[0], snap_start_a, snap_start_b);
        let end = choose_better(chain[chain.len() - 1], snap_end_a, snap_end_b);

        // Only snap if it is a local correction rather than a gross relocation.
        let local_scale = chain
            .windows(2)
            .map(|w| (w[1] - w[0]).length())
            .filter(|d| d.is_finite() && *d > 0.0)
            .fold(f64::INFINITY, f64::min)
            .min(1.0);
        let snap_tol = (local_scale * 4.0)
            .max(TOLERANCE_RETRY_LADDER_COARSE)
            .max(self.ff_tol(f1, f2));

        if (start - chain[0]).length() <= snap_tol {
            chain[0] = start;
        }
        if (end - chain[chain.len() - 1]).length() <= snap_tol {
            let last = chain.len() - 1;
            chain[last] = end;
        }
    }

    pub fn perform(&mut self) {
        // �?OCCT-aligned: no extreme-geometry pre-analysis �?OCCT Prepare builds
        //    pcurves on planar faces; rcad does this in DS::build_face_reps.
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

        // OCCT L145: BOPDS_Iterator builds BVH trees for all shape types.
        //   rcad: build DS-level BVHs for VE/EE/VF/EF pair culling.
        let bvh_verts_a = self.build_ds_bvh(true, false);
        let bvh_verts_b = self.build_ds_bvh(false, false);
        let bvh_edges_a = self.build_ds_bvh(true, true);
        let bvh_edges_b = self.build_ds_bvh(false, true);
        // Face BVHs are built from the same index ranges as the face BVH trees.
        let bvh_faces_a = self.build_ds_bvh_face(true);
        let bvh_faces_b = self.build_ds_bvh_face(false);

        if !skip_ve {
            self.perform_ve_bvh(&bvh_verts_a, &bvh_edges_b);
        }
        // �?OCCT-aligned: UpdatePaveBlocksWithSDVertices (PerformInternal L266)
        self.ds.update_pave_blocks_with_sd_vertices();

        let ee_survivors: Vec<usize> = if !skip_ee {
            self.perform_ee_bvh(&bvh_edges_a, &bvh_edges_b);
            // �?OCCT-aligned: TreatNewVertices �?merge new vertices created by EE intersection.
            //    OCCT PaveFiller_5.cxx L570: PerformNewVertices(aMVCPB, ..., false)
            let survivors = self.treat_new_vertices();
            // �?OCCT-aligned: UpdatePaveBlocksWithSDVertices (PerformInternal L273)
            self.ds.update_pave_blocks_with_sd_vertices();
            survivors
        } else { vec![] };

        if !skip_vf {
            self.perform_vf_bvh(&bvh_verts_a, &bvh_faces_b);
            self.perform_vf_bvh(&bvh_verts_b, &bvh_faces_a);
        }
        // �?OCCT-aligned: UpdatePaveBlocksWithSDVertices (PerformInternal L280)
        self.ds.update_pave_blocks_with_sd_vertices();

        if !skip_ef {
            self.perform_ef_bvh(&bvh_edges_a, &bvh_faces_b);
            self.perform_ef_bvh(&bvh_edges_b, &bvh_faces_a);
            // �?OCCT-aligned: TreatNewVertices �?merge new vertices created by EF intersection.
            //    OCCT PaveFiller_5.cxx L570: PerformNewVertices(aMVCPB, ..., false)
            let ef_survivors = self.treat_new_vertices();

            // �?OCCT-aligned: RepeatIntersection (PaveFiller.cxx L296-299, L359-420).
            //    After EF, before FF, re-run VV/VE/VF for vertices with increased tolerance.
            //    OCCT reads from myIncreasedSS (populated by TreatNewVertices).
            //    rcad: ds.increased_ss is populated by treat_new_vertices above.
            self.ds.update_pave_blocks_with_sd_vertices();
            self.update_interfs_with_sd_vertices();
            self.repeat_intersection();
        }

        // �?OCCT-aligned: ForceInterfEE (PaveFiller_3.cxx L978-1276)
        //    OCCT L302: ForceInterfEE �?after RepeatIntersection, force intersection
        //    of edge pairs sharing a vertex with increased tolerance, detecting
        //    collinear/coincident edges (common block).
        //    �?rcad: simplified, only checks line-line edge pairs sharing a pave vertex.
        if !skip_ee {
            self.force_interf_ee();
        }

        // �?OCCT-aligned: ForceInterfEF (PaveFiller_5.cxx L764-1099+)
        //    OCCT L309: ForceInterfEF �?after ForceInterfEE, force intersection of
        //    edges whose both endpoints are on a face with increased tolerance.
        //    �?rcad: simplified, only checks edge-face pairs where both endpoints are on the face.
        if !skip_ef {
            self.force_interf_ef();
        }

        if !skip_ff {
            self.perform_ff();

            // �?OCCT-aligned: MakeSDVerticesFF (PaveFiller_6.cxx L1113)
            //    After FF, create shared SD vertices for same-domain (coplanar) face
            //    overlap boundaries so that overlap polygon vertices are shared between
            //    both faces and registered in face_info.vertices_in.
            self.make_sd_vertices_ff();
        }

        // �?OCCT-aligned: PostTreatFF (PaveFiller_6.cxx)
        //    Reconcile FF interference data with face info. Iterates all FF interferences
        //    and updates face_info.curves_sc + vertices_in from curve endpoints.
        self.post_treat_ff();

        // �?OCCT-aligned: UpdateBlocksWithSharedVertices (PerformInternal L318)
        self.update_blocks_with_shared_vertices();

        // �?OCCT-aligned: RefineFaceInfoIn �?before MakeSplitEdges, remove
        //    On-overlapping In pave blocks (PerformInternal L320, BOPDS_DS::RefineFaceInfoIn).
        for fi in 0..self.ds.faces.len() {
            self.ds.refine_face_info_in(fi);
        }

        // �?OCCT-aligned: MakeSplitEdges �?create split edges from PaveBlocks (PerformInternal L322).
        //   rcad: build_split_edges() = MakeSplitEdges under OCCT name.
        self.build_split_edges();

        // �?OCCT-aligned: UpdatePaveBlocksWithSDVertices (PerformInternal L328)
        self.ds.update_pave_blocks_with_sd_vertices();

        // �?OCCT-aligned: MakeBlocks �?inject EF/EE vertices onto FF curves (PerformInternal L330)
        self.make_blocks();

        // �?OCCT-aligned: CheckSelfInterference (PerformInternal L336, BOPAlgo_PaveFiller_11.cxx L28-221)
        //    OCCT uses AddWarning �?non-fatal, the operation continues.
        if let Err(msg) = self.check_self_interference() {
            eprintln!("[PAVEFILLER] {}", msg);
        }

        // �?OCCT-aligned: UpdateInterfsWithSDVertices (PerformInternal L338)
        self.update_interfs_with_sd_vertices();

        // �?OCCT-aligned: ReleasePaveBlocks �?free unused pave block memory (PerformInternal L339).
        self.ds.pave_blocks.clear();

        // �?OCCT-aligned: RefineFaceInfoOn �?after ReleasePaveBlocks, remove
        //    zero-length On pave blocks (PerformInternal L340, BOPDS_DS::RefineFaceInfoOn).
        for fi in 0..self.ds.faces.len() {
            self.ds.refine_face_info_on(fi);
        }

        // �?OCCT-aligned: RemoveMicroEdges �?after MakeBlocks, before MakePCurves
        //    (PerformInternal L342, PaveFiller_6.cxx L4229-4270).
        self.remove_micro_edges();

        // �?OCCT-aligned: MakePCurves �?after RemoveMicroEdges (PerformInternal L344)
        self.make_pcurves();

        // �?OCCT-aligned: ProcessDE �?after MakePCurves (PerformInternal L350)
        self.process_de();
    }

    // ===== BVH-based pair enumeration (OCCT BOPDS_Iterator) =====













    // 鈹€鈹€鈹€ Pass 1: Vertex-Vertex 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€


    // 鈹€鈹€鈹€ Pass 2: Vertex-Edge 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€



    // 鈹€鈹€鈹€ Pass 3: Edge-Edge 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€




    // 鈹€鈹€鈹€ Pass 4: Vertex-Face 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€





// 鈹€鈹€鈹€ Pass 5: Edge-Face 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€






    // 鈹€鈹€鈹€ Pass 6: Face-Face 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€


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
            //   rcad: not needed �?FaceInfo data is already populated.

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
        // Build vertex �?SD vertex lookup (OCCT HasShapeSD equivalent)
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
    fn make_section_edges_from_curve_pbs(&mut self) {
        let n_edges_before = self.ds.edges.len();
        // Collect section edge data per curve to avoid borrow conflicts
        struct SECurve { curve_idx: usize, sv: usize, ev: usize, curve: Curve3, geom_tol: f64, t_range: [f64; 2], pbs: Vec<PaveBlock> }
        let mut se_data: Vec<SECurve> = Vec::new();

        for ci in 0..self.ds.intersection_curves.len() {
            let ic = &self.ds.intersection_curves[ci];
            // Must have exactly one PB (init_pave_block1 was called)
            if ic.pave_blocks.len() != 1 { continue; }
            let pb = &ic.pave_blocks[0];

            // Clone all data before mutable access
            let mut pb_clone = pb.clone();
            let sub_pbs = if pb_clone.is_to_update() {
                pb_clone.update(false) // flag=false: ext_paves only, no boundary paves
            } else {
                // OCCT-aligned: curves without ext_paves produce a single section edge
                // spanning the entire curve (no split points).
                vec![PaveBlock::new(
                    crate::bopds::pave::NO_EDGE,
                    Pave { vertex_idx: ic.start_vertex, param: ic.t_range[0] },
                    Pave { vertex_idx: ic.end_vertex, param: ic.t_range[1] },
                )]
            };

            // Find the two faces for IsValidBlockForFaces check (OCCT L906-918)
            let face_ids = find_face_idxs_for_curve(&self.ds, ci);
            let ff_tol = if face_ids[0] != usize::MAX && face_ids[1] != usize::MAX {
                self.ff_tol(face_ids[0], face_ids[1])
            } else { ic.geom_tol };
            // Pre-extract surface references for borrow-free comparison
            let surf0 = if face_ids[0] != usize::MAX { Some(self.ds.faces[face_ids[0]].surface.clone()) } else { None };
            let surf1 = if face_ids[1] != usize::MAX { Some(self.ds.faces[face_ids[1]].surface.clone()) } else { None };

            let mut sub_with_edge: Vec<PaveBlock> = Vec::new();
            for mut sub_pb in sub_pbs {
                let (nV1, nV2) = sub_pb.indices();
                let (aT1, aT2) = sub_pb.range();
                if (aT2 - aT1).abs() < crate::tolerance::TOLERANCE_ABS {
                    continue;
                }
                // OCCT L906-918: IsValidBlockForFaces �?check midpoint of sub-PB against both faces
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
                            if !self.context.is_point_in_on_face(self.ds, fi, uv) {
                                b_flag = false; break;
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
                // OCCT L936-947: FindValidRange check �?skip micro-edges where vertex tolerance
                //   spheres cover the entire parameter range.
                if nV1 < self.ds.vertices.len() && nV2 < self.ds.vertices.len() {
                    let v1_pt = self.ds.vertices[nV1].point;
                    let v2_pt = self.ds.vertices[nV2].point;
                    let v1_tol = ff_tol.max(self.ds.vertices[nV1].geom_tol);
                    let v2_tol = ff_tol.max(self.ds.vertices[nV2].geom_tol);
                    if find_valid_range(&ic.curve, aT1, aT2, v1_pt, v1_tol, v2_pt, v2_tol).is_none() {
                        continue;
                    }
                }
                // Create new DSEdge for this sub-PB
                let new_ei = self.ds.edges.len();
                self.ds.edges.push(DSEdge {
                    start_vertex: nV1, end_vertex: nV2,
                    curve: ic.curve.clone(),
                    t_range: [aT1, aT2],
                    origin: ShapeOrigin::ShapeA,
                    geom_tol: ic.geom_tol,
                    paves: Vec::new(),
                    pave_blocks: Vec::new(),
                    face_reps: Vec::new(),
                    is_internal: false,
                });
                sub_pb.new_edge = Some(new_ei);
                sub_with_edge.push(sub_pb);
            }

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
                    // �?rcad: shrunk data not available �?skip valid-range check
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

    fn force_interf_ee(&mut self) {
        // OCCT L1008-1023: initialize PBs for interfered vertices
        // rcad: build vertex �?edge mapping from edge.paves
        // OCCT L1047-1051: skip degenerated edges (HasFlag)
        // OCCT L1041-1045: HasReference �?non-empty pave_blocks
        // OCCT L1047-1051: HasFlag �?skip degenerated edges
        let mut vert_edges: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
        for (ei, edge) in self.ds.edges.iter().enumerate() {
            if edge.paves.is_empty() { continue; }
            for pave in &edge.paves {
                vert_edges.entry(pave.vertex_idx).or_default().push(ei);
            }
        }

        // OCCT L1060-1177: intersect PaveBlock pairs sharing a vertex
        for (&vi, edges) in &vert_edges {
            if edges.len() < 2 { continue; }
            for i in 0..edges.len() {
                let e1 = edges[i];
                for j in (i + 1)..edges.len() {
                    let e2 = edges[j];
                    // OCCT L1113-1121: check edges are from different operands
                    let o1 = self.ds.edges[e1].origin;
                    let o2 = self.ds.edges[e2].origin;
                    if o1 == o2 { continue; }

                    // OCCT L1123-1128: skip edge pairs already forming a CommonBlock
                    if self.has_ee_interf(e1, e2) { continue; }

                    // OCCT L1077-1083: aTolAdd = 2 * max(tol(V1), tol(V2))
                    let v_tol = self.ds.vertices[vi].geom_tol;
                    let tol_add = 2.0 * v_tol;

                    // Only check line-line edges (OCCT L1138-1157: angle check)
                    let c1 = &self.ds.edges[e1].curve;
                    let c2 = &self.ds.edges[e2].curve;
                    match (c1, c2) {
                        (Curve3::Line(l1), Curve3::Line(l2)) => {
                            // OCCT L1097-1102: midpoint direction vector, check angle
                            let d1 = l1.direction;
                            let d2 = l2.direction;
                            let cos_angle = d1.dot(d2).abs();
                            // OCCT L1155: angle > 25�?�?cos < 0.9063 �?skip addTol
                            let fuzzy = if cos_angle >= 0.9063 {
                                self.ds.fuzzy_tol + tol_add
                            } else {
                                self.ds.fuzzy_tol
                            };

                            // Check collinearity (OCCT EdgeEdge intersection)
                            // intersect_line_line returns Option<(f64,f64,DVec3)>
                            if let Some((t1, t2, pt)) = intersect_line_line(
                                l1, self.ds.edges[e1].t_range,
                                l2, self.ds.edges[e2].t_range, fuzzy)
                            {
                                self.ds.interferences.push(Interference::EdgeEdge {
                                    e1, e2, point: pt, param1: t1, param2: t2, new_vertex: vi,
                                });
                            }
                        }
                        (Curve3::Circle(circ), Curve3::Circle(_)) => {
                            // �?circle-circle coincidence detection simplified: uses normal tolerance
                            // intersect_circle_circle returns Vec<DVec3>
                            let fuzzy = self.ds.fuzzy_tol + tol_add;
                            let cp_hits = intersect_circle_circle(circ, circ, fuzzy);
                            if let Some(&pt) = cp_hits.first() {
                                self.ds.interferences.push(Interference::EdgeEdge {
                                    e1, e2, point: pt, param1: 0.0, param2: 0.0, new_vertex: vi,
                                });
                            }
                        }
                        _ => {
                            // OCCT L1138-1157: IntTools_EdgeEdge numerical intersection
                            // with recursive adaptive subdivision (OCCT IntTools_CurveRange).
                            let tr1 = self.ds.edges[e1].t_range;
                            let tr2 = self.ds.edges[e2].t_range;
                            let mid_t1 = (tr1[0] + tr1[1]) * 0.5;
                            let mid_t2 = (tr2[0] + tr2[1]) * 0.5;
                            let tgt1 = c1.tangent_at(mid_t1);
                            let tgt2 = c2.tangent_at(mid_t2);
                            let cos_angle = if tgt1.length_squared() > 1e-30 && tgt2.length_squared() > 1e-30 {
                                tgt1.normalize().dot(tgt2.normalize()).abs()
                            } else { 0.0 };
                            let fuzzy = if cos_angle >= 0.9063 {
                                self.ds.fuzzy_tol + tol_add
                            } else {
                                self.ds.fuzzy_tol
                            };
                            // OCCT IntTools_EdgeEdge: coarse �?adaptive �?Newton
                            // (1) Coarse 21�?1 grid �?find best (t1,t2)
                            // (2) Recursive subdivision around best: 2�?denser per level
                            // (3) Converge when distance < fuzzy OR subrange < 1e-6
                            let mut best_t1 = mid_t1;
                            let mut best_t2 = mid_t2;
                            let mut best_d = f64::MAX;
                            // OCCT N=20 �?21 samples per curve
                            for si in 0..21 {
                                let t1 = tr1[0] + (tr1[1] - tr1[0]) * (si as f64 / 20.0);
                                let p1 = c1.point_at(t1);
                                for sj in 0..21 {
                                    let t2 = tr2[0] + (tr2[1] - tr2[0]) * (sj as f64 / 20.0);
                                    let d = p1.distance(c2.point_at(t2));
                                    if d < best_d { best_d = d; best_t1 = t1; best_t2 = t2; }
                                }
                            }
                            // (2) Adaptive refinement: subdivide around min point
                            let mut r1_lo = (best_t1 - (tr1[1] - tr1[0]) / 20.0).max(tr1[0]);
                            let mut r1_hi = (best_t1 + (tr1[1] - tr1[0]) / 20.0).min(tr1[1]);
                            let mut r2_lo = (best_t2 - (tr2[1] - tr2[0]) / 20.0).max(tr2[0]);
                            let mut r2_hi = (best_t2 + (tr2[1] - tr2[0]) / 20.0).min(tr2[1]);
                            for _ in 0..4 {
                                let mid1 = (r1_lo + r1_hi) * 0.5;
                                let mid2 = (r2_lo + r2_hi) * 0.5;
                                let test_t1 = [r1_lo, mid1, r1_hi];
                                let test_t2 = [r2_lo, mid2, r2_hi];
                                for &t1 in &test_t1 {
                                    let pt1 = c1.point_at(t1);
                                    for &t2 in &test_t2 {
                                        let d = pt1.distance(c2.point_at(t2));
                                        if d < best_d { best_d = d; best_t1 = t1; best_t2 = t2; }
                                    }
                                }
                                let span = (r1_hi - r1_lo) * 0.5;
                                r1_lo = (best_t1 - span).max(tr1[0]);
                                r1_hi = (best_t1 + span).min(tr1[1]);
                                r2_lo = (best_t2 - span).max(tr2[0]);
                                r2_hi = (best_t2 + span).min(tr2[1]);
                            }
                            // (3) OCCT IntTools_CurveRange L230-260: Newton-Raphson iteration
                            // Minimize F(t1,t2) = ||C1(t1)-C2(t2)||�?using gradient+Hessian.
                            let mut nr_t1 = best_t1;
                            let mut nr_t2 = best_t2;
                            for _ in 0..8 {
                                let p1 = c1.point_at(nr_t1);
                                let p2 = c2.point_at(nr_t2);
                                let diff = p1 - p2;
                                if diff.length_squared() < 1e-30 { break; }
                                let t1 = c1.tangent_at(nr_t1);
                                let t2 = c2.tangent_at(nr_t2);
                                if t1.length_squared() < 1e-30 || t2.length_squared() < 1e-30 { break; }
                                let d1 = t1.normalize();
                                let d2 = t2.normalize();
                                // Hessian H and gradient 鈭嘑 of F(t1,t2) = ||C1-C2||�?
                                let h00 = 2.0;  // H = 2*M, M = [[d1路d1, -d1路d2], [-d2路d1, d2路d2]]
                                let h01 = -2.0 * d1.dot(d2);
                                let h10 = h01;  // symmetric
                                let h11 = 2.0;
                                // OCCT IntTools_CurveRange: R[0]=-(C1-C2)路C1', R[1]=(C1-C2)路C2'
                                // H = 2*M where M = [[1, -cos], [-cos, 1]], RHS = [2*R[0], 2*R[1]]
                                let g0 = 2.0 * diff.dot(d1);   // = -2*R[0]
                                let g1 = 2.0 * diff.dot(d2);   // = 2*R[1]
                                let det = h00 * h11 - h01 * h01;
                                if det.abs() < 1e-30 { break; }
                                // H路螖t = [-g0, g1] �?M路螖t = [R[0], R[1]] (OCCT L245-250)
                                let dt1 = (-g0 * h11 - g1 * h01) / det;
                                let dt2 = (g1 * h00 + g0 * h10) / det;
                                let new_t1 = (nr_t1 + dt1).clamp(tr1[0], tr1[1]);
                                let new_t2 = (nr_t2 + dt2).clamp(tr2[0], tr2[1]);
                                if (new_t1 - nr_t1).abs() < 1e-12 && (new_t2 - nr_t2).abs() < 1e-12 { break; }
                                nr_t1 = new_t1; nr_t2 = new_t2;
                            }
                            let nr_d = c1.point_at(nr_t1).distance(c2.point_at(nr_t2));
                            if nr_d < best_d { best_d = nr_d; best_t1 = nr_t1; best_t2 = nr_t2; }
                            if best_d <= fuzzy {
                                let best_pt = c1.point_at(best_t1);
                                self.ds.interferences.push(Interference::EdgeEdge {
                                    e1, e2, point: best_pt,
                                    param1: best_t1, param2: best_t2, new_vertex: vi,
                                });
                            }
                        }
                    }
                }
            }
        }
    }

    fn force_interf_ef(&mut self) {
        // OCCT L787-821: collect all PBs (skip edges without PBs or degenerated)
        for ei in 0..self.ds.edges.len() {
            let edge = &self.ds.edges[ei];
            if edge.pave_blocks.is_empty() { continue; }
            // OCCT L804-808: skip degenerated edges
            if self.ds.is_edge_degenerated(ei) { continue; }

            for fi in 0..self.ds.faces.len() {
                // OCCT L1150: skip same-origin pairs
                if edge.origin == self.ds.faces[fi].origin { continue; }

                // OCCT L953-955: skip PBs already in face's On/In/Sc sets
                if self.ds.interferences.iter().any(|inf| {
                    matches!(inf, Interference::EdgeFace { edge: e, face: f, .. } if *e == ei && *f == fi)
                }) { continue; }

                let face = &self.ds.faces[fi];

                // OCCT L915-924: collect ALL face vertices (VerticesOn + VerticesIn +
                //   VerticesSc).  rcad: vertices_on + vertices_in + curves_sc endpoints.
                let mut face_verts = face.face_info.vertices_on.clone();
                face_verts.extend(&face.face_info.vertices_in);
                // OCCT VerticesSc: section-curve vertices from FF intersection.
                for &ci in &face.face_info.curves_sc {
                    if ci < self.ds.intersection_curves.len() {
                        let ic = &self.ds.intersection_curves[ci];
                        face_verts.insert(ic.start_vertex);
                        face_verts.insert(ic.end_vertex);
                    }
                }

                // OCCT L958-964: check if PB's vertices are in the face's vertex set
                for pb in &edge.pave_blocks {
                    if !face_verts.contains(&pb.pave1.vertex_idx)
                        || !face_verts.contains(&pb.pave2.vertex_idx)
                    {
                        continue;
                    }

                    // OCCT L970-976: aTolAdd = max(tol(V1), tol(V2))
                    let gt1 = self.ds.vertices.get(pb.pave1.vertex_idx).map(|v| v.geom_tol).unwrap_or(0.0);
                    let gt2 = self.ds.vertices.get(pb.pave2.vertex_idx).map(|v| v.geom_tol).unwrap_or(0.0);
                    let a_tol_add = gt1.max(gt2);
                    // OCCT L1053: SetFuzzyValue(myFuzzyValue + aTolAdd)
                    let fuzzy = self.ds.fuzzy_tol + a_tol_add;

                    // OCCT L982-1000: project PB midpoint onto face surface
                    let mid_t = (pb.pave1.param + pb.pave2.param) * 0.5;
                    let mid_pt = edge.curve.point_at(mid_t);
                    let on_face = match &face.surface {
                        Surface3::Plane(pl) => {
                            let d = mid_pt - pl.origin;
                            let proj = d - d.dot(pl.normal) * pl.normal;
                            proj.length() <= fuzzy
                        }
                        _ => {
                            let (_, proj_pt) = crate::extrema::closest_point_on_surface(&face.surface, mid_pt);
                            mid_pt.distance(proj_pt) <= fuzzy
                        }
                    };

                    if on_face {
                        // OCCT L1007-1025: create EdgeFace interference
                        self.ds.interferences.push(Interference::EdgeFace {
                            edge: ei, face: fi,
                            point: mid_pt, edge_param: mid_t,
                            new_vertex: pb.pave1.vertex_idx,
                        });
                        break;
                    }
                }
            }
        }
    }

    fn force_interf_ve(&mut self) {
        // Build set of existing VE interferences for dedup
        let mut ve_done: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
        for inf in &self.ds.interferences {
            if let Interference::VertexEdge { vertex, edge, .. } = inf {
                ve_done.insert((*vertex, *edge));
            }
        }

        for fi in 0..self.ds.faces.len() {
            // Collect all boundary edges of this face (outer + inner wires)
            let face = &self.ds.faces[fi];
            let face_vertices: Vec<usize> = face.face_info.vertices_on
                .iter()
                .chain(face.face_info.vertices_in.iter())
                .copied()
                .collect();
            if face_vertices.is_empty() {
                continue;
            }

            let boundary_edges: Vec<usize> = {
                let f = &self.ds.faces[fi];
                let mut edges = f.boundary_edges.clone();
                for inner in &f.inner_boundary_edges {
                    for &(ei, _) in inner {
                        edges.push(ei);
                    }
                }
                edges
            };

            for &vi in &face_vertices {
                let v_origin = self.ds.vertices[vi].origin;
                if v_origin.is_none() { continue; }
                for &ei in &boundary_edges {
                    let e_origin = self.ds.edges[ei].origin;
                    if e_origin == v_origin.unwrap() { continue; }
                    if self.ds.is_edge_degenerated(ei) { continue; }
                    if ve_done.contains(&(vi, ei)) { continue; }

                    self.check_vertex_edge(vi, ei);
                }
            }
        }
    }

    fn force_interf_vf(&mut self) {
        // Build set of existing VF interferences for dedup
        let mut vf_done: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
        for inf in &self.ds.interferences {
            if let Interference::VertexFace { vertex, face } = inf {
                vf_done.insert((*vertex, *face));
            }
        }

        for vi in 0..self.ds.vertices.len() {
            let v_origin = self.ds.vertices[vi].origin;
            let opposite_faces: Vec<usize> = match v_origin {
                Some(ShapeOrigin::ShapeA) => self.faces_of(ShapeOrigin::ShapeB),
                Some(ShapeOrigin::ShapeB) => self.faces_of(ShapeOrigin::ShapeA),
                _ => continue,
            };

            for &fi in &opposite_faces {
                if vf_done.contains(&(vi, fi)) { continue; }
                self.check_vertex_face(vi, fi);
            }
        }
    }

    fn put_se_in_other_faces(&mut self) {
        let n_faces = self.ds.faces.len();
        let ics = self.ds.intersection_curves.clone();
        let mut ic_creators: Vec<Vec<usize>> = vec![Vec::new(); ics.len()];
        for inf in &self.ds.interferences {
            if let Interference::FaceFace { f1, f2, curves, .. } = inf {
                for &ci in curves { if ci < ic_creators.len() { ic_creators[ci].push(*f1); ic_creators[ci].push(*f2); } }
            }
        }
        for (ci, ic) in ics.iter().enumerate() {
            let creators = &ic_creators[ci];
            if creators.is_empty() { continue; }
            let mid_t = (ic.t_range[0] + ic.t_range[1]) * 0.5;
            let params = if (ic.t_range[1] - ic.t_range[0]).abs() < TOLERANCE_ABS { vec![mid_t] }
            else { vec![ic.t_range[0]*0.9+ic.t_range[1]*0.1, mid_t, ic.t_range[0]*0.1+ic.t_range[1]*0.9] };
            for fi in 0..n_faces {
                if creators.contains(&fi) { continue; }
                if !self.ds.faces[fi].face_info.has_any_interference() { continue; }
                let on_face = params.iter().any(|&t| {
                    use rcad_kernel::geom::CurveEval;
                    let pt = ic.curve.point_at(t);
                    let tol = self.ds.faces[fi].geom_tol.max(TOLERANCE_ABS);
                    self.point_on_face(pt, fi, tol)
                });
                if on_face {
                    self.ds.faces[fi].face_info.curves_sc.insert(ci);
                    self.ds.faces[fi].face_info.vertices_in.insert(ic.start_vertex);
                    self.ds.faces[fi].face_info.vertices_in.insert(ic.end_vertex);
                }
            }
        }
    }

    fn point_on_face(&self, pt: DVec3, fi: usize, tol: f64) -> bool {
        use rcad_kernel::geom::{SurfaceEval, Surface3};
        let face = &self.ds.faces[fi];
        match &face.surface {
            Surface3::Plane(p) => (pt - p.origin).dot(p.normal).abs() <= tol,
            Surface3::Sphere(s) => ((pt - s.center).length() - s.radius).abs() <= tol,
            Surface3::Cylinder(c) => { let v = pt - c.origin; let radial = v - c.axis.normalize() * v.dot(c.axis.normalize()); (radial.length() - c.radius).abs() <= tol }
            Surface3::Cone(c) => { let v = pt - c.apex; let a = c.axis_dir(); let pj = v.dot(a); (v - a*pj).length(); false }
            _ => false,
        }
    }

    fn process_de(&mut self) {
        let mut fv: Vec<Vec<usize>> = vec![Vec::new(); self.ds.faces.len()];
        let mut degen_flags: Vec<(usize, usize)> = Vec::new();
        for fi in 0..self.ds.faces.len() {
            let is_periodic = matches!(self.ds.faces[fi].surface,
                Surface3::Sphere(_) | Surface3::Cylinder(_) | Surface3::Cone(_));
            if !is_periodic { continue; }
            let mut vs = Vec::new();
            let boundary_edges: Vec<usize> = self.ds.faces[fi].boundary_edges.clone();
            for &ei in &boundary_edges {
                if ei >= self.ds.edges.len() { continue; }
                let e = &self.ds.edges[ei];
                if e.start_vertex == e.end_vertex {
                    degen_flags.push((ei, fi + 1));
                    vs.push(e.start_vertex);
                } else {
                    let d = (self.ds.vertices[e.start_vertex].point - self.ds.vertices[e.end_vertex].point).length();
                    if d < TOLERANCE_ABS*100.0 {
                        degen_flags.push((ei, fi + 1));
                        vs.push(e.start_vertex); vs.push(e.end_vertex);
                    }
                }
            }
            if let Surface3::Sphere(sp) = &self.ds.faces[fi].surface.clone() {
                let tl = self.ds.faces[fi].geom_tol.max(TOLERANCE_ABS);
                for pp in [sp.center + sp.axis * sp.radius, sp.center - sp.axis * sp.radius] {
                    if let Some(vi) = self.ds.find_vertex_near(pp, tl) { vs.push(vi); }
                }
            }
            fv[fi] = vs;
        }
        for (fi, vs) in fv.iter().enumerate() { for &vi in vs { self.ds.faces[fi].face_info.vertices_in.insert(vi); } }
        for (ei, flag_val) in degen_flags { self.ds.set_edge_flag(ei, flag_val); }
    }

    fn fill_shrunk_data(&mut self) {
        let ec: Vec<Curve3> = self.ds.edges.iter().map(|e| e.curve.clone()).collect();
        let et: Vec<f64> = self.ds.edges.iter().map(|e| e.geom_tol).collect();
        let cf = TOLERANCE_ABS;
        for ei in 0..self.ds.edges.len() {
            for pb in &mut self.ds.edges[ei].pave_blocks {
                let v1 = self.ds.vertices[pb.pave1.vertex_idx].geom_tol;
                let v2 = self.ds.vertices[pb.pave2.vertex_idx].geom_tol;
                if let Some(sr) = crate::inttools::curve_range::shrunk_range(&ec[ei], [pb.pave1.param, pb.pave2.param], v1, v2, et[ei]) {
                    pb.shrunk_range = Some(sr); pb.is_splittable = (sr[1]-sr[0]) > 2.0*et[ei] + 2.0*cf;
                } else { pb.shrunk_range = None; pb.is_splittable = false; }
            }
        }
        for pb in &mut self.ds.pave_blocks {
            if pb.original_edge >= self.ds.edges.len() { continue; }
            let v1 = self.ds.vertices[pb.pave1.vertex_idx].geom_tol;
            let v2 = self.ds.vertices[pb.pave2.vertex_idx].geom_tol;
            if let Some(sr) = crate::inttools::curve_range::shrunk_range(&ec[pb.original_edge], [pb.pave1.param, pb.pave2.param], v1, v2, et[pb.original_edge]) {
                pb.shrunk_range = Some(sr); pb.is_splittable = (sr[1]-sr[0]) > 2.0*et[pb.original_edge] + 2.0*cf;
            } else { pb.shrunk_range = None; pb.is_splittable = false; }
        }
    }

    fn existing_pave_block(&self, ei: usize, vi: usize) -> bool {
        for pb in &self.ds.edges[ei].pave_blocks {
            if pb.pave1.vertex_idx == vi || pb.pave2.vertex_idx == vi { return true; }
        }
        false
    }

    fn split_ics_at_periodic_boundary(&mut self) {
        let n_curves = self.ds.intersection_curves.len();
        let mut curve_faces = std::collections::HashMap::new();
        for inf in &self.ds.interferences {
            if let Interference::FaceFace { f1, f2, curves, .. } = inf {
                for &ci in curves {
                    curve_faces.entry(ci).or_insert((*f1, *f2));
                }
            }
        }
        for ci in 0..n_curves {
            let Some(&(fa, fb)) = curve_faces.get(&ci) else { continue };
            let ic = &self.ds.intersection_curves[ci];
            let t0 = ic.t_range[0]; let t1 = ic.t_range[1];
            if (t1 - t0).abs() < TOLERANCE_ABS { continue; }
            for &(fi, use_a) in &[(fa, true), (fb, false)] {
                let (period, _, _) = match &self.ds.faces[fi].surface {
                    Surface3::Sphere(_) => (std::f64::consts::TAU, 0.0, std::f64::consts::TAU),
                    Surface3::Cylinder(_) => (std::f64::consts::TAU, 0.0, std::f64::consts::TAU),
                    _ => continue,
                };
                let pcurve = if use_a { ic.pcurve_on_a.as_ref() } else { ic.pcurve_on_b.as_ref() };
                let Some(pc) = pcurve else { continue };
                const N_SAMP: usize = 64;
                let mut prev_u: Option<f64> = None;
                for i in 0..=N_SAMP {
                    let frac = i as f64 / N_SAMP as f64;
                    let u = pc.point_at(t0 + (t1 - t0) * frac).x;
                    if let Some(pu) = prev_u {
                        if (u - pu).abs() > period * 0.5 {
                            let mut lo = (i as f64 - 1.0) / N_SAMP as f64;
                            let mut hi = i as f64 / N_SAMP as f64;
                            for _ in 0..32 {
                                let mid = (lo + hi) * 0.5;
                                let u_mid = pc.point_at(t0 + (t1 - t0) * mid).x;
                                if u_mid < (pu / period).round() * period { lo = mid; }
                                else { hi = mid; }
                            }
                            let t_cross = t0 + (t1 - t0) * (lo + hi) * 0.5;
                            let pt = ic.curve.point_at(t_cross);
                            let v_new = self.ds.find_vertex_near(pt, TOLERANCE_ABS * 1000.0)
                                .unwrap_or_else(|| {
                                    let vi = self.ds.vertices.len();
                                    self.ds.vertices.push(crate::bopds::ds::DSVertex {
                                        point: pt, geom_tol: TOLERANCE_ABS, origin: None,
                                        is_internal: false,
                                    });
                                    vi
                                });
                            self.ds.faces[fa].face_info.vertices_in.insert(v_new);
                            self.ds.faces[fb].face_info.vertices_in.insert(v_new);
                        }
                    }
                    prev_u = Some(u);
                }
            }
        }
    }

    fn put_pave_on_curve(&mut self, nV: usize, curve_idx: usize) -> Option<f64> {
        let t = self.project_vertex_on_curve(nV, &self.ds.intersection_curves[curve_idx])?;
        let a_tol_r3d = self.ds.intersection_curves[curve_idx].geom_tol;
        let v_tol = self.ds.vertices[nV].geom_tol;
        let a_ptol = a_tol_r3d.max(v_tol); // OCCT L3002: Resolution(max(aTolR3D, aTolV))
        if let Some(pb) = self.ds.intersection_curves[curve_idx].change_pave_block1() {
            let mut n_v_used = 0usize;
            // OCCT L3004: bExist = aPB->ContainsParameter(aT, aPTol, nVUsed)
            if pb.contains_parameter(t, a_ptol, &mut n_v_used) {
                return Some(t);
            }
            pb.append_ext_pave(Pave { vertex_idx: nV, param: t });
        }
        Some(t)
    }

    fn put_paves_on_curve(&mut self, curve_idx: usize) {
        let ic = self.ds.intersection_curves[curve_idx].clone();
        let aBoxC = curve_bounding_box_simple(&ic.curve, ic.geom_tol.max(TOLERANCE_ABS));
        let aTolR3D = ic.geom_tol;
        let ef_vertices: Vec<usize> = self.ds.interferences.iter()
            .filter_map(|inf| {
                if let Interference::EdgeFace { new_vertex, .. } = inf { Some(*new_vertex) } else { None }
            }).collect();
        let ef_set: HashSet<usize> = ef_vertices.iter().copied().collect();
        let in_vertices: Vec<usize> = (0..self.ds.faces.len())
            .flat_map(|fi| self.ds.faces[fi].face_info.vertices_in.iter().copied()).collect();
        for &vi in &ef_vertices { self.put_pave_on_curve(vi, curve_idx); }
        for &vi in &in_vertices {
            if ef_set.contains(&vi) { continue; }
            if let Some([c_min, c_max]) = aBoxC {
                let v_pt = self.ds.vertices[vi].point;
                let v_tol = self.ds.vertices[vi].geom_tol.max(aTolR3D);
                let v_min = v_pt - DVec3::splat(v_tol);
                let v_max = v_pt + DVec3::splat(v_tol);
                if v_max.x < c_min.x || v_min.x > c_max.x || v_max.y < c_min.y || v_min.y > c_max.y || v_max.z < c_min.z || v_min.z > c_max.z { continue; }
            }
            if !self.ds.is_new_vertex(vi) { continue; }
            self.put_pave_on_curve(vi, curve_idx);
        }
    }

    fn project_vertex_on_curve(&self, vi: usize, ic: &IntersectionCurve) -> Option<f64> {
        use rcad_kernel::geom::CurveEval;
        let v_tol = self.ds.vertices[vi].geom_tol;
        let c_tol = ic.geom_tol;
        let f_tol = self.ds.fuzzy_tol;
        // OCCT IsVertexOnLine L800: aTolSum = aTolV + aTolC
        //   where aTolC = aTolR3D + myFuzzyValue (PutPaveOnCurve L2976)
        let raw_sum = v_tol + c_tol + f_tol;
        // OCCT L806-819: aTolSum = 2 * aTolSum, clamped to >= 1e-6
        let tl = (2.0 * raw_sum).max(1e-6);
        let result = self.project_vertex_on_curve_with_tol(vi, ic, tl);
        if result.is_some() { return result; }
        let ext_tol = self.extended_tolerance_occt(vi);
        if ext_tol > tl { self.project_vertex_on_curve_with_tol(vi, ic, ext_tol) } else { None }
    }

    fn project_vertex_on_curve_with_tol(&self, vi: usize, ic: &IntersectionCurve, tl: f64) -> Option<f64> {
        use rcad_kernel::geom::CurveEval;
        let vp = self.ds.vertices[vi].point;
        match &ic.curve {
            Curve3::Line(l) => { let v = vp - l.origin; let t = v.dot(l.direction);
                if (v - l.direction*t).length() <= tl { Some(t.clamp(ic.t_range[0], ic.t_range[1])) } else { None } }
            Curve3::Circle(c) => {
                let v = vp - c.center;
                let tl_cap = tl.min(c.radius.max(1e-7) * 1e-4);
                if (v.length() - c.radius).abs() > tl_cap { return None; }
                let nm = c.normal.normalize();
                if v.dot(nm).abs() > tl_cap { return None; }
                let xa = any_perpendicular(nm).normalize();
                let ya = nm.cross(xa);
                Some(v.dot(xa).atan2(v.dot(ya)).rem_euclid(std::f64::consts::TAU).clamp(ic.t_range[0], ic.t_range[1]))
            }
            _ => {
                for i in 0..20 {
                    let t = ic.t_range[0] + (ic.t_range[1]-ic.t_range[0])*i as f64/20.0;
                    if (ic.curve.point_at(t)-vp).length() <= tl { return Some(t); }
                }
                None
            }
        }
    }

    fn reduce_intersection_range(&self, ei: usize, param: f64, tol: f64) -> [f64; 2] {
        let edge = &self.ds.edges[ei];
        let mut t_range = edge.t_range;
        for pave in &edge.paves {
            let d = (pave.param - param).abs();
            if d < tol * 10.0 {
                if pave.param < param && pave.param > t_range[0] { t_range[0] = pave.param; }
                if pave.param > param && pave.param < t_range[1] { t_range[1] = pave.param; }
            }
        }
        t_range
    }

    fn make_pcurves(&mut self) {
        use crate::bopds::ds::DSRepOnFace;

        // Pre-collect (pb_idx, ei) pairs for all faces to avoid borrow conflicts
        let mut pcurves_to_add: Vec<(usize, usize, Curve2d, f64, f64)> = Vec::new();

        // OCCT L606: for each face in the FaceInfo pool
        for fi in 0..self.ds.faces.len() {
            let face_surface = self.ds.faces[fi].surface.clone();

            // OCCT L618-631: PaveBlocksIn + L633-699: PaveBlocksOn
            let pb_indices: Vec<usize> = self.ds.faces[fi].face_info.pave_blocks_in
                .iter()
                .chain(self.ds.faces[fi].face_info.pave_blocks_on.iter())
                .copied()
                .collect();

            for &pb_idx in &pb_indices {
                // OCCT L624: nE = aPB->Edge()
                if pb_idx >= self.ds.pave_blocks.len() { continue; }
                let ei = self.ds.pave_blocks[pb_idx].original_edge;
                if ei >= self.ds.edges.len() { continue; }

                // OCCT L641: check if pcurve already exists (HasCurveOnSurface)
                if self.ds.edges[ei].face_reps.iter().any(|r| r.face_idx == fi) { continue; }

                // Compute pcurve: project edge 3D curve onto face surface
                let edge_curve = &self.ds.edges[ei].curve;
                if let Some((pcurve, len)) = DS::compute_edge_pcurve(edge_curve, &face_surface) {
                    pcurves_to_add.push((ei, fi, pcurve, self.ds.edges[ei].t_range[0], self.ds.edges[ei].t_range[1]));
                    let _ = len; // length may be unused
                }
            }
        }

        // Apply all collected pcurves
        for (ei, fi, pcurve, t0, t1) in pcurves_to_add {
            self.ds.edges[ei].face_reps.push(DSRepOnFace {
                face_idx: fi,
                pcurve,
                pcurve2: None,
                pcurve_range: [t0, t1],
                start_param: t0,
                end_param: t1,
            });
        }
    }
    fn update_face_info(&mut self, fi: usize) {
        let face = &self.ds.faces[fi].clone();
        let mut sc_curves: Vec<usize> = face.face_info.curves_sc.iter().copied().collect();
        // Recompute vertices_in from all IC endpoints on this face
        for &ci in &sc_curves {
            if ci < self.ds.intersection_curves.len() {
                let ic = &self.ds.intersection_curves[ci];
                self.ds.faces[fi].face_info.vertices_in.insert(ic.start_vertex);
                self.ds.faces[fi].face_info.vertices_in.insert(ic.end_vertex);
            }
        }
        // Update pave_blocks_in/pave_blocks_on from edge pave blocks
        for &ei in &face.boundary_edges {
            if ei < self.ds.edges.len() {
                for pb in &self.ds.edges[ei].pave_blocks {
                    let pb_idx = self.ds.pave_blocks.len();
                    self.ds.pave_blocks.push(pb.clone());
                    self.ds.faces[fi].face_info.pave_blocks_in.insert(pb_idx);
                }
            }
        }
    }

    fn check_planes(&self, f1: usize, f2: usize) -> bool {
        let s1 = &self.ds.faces[f1].surface;
        let s2 = &self.ds.faces[f2].surface;
        match (s1, s2) {
            (Surface3::Plane(p1), Surface3::Plane(p2)) => {
                let dot = p1.normal.dot(p2.normal).abs();
                if dot > 0.9999 {
                    let d = (p2.origin - p1.origin).dot(p1.normal).abs();
                    d < TOLERANCE_ABS * 100.0
                } else { false }
            }
            _ => false,
        }
    }

    fn make_sd_vertices(&mut self, verts: &[usize]) {
        if verts.len() < 2 { return; }
        let tol = TOLERANCE_ABS * 100.0;
        let mut groups: Vec<Vec<usize>> = Vec::new();
        let mut used = vec![false; verts.len()];
        for i in 0..verts.len() {
            if used[i] { continue; }
            let mut group = vec![verts[i]];
            used[i] = true;
            for j in (i+1)..verts.len() {
                if used[j] { continue; }
                let d = (self.ds.vertices[verts[i]].point - self.ds.vertices[verts[j]].point).length();
                if d < tol { group.push(verts[j]); used[j] = true; }
            }
            if group.len() > 1 { groups.push(group); }
        }
        // For each group, keep the first vertex as SD, remap others
        for group in &groups {
            let sd_vi = group[0];
            for &vi in &group[1..] {
                for ei in 0..self.ds.edges.len() {
                    for pb in &mut self.ds.edges[ei].pave_blocks {
                        if pb.pave1.vertex_idx == vi { pb.pave1.vertex_idx = sd_vi; }
                        if pb.pave2.vertex_idx == vi { pb.pave2.vertex_idx = sd_vi; }
                    }
                }
                for fi in 0..self.ds.faces.len() {
                    if self.ds.faces[fi].face_info.vertices_in.contains(&vi) {
                        self.ds.faces[fi].face_info.vertices_in.remove(&vi);
                        self.ds.faces[fi].face_info.vertices_in.insert(sd_vi);
                    }
                }
            }
        }
    }
    fn get_stick_vertices(&self, fi: usize) -> Vec<usize> {
        let mut verts: Vec<usize> = Vec::new();
        for inf in &self.ds.interferences {
            match inf {
                Interference::VertexVertex { v1, v2, .. } => {
                    if self.vertex_on_face(*v1, fi) { verts.push(*v1); }
                    if self.vertex_on_face(*v2, fi) { verts.push(*v2); }
                }
                Interference::VertexEdge { vertex, .. } => {
                    if self.vertex_on_face(*vertex, fi) { verts.push(*vertex); }
                }
                Interference::VertexFace { vertex, face } => {
                    if *face == fi { verts.push(*vertex); }
                }
                _ => {}
            }
        }
        verts.sort(); verts.dedup();
        verts
    }

    fn vertex_on_face(&self, vi: usize, fi: usize) -> bool {
        if vi >= self.ds.vertices.len() || fi >= self.ds.faces.len() { return false; }
        let face = &self.ds.faces[fi];
        let tol = face.geom_tol.max(TOLERANCE_ABS);
        for &bvi in &face.boundary_verts {
            if bvi == vi { return true; }
        }
        self.ds.faces[fi].face_info.vertices_in.contains(&vi)
    }

    fn is_existing_vertex(&self, ci: usize, param: f64) -> bool {
        let ic = &self.ds.intersection_curves[ci];
        let tol = ic.geom_tol.max(TOLERANCE_ABS) * 100.0;
        for inf in &self.ds.interferences {
            if let Interference::FaceFace { curves, .. } = inf {
                if curves.contains(&ci) { continue; }
            }
        }
        // Check if any vertex already exists at this parameter
        for fi in 0..self.ds.faces.len() {
            for &vi in &self.ds.faces[fi].face_info.vertices_in {
                if let Some(t) = self.project_vertex_on_curve(vi, ic) {
                    if (t - param).abs() < tol { return true; }
                }
            }
        }
        false
    }

    fn estimate_pave_on_curve(&self, ci: usize, vi: usize) -> Option<f64> {
        let ic = &self.ds.intersection_curves[ci];
        let v_tol = self.ds.vertices[vi].geom_tol;
        let c_tol = ic.geom_tol;
        // OCCT IsVertexOnLine 3-param overload (L4066):
        //   aTolV = BRep_Tool::Tolerance(aV)        �?v_tol
        //   aTolC = aTolR3D                         �?c_tol (no fuzzy)
        //   5-param: aTolSum = aTolV + aTolC; ×2
        let tl = (2.0 * (v_tol + c_tol)).max(1e-6);
        self.project_vertex_on_curve_with_tol(vi, ic, tl)
    }

    fn extended_tolerance_occt(&self, vi: usize) -> f64 {
        let base_tol = if vi < self.ds.vertices.len() {
            self.ds.vertices[vi].geom_tol
        } else { TOLERANCE_ABS };
        let mut max_tol = base_tol;

        // For newly created vertices, check EE/EF interferences
        for inf in &self.ds.interferences {
            let (ei, param) = match inf {
                Interference::EdgeEdge { e1, param1, new_vertex, .. } if *new_vertex == vi => {
                    (*e1, *param1)
                }
                Interference::EdgeFace { edge, edge_param, new_vertex, .. } if *new_vertex == vi => {
                    (*edge, *edge_param)
                }
                _ => continue,
            };
            // Compute distance from vertex to edge endpoints at the parameter
            if ei < self.ds.edges.len() {
                if let Some(pt) = {
                    use rcad_kernel::geom::CurveEval;
                    Some(self.ds.edges[ei].curve.point_at(param))
                } {
                    let v_pt = self.ds.vertices[vi].point;
                    let d = (v_pt - pt).length();
                    if d > max_tol { max_tol = d; }
                }
            }
        }
        max_tol * 2.0  // Safety factor matching OCCT's approach
    }

    fn get_ef_pnts(&self, ei: usize) -> Vec<(usize, f64)> {
        let mut pnts: Vec<(usize, f64)> = Vec::new();
        for inf in &self.ds.interferences {
            if let Interference::EdgeFace { edge, edge_param, new_vertex, .. } = inf {
                if *edge == ei { pnts.push((*new_vertex, *edge_param)); }
            }
        }
        pnts
    }
    fn treat_vertices_ee(&mut self) {
        let mut to_merge: Vec<(usize, usize)> = Vec::new();
        for inf in &self.ds.interferences {
            if let Interference::EdgeEdge { new_vertex, e1, e2, .. } = inf {
                let vi = *new_vertex;
                if vi >= self.ds.vertices.len() { continue; }
                let v_pt = self.ds.vertices[vi].point;
                for inf2 in &self.ds.interferences {
                    if let Interference::EdgeEdge { new_vertex: vi2, .. } = inf2 {
                        if *vi2 != vi && *vi2 < self.ds.vertices.len() {
                            let d = (v_pt - self.ds.vertices[*vi2].point).length();
                            if d < TOLERANCE_ABS * 100.0 {
                                to_merge.push((vi, *vi2));
                            }
                        }
                    }
                }
            }
        }
        for (keep, remove) in &to_merge {
            for ei in 0..self.ds.edges.len() {
                for pb in &mut self.ds.edges[ei].pave_blocks {
                    if pb.pave1.vertex_idx == *remove { pb.pave1.vertex_idx = *keep; }
                    if pb.pave2.vertex_idx == *remove { pb.pave2.vertex_idx = *keep; }
                }
            }
        }
    }

    fn check_face_paves(&self, fi: usize, vi: usize) -> bool {
        let face = &self.ds.faces[fi];
        let tol = face.geom_tol.max(TOLERANCE_ABS);
        for &ei in &face.boundary_edges {
            if ei >= self.ds.edges.len() { continue; }
            for pb in &self.ds.edges[ei].pave_blocks {
                if pb.pave1.vertex_idx == vi || pb.pave2.vertex_idx == vi {
                    return true;
                }
            }
        }
        false
    }

    fn get_full_shape_map(&self, fi: usize) -> Vec<usize> {
        let mut indices: Vec<usize> = Vec::new();
        indices.push(fi);
        for &ei in &self.ds.faces[fi].boundary_edges { indices.push(ei); }
        for &vi in &self.ds.faces[fi].boundary_verts { indices.push(vi); }
        indices
    }

    fn remove_used_vertices(&self, verts: &mut Vec<usize>, used: &std::collections::BTreeSet<usize>) {
        verts.retain(|v| !used.contains(v));
    }

    fn correct_t_range(&self, ei: usize, t_start: f64, t_end: f64) -> [f64; 2] {
        let edge = &self.ds.edges[ei];
        let mut ts = t_start.max(edge.t_range[0]);
        let mut te = t_end.min(edge.t_range[1]);
        if te < ts { std::mem::swap(&mut ts, &mut te); }
        [ts, te]
    }

    fn is_block_in_on_face(&self, ei: usize, pbi: usize, fi: usize) -> bool {
        if ei >= self.ds.edges.len() { return false; }
        if fi >= self.ds.faces.len() { return false; }
        let pb = &self.ds.edges[ei].pave_blocks[pbi];
        let mid_param = (pb.pave1.param + pb.pave2.param) * 0.5;
        let mid_pt = self.ds.edges[ei].curve.point_at(mid_param);
        let face = &self.ds.faces[fi];
        let tol = face.geom_tol.max(TOLERANCE_ABS);
        match &face.surface {
            Surface3::Plane(p) => (mid_pt - p.origin).dot(p.normal).abs() <= tol,
            _ => false,
        }
    }

    fn post_treat_ff(&mut self) {
        // OCCT PaveFiller_6.cxx L1165-1397: PostTreatFF handles VertsUnused,
        // MSCPB processing, missing-curve recomputation via sub-PaveFiller, and
        // PreparePostTreatFF.  rcad distributes these across:
        //   - post_treat_ff (here): register curves_sc + vertices_in from FF curves
        //   - make_sd_vertices_ff (below): SD vertex creation
        //   - refine_face_info_in/on (DS): face info refinement
        //   - make_blocks / remove_micro_edges (elsewhere): PB/micro-edge handling
        let n_faces = self.ds.faces.len();
        let mut face_boundary_verts: Vec<Vec<usize>> = Vec::with_capacity(n_faces);
        for fi in 0..n_faces {
            face_boundary_verts.push(self.ds.faces[fi].boundary_verts.clone());
        }

        for inf in &self.ds.interferences.clone() {
            if let Interference::FaceFace { f1, f2, curves, .. } = inf {
                if curves.is_empty() { continue; }

                for &ci in curves {
                    self.ds.faces[*f1].face_info.curves_sc.insert(ci);
                    self.ds.faces[*f2].face_info.curves_sc.insert(ci);

                    // �?OCCT-aligned: register section edge PaveBlocks into FaceInfo::PaveBlocksSc
                    //   (PaveFiller_6.cxx L1700-1734).
                    if ci < self.ds.intersection_curves.len() {
                        // Extract all data before mutating self.ds (avoid borrow conflicts)
                        let ic = &self.ds.intersection_curves[ci];
                        let sv = ic.start_vertex;
                        let ev = ic.end_vertex;
                        let sc_pbs: Vec<PaveBlock> = ic.pave_blocks.iter()
                            .filter(|pb| pb.new_edge.is_some())
                            .cloned().collect();
                        let pb_count = sc_pbs.len();
                        drop(ic); // release immutable borrow before mutable access

                        for pb in sc_pbs {
                            let g_pb_idx = self.ds.allocate_pave_block(pb);
                            self.ds.faces[*f1].face_info.pave_blocks_sc.insert(g_pb_idx);
                            self.ds.faces[*f2].face_info.pave_blocks_sc.insert(g_pb_idx);
                        }

                        self.ds.faces[*f1].face_info.vertices_in.insert(sv);
                        self.ds.faces[*f1].face_info.vertices_in.insert(ev);
                        self.ds.faces[*f2].face_info.vertices_in.insert(sv);
                        self.ds.faces[*f2].face_info.vertices_in.insert(ev);

                        // Also register curve endpoints as vertices_on if they match
                        // boundary vertices of either face
                        for &fi in &[*f1, *f2] {
                            if fi < face_boundary_verts.len() {
                                for &bvi in &face_boundary_verts[fi] {
                                    if bvi == sv || bvi == ev {
                                        self.ds.faces[fi].face_info.vertices_on.insert(bvi);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    fn make_sd_vertices_ff(&mut self) {
        let overlaps = self.ds.same_domain_overlaps.clone();
        for (f1, f2, polygon) in &overlaps {
            let tol = self.ff_tol(*f1, *f2);
            for &pt in polygon {
                let v_idx = if let Some(idx) = self.ds.find_vertex_near(pt, tol) {
                    idx
                } else {
                    self.ds.add_vertex(pt)
                };
                self.ds.faces[*f1].face_info.vertices_in.insert(v_idx);
                self.ds.faces[*f2].face_info.vertices_in.insert(v_idx);
            }
        }
    }

    fn perform_ff(&mut self) {
        // OCCT PaveFiller_6.cxx: FillShrunkData + BVH pair iteration
        self.fill_shrunk_data(); // OCCT: FillShrunkData(FACE, FACE)
        let a_faces = self.faces_of(ShapeOrigin::ShapeA);
        let b_faces = self.faces_of(ShapeOrigin::ShapeB);

        if a_faces.is_empty() || b_faces.is_empty() {
            return;
        }

        if let (Some(bvh_a), Some(bvh_b)) = (self.bvh_a, self.bvh_b) {
            // Build reverse maps: BRep face index �?position in a_faces/b_faces
            let a_max_idx = a_faces.iter().map(|&dsi| self.ds.faces[dsi].source_face_idx).max().unwrap_or(0);
            let b_max_idx = b_faces.iter().map(|&dsi| self.ds.faces[dsi].source_face_idx).max().unwrap_or(0);
            let mut a_rev = vec![usize::MAX; a_max_idx + 1];
            for (pos, &dsi) in a_faces.iter().enumerate() {
                a_rev[self.ds.faces[dsi].source_face_idx] = pos;
            }
            let mut b_rev = vec![usize::MAX; b_max_idx + 1];
            for (pos, &dsi) in b_faces.iter().enumerate() {
                b_rev[self.ds.faces[dsi].source_face_idx] = pos;
            }

            let candidates = Bvh::candidate_pairs(bvh_a, bvh_b);
            let mut processed_pairs = std::collections::HashSet::new();
            for (fa_brep, fb_brep) in candidates {
                if let (Some(&ai), Some(&bi)) = (a_rev.get(fa_brep), b_rev.get(fb_brep))
                    && ai != usize::MAX && bi != usize::MAX {
                        // �?OCCT-aligned: BVH may produce duplicate candidate pairs when a face appears
                        //    in multiple intersecting leaf nodes, causing duplicate intersection curves.
                        //    OCCT PaveFiller processes each face pair once (FF matrix uses BOPDS_IndexRange
                        //    to mark pairs as already processed).
                        if !processed_pairs.insert((ai, bi)) { continue; }
                        let af = a_faces[ai];
                        let bf = b_faces[bi];
                        if self.should_skip_glued_face_pair(af, bf) {
                            continue;
                        }
                        self.intersect_face_face(af, bf);
                    }
            }
        } else {
            // �?OCCT-aligned: BOPDS_Iterator cross-group face pair iteration.
            let a_fcount = self.ds.a_face_count;
            let mut fit = crate::bopds::ds::PairIterator::prepare_ab(a_fcount, self.ds.faces.len());
            while fit.more() {
                let pk = fit.value();
                let af = pk.i1; let bf = pk.i2;
                // OCCT: myDS->HasInterf(nF1, nF2) �?skip if already interfered
                if self.ds.has_interf_ff(af, bf) { fit.next(); continue; }
                if !self.should_skip_glued_face_pair(af, bf) {
                    self.intersect_face_face(af, bf);
                }
                fit.next();
            }
        }
    }

    fn should_skip_glued_face_pair(&self, f1: usize, f2: usize) -> bool {
        if !self.use_glue {
            return false;
        }

        // Use pre-detected fully-glued faces if available
        if self.ds.is_fully_glued_face_pair(f1, f2) {
            return true;
        }

        let face1 = &self.ds.faces[f1];
        let face2 = &self.ds.faces[f2];
        if face1.origin == face2.origin {
            return false;
        }
        if !self.surfaces_glue_compatible(&face1.surface, &face2.surface) {
            return false;
        }

        let n1_len2 = face1.normal.length_squared();
        let n2_len2 = face2.normal.length_squared();
        if n1_len2 <= TOLERANCE_ABS || n2_len2 <= TOLERANCE_ABS {
            return false;
        }
        let n1 = face1.normal / n1_len2.sqrt();
        let n2 = face2.normal / n2_len2.sqrt();
        if n1.dot(n2) > -0.99 {
            return false;
        }

        self.boundaries_fully_overlap(f1, f2)
    }

    fn is_seam_edge(&self, edge_idx: usize, face_idx: usize) -> bool {
        let face = &self.ds.faces[face_idx];
        let edge = &self.ds.edges[edge_idx];

        match &face.surface {
            Surface3::Cylinder(cyl) => {
                // Cylinder seam edge: Line3 parallel to axis
                if let Curve3::Line(line) = &edge.curve {
                    let dir = line.direction.normalize();
                    let axis = cyl.axis.normalize();
                    dir.dot(axis).abs() > 1.0 - TOLERANCE_ABS
                } else {
                    false
                }
            }
            Surface3::Sphere(sph) => {
                // OCCT BOPAlgo_PaveFiller_6.cxx L106-134 IsClosedFF:
                // Sphere seam edge = great circle arc in meridian plane (U=0 boundary).
                // Checks mirror OCCT exactly:
                //   (1) Curve is Geom_Circle  �? Curve3::Circle
                //   (2) |center - S.Location()| < Precision::Confusion()  �? TOLERANCE_ABS_SQ
                //   (3) |radius - S.Radius| < Precision::Confusion()     �? TOLERANCE_ABS
                //   (4) |circle_normal �?sphere_axis| < Precision::Angular()  �? perp check
                match &edge.curve {
                    Curve3::Circle(c) => {
                        (c.center - sph.center).length_squared() < TOLERANCE_ABS_SQ
                        && (c.radius - sph.radius).abs() < TOLERANCE_ABS
                        && c.normal.normalize().dot(sph.axis.normalize()).abs() < 1e-12
                    }
                    _ => false,
                }
            }
            Surface3::Torus(tor) => {
                // OCCT IsClosedFF: torus has TWO periodic boundaries.
                // U-seam: major circle, center = torus center, radius = major_radius,
                //         normal �?torus axis.
                // V-seam: minor circle, center on major circle, radius = minor_radius,
                //         normal �?torus axis.
                // All tolerances match OCCT Precision::Confusion/Angular.
                match &edge.curve {
                    Curve3::Circle(c) => {
                        let axis = tor.axis.normalize();
                        let c_normal = c.normal.normalize();
                        let center_dist = (c.center - tor.center).length();
                        // U-seam: center at torus center, normal �?axis, radius = major
                        let is_u_seam = center_dist < TOLERANCE_ABS
                            && c_normal.dot(axis).abs() > 1.0 - 1e-12
                            && (c.radius - tor.major_radius).abs() < TOLERANCE_ABS;
                        // V-seam: center on major circle, normal �?axis, radius = minor
                        let on_major = (center_dist - tor.major_radius).abs() < TOLERANCE_ABS;
                        let is_v_seam = on_major
                            && c_normal.dot(axis).abs() < 1e-12
                            && (c.radius - tor.minor_radius).abs() < TOLERANCE_ABS;
                        is_u_seam || is_v_seam
                    }
                    _ => false,
                }
            }
            _ => false,
        }
    }

    fn check_seam_edge_shift(&self, f1: usize, f2: usize) -> Option<SeamEdgeShift> {
        let s1 = &self.ds.faces[f1].surface;
        let s2 = &self.ds.faces[f2].surface;

        // Skip if both faces are Planes (seam edges only exist on periodic surfaces)
        if matches!(s1, Surface3::Plane(_)) && matches!(s2, Surface3::Plane(_)) {
            return None;
        }

        for &e1 in &self.ds.faces[f1].boundary_edges {
            let is_closed1 = self.is_seam_edge(e1, f1);
            for &e2 in &self.ds.faces[f2].boundary_edges {
                let is_closed2 = self.is_seam_edge(e2, f2);
                if !is_closed1 && !is_closed2 {
                    continue;
                }

                // Look for EE interference between this edge pair
                for inf in &self.ds.interferences {
                    if let Interference::EdgeEdge {
                        e1: ee1,
                        e2: ee2,
                        point,
                        new_vertex,
                        ..
                    } = inf
                    {
                        if !((*ee1 == e1 && *ee2 == e2) || (*ee1 == e2 && *ee2 == e1)) {
                            continue;
                        }

                        // Project the EE vertex point onto both edges' 3D curves
                        // (OCCT: GeomAPI_ProjectPointOnCurve)
                        let curve1 = &self.ds.edges[e1].curve;
                        let curve2 = &self.ds.edges[e2].curve;
                        let proj1 = closest_point_on_curve(curve1, *point, 64);
                        let proj2 = closest_point_on_curve(curve2, *point, 64);

                        let a_p1 = proj1.point;
                        let a_p2 = proj2.point;
                        let shift_dist = a_p1.distance(a_p2);

                        // OCCT-aligned: the seam edge shift is a SMALL tolerance
                        // correction, not a geometric transformation.  Verify both
                        // projections are close to the EE vertex �?if either is
                        // far, the vertex is not near both edges and shifting would
                        // be invalid (e.g. sphere center jumps by 1 unit).
                        let vtx_pt = *point;
                        let d1 = a_p1.distance(vtx_pt);
                        let d2 = a_p2.distance(vtx_pt);
                        // OCCT's shift is a sub-tolerance adjustment.  A projection
                        // error exceeding 1e-4 means the vertex is not on this edge.
                        let sanity_tol = TOLERANCE_ABS * 1000.0;
                        if d1 > sanity_tol || d2 > sanity_tol {
                            continue;
                        }

                        // Check if the shift exceeds vertex tolerance
                        let vtx_tol = self.ds.vertices[*new_vertex].geom_tol;
                        if shift_dist > vtx_tol {
                            // OCCT: shift the face with the closed/seam edge
                            let shift_vector = if is_closed1 {
                                a_p2 - a_p1 // Shift f1: move aP1 toward aP2
                            } else {
                                a_p1 - a_p2 // Shift f2: move aP2 toward aP1
                            };

                            return Some(SeamEdgeShift {
                                shift_vector,
                                shift_value: shift_dist,
                                shifted_face: if is_closed1 { 1 } else { 2 },
                            });
                        }
                    }
                }
            }
        }
        None
    }

    fn reverse_seam_edge_shift(&mut self, f1: usize, f2: usize, shift: &SeamEdgeShift) {
        let inv_vec = if shift.shifted_face == 1 {
            -shift.shift_vector
        } else {
            shift.shift_vector
        };

        // Collect curve indices from the FaceFace interference for this pair
        let mut curve_indices: Vec<usize> = Vec::new();
        for inf in &self.ds.interferences {
            if let Interference::FaceFace {
                f1: a,
                f2: b,
                curves,
                ..
            } = inf
            {
                if (*a == f1 && *b == f2) || (*a == f2 && *b == f1) {
                    curve_indices = curves.clone();
                    break;
                }
            }
        }

        // Reverse shift on each curve
        for &ci in &curve_indices {
            if ci >= self.ds.intersection_curves.len() {
                continue;
            }
            let ic = &mut self.ds.intersection_curves[ci];

            // Translate 3D curve back by inverse shift
            ic.curve = translate_curve3(&ic.curve, inv_vec);

            // Translate polyline points if any
            for p in &mut ic.polyline {
                *p += inv_vec;
            }

            // Translate vertex positions back
            let sv = ic.start_vertex;
            let ev = ic.end_vertex;
            if sv < self.ds.vertices.len() {
                self.ds.vertices[sv].point += inv_vec;
            }
            if ev < self.ds.vertices.len() {
                self.ds.vertices[ev].point += inv_vec;
            }
        }
    }

    fn intersect_face_face(&mut self, f1: usize, f2: usize) {
        // 鈹€鈹€ Seam Edge Shift (OCCT PaveFiller_6.cxx L393-479) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        let shift_info = self.check_seam_edge_shift(f1, f2);
        let old_shift_tol = self.seam_shift_tol;
        if let Some(ref info) = shift_info {
            self.seam_shift_tol = info.shift_value;
        }

        let s1_orig = self.ds.faces[f1].surface.clone();
        let s2_orig = self.ds.faces[f2].surface.clone();

        // Apply seam edge shift to surface clones if needed
        let s1 = match &shift_info {
            Some(info) if info.shifted_face == 1 => {
                apply_shift_to_surface(&s1_orig, info.shift_vector)
            }
            _ => s1_orig,
        };
        let s2 = match &shift_info {
            Some(info) if info.shifted_face == 2 => {
                apply_shift_to_surface(&s2_orig, info.shift_vector)
            }
            _ => s2_orig,
        };

        // OCCT IntPatch_Intersection does NOT demote BSpline surfaces.
        // BSpline stays as Parametric (ts=0); Plane stays as Geom (ts=1).
        // The (ts1 != ts2) condition triggers ImpPrmIntersection path (marching).

        // 鈺愨晲鈺?OCCT IntPatch_Intersection 3-category dispatch 鈺愨晲鈺?
        //   OCCT IntPatch_Intersection.cxx L1298-1339 classifies surface pairs:
        //   - ts1 == ts2 == 1 : Geom-Geom (both analytic) �?ImpImpIntersection
        //   - ts1 != ts2      : Geom-Param (one analytic, one parametric) �?ImpPrmIntersection
        //   - ts1 == ts2 == 0 : Param-Param (both parametric) �?PrmPrmIntersection
        let (cat1, cat2) = (classify_surface_type(&s1), classify_surface_type(&s2));
        match (cat1, cat2) {
            // 鈹€鈹€ Geom-Geom: both analytic surfaces 鈹€鈹€
            //   OCCT ImpImpIntersection handles all analytic-analytic pairs.
            //   rcad dispatches to specialized functions per combination.
            (SurfaceCategory::GeomGeom, SurfaceCategory::GeomGeom) => {
                match (&s1, &s2) {
                    (Surface3::Plane(p1), Surface3::Plane(p2)) => {
                        self.intersect_plane_plane_faces(f1, f2, p1, p2);
                    }
                    (Surface3::Plane(pl), Surface3::Sphere(sph))
                    | (Surface3::Sphere(sph), Surface3::Plane(pl)) => {
                        self.intersect_plane_sphere_faces(f1, f2, pl, sph);
                    }
                    (Surface3::Plane(pl), Surface3::Cylinder(cyl))
                    | (Surface3::Cylinder(cyl), Surface3::Plane(pl)) => {
                        self.intersect_plane_cylinder_faces(f1, f2, pl, cyl);
                    }
                    (Surface3::Sphere(sph1), Surface3::Sphere(sph2)) => {
                        let (sph1, sph2) = (*sph1, *sph2);
                        self.intersect_sphere_sphere_faces(f1, f2, &sph1, &sph2);
                    }
                    (Surface3::Sphere(sph), Surface3::Cylinder(cyl))
                    | (Surface3::Cylinder(cyl), Surface3::Sphere(sph)) => {
                        let (sph, cyl) = (*sph, *cyl);
                        self.intersect_sphere_cylinder_faces(f1, f2, &sph, &cyl);
                    }
                    (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => {
                        let (c1, c2) = (*c1, *c2);
                        self.intersect_cylinder_cylinder_faces(f1, f2, &c1, &c2);
                    }
                    (Surface3::Plane(pl), Surface3::Cone(cone))
                    | (Surface3::Cone(cone), Surface3::Plane(pl)) => {
                        self.intersect_plane_cone_faces(f1, f2, pl, cone);
                    }
                    (Surface3::Cylinder(cyl), Surface3::Cone(cone))
                    | (Surface3::Cone(cone), Surface3::Cylinder(cyl)) => {
                        let (cyl, cone) = (*cyl, *cone);
                        self.intersect_cylinder_cone_faces(f1, f2, &cyl, &cone);
                    }
                    (Surface3::Cone(cone1), Surface3::Cone(cone2)) => {
                        let (cone1, cone2) = (*cone1, *cone2);
                        self.intersect_cone_cone_faces(f1, f2, &cone1, &cone2);
                    }
                    (Surface3::Plane(pl), Surface3::Torus(tor))
                    | (Surface3::Torus(tor), Surface3::Plane(pl)) => {
                        self.intersect_torus_plane_faces(f1, f2, tor, pl);
                    }
                    (Surface3::Sphere(sph), Surface3::Torus(tor))
                    | (Surface3::Torus(tor), Surface3::Sphere(sph)) => {
                        self.intersect_torus_sphere_faces(f1, f2, tor, sph);
                    }
                    (Surface3::Cylinder(cyl), Surface3::Torus(tor))
                    | (Surface3::Torus(tor), Surface3::Cylinder(cyl)) => {
                        self.intersect_torus_cylinder_faces(f1, f2, tor, cyl);
                    }
                    (Surface3::Cone(cone), Surface3::Torus(tor))
                    | (Surface3::Torus(tor), Surface3::Cone(cone)) => {
                        self.intersect_torus_cone_faces(f1, f2, tor, cone);
                    }
                    (Surface3::Torus(tor1), Surface3::Torus(tor2)) => {
                        self.intersect_torus_torus_faces(f1, f2, tor1, tor2);
                    }
                    (Surface3::Sphere(sph), Surface3::Cone(cone))
                    | (Surface3::Cone(cone), Surface3::Sphere(sph)) => {
                        let (sph, cone) = (*sph, *cone);
                        self.intersect_sphere_cone_faces(f1, f2, &sph, &cone);
                    }
                    _ => {}
                }
            }
            // 鈹€鈹€ Geom-Param: one analytic, one parametric 鈹€鈹€
            //   OCCT ImpPrmIntersection handles this category (IntPatch_Intersection.cxx L1326-1330).
            //   ts1 != ts2: one analytic (GeomGeom) + one parametric (ParamParam).
            //   rcad: directly use marching �?no demotion, no plane-plane redirect.
            (SurfaceCategory::GeomGeom, SurfaceCategory::ParamParam)
            | (SurfaceCategory::ParamParam, SurfaceCategory::GeomGeom) => {
                self.intersect_ff_by_marching(f1, f2);
            }
            _ => {
                // General case: numerical marching
                self.intersect_ff_by_marching(f1, f2);
            }
        }

        // 鈹€鈹€ Reverse Seam Edge Shift (OCCT ApplyTrsf L560) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        if let Some(ref info) = shift_info {
            self.reverse_seam_edge_shift(f1, f2, info);
        }
        // 鈹€鈹€ Restore seam shift tol 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        self.seam_shift_tol = old_shift_tol;

        // �?OCCT-aligned:ComputeTolReached3d + PrepareLines3D �?post-process all
        // intersection curves for this face pair.  Runs for every path (analytic,
        // numeric_intss, marching) to ensure consistent curve tolerance and
        // closed-curve splitting.
        if let Some(ff_curves) = self.find_face_face_curve_indices(f1, f2) {
            let t_a = self.ff_tol(f1, f1);
            let t_b = self.ff_tol(f2, f2);
            for &ci in &ff_curves {
                let (curve, pca, pcb, sv, ev, tr) = {
                    let ic = &self.ds.intersection_curves[ci];
                    (ic.curve.clone(), ic.pcurve_on_a.clone(), ic.pcurve_on_b.clone(),
                     ic.start_vertex, ic.end_vertex, ic.t_range)
                };
                let (new_tol, _) = inttools::pcurve_derive::compute_intersection_curve_tolerance(
                    &curve, pca.as_ref(), pcb.as_ref(),
                    &self.ds.faces[f1].surface, &self.ds.faces[f2].surface, tr,
                    t_a, t_b, 0.0,
                );
                if new_tol > TOLERANCE_ABS {
                    let vt = new_tol.min(TOLERANCE_MESH_LEGACY);
                    self.ds.vertices[sv].geom_tol = self.ds.vertices[sv].geom_tol.max(vt);
                    self.ds.vertices[ev].geom_tol = self.ds.vertices[ev].geom_tol.max(vt);
                }
            }
            // PrepareLines3D �?split closed curves
            inttools::pcurve_derive::prepare_lines_3d(&mut self.ds.intersection_curves);
            // �?OCCT-aligned: After PrepareLines3D splits closed curves, new curve endpoints
            //    must be updated to the split points. OCCT's BRepBuilderAPI_MakeEdge auto-sets
            //    endpoints when creating edges. rcad's IntersectionCurve requires explicit update:
            //    for start==end but non-full-period t_range (i.e. split half-circle), compute
            //    correct endpoint positions via point_at and create new DS vertices.
            for ci in 0..self.ds.intersection_curves.len() {
                let needs_fix = {
                    let ic = &self.ds.intersection_curves[ci];
                    let half_circle = match &ic.curve {
                        rcad_kernel::geom::Curve3::Circle(_) | rcad_kernel::geom::Curve3::Ellipse(_) => {
                            (ic.t_range[1] - ic.t_range[0] - std::f64::consts::TAU).abs() >= TOLERANCE_ANG
                        }
                        _ => false,
                    };
                    half_circle && ic.start_vertex == ic.end_vertex
                };
                if needs_fix {
                    let t0 = self.ds.intersection_curves[ci].t_range[0];
                    let t1 = self.ds.intersection_curves[ci].t_range[1];
                    let p_start = self.ds.intersection_curves[ci].curve.point_at(t0);
                    let p_end = self.ds.intersection_curves[ci].curve.point_at(t1);
                    let v_start = self.ds.add_vertex(p_start);
                    let v_end = self.ds.add_vertex(p_end);
                    self.ds.intersection_curves[ci].start_vertex = v_start;
                    self.ds.intersection_curves[ci].end_vertex = v_end;
                }
            }
            // �?OCCT-aligned: InitPaveBlock1 for all curves (PaveFiller_6.cxx L800).
            //   Creates an initial PaveBlock on each curve for ext_pave tracking.
            for ci in 0..self.ds.intersection_curves.len() {
                self.ds.intersection_curves[ci].init_pave_block1();
            }
        }
    }

    fn intersect_plane_plane_faces(&mut self, f1: usize, f2: usize, p1: &Plane, p2: &Plane) {
        use inttools::pcurve_derive::line_pcurve_on_plane;

        match inttools::plane_plane::intersect_plane_plane(p1, p2) {
            inttools::plane_plane::PlanePlaneResult::Parallel => {
            }
            inttools::plane_plane::PlanePlaneResult::Coincident => {
                self.handle_coplanar_faces(f1, f2, p1);
            }
            inttools::plane_plane::PlanePlaneResult::Line(line) => {
                let verts1 = self.ds.face_boundary_points(f1);
                let verts2 = self.ds.face_boundary_points(f2);
                let clip_tol = self.ff_tol(f1, f2);

                let ranges1 =
                    inttools::edge_face::clip_line_to_polygon_with_tol(&line, p1, &verts1, clip_tol);
                let ranges2 =
                    inttools::edge_face::clip_line_to_polygon_with_tol(&line, p2, &verts2, clip_tol);

                for &(t1_min, t1_max) in &ranges1 {
                    for &(t2_min, t2_max) in &ranges2 {
                        let t_min = t1_min.max(t2_min);
                        let t_max = t1_max.min(t2_max);
                        // Keep strict: overlap length along the intersection line is parametric, not V鈥揤
                        // coincidence �?tying this to `fuzzy_tol` can change sphere鈥揵ox trims and area.
                        if t_max - t_min < TOLERANCE_ABS {
                            continue;
                        }

                        let p_start = line.origin + line.direction * t_min;
                        let p_end = line.origin + line.direction * t_max;

                        let v_start = self.ds.add_vertex(p_start);
                        let v_end = self.ds.add_vertex(p_end);

                        let curve_idx = self.ds.intersection_curves.len();
                        let pca = line_pcurve_on_plane(&line, p1);
                        let pcb = line_pcurve_on_plane(&line, p2);
                        self.ds.intersection_curves.push(IntersectionCurve {
                            curve: Curve3::Line(line),
                            polyline: vec![],
                            start_vertex: v_start,
                            end_vertex: v_end,
                            t_range: [t_min, t_max],
                            pcurve_on_a: Some(pca),
                            pcurve_on_b: Some(pcb),
                            geom_tol: crate::tolerance::TOLERANCE_ABS,
                        pave_blocks: Vec::new(),
                        });

                        self.ds.interferences.push(Interference::FaceFace {
                            f1,
                            f2,
                            curves: vec![curve_idx],
                            points: vec![],
                        });

                        self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
                        self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
                        self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                        self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                        self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                        self.ds.faces[f2].face_info.vertices_in.insert(v_end);
                    }
                }
            }
        }
    }

    fn handle_coplanar_faces(&mut self, f1: usize, f2: usize, plane: &Plane) {
        let verts1 = self.ds.face_boundary_points(f1);
        let verts2 = self.ds.face_boundary_points(f2);

        let result = inttools::coplanar::analyze_coplanar_faces(&verts1, &verts2, plane);

        if !result.overlap.is_empty() {
            // �?OCCT-aligned: create IC for each overlap edge (BOPAlgo_PaveFiller_6.cxx:285-622)
            let plane1 = self.ds.face_plane(f1);
            let plane2 = self.ds.face_plane(f2);

            for overlap_poly in &result.overlap {
                if overlap_poly.len() < 3 { continue; }
                for i in 0..overlap_poly.len() {
                    let j = (i + 1) % overlap_poly.len();
                    let p_start = overlap_poly[i];
                    let p_end = overlap_poly[j];
                    if (p_end - p_start).length_squared() < TOLERANCE_ABS_SQ { continue; }

                    let v_start = self.ds.add_vertex(p_start);
                    let v_end = self.ds.add_vertex(p_end);
                    let dir = (p_end - p_start).normalize();
                    let len = (p_end - p_start).length();
                    let line = Line3 { origin: p_start, direction: dir };

                    let pca = inttools::pcurve_derive::line_pcurve_on_plane(&line, &plane1);
                    let pcb = inttools::pcurve_derive::line_pcurve_on_plane(&line, &plane2);

                    let curve_idx = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Line(line),
                        polyline: vec![],
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, len],
                        pcurve_on_a: Some(pca),
                        pcurve_on_b: Some(pcb),
                        geom_tol: crate::tolerance::TOLERANCE_ABS,
                    pave_blocks: Vec::new(),
                    });

                    self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
                    self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);
                }
            }

            // Keep existing same_domain_overlaps for backward compatibility
            self.ds.interferences.push(Interference::FaceFace {
                f1,
                f2,
                curves: vec![],
                points: vec![],
            });
            if let Some(overlap) = result.overlap.into_iter().max_by_key(|poly| poly.len()) {
                self.ds.same_domain_overlaps.push((f1, f2, overlap));
            }
        }
    }

    // 鈹€鈹€ Plane �?Sphere analytic face-face intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    fn clip_circle_to_faces(
        &self, circle: &rcad_kernel::geom::Circle3,
        f1: usize, f2: usize,
    ) -> Option<[f64; 2]> {
        use crate::inttools::edge_face::plane_local_basis;
        let tol = self.ff_tol(f1, f2);
        let mut result: Option<[f64; 2]> = None;
        let [t0, t1] = [0.0, std::f64::consts::TAU];

        for &fi in &[f1, f2] {
            let Surface3::Plane(plane) = &self.ds.faces[fi].surface else { continue };
            let bnd = self.ds.face_boundary_points(fi);
            if bnd.len() < 3 { continue; }
            let (u_ax, v_ax) = plane_local_basis(plane);
            let to2 = |pt: DVec3| -> DVec2 {
                let d = pt - plane.origin; DVec2::new(d.dot(u_ax), d.dot(v_ax))
            };
            let c2d = to2(DVec3::from(circle.center));
            let b2d: Vec<DVec2> = bnd.iter().map(|&pt| to2(pt)).collect();
            let mut tc: Vec<f64> = Vec::new();
            for k in 0..b2d.len() {
                let l = (k + 1) % b2d.len();
                let a = b2d[k]; let b = b2d[l];
                let ab = b - a; let ac = a - c2d;
                let qa = ab.dot(ab); let qb = 2.0 * ab.dot(ac);
                let qc = ac.dot(ac) - circle.radius * circle.radius;
                let disc = qb * qb - 4.0 * qa * qc;
                if disc < 0.0 { continue; }
                for &sign in &[-1.0_f64, 1.0_f64] {
                    let t = (-qb + sign * disc.sqrt()) / (2.0 * qa);
                    if t >= -1e-12 && t <= 1.0 + 1e-12 {
                        let pt2 = a + t.clamp(0.0, 1.0) * ab;
                        let ang = (pt2 - c2d).to_angle();
                        let mut a2 = ang;
                        if a2 < t0 { a2 += std::f64::consts::TAU; }
                        while a2 > t0 + std::f64::consts::TAU - 1e-12 { a2 -= std::f64::consts::TAU; }
                        if a2 >= t0 && a2 <= t1 { tc.push(a2); }
                    }
                }
            }
            if tc.is_empty() {
                let mut inside = false; let mut j = b2d.len() - 1;
                for i in 0..b2d.len() {
                    if ((b2d[i].y > c2d.y) != (b2d[j].y > c2d.y))
                        && (c2d.x < (b2d[j].x - b2d[i].x) * (c2d.y - b2d[i].y) / (b2d[j].y - b2d[i].y) + b2d[i].x)
                    { inside = !inside; }
                    j = i;
                }
                if !inside { return Some([0.0, 0.0]); }
                continue;
            }
            tc.sort_by(|a, b| a.partial_cmp(b).unwrap());
            // �?OCCT-aligned: deduplicate nearby angles (same intersection detected by adjacent edges)
            tc.dedup_by(|a, b| (*a - *b).abs() < TOLERANCE_ABS * 1000.0);
            // �?OCCT-aligned: select candidate arc via midpoint face-in test (IntTools_FaceFace.cxx L1084-1101)
            //    OCCT splits the full circle into 18 samples and uses dom->Classify() to test
            //    whether each UV is inside the face. rcad tests candidate arc midpoints in 2D polygon.
            let mut best: Option<[f64;2]> = None;
            for i in 0..tc.len() {
                let j = (i + 1) % tc.len();
                let a_start = tc[i];
                let a_end = if j == 0 { tc[j] + std::f64::consts::TAU } else { tc[j] };
                let a_len = a_end - a_start;
                if a_len <= tol { continue; }
                let mid_angle = (a_start + a_end) * 0.5;
                let mid_3d = circle.center + circle.radius * (mid_angle.cos() * u_ax + mid_angle.sin() * v_ax);
                let mid_2d = to2(mid_3d);
                let mut inside = false;
                {
                    let mut k = b2d.len() - 1;
                    for i2 in 0..b2d.len() {
                        if ((b2d[i2].y > mid_2d.y) != (b2d[k].y > mid_2d.y))
                            && (mid_2d.x < (b2d[k].x - b2d[i2].x) * (mid_2d.y - b2d[i2].y) / (b2d[k].y - b2d[i2].y) + b2d[i2].x)
                        { inside = !inside; }
                        k = i2;
                    }
                }
                if !inside { continue; }
                let ws = a_start % std::f64::consts::TAU;
                let we = tc[j];
                best = Some(if ws <= we { [ws, we] } else { [ws, std::f64::consts::TAU] });
                break;
            }
            let nr = match best { Some(r) => r, None => continue, };
            if nr[1] - nr[0] > tol {
                result = Some(match result {
                    Some(prev) => [prev[0].max(nr[0]), prev[1].min(nr[1])],
                    None => nr,
                });
            }
        }
        result
    }

    fn intersect_plane_sphere_faces(
        &mut self,
        f1: usize,
        f2: usize,
        plane: &Plane,
        sphere: &SphericalSurface,
    ) {
        use inttools::pcurve_derive::{
            circle_pcurve_on_plane, circle_pcurve_on_sphere, fallback_pcurve_by_projection,
        };
        use inttools::plane_sphere::{PlaneSphereResult, intersect_plane_sphere};

        // Determine which face carries the plane (for correct pcurve_on_a/b assignment)
        let plane_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Plane(_));

        let ps_result = intersect_plane_sphere(plane, sphere);
        match ps_result {
            PlaneSphereResult::NoIntersection => {}
            PlaneSphereResult::TangentPoint(pt) => {
                let verts1 = self.ds.face_boundary_points(f1);
                let verts2 = self.ds.face_boundary_points(f2);
                let tff = self.ff_tol(f1, f2);
                if inttools::edge_face::point_in_planar_face_with_tol(pt, plane, &verts1, tff)
                    && point_in_sphere_face(pt, &verts2, self.ds)
                {
                    let v = self.ds.add_vertex(pt);
                    self.ds.interferences.push(Interference::FaceFace {
                        f1,
                        f2,
                        curves: vec![],
                        points: vec![v],
                    });
                }
            }
            PlaneSphereResult::Circle(circle) => {
                // �?OCCT-aligned: clip Circle3 to planar face polygon boundaries
                // map to TWO vertical meridians in sphere UV space. We must create
                // two separate IntersectionCurves �?one per UV branch �?because
                // a single BSpline pcurve cannot span the atan2-wrap discontinuity.
                let is_great = (circle.center - sphere.center).length_squared() < TOLERANCE_ABS_SQ;
                let axis_dot_normal = sphere
                    .axis
                    .normalize()
                    .dot(plane.normal.normalize())
                    .abs();
                let _passes_poles = is_great && axis_dot_normal < TOLERANCE_ABS;

                // �?OCCT-aligned: all plane-sphere ICs go through clip_circle_to_faces unified path.
                //    add_great_circle_curves is disabled �?double-half-arc branches are rcad's own design,
                //    OCCT IntTools_FaceFace clips ICs directly using PutBoundPaveOnCurve.

                // �?OCCT-aligned: clip Circle3 to planar face polygon boundaries
                //    OCCT IntTools_Curve limits range to face boundary at creation time.
                //    rcad: project Circle3 onto plane face 2D polygon,
                //    intersect to get valid parameter range within the face, use its endpoints
                //    as start/end vertices of the curve.
                let clipped_range = self.clip_circle_to_faces(&circle, f1, f2);
                let clipped = match clipped_range {
                    Some(r) => r,
                    None => { return; }
                };
                if clipped[1] - clipped[0] <= TOLERANCE_ABS { return; }
                let (effective_t0, effective_t1) = (clipped[0], clipped[1]);

                if circle.radius <= TOLERANCE_MESH_LEGACY + TOLERANCE_ABS { return; }
                let valid_arc = clipped_range.map(|r| r[1] - r[0]).unwrap_or(0.0);
                if valid_arc <= TOLERANCE_ABS { return; }

                let pcurve_plane = circle_pcurve_on_plane(&circle, plane);
                let pcurve_sphere = fallback_pcurve_by_projection(
                    &Curve3::Circle(circle),
                    &[effective_t0, effective_t1],
                    &Surface3::Sphere(*sphere),
                );
                let (pcurve_on_a, pcurve_on_b) = if plane_is_f1 {
                    (Some(pcurve_plane), Some(pcurve_sphere))
                } else {
                    (Some(pcurve_sphere), Some(pcurve_plane))
                };

                // �?OCCT-aligned: IC endpoints use plane_local_basis (consistent with clip_circle_to_faces)
                //    circle.point_at uses Circle3.normal's any_perpendicular axis,
                //    which may be opposite to plane_local_basis direction, flipping endpoint positions.
                let (u_ax_p, v_ax_p) = crate::inttools::edge_face::plane_local_basis(plane);
                let p_start = circle.center + circle.radius * (effective_t0.cos() * u_ax_p + effective_t0.sin() * v_ax_p);
                let p_end = circle.center + circle.radius * (effective_t1.cos() * u_ax_p + effective_t1.sin() * v_ax_p);
                if p_start.distance_squared(p_end) < TOLERANCE_ABS_SQ { return; }
                // �?OCCT-aligned: try to reuse existing DS vertex (PutPaveOnCurve).
                //    OCCT's IsVertexOnLine detects boundary vertices ON the curve and
                //    places their DS index into the pave block, so the section edge
                //    shares the same TopoDS_Vertex as the boundary edge.  rcad: find
                //    existing vertex within tolerance; only create new if none found.
                //    The tolerance TOLERANCE_ABS*1000 (1e-4) covers intersection noise.
                const IC_VERTEX_MERGE_TOL: f64 = crate::tolerance::TOLERANCE_ABS * 1000.0;
                let v_start = self.ds.find_vertex_near(p_start, IC_VERTEX_MERGE_TOL)
                    .unwrap_or_else(|| self.ds.add_vertex(p_start));
                let v_end = self.ds.find_vertex_near(p_end, IC_VERTEX_MERGE_TOL)
                    .unwrap_or_else(|| self.ds.add_vertex(p_end));
                // OCCT-aligned: inherit tolerance from parent faces (BRep_Tool::Tolerance).
                // OCCT vertices on edges carry source edge/face tolerances (typically 1e-4 to 1e-6);
                // rcad defaults to TOLERANCE_ABS (1e-7) which is too tight for pcurve comparison.
                let parent_tol = self.ds.faces[f1].geom_tol
                    .max(self.ds.faces[f2].geom_tol)
                    .max(self.seam_shift_tol);
                if v_start < self.ds.vertices.len() {
                    self.ds.vertices[v_start].geom_tol = self.ds.vertices[v_start].geom_tol.max(parent_tol);
                }
                if v_end < self.ds.vertices.len() {
                    self.ds.vertices[v_end].geom_tol = self.ds.vertices[v_end].geom_tol.max(parent_tol);
                }

                let curve_idx = self.ds.intersection_curves.len();
                self.ds.intersection_curves.push(IntersectionCurve {
                    curve: Curve3::Circle(circle),
                    polyline: vec![],
                    start_vertex: v_start,
                    end_vertex: v_end,
                    t_range: [effective_t0, effective_t1],
                    pcurve_on_a,
                    pcurve_on_b,
                    geom_tol: crate::tolerance::TOLERANCE_ABS,
                pave_blocks: Vec::new(),
                });

                self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
                self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
                self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                self.ds.faces[f2].face_info.vertices_in.insert(v_end);

                self.ds.interferences.push(Interference::FaceFace {
                    f1,
                    f2,
                    curves: vec![curve_idx],
                    points: vec![],
                });
            }
        }
    }

    // 鈹€鈹€ Sphere �?Sphere analytic face-face intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    fn intersect_sphere_sphere_faces(
        &mut self,
        f1: usize,
        f2: usize,
        sph1: &SphericalSurface,
        sph2: &SphericalSurface,
    ) {
        use inttools::pcurve_derive::fallback_pcurve_by_projection;
        use std::f64::consts::TAU;

        let d_vec = sph2.center - sph1.center;
        let d = d_vec.length();

        // No intersection if disjoint or one contains the other
        if d < TOLERANCE_FLOAT_LOOSE {
            // Concentric spheres: same-domain (same center). Record empty FaceFace.
            self.ds.interferences.push(Interference::FaceFace {
                f1, f2, curves: vec![], points: vec![],
            });
            return;
        }
        if d >= sph1.radius + sph2.radius || d <= (sph1.radius - sph2.radius).abs() {
            return;
        }

        // Distance from sph1 center to the radical plane
        let h = (d * d + sph1.radius * sph1.radius - sph2.radius * sph2.radius) / (2.0 * d);
        let r_circ_sq = sph1.radius * sph1.radius - h * h;
        if r_circ_sq <= 0.0 {
            return; // Tangent or near-tangent
        }
        let r_circ = r_circ_sq.sqrt();

        // Normal of the intersection circle (axis of the radical plane)
        let normal = d_vec.normalize();
        // Center of the intersection circle
        let center = sph1.center + normal * h;

        let circle = Circle3 {
            center,
            normal,
            radius: r_circ,
        };

        let curve3 = Curve3::Circle(circle);
        let t_range = [0.0_f64, TAU];
        let pcurve_a = fallback_pcurve_by_projection(&curve3, &t_range, &Surface3::Sphere(*sph1));
        let pcurve_b = fallback_pcurve_by_projection(&curve3, &t_range, &Surface3::Sphere(*sph2));

        let pts = sample_circle_arc(&circle, 0.0, TAU, 32);
        if pts.len() < 2 {
            return;
        }

        let v_start = self.ds.add_vertex(pts[0]);
        let v_end = self.ds.add_vertex(pts[pts.len() - 1]);

        let curve_idx = self.ds.intersection_curves.len();
        self.ds.intersection_curves.push(IntersectionCurve {
            curve: curve3,
            polyline: vec![],
            start_vertex: v_start,
            end_vertex: v_end,
            t_range: [0.0, TAU],
            pcurve_on_a: Some(pcurve_a),
            pcurve_on_b: Some(pcurve_b),
            geom_tol: crate::tolerance::TOLERANCE_ABS,
        pave_blocks: Vec::new(),
        });

        self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
        self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
        self.ds.faces[f1].face_info.vertices_in.insert(v_start);
        self.ds.faces[f1].face_info.vertices_in.insert(v_end);
        self.ds.faces[f2].face_info.vertices_in.insert(v_start);
        self.ds.faces[f2].face_info.vertices_in.insert(v_end);

        self.ds.interferences.push(Interference::FaceFace {
            f1,
            f2,
            curves: vec![curve_idx],
            points: vec![],
        });
    }

    // 鈹€鈹€ Sphere �?Cylinder analytic face-face intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    fn intersect_sphere_cylinder_faces(
        &mut self,
        f1: usize,
        f2: usize,
        sphere: &SphericalSurface,
        cyl: &CylindricalSurface,
    ) {
        use inttools::pcurve_derive::{
            circle_pcurve_on_cylinder, circle_pcurve_on_sphere, polyline_pcurve_by_projection,
        };
        use inttools::sphere_cylinder::{SphereCylinderResult, intersect_sphere_cylinder};
        use std::f64::consts::TAU;

        // Determine which face is the sphere face (for pcurve_on_a/b ordering)
        let sphere_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Sphere(_));

        let make_pcurves = |pca: Curve2d, pcb: Curve2d| -> (Option<Curve2d>, Option<Curve2d>) {
            if sphere_is_f1 {
                (Some(pca), Some(pcb))
            } else {
                (Some(pcb), Some(pca))
            }
        };

        // Helper: add one intersection circle to the DS and return its index.
        let add_circle =
            |ds: &mut DS,
             circle: &Circle3,
             pcurve_on_a: Option<Curve2d>,
             pcurve_on_b: Option<Curve2d>,
             f1: usize,
             f2: usize|
             -> usize {
                let pts = sample_circle_arc(circle, 0.0, TAU, 32);
                let v_start = ds.add_vertex(pts[0]);
                let v_end = ds.add_vertex(pts[pts.len() - 1]);
                let curve_idx = ds.intersection_curves.len();
                ds.intersection_curves.push(IntersectionCurve {
                    curve: Curve3::Circle(*circle),
                    polyline: vec![],
                    start_vertex: v_start,
                    end_vertex: v_end,
                    t_range: [0.0, TAU],
                    pcurve_on_a,
                    pcurve_on_b,
                    geom_tol: crate::tolerance::TOLERANCE_ABS,
                pave_blocks: Vec::new(),
                });
                ds.faces[f1].face_info.curves_sc.insert(curve_idx);
                ds.faces[f2].face_info.curves_sc.insert(curve_idx);
                ds.faces[f1].face_info.vertices_in.insert(v_start);
                ds.faces[f1].face_info.vertices_in.insert(v_end);
                ds.faces[f2].face_info.vertices_in.insert(v_start);
                ds.faces[f2].face_info.vertices_in.insert(v_end);
                curve_idx
            };

        // Closure to compute pcurves for one intersection circle.
        // The intersection circle is always a latitude line on the sphere
        // (�?= acos((h �?h_c) / R)), so `circle_pcurve_on_sphere` is exact
        // here regardless of whether the sphere and cylinder axes are parallel.
        let make_circle_pcurves = |circle: &Circle3| -> (Option<Curve2d>, Option<Curve2d>) {
            let pcurve_sph = circle_pcurve_on_sphere(circle, sphere);
            let pcurve_cyl = circle_pcurve_on_cylinder(circle, cyl);
            make_pcurves(pcurve_sph, pcurve_cyl)
        };

        match intersect_sphere_cylinder(sphere, cyl) {
            SphereCylinderResult::NoIntersection => (),
            SphereCylinderResult::General => {
                // Fall back to numeric marching for the quartic case.
                self.intersect_ff_by_marching(f1, f2);
            }
            SphereCylinderResult::TangentCircle(circle) => {
                let (pca, pcb) = make_circle_pcurves(&circle);
                let ci = add_circle(self.ds, &circle, pca, pcb, f1, f2);
                self.ds.interferences.push(Interference::FaceFace {
                    f1,
                    f2,
                    curves: vec![ci],
                    points: vec![],
                });
            }
            SphereCylinderResult::TwoCircles(c1, c2) => {
                let (pca1, pcb1) = make_circle_pcurves(&c1);
                let ci1 = add_circle(self.ds, &c1, pca1, pcb1, f1, f2);
                let (pca2, pcb2) = make_circle_pcurves(&c2);
                let ci2 = add_circle(self.ds, &c2, pca2, pcb2, f1, f2);
                self.ds.interferences.push(Interference::FaceFace {
                    f1,
                    f2,
                    curves: vec![ci1, ci2],
                    points: vec![],
                });
            }

            SphereCylinderResult::SkewQuartic(branches) => {
                let s1 = Surface3::Sphere(*sphere);
                let s2 = Surface3::Cylinder(*cyl);
                let mut curve_indices = Vec::new();
                for branch in branches {
                    if branch.len() < 2 {
                        continue;
                    }
                    let v_start = self.ds.add_vertex(branch[0]);
                    let v_end = self.ds.add_vertex(branch[branch.len() - 1]);
                    let dir = (branch[branch.len() - 1] - branch[0])
                        .normalize_or_zero();
                    let ci = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Line(Line3 {
                            origin: branch[0],
                            direction: if dir.length_squared() > 0.5 {
                                dir
                            } else {
                                DVec3::X
                            },
                        }),
                        polyline: branch.clone(),
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, 1.0],
                        pcurve_on_a: polyline_pcurve_by_projection(&branch, &s1),
                        pcurve_on_b: polyline_pcurve_by_projection(&branch, &s2),
                        geom_tol: crate::tolerance::TOLERANCE_ABS,
                    pave_blocks: Vec::new(),
                    });
                    curve_indices.push(ci);
                    self.ds.faces[f1].face_info.curves_sc.insert(ci);
                    self.ds.faces[f2].face_info.curves_sc.insert(ci);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);
                }
                if !curve_indices.is_empty() {
                    self.ds.interferences.push(Interference::FaceFace {
                        f1,
                        f2,
                        curves: curve_indices,
                        points: vec![],
                    });
                }
            }
        }
    }

    // 鈹€鈹€ Cylinder �?Cylinder analytic face-face intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    fn intersect_cylinder_cylinder_faces(
        &mut self,
        f1: usize,
        f2: usize,
        cyl1: &CylindricalSurface,
        cyl2: &CylindricalSurface,
    ) {
        use inttools::cylinder_cylinder::{CylinderCylinderResult, intersect_cylinder_cylinder};
        use inttools::pcurve_derive::{
            circle_pcurve_on_cylinder, ellipse_pcurve_on_cylinder, line_pcurve_on_cylinder, polyline_pcurve_by_projection,
        };
        use std::f64::consts::TAU;

        // Determine which face is cyl1 (for pcurve_on_a/b ordering)
        let cyl1_is_f1 = {
            if let Surface3::Cylinder(c) = &self.ds.faces[f1].surface {
                (c.origin - cyl1.origin).length_squared() < TOLERANCE_LINEAR_ULTRA_STRICT * TOLERANCE_LINEAR_ULTRA_STRICT
                    && (c.axis - cyl1.axis).length_squared() < TOLERANCE_LINEAR_ULTRA_STRICT * TOLERANCE_LINEAR_ULTRA_STRICT
            } else {
                false
            }
        };

        let make_pcurves = |pca: Curve2d, pcb: Curve2d| -> (Option<Curve2d>, Option<Curve2d>) {
            if cyl1_is_f1 {
                (Some(pca), Some(pcb))
            } else {
                (Some(pcb), Some(pca))
            }
        };

        // Helper: push a circle intersection curve and register it with both faces.
        let add_circle =
            |ds: &mut DS,
             circle: &Circle3,
             pcurve_on_a: Option<Curve2d>,
             pcurve_on_b: Option<Curve2d>,
             f1: usize,
             f2: usize|
             -> usize {
                let pts = sample_circle_arc(circle, 0.0, TAU, 32);
                let v_start = ds.add_vertex(pts[0]);
                let v_end = ds.add_vertex(pts[pts.len() - 1]);
                let ci = ds.intersection_curves.len();
                ds.intersection_curves.push(IntersectionCurve {
                    curve: Curve3::Circle(*circle),
                    polyline: vec![],
                    start_vertex: v_start,
                    end_vertex: v_end,
                    t_range: [0.0, TAU],
                    pcurve_on_a,
                    pcurve_on_b,
                    geom_tol: crate::tolerance::TOLERANCE_ABS,
                pave_blocks: Vec::new(),
                });
                ds.faces[f1].face_info.curves_sc.insert(ci);
                ds.faces[f2].face_info.curves_sc.insert(ci);
                ds.faces[f1].face_info.vertices_in.insert(v_start);
                ds.faces[f1].face_info.vertices_in.insert(v_end);
                ds.faces[f2].face_info.vertices_in.insert(v_start);
                ds.faces[f2].face_info.vertices_in.insert(v_end);
                ci
            };

        // Helper: push a line generator intersection and register it.
        let add_line = |ds: &mut DS,
                        line: &Line3,
                        t_range: [f64; 2],
                        pcurve_on_a: Option<Curve2d>,
                        pcurve_on_b: Option<Curve2d>,
                        f1: usize,
                        f2: usize|
         -> usize {
            use rcad_kernel::CurveEval;
            let v_start = ds.add_vertex(Curve3::Line(*line).point_at(t_range[0]));
            let v_end = ds.add_vertex(Curve3::Line(*line).point_at(t_range[1]));
            let ci = ds.intersection_curves.len();
            ds.intersection_curves.push(IntersectionCurve {
                curve: Curve3::Line(*line),
                polyline: vec![],
                start_vertex: v_start,
                end_vertex: v_end,
                t_range,
                pcurve_on_a,
                pcurve_on_b,
                geom_tol: crate::tolerance::TOLERANCE_ABS,
            pave_blocks: Vec::new(),
            });
            ds.faces[f1].face_info.curves_sc.insert(ci);
            ds.faces[f2].face_info.curves_sc.insert(ci);
            ds.faces[f1].face_info.vertices_in.insert(v_start);
            ds.faces[f1].face_info.vertices_in.insert(v_end);
            ds.faces[f2].face_info.vertices_in.insert(v_start);
            ds.faces[f2].face_info.vertices_in.insert(v_end);
            ci
        };

        // Helper: push an ellipse intersection and register it.
        let add_ellipse = |ds: &mut DS,
                           ellipse: &Ellipse3,
                           pcurve_on_a: Option<Curve2d>,
                           pcurve_on_b: Option<Curve2d>,
                           f1: usize,
                           f2: usize|
         -> usize {
            let pts = sample_circle_arc(
                &Circle3 {
                    center: ellipse.center,
                    normal: ellipse.normal,
                    radius: ellipse.major_radius.max(ellipse.minor_radius),
                },
                0.0,
                TAU,
                32,
            );
            let v_start = ds.add_vertex(pts[0]);
            let v_end = ds.add_vertex(pts[pts.len() - 1]);
            let ci = ds.intersection_curves.len();
            ds.intersection_curves.push(IntersectionCurve {
                curve: Curve3::Ellipse(*ellipse),
                polyline: vec![],
                start_vertex: v_start,
                end_vertex: v_end,
                t_range: [0.0, TAU],
                pcurve_on_a,
                pcurve_on_b,
                geom_tol: crate::tolerance::TOLERANCE_ABS,
            pave_blocks: Vec::new(),
            });
            ds.faces[f1].face_info.curves_sc.insert(ci);
            ds.faces[f2].face_info.curves_sc.insert(ci);
            ds.faces[f1].face_info.vertices_in.insert(v_start);
            ds.faces[f1].face_info.vertices_in.insert(v_end);
            ds.faces[f2].face_info.vertices_in.insert(v_start);
            ds.faces[f2].face_info.vertices_in.insert(v_end);
            ci
        };

        let extent = 20.0_f64;
        let mut curve_indices = Vec::new();

        match intersect_cylinder_cylinder(cyl1, cyl2) {
            CylinderCylinderResult::NoIntersection => return,
            CylinderCylinderResult::Coaxial => {
                // Same-domain coaxial cylinders: record empty-curves FaceFace so
                // the Builder treats this pair as coincident (no intersection to split).
                self.ds.interferences.push(Interference::FaceFace {
                    f1, f2, curves: vec![], points: vec![],
                });
                return;
            }

            CylinderCylinderResult::PerpendicularOffsetCurves {
                cyl1: off_cyl1,
                cyl2: off_cyl2,
                ..
            } => {
                // Perpendicular cylinders with offset (non-intersecting) axes.
                // Parametrization on cyl1's surface:
                //   P(�? = O1 + v(�?*a1 + R1*(cos(�?*U1 + sin(�?*V1)
                //   v(�? = dz �?�?R2�?- (R1路cos(�? - dx)�?
                //
                // Two closed-loop intersection curves per face, one per �?interval:
                //   Loop 1 (胃鈭圼t_low, t_high]): forward branch+  back branch-
                //   Loop 2 (胃鈭圼蟿-t_high, �?t_low]): forward branch+  back branch-
                // Each loop is a single IntersectionCurve whose start/end vertex
                // coincide (same 3D tangent point) �?the boolean builder sees a
                // single closed trim boundary per loop.
                let a1 = off_cyl1.axis.normalize();
                let a2 = off_cyl2.axis.normalize();
                let r1 = off_cyl1.radius;
                let r2 = off_cyl2.radius;
                let r2_sq = r2 * r2;
                let w = off_cyl1.origin - off_cyl2.origin;
                let denom = 1.0 - a1.dot(a2) * a1.dot(a2);
                if denom.abs() < 1e-12 { return; }
                let d1 = a1.dot(w); let d2 = a2.dot(w);
                let t = (a1.dot(a2) * d2 - d1) / denom;
                let s = (d2 - a1.dot(a2) * d1) / denom;
                let conn = (off_cyl1.origin + a1 * t) - (off_cyl2.origin + a2 * s);
                let conn_len = conn.length();
                let u1 = if conn_len < 1e-12 {
                    // Axes intersect (zero or near-zero offset along the
                    // connecting vector).  Pick any direction perpendicular to
                    // a1 �?a1 �?a2 works since the axes are perpendicular.
                    a1.cross(a2).normalize()
                } else {
                    conn / conn_len
                };
                let v1 = a1.cross(u1).normalize();
                let delta = off_cyl2.origin - off_cyl1.origin;
                let dx = delta.dot(u1);
                let dz = delta.dot(a1);
                let cos_min = ((dx - r2) / r1).clamp(-1.0, 1.0);
                let cos_max = ((dx + r2) / r1).clamp(-1.0, 1.0);
                if cos_min > cos_max { return; }
                let t_low = cos_max.acos();
                let t_high = cos_min.acos();

                let surface1 = Surface3::Cylinder(off_cyl1);
                let surface2 = Surface3::Cylinder(off_cyl2);
                let n_per = 9;

                for (t_start, t_end) in [(t_low, t_high), (TAU - t_high, TAU - t_low)] {
                    let n_pts = n_per * 2 + 1; // forward + backward (share the turn-around point)
                    let mut pts: Vec<DVec3> = Vec::with_capacity(n_pts);

                    // Forward: branch = +1, �?= t_start �?t_end
                    for i in 0..=n_per {
                        let theta = t_start + (t_end - t_start) * i as f64 / n_per as f64;
                        let (ct, st) = (theta.cos(), theta.sin());
                        let diff = r1 * ct - dx;
                        let disc = (r2_sq - diff * diff).max(0.0).sqrt();
                        let v_z = dz + disc; // branch sign +1
                        pts.push(off_cyl1.origin + v_z * a1 + r1 * (ct * u1 + st * v1));
                    }
                    // Backward: branch = -1, �?= t_end �?t_start (reversed)
                    for i in 1..=n_per {
                        let theta = t_end - (t_end - t_start) * i as f64 / n_per as f64;
                        let (ct, st) = (theta.cos(), theta.sin());
                        let diff = r1 * ct - dx;
                        let disc = (r2_sq - diff * diff).max(0.0).sqrt();
                        let v_z = dz - disc; // branch sign -1
                        pts.push(off_cyl1.origin + v_z * a1 + r1 * (ct * u1 + st * v1));
                    }
                    if pts.len() < 2 { continue; }

                    let pca = polyline_pcurve_by_projection(&pts, &surface1);
                    let pcb = polyline_pcurve_by_projection(&pts, &surface2);
                    let (pca, pcb) = match (pca, pcb) {
                        (Some(a), Some(b)) => make_pcurves(a, b),
                        _ => continue,
                    };

                    // add_vertex dedup: pts[0] and pts[pts.len()-1] are the same
                    // tangent point (�?t_start, disc=0, both branches coincide).
                    let v_start = self.ds.add_vertex(pts[0]);
                    let v_end = self.ds.add_vertex(pts[pts.len() - 1]);
                    let ci = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Line(Line3 {
                            origin: pts[0],
                            direction: (pts[pts.len() - 1] - pts[0]).normalize_or(DVec3::X),
                        }),
                        polyline: pts, start_vertex: v_start, end_vertex: v_end,
                        t_range: [0.0, 1.0], pcurve_on_a: pca, pcurve_on_b: pcb,
                        geom_tol: crate::tolerance::TOLERANCE_ABS,
                    pave_blocks: Vec::new(),
                    });
                    self.ds.faces[f1].face_info.curves_sc.insert(ci);
                    self.ds.faces[f2].face_info.curves_sc.insert(ci);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);
                    curve_indices.push(ci);
                }
            }

            CylinderCylinderResult::General => {
                // Fall back to numeric marching for skew/oblique axes.
                self.intersect_ff_by_marching(f1, f2);
                return;
            }

            CylinderCylinderResult::OneGeneratorLine(line) => {
                let pca = line_pcurve_on_cylinder(&line, cyl1);
                let pcb = line_pcurve_on_cylinder(&line, cyl2);
                let (pca, pcb) = make_pcurves(pca, pcb);
                let ci = add_line(self.ds, &line, [-extent, extent], pca, pcb, f1, f2);
                curve_indices.push(ci);
            }

            CylinderCylinderResult::TwoGeneratorLines(l1, l2) => {
                let pca1 = line_pcurve_on_cylinder(&l1, cyl1);
                let pcb1 = line_pcurve_on_cylinder(&l1, cyl2);
                let (pca1, pcb1) = make_pcurves(pca1, pcb1);
                let ci1 = add_line(self.ds, &l1, [-extent, extent], pca1, pcb1, f1, f2);

                let pca2 = line_pcurve_on_cylinder(&l2, cyl1);
                let pcb2 = line_pcurve_on_cylinder(&l2, cyl2);
                let (pca2, pcb2) = make_pcurves(pca2, pcb2);
                let ci2 = add_line(self.ds, &l2, [-extent, extent], pca2, pcb2, f1, f2);

                curve_indices.push(ci1);
                curve_indices.push(ci2);
            }

            CylinderCylinderResult::TwoCircles(c1, c2) => {
                // Perpendicular Steinmetz equal-radii: circles in diagonal planes.
                let pca1 = circle_pcurve_on_cylinder(&c1, cyl1);
                let pcb1 = circle_pcurve_on_cylinder(&c1, cyl2);
                let (pca1, pcb1) = make_pcurves(pca1, pcb1);
                let ci1 = add_circle(self.ds, &c1, pca1, pcb1, f1, f2);

                let pca2 = circle_pcurve_on_cylinder(&c2, cyl1);
                let pcb2 = circle_pcurve_on_cylinder(&c2, cyl2);
                let (pca2, pcb2) = make_pcurves(pca2, pcb2);
                let ci2 = add_circle(self.ds, &c2, pca2, pcb2, f1, f2);

                curve_indices.push(ci1);
                curve_indices.push(ci2);
            }

            CylinderCylinderResult::TwoEllipses(e1, e2) => {
                // Perpendicular Steinmetz unequal-radii.
                // Use analytic pcurves on both cylinders (no iterative projection).
                let pca1 = ellipse_pcurve_on_cylinder(&e1, cyl1);
                let pcb1 = ellipse_pcurve_on_cylinder(&e1, cyl2);
                let (pca1, pcb1) = make_pcurves(pca1, pcb1);
                let ci1 = add_ellipse(self.ds, &e1, pca1, pcb1, f1, f2);

                let pca2 = ellipse_pcurve_on_cylinder(&e2, cyl1);
                let pcb2 = ellipse_pcurve_on_cylinder(&e2, cyl2);
                let (pca2, pcb2) = make_pcurves(pca2, pcb2);
                let ci2 = add_ellipse(self.ds, &e2, pca2, pcb2, f1, f2);

                curve_indices.push(ci1);
                curve_indices.push(ci2);
            }

            CylinderCylinderResult::SkewQuartic(branches) => {
                let s1 = Surface3::Cylinder(*cyl1);
                let s2 = Surface3::Cylinder(*cyl2);
                for branch in branches {
                    if branch.len() < 2 {
                        continue;
                    }
                    let v_start = self.ds.add_vertex(branch[0]);
                    let v_end = self.ds.add_vertex(branch[branch.len() - 1]);
                    let dir = (branch[branch.len() - 1] - branch[0])
                        .normalize_or_zero();
                    let ci = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Line(Line3 {
                            origin: branch[0],
                            direction: if dir.length_squared() > 0.5 {
                                dir
                            } else {
                                DVec3::X
                            },
                        }),
                        polyline: branch.clone(),
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, 1.0],
                        pcurve_on_a: polyline_pcurve_by_projection(&branch, &s1),
                        pcurve_on_b: polyline_pcurve_by_projection(&branch, &s2),
                        geom_tol: crate::tolerance::TOLERANCE_ABS,
                    pave_blocks: Vec::new(),
                    });
                    curve_indices.push(ci);
                    self.ds.faces[f1].face_info.curves_sc.insert(ci);
                    self.ds.faces[f2].face_info.curves_sc.insert(ci);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);
                }
            }
        }

        if !curve_indices.is_empty() {
            self.ds.interferences.push(Interference::FaceFace {
                f1,
                f2,
                curves: curve_indices,
                points: vec![],
            });
        }
    }

    fn cylinder_face_v_range(&self, face_idx: usize, cyl: &CylindricalSurface) -> [f64; 2] {
        let axis = cyl.axis.normalize();
        let mut v_min = f64::INFINITY;
        let mut v_max = f64::NEG_INFINITY;
        for &ei in &self.ds.faces[face_idx].boundary_edges {
            if let Some(edge) = self.ds.edges.get(ei) {
                if let Some(v) = self.ds.vertices.get(edge.start_vertex) {
                    let proj = (v.point - cyl.origin).dot(axis);
                    v_min = v_min.min(proj);
                    v_max = v_max.max(proj);
                }
                if let Some(v) = self.ds.vertices.get(edge.end_vertex) {
                    let proj = (v.point - cyl.origin).dot(axis);
                    v_min = v_min.min(proj);
                    v_max = v_max.max(proj);
                }
            }
        }
        let r = if v_min.is_infinite() {
            // Fallback: a generous default range
            [-20.0, 20.0]
        } else {
            [v_min, v_max]
        };
        r
    }

    // 鈹€鈹€ Plane �?Cylinder analytic face-face intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    fn intersect_plane_cylinder_faces(
        &mut self,
        f1: usize,
        f2: usize,
        plane: &Plane,
        cyl: &CylindricalSurface,
    ) {
        use inttools::pcurve_derive::{
            circle_pcurve_on_cylinder, circle_pcurve_on_plane, ellipse_pcurve_on_cylinder,
            ellipse_pcurve_on_plane, line_pcurve_on_cylinder,
            line_pcurve_on_plane,
        };
        use inttools::plane_cylinder::{PlaneCylinderResult, intersect_plane_cylinder};
        use rcad_kernel::CurveEval;
        use std::f64::consts::TAU;

        let result = intersect_plane_cylinder(plane, cyl);

        // Determine which face carries the plane (for correct pcurve_on_a/b assignment)
        let plane_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Plane(_));

        let make_pcurves = |pca: Curve2d, pcb: Curve2d| -> (Option<Curve2d>, Option<Curve2d>) {
            if plane_is_f1 {
                (Some(pca), Some(pcb))
            } else {
                (Some(pcb), Some(pca))
            }
        };

        let add_curve = |ds: &mut DS,
                         curve: Curve3,
                         t_range: [f64; 2],
                         pcurve_on_a: Option<Curve2d>,
                         pcurve_on_b: Option<Curve2d>,
                         f1: usize,
                         f2: usize|
         -> usize {
            let p_start = curve.point_at(t_range[0]);
            let p_end = curve.point_at(t_range[1]);
            let v_start = ds.add_vertex(p_start);
            let v_end = ds.add_vertex(p_end);
            let curve_idx = ds.intersection_curves.len();
            ds.intersection_curves.push(IntersectionCurve {
                curve,
                polyline: vec![],
                start_vertex: v_start,
                end_vertex: v_end,
                t_range,
                pcurve_on_a,
                pcurve_on_b,
                geom_tol: crate::tolerance::TOLERANCE_ABS,
            pave_blocks: Vec::new(),
            });
            ds.faces[f1].face_info.curves_sc.insert(curve_idx);
            ds.faces[f2].face_info.curves_sc.insert(curve_idx);
            ds.faces[f1].face_info.vertices_in.insert(v_start);
            ds.faces[f1].face_info.vertices_in.insert(v_end);
            ds.faces[f2].face_info.vertices_in.insert(v_start);
            ds.faces[f2].face_info.vertices_in.insert(v_end);
            // �?OCCT-aligned:Propagate IC vertices to all faces sharing boundary edges
            //    (BOPDS_FaceInfo::AppendBlock equivalent).
            propagate_ic_vertices_to_shared_faces(ds, &[v_start, v_end], &[f1, f2]);
            curve_idx
        };

        let mut curve_indices = Vec::new();

        match result {
            PlaneCylinderResult::NoIntersection => return,
            PlaneCylinderResult::TangentLine(line) => {
                // �?OCCT aligned: tangent lines are also valid intersection curves,
                //    used to split the cylinder face. OCCT IntTools_FaceFace::MakeCurve
                //    creates BRep edges for tangent lines too.
                // Clip to the cylinder face's parametric V range along the axis
                // so the intersection curve doesn't extend beyond the actual face.
                let cyl_fi = if plane_is_f1 { f2 } else { f1 };
                let v_range = self.cylinder_face_v_range(cyl_fi, cyl);
                let (pca, pcb) = make_pcurves(
                    line_pcurve_on_plane(&line, plane),
                    line_pcurve_on_cylinder(&line, cyl),
                );
                let ci = add_curve(self.ds, Curve3::Line(line), v_range, pca, pcb, f1, f2);
                curve_indices.push(ci);
            }
            PlaneCylinderResult::TwoLines(l1, l2) => {
                // Clip each line to the cylinder face's parametric V range
                let cyl_fi = if plane_is_f1 { f2 } else { f1 };
                let v_range = self.cylinder_face_v_range(cyl_fi, cyl);
                let (pca1, pcb1) = make_pcurves(
                    line_pcurve_on_plane(&l1, plane),
                    line_pcurve_on_cylinder(&l1, cyl),
                );
                let ci1 = add_curve(
                    self.ds,
                    Curve3::Line(l1),
                    v_range,
                    pca1,
                    pcb1,
                    f1,
                    f2,
                );
                let (pca2, pcb2) = make_pcurves(
                    line_pcurve_on_plane(&l2, plane),
                    line_pcurve_on_cylinder(&l2, cyl),
                );
                let ci2 = add_curve(
                    self.ds,
                    Curve3::Line(l2),
                    v_range,
                    pca2,
                    pcb2,
                    f1,
                    f2,
                );
                curve_indices.push(ci1);
                curve_indices.push(ci2);
            }
            PlaneCylinderResult::Circle(circle) => {
                let (pca, pcb) = make_pcurves(
                    circle_pcurve_on_plane(&circle, plane),
                    circle_pcurve_on_cylinder(&circle, cyl),
                );
                let ci = add_curve(
                    self.ds,
                    Curve3::Circle(circle),
                    [0.0, TAU],
                    pca,
                    pcb,
                    f1,
                    f2,
                );
                curve_indices.push(ci);
                // �?OCCT-aligned: when circle lies on cylinder V-boundary, remove IC from cylinder face.
                //    OCCT PerformLoops does not create sub-faces outside face domain, but rcad's BooleanBuilder does.
                //    Plane face retains the IC (split arc segments) for correct box face splitting.
                let cyl_fi = if plane_is_f1 { f2 } else { f1 };
                let v_range = self.cylinder_face_v_range(cyl_fi, cyl);
                let v = (circle.center - cyl.origin).dot(cyl.axis.normalize());
                let boundary_tol = TOLERANCE_ABS * 1000.0;
                if (v - v_range[0]).abs() < boundary_tol
                    || (v - v_range[1]).abs() < boundary_tol
                {
                    self.ds.faces[cyl_fi].face_info.curves_sc.remove(&ci);
                }
            }
            PlaneCylinderResult::Ellipse(ellipse) => {
                let pca_plane = ellipse_pcurve_on_plane(&ellipse, plane);
                let pcb_cyl = ellipse_pcurve_on_cylinder(&ellipse, cyl);
                let (pca, pcb) = make_pcurves(pca_plane, pcb_cyl);
                let ci = add_curve(
                    self.ds,
                    Curve3::Ellipse(ellipse),
                    [0.0, TAU],
                    pca,
                    pcb,
                    f1,
                    f2,
                );
                curve_indices.push(ci);
            }
        }

        if !curve_indices.is_empty() {
            self.ds.interferences.push(Interference::FaceFace {
                f1,
                f2,
                curves: curve_indices,
                points: vec![],
            });
        }
    }

    fn trim_curve_to_faces(
        ds: &DS,
        curve: &Curve3,
        search_range: [f64; 2],
        f1: usize,
        f2: usize,
    ) -> Option<[f64; 2]> {
        use crate::medial_axis::point_in_polygon_2d;
        use rcad_kernel::projection::closest_point_on_surface;
        use std::f64::consts::TAU;

        const N: usize = 256;

        let face1 = &ds.faces[f1];
        let face2 = &ds.faces[f2];
        let uv_bnd1 = face1.uv_boundary.as_ref()?;
        let uv_bnd2 = face2.uv_boundary.as_ref()?;
        let s1 = &face1.surface;
        let s2 = &face2.surface;

        // UV from 3D point on a surface, normalising u �?[0, 2蟺].
        let uv_on_surface = |surface: &Surface3, p: DVec3| -> DVec2 {
            match surface {
                Surface3::Cone(cone) => {
                    let uv = cone.world_to_uv(p);
                    DVec2::new(if uv.x < 0.0 { uv.x + TAU } else { uv.x }, uv.y)
                }
                Surface3::Sphere(sph) => sph.world_to_uv(p),
                Surface3::Cylinder(cyl) => {
                    let x_ax = cyl.ref_dir.normalize();
                    let y_ax = cyl.axis.cross(x_ax).normalize();
                    let local = p - cyl.origin;
                    let u = local.dot(y_ax).atan2(local.dot(x_ax));
                    DVec2::new(if u < 0.0 { u + TAU } else { u }, local.dot(cyl.axis))
                }
                _ => {
                    let proj = closest_point_on_surface(surface, p, 16);
                    DVec2::new(proj.params.0, proj.params.1)
                }
            }
        };

        // True when the curve point at t is inside *both* faces' UV boundaries.
        // For planar faces the 3D point must actually lie on the plane (not just
        // project there), else an off-surface point would be a false positive.
        let point_in_both = |t: f64| -> bool {
            let pt = curve.point_at(t);
            for (sf, bnd) in &[(s1, uv_bnd1), (s2, uv_bnd2)] {
                if let Surface3::Plane(pl) = sf {
                    if (pt - pl.origin).dot(pl.normal).abs() > TOLERANCE_COORD_SUB {
                        return false;
                    }
                }
                let uv = uv_on_surface(sf, pt);
                if !point_in_polygon_2d(uv, bnd) {
                    return false;
                }
            }
            true
        };

        let [t0, t1] = search_range;
        let step = (t1 - t0) / N as f64;
        let mut seg_start: Option<(usize, f64)> = None;
        let mut segments: Vec<(usize, usize, f64, f64)> = Vec::new();

        for i in 0..=N {
            let t = t0 + step * i as f64;
            let inside = point_in_both(t);

            if inside {
                if seg_start.is_none() {
                    seg_start = Some((i, t));
                }
            } else if let Some((si, st)) = seg_start.take() {
                if t - st > TOLERANCE_LINEAR_ULTRA_STRICT {
                    segments.push((si, i, st, t));
                }
            }
        }
        if let Some((si, st)) = seg_start.take() {
            if t1 - st > TOLERANCE_LINEAR_ULTRA_STRICT {
                segments.push((si, N, st, t1));
            }
        }

        // Longest segment
        let (si, ei, rough_start, rough_end) = segments.into_iter().max_by(|a, b| {
            (a.3 - a.2)
                .partial_cmp(&(b.3 - b.2))
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;

        // 鈹€鈹€ binary-search refinement of both endpoints 鈹€鈹€
        // Start: between sample (si-1, outside) and sample (si, inside).
        let refined_start = if si > 0 {
            let t_out = t0 + step * (si - 1) as f64;
            let mut lo = t_out;   // outside
            let mut hi = rough_start; // inside
            for _ in 0..48 {
                let mid = 0.5 * (lo + hi);
                if point_in_both(mid) {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            hi
        } else {
            rough_start
        };

        // End: between sample (ei-1, inside) and sample (ei, outside).
        let refined_end = if ei < N {
            let t_in = t0 + step * (ei - 1) as f64;
            let t_out = t0 + step * ei as f64;
            let mut lo = t_in;  // inside
            let mut hi = t_out; // outside
            for _ in 0..48 {
                let mid = 0.5 * (lo + hi);
                if point_in_both(mid) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            lo
        } else {
            rough_end
        };

        if refined_end - refined_start > TOLERANCE_LINEAR_ULTRA_STRICT {
            Some([refined_start, refined_end])
        } else {
            None
        }
    }

    // 鈹€鈹€ Plane �?Cone analytic face-face intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    fn intersect_plane_cone_faces(
        &mut self,
        f1: usize,
        f2: usize,
        plane: &Plane,
        cone: &ConicalSurface,
    ) {
        use inttools::pcurve_derive::{
            circle_pcurve_on_cone, circle_pcurve_on_plane, ellipse_pcurve_on_cone,
            ellipse_pcurve_on_plane, fallback_pcurve_by_projection,
            line_pcurve_on_cone, line_pcurve_on_plane, sampled_pcurve_on_cone,
        };
        use inttools::plane_cone::{PlaneConicalResult, intersect_plane_cone};
        use std::f64::consts::TAU;

        // Determine which face carries the plane
        let plane_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Plane(_));

        let make_pcurves = |pca: Curve2d, pcb: Curve2d| -> (Option<Curve2d>, Option<Curve2d>) {
            if plane_is_f1 {
                (Some(pca), Some(pcb))
            } else {
                (Some(pcb), Some(pca))
            }
        };

        // Helper: push a generic curve and register it with both faces.
        let add_curve = |ds: &mut DS,
                         curve: Curve3,
                         t_range: [f64; 2],
                         pcurve_on_a: Option<Curve2d>,
                         pcurve_on_b: Option<Curve2d>,
                         f1: usize,
                         f2: usize|
         -> usize {
            let p_start = curve.point_at(t_range[0]);
            let p_end = curve.point_at(t_range[1]);
            let v_start = ds.add_vertex(p_start);
            let v_end = ds.add_vertex(p_end);
            let ci = ds.intersection_curves.len();
            ds.intersection_curves.push(IntersectionCurve {
                curve,
                polyline: vec![],
                start_vertex: v_start,
                end_vertex: v_end,
                t_range,
                pcurve_on_a,
                pcurve_on_b,
                geom_tol: crate::tolerance::TOLERANCE_ABS,
            pave_blocks: Vec::new(),
            });
            ds.faces[f1].face_info.curves_sc.insert(ci);
            ds.faces[f2].face_info.curves_sc.insert(ci);
            ds.faces[f1].face_info.vertices_in.insert(v_start);
            ds.faces[f1].face_info.vertices_in.insert(v_end);
            ds.faces[f2].face_info.vertices_in.insert(v_start);
            ds.faces[f2].face_info.vertices_in.insert(v_end);
            // �?OCCT-aligned:Propagate IC vertices to all faces sharing boundary edges
            //    (BOPDS_FaceInfo::AppendBlock equivalent).
            propagate_ic_vertices_to_shared_faces(ds, &[v_start, v_end], &[f1, f2]);
            ci
        };

        let mut curve_indices = Vec::new();

        match intersect_plane_cone(plane, cone) {
            PlaneConicalResult::NoIntersection | PlaneConicalResult::Point(_) => {
                return;
            }

            PlaneConicalResult::SingleLine(line) => {
                if let Some(trimmed) = Self::trim_curve_to_faces(
                    self.ds,
                    &Curve3::Line(line),
                    [-30.0, 30.0],
                    f1,
                    f2,
                ) {
                    let pca_plane = line_pcurve_on_plane(&line, plane);
                    let pcb_cone = line_pcurve_on_cone(&line, cone);
                    let (pca, pcb) = make_pcurves(pca_plane, pcb_cone);
                    let ci = add_curve(
                        self.ds,
                        Curve3::Line(line),
                        trimmed,
                        pca,
                        pcb,
                        f1,
                        f2,
                    );
                    curve_indices.push(ci);
                }
            }

            PlaneConicalResult::TwoLines(l1, l2) => {
                if let Some(t1) = Self::trim_curve_to_faces(
                    self.ds,
                    &Curve3::Line(l1),
                    [-30.0, 30.0],
                    f1,
                    f2,
                ) {
                    let pca1 = line_pcurve_on_plane(&l1, plane);
                    let pcb1 = line_pcurve_on_cone(&l1, cone);
                    let (pca1, pcb1) = make_pcurves(pca1, pcb1);
                    let ci1 = add_curve(
                        self.ds,
                        Curve3::Line(l1),
                        t1,
                        pca1,
                        pcb1,
                        f1,
                        f2,
                    );
                    curve_indices.push(ci1);
                }

                if let Some(t2) = Self::trim_curve_to_faces(
                    self.ds,
                    &Curve3::Line(l2),
                    [-30.0, 30.0],
                    f1,
                    f2,
                ) {
                    let pca2 = line_pcurve_on_plane(&l2, plane);
                    let pcb2 = line_pcurve_on_cone(&l2, cone);
                    let (pca2, pcb2) = make_pcurves(pca2, pcb2);
                    let ci2 = add_curve(
                        self.ds,
                        Curve3::Line(l2),
                        t2,
                        pca2,
                        pcb2,
                        f1,
                        f2,
                    );
                    curve_indices.push(ci2);
                }
            }

            PlaneConicalResult::Circle(circle) => {
                let pca_plane = circle_pcurve_on_plane(&circle, plane);
                let pcb_cone = circle_pcurve_on_cone(&circle, cone);
                let (pca, pcb) = make_pcurves(pca_plane, pcb_cone);
                let ci = add_curve(self.ds, Curve3::Circle(circle), [0.0, TAU], pca, pcb, f1, f2);
                curve_indices.push(ci);
            }

            PlaneConicalResult::Ellipse(ellipse) => {
                let pca_plane = ellipse_pcurve_on_plane(&ellipse, plane);
                let pcb_cone = ellipse_pcurve_on_cone(&ellipse, cone);
                let (pca, pcb) = make_pcurves(pca_plane, pcb_cone);
                let ci = add_curve(self.ds, Curve3::Ellipse(ellipse), [0.0, TAU], pca, pcb, f1, f2);
                curve_indices.push(ci);
            }

            PlaneConicalResult::Parabola(parabola) => {
                if let Some(trimmed) = Self::trim_curve_to_faces(
                    self.ds,
                    &Curve3::Parabola(parabola),
                    [-30.0, 30.0],
                    f1,
                    f2,
                ) {
                    let pca_plane = fallback_pcurve_by_projection(
                        &Curve3::Parabola(parabola),
                        &trimmed,
                        &Surface3::Plane(*plane),
                    );
                    let pcb_cone = sampled_pcurve_on_cone(
                        &Curve3::Parabola(parabola),
                        &trimmed,
                        cone,
                    );
                    let (pca, pcb) = make_pcurves(pca_plane, pcb_cone);
                    let ci = add_curve(
                        self.ds,
                        Curve3::Parabola(parabola),
                        trimmed,
                        pca,
                        pcb,
                        f1,
                        f2,
                    );
                    curve_indices.push(ci);
                }
            }

            PlaneConicalResult::Hyperbola(hyperbola) => {
                if let Some(trimmed) = Self::trim_curve_to_faces(
                    self.ds,
                    &Curve3::Hyperbola(hyperbola),
                    [-30.0, 30.0],
                    f1,
                    f2,
                ) {
                    let pca_plane = fallback_pcurve_by_projection(
                        &Curve3::Hyperbola(hyperbola),
                        &trimmed,
                        &Surface3::Plane(*plane),
                    );
                    let pcb_cone = sampled_pcurve_on_cone(
                        &Curve3::Hyperbola(hyperbola),
                        &trimmed,
                        cone,
                    );
                    let (pca, pcb) = make_pcurves(pca_plane, pcb_cone);
                    let ci = add_curve(
                        self.ds,
                        Curve3::Hyperbola(hyperbola),
                        trimmed,
                        pca,
                        pcb,
                        f1,
                        f2,
                    );
                    curve_indices.push(ci);
                }
            }
        }

        if !curve_indices.is_empty() {
            self.ds.interferences.push(Interference::FaceFace {
                f1,
                f2,
                curves: curve_indices,
                points: vec![],
            });
        }
    }

    // 鈹€鈹€ Cylinder �?Cone analytic face-face intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    fn intersect_cylinder_cone_faces(
        &mut self,
        f1: usize,
        f2: usize,
        cyl: &CylindricalSurface,
        cone: &ConicalSurface,
    ) {
        use inttools::cylinder_cone::{CylinderConeResult, intersect_cylinder_cone};
        use inttools::pcurve_derive::{
            circle_pcurve_on_cone, circle_pcurve_on_cylinder,
            polyline_pcurve_by_projection,
        };
        use std::f64::consts::TAU;

        // Determine which face carries the cylinder (for pcurve_on_a/b ordering).
        let cyl_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Cylinder(_));

        let make_pcurves = |pca: Curve2d, pcb: Curve2d| -> (Option<Curve2d>, Option<Curve2d>) {
            if cyl_is_f1 { (Some(pca), Some(pcb)) } else { (Some(pcb), Some(pca)) }
        };

        match intersect_cylinder_cone(cyl, cone) {
            CylinderConeResult::NoIntersection => (),

            CylinderConeResult::General => {
                self.intersect_ff_by_marching(f1, f2);
            }

            CylinderConeResult::SkewQuartic(branches) => {
                let s1 = Surface3::Cylinder(*cyl);
                let s2 = Surface3::Cone(*cone);
                let mut curve_indices = Vec::new();
                for branch in branches {
                    if branch.len() < 2 {
                        continue;
                    }
                    let v_start = self.ds.add_vertex(branch[0]);
                    let v_end = self.ds.add_vertex(branch[branch.len() - 1]);
                    let dir = (branch[branch.len() - 1] - branch[0])
                        .normalize_or_zero();
                    let ci = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Line(Line3 {
                            origin: branch[0],
                            direction: if dir.length_squared() > 0.5 {
                                dir
                            } else {
                                DVec3::X
                            },
                        }),
                        polyline: branch.clone(),
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, 1.0],
                        pcurve_on_a: polyline_pcurve_by_projection(&branch, &s1),
                        pcurve_on_b: polyline_pcurve_by_projection(&branch, &s2),
                        geom_tol: crate::tolerance::TOLERANCE_ABS,
                    pave_blocks: Vec::new(),
                    });
                    curve_indices.push(ci);
                    self.ds.faces[f1].face_info.curves_sc.insert(ci);
                    self.ds.faces[f2].face_info.curves_sc.insert(ci);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);
                }
                if !curve_indices.is_empty() {
                    self.ds.interferences.push(Interference::FaceFace {
                        f1,
                        f2,
                        curves: curve_indices,
                        points: vec![],
                    });
                }
            }

            CylinderConeResult::ParallelOffsetPolyline(branches) => {
                let s1 = Surface3::Cylinder(*cyl);
                let s2 = Surface3::Cone(*cone);
                let mut curve_indices = Vec::new();
                for branch in branches {
                    if branch.len() < 2 {
                        continue;
                    }
                    let v_start = self.ds.add_vertex(branch[0]);
                    let v_end = self.ds.add_vertex(branch[branch.len() - 1]);
                    let dir = (branch[branch.len() - 1] - branch[0])
                        .normalize_or_zero();
                    let ci = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Line(Line3 {
                            origin: branch[0],
                            direction: if dir.length_squared() > 0.5 {
                                dir
                            } else {
                                DVec3::X
                            },
                        }),
                        polyline: branch.clone(),
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, 1.0],
                        pcurve_on_a: polyline_pcurve_by_projection(&branch, &s1),
                        pcurve_on_b: polyline_pcurve_by_projection(&branch, &s2),
                        geom_tol: crate::tolerance::TOLERANCE_ABS,
                    pave_blocks: Vec::new(),
                    });
                    curve_indices.push(ci);
                    self.ds.faces[f1].face_info.curves_sc.insert(ci);
                    self.ds.faces[f2].face_info.curves_sc.insert(ci);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);
                }
                if !curve_indices.is_empty() {
                    self.ds.interferences.push(Interference::FaceFace {
                        f1,
                        f2,
                        curves: curve_indices,
                        points: vec![],
                    });
                }
            }

            CylinderConeResult::CoaxialCircle(circle) => {
                let pca_cyl = circle_pcurve_on_cylinder(&circle, cyl);
                let pcb_cone = circle_pcurve_on_cone(&circle, cone);
                let (pca, pcb) = make_pcurves(pca_cyl, pcb_cone);

                let pts = sample_circle_arc(&circle, 0.0, TAU, 32);
                let v_start = self.ds.add_vertex(pts[0]);
                let v_end = self.ds.add_vertex(pts[pts.len() - 1]);
                let ci = self.ds.intersection_curves.len();
                self.ds.intersection_curves.push(IntersectionCurve {
                    curve: Curve3::Circle(circle),
                    polyline: vec![],
                    start_vertex: v_start,
                    end_vertex: v_end,
                    t_range: [0.0, TAU],
                    pcurve_on_a: pca,
                    pcurve_on_b: pcb,
                    geom_tol: crate::tolerance::TOLERANCE_ABS,
                pave_blocks: Vec::new(),
                });
                self.ds.faces[f1].face_info.curves_sc.insert(ci);
                self.ds.faces[f2].face_info.curves_sc.insert(ci);
                self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                self.ds.faces[f2].face_info.vertices_in.insert(v_end);

                self.ds.interferences.push(Interference::FaceFace {
                    f1,
                    f2,
                    curves: vec![ci],
                    points: vec![],
                });
            }

            CylinderConeResult::CoaxialTwoCircles(c1, c2) => {
                let mut curve_indices = Vec::new();
                for circle in [c1, c2] {
                    let pca_cyl = circle_pcurve_on_cylinder(&circle, cyl);
                    let pcb_cone = circle_pcurve_on_cone(&circle, cone);
                    let (pca, pcb) = make_pcurves(pca_cyl, pcb_cone);

                    let pts = sample_circle_arc(&circle, 0.0, TAU, 32);
                    let v_start = self.ds.add_vertex(pts[0]);
                    let v_end = self.ds.add_vertex(pts[pts.len() - 1]);
                    let ci = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Circle(circle),
                        polyline: vec![],
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, TAU],
                        pcurve_on_a: pca,
                        pcurve_on_b: pcb,
                        geom_tol: crate::tolerance::TOLERANCE_ABS,
                    pave_blocks: Vec::new(),
                    });
                    self.ds.faces[f1].face_info.curves_sc.insert(ci);
                    self.ds.faces[f2].face_info.curves_sc.insert(ci);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);
                    curve_indices.push(ci);
                }
                self.ds.interferences.push(Interference::FaceFace {
                    f1,
                    f2,
                    curves: curve_indices,
                    points: vec![],
                });
            }
        }
    }

    // 鈹€鈹€ Cone �?Cone analytic face-face intersection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    fn intersect_cone_cone_faces(
        &mut self,
        f1: usize,
        f2: usize,
        cone1: &ConicalSurface,
        cone2: &ConicalSurface,
    ) {
        use inttools::cone_cone::{ConeConeResult, intersect_cone_cone};
        use inttools::pcurve_derive::{circle_pcurve_on_cone, polyline_pcurve_by_projection};
        use std::f64::consts::TAU;

        // Determine which face is cone1 (for pcurve_on_a/b ordering).
        let cone1_is_f1 = {
            if let Surface3::Cone(c) = &self.ds.faces[f1].surface {
                (c.apex - cone1.apex).length_squared() < TOLERANCE_LINEAR_ULTRA_STRICT * TOLERANCE_LINEAR_ULTRA_STRICT
                    && (c.axis - cone1.axis).length_squared() < TOLERANCE_LINEAR_ULTRA_STRICT * TOLERANCE_LINEAR_ULTRA_STRICT
            } else {
                false
            }
        };

        let make_pcurves = |pca: Curve2d, pcb: Curve2d| -> (Option<Curve2d>, Option<Curve2d>) {
            if cone1_is_f1 { (Some(pca), Some(pcb)) } else { (Some(pcb), Some(pca)) }
        };

        match intersect_cone_cone(cone1, cone2) {
            ConeConeResult::NoIntersection => (),
            ConeConeResult::Coaxial => {
                self.ds.interferences.push(Interference::FaceFace {
                    f1, f2, curves: vec![], points: vec![],
                });
            }

            ConeConeResult::CoaxialPoint(_pt) => {
                // Single shared apex �?a point contact, not a curve.
            }

            ConeConeResult::General => {
                self.intersect_ff_by_marching(f1, f2);
            }

            ConeConeResult::SkewQuartic(branches) => {
                let s1 = Surface3::Cone(*cone1);
                let s2 = Surface3::Cone(*cone2);
                let mut curve_indices = Vec::new();
                for branch in branches {
                    if branch.len() < 2 {
                        continue;
                    }
                    let v_start = self.ds.add_vertex(branch[0]);
                    let v_end = self.ds.add_vertex(branch[branch.len() - 1]);
                    let dir = (branch[branch.len() - 1] - branch[0])
                        .normalize_or_zero();
                    let ci = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Line(Line3 {
                            origin: branch[0],
                            direction: if dir.length_squared() > 0.5 {
                                dir
                            } else {
                                DVec3::X
                            },
                        }),
                        polyline: branch.clone(),
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, 1.0],
                        pcurve_on_a: polyline_pcurve_by_projection(&branch, &s1),
                        pcurve_on_b: polyline_pcurve_by_projection(&branch, &s2),
                        geom_tol: crate::tolerance::TOLERANCE_ABS,
                    pave_blocks: Vec::new(),
                    });
                    curve_indices.push(ci);
                    self.ds.faces[f1].face_info.curves_sc.insert(ci);
                    self.ds.faces[f2].face_info.curves_sc.insert(ci);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);
                }
                if !curve_indices.is_empty() {
                    self.ds.interferences.push(Interference::FaceFace {
                        f1,
                        f2,
                        curves: curve_indices,
                        points: vec![],
                    });
                }
            }

            ConeConeResult::CoaxialCircle(circle) => {
                let pca_cone1 = circle_pcurve_on_cone(&circle, cone1);
                let pcb_cone2 = circle_pcurve_on_cone(&circle, cone2);
                let (pca, pcb) = make_pcurves(pca_cone1, pcb_cone2);

                let pts = sample_circle_arc(&circle, 0.0, TAU, 32);
                let v_start = self.ds.add_vertex(pts[0]);
                let v_end = self.ds.add_vertex(pts[pts.len() - 1]);
                let ci = self.ds.intersection_curves.len();
                self.ds.intersection_curves.push(IntersectionCurve {
                    curve: Curve3::Circle(circle),
                    polyline: vec![],
                    start_vertex: v_start,
                    end_vertex: v_end,
                    t_range: [0.0, TAU],
                    pcurve_on_a: pca,
                    pcurve_on_b: pcb,
                    geom_tol: crate::tolerance::TOLERANCE_ABS,
                pave_blocks: Vec::new(),
                });
                self.ds.faces[f1].face_info.curves_sc.insert(ci);
                self.ds.faces[f2].face_info.curves_sc.insert(ci);
                self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                self.ds.faces[f2].face_info.vertices_in.insert(v_end);

                self.ds.interferences.push(Interference::FaceFace {
                    f1,
                    f2,
                    curves: vec![ci],
                    points: vec![],
                });
            }
        }
    }

    // 鈹€鈹€ Torus intersection helpers 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    fn register_torus_intersection(
        &mut self,
        f1: usize,
        f2: usize,
        s1: &Surface3,
        s2: &Surface3,
        torus_is_f1: bool,
    ) {
        use inttools::intss::{intersect_surfaces_with_density_tol, SurfaceCurve};
        use inttools::pcurve_derive::polyline_pcurve_by_projection;

        let result = intersect_surfaces_with_density_tol(s1, s2, 48, self.ff_tol(f1, f2));
        if result.is_empty() {
            return;
        }

        for sir in &result.curves {
            match &sir.curve_3d {
                SurfaceCurve::Circle(circle) => {
                    // Only split into half-circles for torus脳cylinder intersections where the
                    // full circle spans 100% of cylinder U (triggers has_full_wrap fallback).
                    // For other surface types the full circle is simpler and more robust.
                    // Note: s1 is always Torus by calling convention, s2 is the other surface.
                    if matches!(s2, Surface3::Cylinder(_)) {
                        let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
                            if torus_is_f1 { (Some(a.clone()), Some(b.clone())) }
                            else { (Some(b.clone()), Some(a.clone())) }
                        } else { (None, None) };

                        let mut curve_indices = Vec::new();
                        for (t0, t1) in [(0.0, std::f64::consts::PI), (std::f64::consts::PI, std::f64::consts::TAU)] {
                            let pts = sample_circle_arc(circle, t0, t1, 16);
                            if pts.len() < 2 { continue; }
                            let v_start = self.ds.add_vertex(pts[0]);
                            let v_end = self.ds.add_vertex(pts[pts.len() - 1]);

                let curve_idx = self.ds.intersection_curves.len();
                eprintln!("[IC] CREATE ci={} f1={} f2={} sv={} ev={}", curve_idx, f1, f2, v_start, v_end);
                            self.ds.intersection_curves.push(IntersectionCurve {
                                curve: Curve3::Circle(*circle),
                                polyline: vec![],
                                start_vertex: v_start,
                                end_vertex: v_end,
                                t_range: [t0, t1],
                                pcurve_on_a: pca.clone(),
                                pcurve_on_b: pcb.clone(),
                                geom_tol: crate::tolerance::TOLERANCE_ABS,
                            pave_blocks: Vec::new(),
                            });

                            self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
                            self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
                            self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                            self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                            self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                            self.ds.faces[f2].face_info.vertices_in.insert(v_end);

                            curve_indices.push(curve_idx);
                        }

                        if !curve_indices.is_empty() {
                            self.ds.interferences.push(Interference::FaceFace {
                                f1, f2, curves: curve_indices, points: vec![],
                            });
                        }
                    } else {
                        let pts = sample_circle_arc(circle, 0.0, std::f64::consts::TAU, 32);
                        if pts.len() < 2 { continue; }
                        let v_start = self.ds.add_vertex(pts[0]);
                        let v_end = self.ds.add_vertex(pts[pts.len() - 1]);
                        let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
                            if torus_is_f1 { (Some(a.clone()), Some(b.clone())) }
                            else { (Some(b.clone()), Some(a.clone())) }
                        } else { (None, None) };

                        let curve_idx = self.ds.intersection_curves.len();
                        self.ds.intersection_curves.push(IntersectionCurve {
                            curve: Curve3::Circle(*circle),
                            polyline: vec![],
                            start_vertex: v_start,
                            end_vertex: v_end,
                            t_range: [0.0, std::f64::consts::TAU],
                            pcurve_on_a: pca,
                            pcurve_on_b: pcb,
                            geom_tol: crate::tolerance::TOLERANCE_ABS,
                        pave_blocks: Vec::new(),
                        });

                        self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
                        self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
                        self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                        self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                        self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                        self.ds.faces[f2].face_info.vertices_in.insert(v_end);

                        self.ds.interferences.push(Interference::FaceFace {
                            f1, f2, curves: vec![curve_idx], points: vec![],
                        });
                    }
                }
                SurfaceCurve::Polyline(pts) => {
                    if pts.len() < 2 {
                        continue;
                    }
                    let v_start = self.ds.add_vertex(pts[0]);
                    let v_end = self.ds.add_vertex(pts[pts.len() - 1]);

                    let arc_len: f64 = pts.windows(2).map(|w| (w[1] - w[0]).length()).sum();
                    let dir = (pts[pts.len() - 1] - pts[0]).normalize_or_zero();

                    let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
                        if torus_is_f1 {
                            (Some(a.clone()), Some(b.clone()))
                        } else {
                            (Some(b.clone()), Some(a.clone()))
                        }
                    } else {
                        (
                            polyline_pcurve_by_projection(pts, s1),
                            polyline_pcurve_by_projection(pts, s2),
                        )
                    };

                    let curve_idx = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Line(Line3 {
                            origin: pts[0],
                            direction: if dir.length_squared() > 0.5 { dir } else { DVec3::X },
                        }),
                        polyline: pts.clone(),
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, arc_len.max(TOLERANCE_LINEAR_ULTRA_STRICT)],
                        pcurve_on_a: pca,
                        pcurve_on_b: pcb,
                        geom_tol: crate::tolerance::TOLERANCE_ABS,
                    pave_blocks: Vec::new(),
                    });

                    self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
                    self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);

                    self.ds.interferences.push(Interference::FaceFace {
                        f1,
                        f2,
                        curves: vec![curve_idx],
                        points: vec![],
                    });
                }
                SurfaceCurve::Ellipse(ellipse) => {
                    let pts = sample_circle_arc(
                        &Circle3 {
                            center: ellipse.center,
                            normal: ellipse.normal,
                            radius: ellipse.major_radius,
                        },
                        0.0,
                        std::f64::consts::TAU,
                        32,
                    );
                    if pts.len() < 2 {
                        continue;
                    }
                    let v_start = self.ds.add_vertex(pts[0]);
                    let v_end = self.ds.add_vertex(pts[pts.len() - 1]);

                    let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
                        if torus_is_f1 {
                            (Some(a.clone()), Some(b.clone()))
                        } else {
                            (Some(b.clone()), Some(a.clone()))
                        }
                    } else {
                        (None, None)
                    };

                    let curve_idx = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Ellipse(*ellipse),
                        polyline: vec![],
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, std::f64::consts::TAU],
                        pcurve_on_a: pca,
                        pcurve_on_b: pcb,
                        geom_tol: crate::tolerance::TOLERANCE_ABS,
                    pave_blocks: Vec::new(),
                    });

                    self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
                    self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);

                    self.ds.interferences.push(Interference::FaceFace {
                        f1,
                        f2,
                        curves: vec![curve_idx],
                        points: vec![],
                    });
                }
                SurfaceCurve::Line(line) => {
                    let pts = self.ds.face_boundary_points(f1);
                    let pts2 = self.ds.face_boundary_points(f2);
                    let bbox1_min = pts.iter().fold(DVec3::INFINITY, |a, &b| a.min(b));
                    let bbox1_max = pts.iter().fold(DVec3::NEG_INFINITY, |a, &b| a.max(b));
                    let bbox2_min = pts2.iter().fold(DVec3::INFINITY, |a, &b| a.min(b));
                    let bbox2_max = pts2.iter().fold(DVec3::NEG_INFINITY, |a, &b| a.max(b));

                    let lo = bbox1_min.min(bbox2_min);
                    let hi = bbox1_max.max(bbox2_max);
                    let extent = (hi - lo).length() * 0.5 + 1.0;

                    let p_start = line.origin + line.direction * (-extent);
                    let p_end = line.origin + line.direction * extent;

                    let v_start = self.ds.add_vertex(p_start);
                    let v_end = self.ds.add_vertex(p_end);

                    let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
                        if torus_is_f1 {
                            (Some(a.clone()), Some(b.clone()))
                        } else {
                            (Some(b.clone()), Some(a.clone()))
                        }
                    } else {
                        (None, None)
                    };

                    let curve_idx = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Line(*line),
                        polyline: vec![p_start, p_end],
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [-extent, extent],
                        pcurve_on_a: pca,
                        pcurve_on_b: pcb,
                        geom_tol: crate::tolerance::TOLERANCE_ABS,
                    pave_blocks: Vec::new(),
                    });

                    self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
                    self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);

                    self.ds.interferences.push(Interference::FaceFace {
                        f1,
                        f2,
                        curves: vec![curve_idx],
                        points: vec![],
                    });
                }
                SurfaceCurve::BSplineCurve(b) => {
                    // Sample the BSpline to produce a polyline for face splitting.
                    use rcad_kernel::geom::CurveEval;
                    let n_samples = 33_usize;
                    let mut pts: Vec<DVec3> = Vec::with_capacity(n_samples);
                    for i in 0..n_samples {
                        let t = i as f64 / (n_samples - 1) as f64;
                        pts.push(b.point_at(t));
                    }
                    if pts.len() < 2 {
                        continue;
                    }
                    let v_start = self.ds.add_vertex(pts[0]);
                    let v_end = self.ds.add_vertex(pts[pts.len() - 1]);

                    let arc_len: f64 = pts.windows(2).map(|w| (w[1] - w[0]).length()).sum();
                    let dir = (pts[pts.len() - 1] - pts[0]).normalize_or_zero();

                    let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
                        if torus_is_f1 {
                            (Some(a.clone()), Some(b.clone()))
                        } else {
                            (Some(b.clone()), Some(a.clone()))
                        }
                    } else {
                        (
                            polyline_pcurve_by_projection(&pts, s1),
                            polyline_pcurve_by_projection(&pts, s2),
                        )
                    };

                    let curve_idx = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::BSpline((**b).clone()),
                        polyline: pts.clone(),
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, arc_len.max(TOLERANCE_LINEAR_ULTRA_STRICT)],
                        pcurve_on_a: pca,
                        pcurve_on_b: pcb,
                        geom_tol: crate::tolerance::TOLERANCE_ABS,
                    pave_blocks: Vec::new(),
                    });

                    self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
                    self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);

                    self.ds.interferences.push(Interference::FaceFace {
                        f1,
                        f2,
                        curves: vec![curve_idx],
                        points: vec![],
                    });
                }
                SurfaceCurve::Point(_) | SurfaceCurve::Parabola(_) | SurfaceCurve::Hyperbola(_) => {
                    // Skip degenerate / unsupported curve types for now
                }
            }
        }
    }

    fn intersect_torus_plane_faces(
        &mut self,
        f1: usize,
        f2: usize,
        torus: &ToroidalSurface,
        plane: &Plane,
    ) {
        let torus_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Torus(_));
        let s1 = Surface3::Torus(*torus);
        let s2 = Surface3::Plane(*plane);
        self.register_torus_intersection(f1, f2, &s1, &s2, torus_is_f1);
    }

    fn intersect_torus_sphere_faces(
        &mut self,
        f1: usize,
        f2: usize,
        torus: &ToroidalSurface,
        sphere: &SphericalSurface,
    ) {
        let torus_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Torus(_));
        let s1 = Surface3::Torus(*torus);
        let s2 = Surface3::Sphere(*sphere);
        self.register_torus_intersection(f1, f2, &s1, &s2, torus_is_f1);
    }

    fn intersect_torus_cylinder_faces(
        &mut self,
        f1: usize,
        f2: usize,
        torus: &ToroidalSurface,
        cylinder: &CylindricalSurface,
    ) {
        let torus_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Torus(_));
        let s1 = Surface3::Torus(*torus);
        let s2 = Surface3::Cylinder(*cylinder);
        self.register_torus_intersection(f1, f2, &s1, &s2, torus_is_f1);
    }

    fn intersect_torus_cone_faces(
        &mut self,
        f1: usize,
        f2: usize,
        torus: &ToroidalSurface,
        cone: &ConicalSurface,
    ) {
        let torus_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Torus(_));
        let s1 = Surface3::Torus(*torus);
        let s2 = Surface3::Cone(*cone);
        self.register_torus_intersection(f1, f2, &s1, &s2, torus_is_f1);
    }

    fn intersect_torus_torus_faces(
        &mut self,
        f1: usize,
        f2: usize,
        torus1: &ToroidalSurface,
        torus2: &ToroidalSurface,
    ) {
        let torus_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Torus(_));
        let s1 = Surface3::Torus(*torus1);
        let s2 = Surface3::Torus(*torus2);
        self.register_torus_intersection(f1, f2, &s1, &s2, torus_is_f1);
    }

    fn intersect_sphere_cone_faces(
        &mut self,
        f1: usize,
        f2: usize,
        sphere: &SphericalSurface,
        cone: &ConicalSurface,
    ) {
        use inttools::sphere_cone::{SphereConeResult, intersect_sphere_cone};
        use inttools::pcurve_derive::{
            circle_pcurve_on_cone, fallback_pcurve_by_projection,
            polyline_pcurve_by_projection,
        };
        use std::f64::consts::TAU;

        let sphere_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Sphere(_));

        let make_pcurves = |pca: Curve2d, pcb: Curve2d| -> (Option<Curve2d>, Option<Curve2d>) {
            if sphere_is_f1 { (Some(pca), Some(pcb)) } else { (Some(pcb), Some(pca)) }
        };

        let s1 = Surface3::Sphere(*sphere);
        let s2 = Surface3::Cone(*cone);

        match intersect_sphere_cone(sphere, cone) {
            SphereConeResult::NoIntersection => (),

            SphereConeResult::General => {
                self.intersect_ff_by_marching(f1, f2);
            }

            SphereConeResult::SingleCircle(circ) => {
                let pca = fallback_pcurve_by_projection(
                    &Curve3::Circle(circ),
                    &[0.0, TAU],
                    &s1,
                );
                let pcb = circle_pcurve_on_cone(&circ, cone);
                let (pca, pcb) = make_pcurves(pca, pcb);
                let pts = sample_circle_arc(&circ, 0.0, TAU, 32);
                let v_start = self.ds.add_vertex(pts[0]);
                let v_end = self.ds.add_vertex(pts[pts.len() - 1]);
                let ci = self.ds.intersection_curves.len();
                self.ds.intersection_curves.push(IntersectionCurve {
                    curve: Curve3::Circle(circ),
                    polyline: vec![],
                    start_vertex: v_start,
                    end_vertex: v_end,
                    t_range: [0.0, TAU],
                    pcurve_on_a: pca,
                    pcurve_on_b: pcb,
                    geom_tol: crate::tolerance::TOLERANCE_ABS,
                pave_blocks: Vec::new(),
                });
                self.ds.faces[f1].face_info.curves_sc.insert(ci);
                self.ds.faces[f2].face_info.curves_sc.insert(ci);
                self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                self.ds.faces[f2].face_info.vertices_in.insert(v_end);
                self.ds.interferences.push(Interference::FaceFace {
                    f1, f2, curves: vec![ci], points: vec![],
                });
            }

            SphereConeResult::TwoCircles(c1, c2) => {
                for circ in [c1, c2] {
                    let pca = fallback_pcurve_by_projection(
                        &Curve3::Circle(circ),
                        &[0.0, TAU],
                        &s1,
                    );
                    let pcb = circle_pcurve_on_cone(&circ, cone);
                    let (pca, pcb) = make_pcurves(pca, pcb);
                    let pts = sample_circle_arc(&circ, 0.0, TAU, 32);
                    let v_start = self.ds.add_vertex(pts[0]);
                    let v_end = self.ds.add_vertex(pts[pts.len() - 1]);
                    let ci = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Circle(circ),
                        polyline: vec![],
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, TAU],
                        pcurve_on_a: pca,
                        pcurve_on_b: pcb,
                        geom_tol: crate::tolerance::TOLERANCE_ABS,
                    pave_blocks: Vec::new(),
                    });
                    self.ds.faces[f1].face_info.curves_sc.insert(ci);
                    self.ds.faces[f2].face_info.curves_sc.insert(ci);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);
                    self.ds.interferences.push(Interference::FaceFace {
                        f1, f2, curves: vec![ci], points: vec![],
                    });
                }
            }

            SphereConeResult::TangentPoint(pt) => {
                let v = self.ds.add_vertex(pt);
                self.ds.faces[f1].face_info.vertices_in.insert(v);
                self.ds.faces[f2].face_info.vertices_in.insert(v);
                self.ds.interferences.push(Interference::FaceFace {
                    f1, f2, curves: vec![], points: vec![v],
                });
            }

            SphereConeResult::Polyline(branches) => {
                let mut curve_indices = Vec::new();
                for branch in branches {
                    if branch.len() < 2 { continue; }
                    let v_start = self.ds.add_vertex(branch[0]);
                    let v_end = self.ds.add_vertex(branch[branch.len() - 1]);
                    let dir = (branch[branch.len() - 1] - branch[0])
                        .normalize_or_zero();
                    let ci = self.ds.intersection_curves.len();
                    self.ds.intersection_curves.push(IntersectionCurve {
                        curve: Curve3::Line(Line3 {
                            origin: branch[0],
                            direction: if dir.length_squared() > 0.5 {
                                dir
                            } else {
                                DVec3::X
                            },
                        }),
                        polyline: branch.clone(),
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, 1.0],
                        pcurve_on_a: polyline_pcurve_by_projection(&branch, &s1),
                        pcurve_on_b: polyline_pcurve_by_projection(&branch, &s2),
                        geom_tol: crate::tolerance::TOLERANCE_ABS,
                    pave_blocks: Vec::new(),
                    });
                    curve_indices.push(ci);
                    self.ds.faces[f1].face_info.curves_sc.insert(ci);
                    self.ds.faces[f2].face_info.curves_sc.insert(ci);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);
                }
                if !curve_indices.is_empty() {
                    self.ds.interferences.push(Interference::FaceFace {
                        f1, f2, curves: curve_indices, points: vec![],
                    });
                }
            }
        }
    }

    fn intersect_ff_by_numeric_intss(
        &mut self,
        f1: usize,
        f2: usize,
        s1: &Surface3,
        s2: &Surface3,
        grid_n: usize,
    ) {
        use inttools::intss::numeric_intss_with_domains;
        use inttools::pcurve_derive::polyline_pcurve_by_projection;

        // Use face-specific UV domains (set up by DS::setup_uv_boundaries)
        // if available.  For cylinders this encodes the actual face height range,
        // ensuring the intersection polyline endpoints fall *inside* the UV
        // boundary rectangle and can be used to split it.
        let dom1 = self.ds.faces[f1]
            .uv_boundary
            .as_ref()
            .and_then(|uv| {
                if uv.len() >= 3 {
                    let u_min = uv.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
                    let u_max = uv.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
                    let v_min = uv.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
                    let v_max = uv.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
                    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
                        return Some([u_min, u_max, v_min, v_max]);
                    }
                }
                None
            });
        let dom2 = self.ds.faces[f2]
            .uv_boundary
            .as_ref()
            .and_then(|uv| {
                if uv.len() >= 3 {
                    let u_min = uv.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
                    let u_max = uv.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
                    let v_min = uv.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
                    let v_max = uv.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
                    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
                        return Some([u_min, u_max, v_min, v_max]);
                    }
                }
                None
            });

        let result = numeric_intss_with_domains(
            s1,
            s2,
            grid_n,
            dom1,
            dom2,
            Some(self.ff_tol(f1, f2)),
        );
        if result.is_empty() {
            return;
        }

        let mut curve_indices = Vec::new();
        for sir in &result.curves {
            let (mut chain, approx_curve) = match &sir.curve_3d {
                crate::inttools::intss::SurfaceCurve::Polyline(pts) => (pts.clone(), None),
                crate::inttools::intss::SurfaceCurve::BSplineCurve(bs) => {
                    // Sample BSpline back to polyline for face-boundary snapping
                    let n = 64usize;
                    let pts: Vec<DVec3> = (0..=n)
                        .map(|i| {
                            let t = i as f64 / n as f64;
                            bs.point_at(t)
                        })
                        .collect();
                    (pts, Some(Curve3::BSpline((**bs).clone())))
                }
                _ => continue,
            };
            if chain.len() < 2 {
                continue;
            }

            self.snap_polyline_endpoints_to_face_boundaries(&mut chain, f1, f2);

            let v_start = self.ds.add_vertex(chain[0]);
            let v_end = self.ds.add_vertex(chain[chain.len() - 1]);

            let arc_len: f64 = chain.windows(2).map(|w| (w[1] - w[0]).length()).sum();
            let dir = (chain[chain.len() - 1] - chain[0]).normalize_or_zero();
            let pcurve_a = sir
                .pcurve_on_a
                .clone()
                .or_else(|| polyline_pcurve_by_projection(&chain, s1));
            let pcurve_b = sir
                .pcurve_on_b
                .clone()
                .or_else(|| polyline_pcurve_by_projection(&chain, s2));

            let curve_idx = self.ds.intersection_curves.len();
            self.ds.intersection_curves.push(IntersectionCurve {
                curve: approx_curve.unwrap_or(Curve3::Line(Line3 {
                    origin: chain[0],
                    direction: if dir.length_squared() > 0.5 {
                        dir
                    } else {
                        DVec3::X
                    },
                })),
                polyline: chain,
                start_vertex: v_start,
                end_vertex: v_end,
                t_range: [0.0, arc_len.max(TOLERANCE_LINEAR_ULTRA_STRICT)],
                pcurve_on_a: pcurve_a,
                pcurve_on_b: pcurve_b,
                geom_tol: crate::tolerance::TOLERANCE_ABS,
            pave_blocks: Vec::new(),
            });

            self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
            self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
            self.ds.faces[f1].face_info.vertices_in.insert(v_start);
            self.ds.faces[f1].face_info.vertices_in.insert(v_end);
            self.ds.faces[f2].face_info.vertices_in.insert(v_start);
            self.ds.faces[f2].face_info.vertices_in.insert(v_end);

            curve_indices.push(curve_idx);
        }

        if !curve_indices.is_empty() {
            self.ds.interferences.push(Interference::FaceFace {
                f1,
                f2,
                curves: curve_indices,
                points: vec![],
            });
        }
    }

    fn intersect_ff_by_marching(&mut self, f1: usize, f2: usize) {
        use inttools::marching::{adaptive_sampling_density, MarchingConfig};

        let s1 = self.ds.faces[f1].surface.clone();
        let s2 = self.ds.faces[f2].surface.clone();

        // OCCT-aligned: use sign-change grid marching (IntTools_FaceFace / IntPatch_ImpPrmIntersection).
        // No BSpline demotion �?BSpline surfaces stay as parametric (ts=0) and use UV grid marching.
        let any_curved = !matches!(&s1, Surface3::Plane(_)) || !matches!(&s2, Surface3::Plane(_));
        if any_curved {
            let char_len = |s: &Surface3| -> f64 {
                match s {
                    Surface3::Sphere(sp) => sp.radius,
                    Surface3::Cylinder(cy) => cy.radius,
                    Surface3::Cone(co) => co.radius.max(0.5),
                    Surface3::Torus(to) => to.major_radius.max(to.minor_radius),
                    Surface3::BSpline(bsp) => {
                        if bsp.control_points.is_empty() {
                            1.0
                        } else {
                            let mut mn = DVec3::splat(f64::INFINITY);
                            let mut mx = DVec3::splat(f64::NEG_INFINITY);
                            for row in &bsp.control_points {
                                for p in row {
                                    mn = mn.min(*p); mx = mx.max(*p);
                                }
                            }
                            (mx - mn).length().max(0.5) * 0.5
                        }
                    }
                    Surface3::Bezier(bez) => {
                        if bez.control_points.is_empty() {
                            1.0
                        } else {
                            let mut mn = DVec3::splat(f64::INFINITY);
                            let mut mx = DVec3::splat(f64::NEG_INFINITY);
                            for row in &bez.control_points {
                                for p in row {
                                    mn = mn.min(*p); mx = mx.max(*p);
                                }
                            }
                            (mx - mn).length().max(0.5) * 0.5
                        }
                    }
                    _ => 1.0,
                }
            };
            let avg_len = (char_len(&s1) + char_len(&s2)) * 0.5;
            let mut grid_n = ((avg_len * 10.0) as usize).max(64).min(256);

            let skew_factor = match (&s1, &s2) {
                (Surface3::Cylinder(c1), Surface3::Cone(c2))
                | (Surface3::Cone(c2), Surface3::Cylinder(c1)) => {
                    let a1 = c1.axis.normalize();
                    let a2 = c2.axis.normalize();
                    let sin_angle = a1.cross(a2).length();
                    (1.0 + sin_angle * 3.0).min(3.0)
                }
                _ => 1.0,
            };
            grid_n = ((grid_n as f64 * skew_factor) as usize).min(256);

            self.intersect_ff_by_numeric_intss(f1, f2, &s1, &s2, grid_n);
            return;
        }

        // Use adaptive sampling density based on surface geometry
        let base_density = 16usize;
        let sampling1 = adaptive_sampling_density(&s1, base_density);
        let sampling2 = adaptive_sampling_density(&s2, base_density);
        // Use the higher density to ensure we don't miss narrow intersections
        let n_u = sampling1.n_u.max(sampling2.n_u);
        let n_v = sampling1.n_v.max(sampling2.n_v);

        let _samples = self.generate_surface_samples_grid(&s1, n_u, n_v);
        // Use multi-scale seed detection for improved robustness
        // Scales: coarse (8x8), medium (16x16), fine (32x32), ultra (64x64)
        let base_step = self.estimate_step_size(&s1, &s2);
        let seed_dedup = (base_step * 2.0).max(self.ff_tol(f1, f2) * 2.0);
        let seeds = inttools::marching::find_seed_points_multiscale(
            &s1,
            &s2,
            |nu, nv| self.generate_surface_samples_grid(&s1, nu, nv),
            &[8, 16, 32, 64],
            seed_dedup,
        );

        if seeds.is_empty() {
            return;
        }

        // Compute a finite bounding box that contains both faces' intersection region.
        // Use boundary vertices (actual face extent) with a generous margin.
        let bounds_from_face = |face_idx: usize| -> (DVec3, DVec3) {
            let mut mn = DVec3::splat(f64::INFINITY);
            let mut mx = DVec3::splat(f64::NEG_INFINITY);
            // Use boundary vertices (from wire edges)
            for &vi in &self.ds.faces[face_idx].boundary_verts {
                let p = self.ds.vertices[vi].point;
                mn = mn.min(p);
                mx = mx.max(p);
            }
            // Also sample boundary edges for curved edges (e.g. circles)
            for &ei in &self.ds.faces[face_idx].boundary_edges {
                if let Some(edge) = self.ds.edges.get(ei) {
                    let [t0, t1] = edge.t_range;
                    for k in 0..=8usize {
                        let t = t0 + (t1 - t0) * k as f64 / 8.0;
                        let p = edge.curve.point_at(t);
                        if p.is_finite() {
                            mn = mn.min(p);
                            mx = mx.max(p);
                        }
                    }
                }
            }
            // If still infinite, use a generous default
            if !mn.is_finite() || !mx.is_finite() {
                mn = DVec3::splat(-10.0);
                mx = DVec3::splat(10.0);
            }
            (mn, mx)
        };

        let (mn1, mx1) = bounds_from_face(f1);
        let (mn2, mx2) = bounds_from_face(f2);
        let margin = 1.0;
        let aabb_min = mn1.min(mn2) - DVec3::splat(margin);
        let aabb_max = mx1.max(mx2) + DVec3::splat(margin);

        // Use adaptive step size based on characteristic lengths
        let char_len = sampling1.characteristic_length.min(sampling2.characteristic_length);
        let step_size = base_step.min(char_len * 0.5).max(TOLERANCE_MESH_LEGACY);

        // Configure marching with convergence monitoring
        let marching_config = MarchingConfig {
            step_size,
            min_step_size: step_size * 0.01,
            max_steps: 500,
            max_oscillations: 3,
            step_reduction_factor: 0.5,
            deflection_tol: step_size * 0.001,
            multiscale_seeds: true,
        };

        let mut curve_indices = Vec::new();
        // Track all points already covered by marched curves, to deduplicate
        // seeds that trace the same intersection curve.
        let mut covered_points: Vec<DVec3> = Vec::new();
        let ff = self.ff_tol(f1, f2);
        let dedup_tol = (step_size * 3.0).max(ff * 2.0);

        for seed in seeds {
            // Skip if this seed is near any point already covered by a previous curve
            if covered_points
                .iter()
                .any(|&cp| (cp - seed).length_squared() < dedup_tol * dedup_tol)
            {
                continue;
            }

            let curve = inttools::marching::march_intersection_with_config(
                &s1,
                &s2,
                seed,
                &marching_config,
                |p: DVec3| p.cmpge(aabb_min).all() && p.cmple(aabb_max).all(),
            );

            if curve.points.len() < 2 {
                continue;
            }

            // Mark all curve points as covered (sample every few for efficiency)
            for (i, &p) in curve.points.iter().enumerate() {
                if i % 5 == 0 {
                    covered_points.push(p);
                }
            }

            let v_start = self.ds.add_vertex(curve.points[0]);
            let v_end = self.ds.add_vertex(curve.points[curve.points.len() - 1]);

            let curve_idx = self.ds.intersection_curves.len();
            // Compute arc-length for t_range
            let arc_len: f64 = curve
                .points
                .windows(2)
                .map(|w| (w[1] - w[0]).length())
                .sum();
            let dir = (curve.points[curve.points.len() - 1] - curve.points[0]).normalize_or_zero();
            let t_range = [0.0, arc_len.max(TOLERANCE_LINEAR_ULTRA_STRICT)];

            // �?OCCT-aligned:reApprox �?validate pcurves; retry with loose tolerance
            // if validation fails.
            let (pcurve_a, pcurve_b) = self.make_marching_pcurves_with_reapprox(
                &curve.points, &s1, &s2, f1, f2, &t_range,
            );

            // OCCT-aligned: approximate marching polyline to BSpline (MakeCurve / GeomInt_IntSS::MakeBSpline)
            let approx_curve = if curve.points.len() >= 4 {
                crate::inttools::intss::polyline_to_bspline(&curve.points, TOLERANCE_TOL_SCALE_MICRO)
                    .filter(|c| matches!(c, Curve3::BSpline(_)))
            } else {
                None
            };

            self.ds.intersection_curves.push(IntersectionCurve {
                curve: approx_curve.unwrap_or(Curve3::Line(Line3 {
                    origin: curve.points[0],
                    direction: if dir.length_squared() > 0.5 {
                        dir
                    } else {
                        DVec3::X
                    },
                })),
                polyline: curve.points.clone(),
                start_vertex: v_start,
                end_vertex: v_end,
                t_range,
                pcurve_on_a: pcurve_a,
                pcurve_on_b: pcurve_b,
                geom_tol: crate::tolerance::TOLERANCE_ABS,
            pave_blocks: Vec::new(),
            });

            self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
            self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
            self.ds.faces[f1].face_info.vertices_in.insert(v_start);
            self.ds.faces[f1].face_info.vertices_in.insert(v_end);
            self.ds.faces[f2].face_info.vertices_in.insert(v_start);
            self.ds.faces[f2].face_info.vertices_in.insert(v_end);

            curve_indices.push(curve_idx);
        }

        if !curve_indices.is_empty() {
            self.ds.interferences.push(Interference::FaceFace {
                f1,
                f2,
                curves: curve_indices,
                points: vec![],
            });
        }
    }

    fn make_marching_pcurves_with_reapprox(
        &self,
        points: &[DVec3],
        s1: &Surface3,
        s2: &Surface3,
        f1: usize,
        f2: usize,
        t_range: &[f64; 2],
    ) -> (Option<Curve2d>, Option<Curve2d>) {
        let uv_bounds1 = s1.default_domain();
        let uv_bounds2 = s2.default_domain();
        let is_u_periodic1 = matches!(s1, Surface3::Cylinder(_) | Surface3::Sphere(_) | Surface3::Torus(_));
        let is_u_periodic2 = matches!(s2, Surface3::Cylinder(_) | Surface3::Sphere(_) | Surface3::Torus(_));
        let u_per1 = if is_u_periodic1 { Some(std::f64::consts::TAU) } else { None };
        let u_per2 = if is_u_periodic2 { Some(std::f64::consts::TAU) } else { None };

        // Attempt 1: default tolerance
        let pca = inttools::pcurve_derive::polyline_pcurve_by_projection(points, s1);
        let pcb = inttools::pcurve_derive::polyline_pcurve_by_projection(points, s2);

        let valid_a = pca.as_ref().map_or(false, |pc| {
            inttools::pcurve_derive::is_curve_valid_2d(pc)
                && inttools::pcurve_derive::check_pcurve_in_face(pc, *t_range, uv_bounds1, u_per1, None)
        });
        let valid_b = pcb.as_ref().map_or(false, |pc| {
            inttools::pcurve_derive::is_curve_valid_2d(pc)
                && inttools::pcurve_derive::check_pcurve_in_face(pc, *t_range, uv_bounds2, u_per2, None)
        });

        if valid_a && valid_b {
            return (pca, pcb);
        }

        // �?OCCT-aligned:reApprox �?fallback with looser validation.
        // Skip the self-intersection check (is_curve_valid_2d) since polyline
        // pcurves from marching can have V-folds that are geometrically correct.
        let valid_a2 = pca.as_ref().map_or(false, |pc| {
            inttools::pcurve_derive::check_pcurve_in_face(pc, *t_range, uv_bounds1, u_per1, None)
        });
        let valid_b2 = pcb.as_ref().map_or(false, |pc| {
            inttools::pcurve_derive::check_pcurve_in_face(pc, *t_range, uv_bounds2, u_per2, None)
        });
        if valid_a2 && valid_b2 {
            return (pca, pcb);
        }

        // Final fallback: return pcurves even if invalid �?the builder handles
        // out-of-face pcurves via its own boundary clipping.
        (pca, pcb)
    }

    fn generate_surface_samples(&self, surface: &Surface3, n1: usize, n2: usize) -> Vec<DVec3> {
        match surface {
            Surface3::Cylinder(cyl) => {
                inttools::marching::sample_cylinder(cyl, [-20.0, 20.0], n1, n2)
            }
            Surface3::Sphere(sph) => inttools::marching::sample_sphere(sph, n1, n2),
            Surface3::Torus(torus) => inttools::marching::sample_torus(torus, n1, n2),
            Surface3::Plane(plane) => sample_plane(plane, 20.0, n1),
            Surface3::Cone(cone) => sample_cone(cone, 0.01, 20.0, n1, n2),
            // Generic fallback: sample via surface.default_domain() UV grid.
            // Works for BSpline, Bezier, Offset, Revolution, Trimmed, LinearExtrusion.
            _ => sample_surface_generic(surface, n1, n2),
        }
    }

    fn generate_surface_samples_grid(
        &self,
        surface: &Surface3,
        n_u: usize,
        n_v: usize,
    ) -> Vec<DVec3> {
        match surface {
            Surface3::Cylinder(cyl) => {
                // u = azimuth index (0..n_u), v = height index (0..n_v)
                // sample_cylinder returns row = height, col = azimuth,
                // so transpose to row = azimuth, col = height for grid indexing.
                // Rebuild in (n_u azimuth) �?(n_v height) order.
                let height_range = [-20.0_f64, 20.0_f64];
                let u_ax = if cyl.axis.x.abs() < 0.9 {
                    cyl.axis.cross(DVec3::X).normalize()
                } else {
                    cyl.axis.cross(DVec3::Y).normalize()
                };
                let v_ax = cyl.axis.cross(u_ax);
                let mut pts = Vec::with_capacity(n_u * n_v);
                for iu in 0..n_u {
                    let theta =
                        2.0 * std::f64::consts::PI * iu as f64 / n_u as f64;
                    for iv in 0..n_v {
                        let h = height_range[0]
                            + (height_range[1] - height_range[0]) * iv as f64
                                / (n_v - 1).max(1) as f64;
                        pts.push(
                            cyl.origin
                                + cyl.axis * h
                                + (u_ax * theta.cos() + v_ax * theta.sin()) * cyl.radius,
                        );
                    }
                }
                pts
            }
            Surface3::Sphere(sph) => {
                let u_ax = if sph.axis.x.abs() < 0.9 {
                    sph.axis.cross(DVec3::X).normalize()
                } else {
                    sph.axis.cross(DVec3::Y).normalize()
                };
                let v_ax = sph.axis.cross(u_ax);
                let mut pts = Vec::with_capacity(n_u * n_v);
                for iu in 0..n_u {
                    let theta =
                        2.0 * std::f64::consts::PI * iu as f64 / n_u as f64;
                    for iv in 0..n_v {
                        let phi = std::f64::consts::PI * iv as f64 / (n_v - 1).max(1) as f64;
                        pts.push(
                            sph.center
                                + sph.radius
                                    * (sph.axis * phi.cos()
                                        + (u_ax * theta.cos() + v_ax * theta.sin()) * phi.sin()),
                        );
                    }
                }
                pts
            }
            _ => {
                // Fallback: generic UV-grid sampling for BSpline, Bezier, Offset, etc.
                sample_surface_generic(surface, n_u, n_v)
            }
        }
    }

    fn estimate_step_size(&self, s1: &Surface3, s2: &Surface3) -> f64 {
        // Use a fraction of the smallest characteristic dimension
        let size1 = match s1 {
            Surface3::Sphere(s) => s.radius,
            Surface3::Cylinder(c) => c.radius,
            Surface3::Cone(c) => c.radius.max(0.5),
            Surface3::Torus(t) => t.minor_radius,
            Surface3::Ellipsoid(e) => e.radius_x.min(e.radius_y).min(e.radius_z),
            Surface3::Pipe(p) => p.radius,
            Surface3::Plane(_)
            | Surface3::Helicoid(_)
            | Surface3::BSpline(_)
            | Surface3::TriBezier(_)
            | Surface3::LinearExtrusion(_)
            | Surface3::Revolution(_)
            | Surface3::Ruled(_)
            | Surface3::Coons(_)
            | Surface3::Bezier(_)
            | Surface3::Offset(_)
            | Surface3::Trimmed(_) => 1.0,
        };
        let size2 = match s2 {
            Surface3::Sphere(s) => s.radius,
            Surface3::Cylinder(c) => c.radius,
            Surface3::Cone(c) => c.radius.max(0.5),
            Surface3::Torus(t) => t.minor_radius,
            Surface3::Ellipsoid(e) => e.radius_x.min(e.radius_y).min(e.radius_z),
            Surface3::Pipe(p) => p.radius,
            Surface3::Plane(_)
            | Surface3::Helicoid(_)
            | Surface3::BSpline(_)
            | Surface3::TriBezier(_)
            | Surface3::LinearExtrusion(_)
            | Surface3::Revolution(_)
            | Surface3::Ruled(_)
            | Surface3::Coons(_)
            | Surface3::Bezier(_)
            | Surface3::Offset(_)
            | Surface3::Trimmed(_) => 1.0,
        };
        size1.min(size2) * 0.1
    }

    // 鈹€鈹€鈹€ MakeBlocks: inject EF/EE vertices onto FF curves (OCCT PaveFiller_6 L647+) 鈹€鈹€

         fn make_blocks(&mut self) {
        // �?OCCT L652-655: GlueOff guard — MakeBlocks should be skipped for GlueFull/GluePartial
        if self.use_glue {
            return;
        }
        // Phase 1: Collect data 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        let n_curves = self.ds.intersection_curves.len();
        let n_faces = self.ds.faces.len();

        // �?OCCT-aligned: collect EF candidate vertices (EdgeFace interferences).
        //    These are the vertices the IC endpoints should share (edge-face intersections).
        let mut all_verts: Vec<(usize, DVec3)> = self.ds.interferences.iter()
            .filter_map(|inf| {
                if let Interference::EdgeFace { new_vertex, point, .. } = inf {
                    Some((*new_vertex, *point))
                } else { None }
            })
            .collect();
        // Add vertices_in/vertices_on (vertices from IC splitting + face interior vertices)
        let mut seen: std::collections::BTreeSet<usize> =
            all_verts.iter().map(|(vi, _)| *vi).collect();
        // OCCT-aligned: Add EE (Edge-Edge), VE (Vertex-Edge), VF (Vertex-Face)
        // interference vertices to all_verts.
        // OCCT PaveFiller_6.cxx: SubShapesOnIn contains all interference types.
        for inf in &self.ds.interferences {
            match inf {
                Interference::EdgeEdge { new_vertex, point, .. } => {
                    if seen.insert(*new_vertex) {
                        all_verts.push((*new_vertex, *point));
                    }
                }
                Interference::VertexEdge { vertex, .. } => {
                    if seen.insert(*vertex) {
                        all_verts.push((*vertex, self.ds.vertices[*vertex].point));
                    }
                }
                Interference::VertexFace { vertex, .. } => {
                    if seen.insert(*vertex) {
                        all_verts.push((*vertex, self.ds.vertices[*vertex].point));
                    }
                }
                _ => {}
            }
        }
        for face in &self.ds.faces {
            for &vi in &face.face_info.vertices_in {
                if seen.insert(vi) {
                    all_verts.push((vi, self.ds.vertices[vi].point));
                }
            }
            for &vi in &face.face_info.vertices_on {
                if seen.insert(vi) {
                    all_verts.push((vi, self.ds.vertices[vi].point));
                }
            }
        }
        if std::env::var("RCAD_DEBUG_SPLIT").is_ok() {
            eprintln!("[SPLIT_DBG] n_curves={} n_faces={} all_verts={}", n_curves, n_faces, all_verts.len());
            for (vi, pt) in &all_verts {
                eprintln!("[SPLIT_DBG]   all_verts[{}] = ({:.9},{:.9},{:.9})", vi, pt.x, pt.y, pt.z);
            }
        }
        // Note: face boundary vertices (e.g. sphere poles) are NOT added to all_verts,
        // to avoid incorrectly assigning boundary vertex indices to IC endpoints.
        // PutBoundPaveOnCurve reads face.boundary_verts directly from DS, bypassing all_verts.

        // Curve snapshots: collect all data upfront to avoid borrow conflicts
        struct CurveSnapshot {
            curve: Curve3,
            t_range: [f64; 2],
            sv: usize, ev: usize,
            sv_pos: DVec3, ev_pos: DVec3,
            pcurve_on_a: Option<Curve2d>,
            pcurve_on_b: Option<Curve2d>,
        }
        let snapshots: Vec<CurveSnapshot> = (0..n_curves).map(|ci| {
            let ic = &self.ds.intersection_curves[ci];
            CurveSnapshot {
                curve: ic.curve.clone(),
                t_range: ic.t_range,
                sv: ic.start_vertex, ev: ic.end_vertex,
                sv_pos: self.ds.vertices[ic.start_vertex].point,
                ev_pos: self.ds.vertices[ic.end_vertex].point,
                pcurve_on_a: ic.pcurve_on_a.clone(),
                pcurve_on_b: ic.pcurve_on_b.clone(),
            }
        }).collect();

        // Face info snapshots: which faces reference which curves.
        // All intersection curves are stored in curves_sc (PaveBlocksSc).
        let face_curves: Vec<Vec<usize>> = (0..n_faces)
            .map(|fi| self.ds.faces[fi].face_info.curves_sc_only())
            .collect();

        // 鈹€鈹€ Phase 2: Compute splits 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        #[derive(Clone)]
        struct SplitAction {
            old_ci: usize,
            split_verts: Vec<(f64, usize)>,
            pca: Option<Curve2d>,
            pcb: Option<Curve2d>,
        }
        let mut actions: Vec<SplitAction> = Vec::new();
        let mut cur_tol = self.tol(); // fallback for post-loop code

        for ci in 0..n_curves {
            // Use per-pair tolerance from the two faces that reference this curve
            let fi_a = (0..n_faces).find(|&fi| face_curves[fi].contains(&ci)).unwrap_or(0);
            let fi_b = (fi_a + 1..n_faces).find(|&fi| face_curves[fi].contains(&ci)).unwrap_or(fi_a);
            cur_tol = self.ff_tol(fi_a, fi_b);
            let tol = cur_tol;
            let snap = &snapshots[ci];
            let [t0, t1] = snap.t_range;
            let mut on_curve: Vec<(f64, usize)> = Vec::new();

            // OCCT-aligned: PutPaveOnCurve (BOPAlgo_PaveFiller_6.cxx L833-900)
            //    Projects ALL vertices face info onto this IC, adding
            //    EF/EE/VE/VF interference vertices to split the curve.
            //    OCCT PutBoundPaveOnCurve runs AFTER FilterPavesOnCurves (below).
            {
                let face_idxs = find_face_idxs_for_curve(&self.ds, ci);
                if let Curve3::Circle(_) = &snap.curve {
                    if std::env::var("RCAD_DEBUG_SPLIT").is_ok() {
                        eprintln!("[SPLIT_DBG] Circle ci={}", ci);
                    }
                }
                let full_paves = put_pave_on_curve_full(&self.ds, ci, &face_idxs);
                on_curve.extend(full_paves);
                on_curve = filter_paves_on_curves(&self.ds, ci, &on_curve);
                if std::env::var("RCAD_DEBUG_SPLIT").is_ok() {
                    eprintln!("[SPLIT_DBG] after filter ci={} on_curve={}", ci, on_curve.len());
                    for (t, vi) in &on_curve {
                        eprintln!("[SPLIT_DBG]   filtered on_curve t={:.9} vi={}", t, vi);
                    }
                }
                let is_closed = matches!(&snap.curve, Curve3::Circle(_));
                put_closing_pave_on_curve(&mut on_curve, is_closed);
            }

            // �?OCCT-aligned: FilterPavesOnCurves post-filter �?remove vertices not in any face boundary or EF
            // OCCT-aligned: sort on_curve by parameter, ensuring vertices in sp are in increasing parameter order.
            //    OCCT PaveBlock::Update gets sorted Pave list from PutPavesOnCurve.
            on_curve.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

            // �?OCCT-aligned: instead of physically splitting the curve (rcad original),
            //   add interior vertices as ext_paves on the curve's PaveBlock via
            //   put_pave_on_curve (= PutPaveOnCurve in OCCT PaveFiller_6.cxx L833-900).
            //   The PaveBlock's update() will later produce sub-PBs from these ext_paves,
            //   and make_section_edges_from_curve_pbs creates DSEdges from each sub-PB.
            let [t0_live, t1_live] = self.ds.intersection_curves[ci].t_range;
            for &(p, vi) in &on_curve {
                // Skip vertices at or near endpoints (already boundary vertices)
                if (p - t0_live).abs() < tol * 0.1 || (p - t1_live).abs() < tol * 0.1 {
                    continue;
                }
                self.put_pave_on_curve(vi, ci);
            }
        }
        // �?OCCT-aligned: ext_paves have been added via put_pave_on_curve above.
        //   make_section_edges_from_curve_pbs will later call update() on each
        //   curve's PB to produce sub-PBs (OCCT PaveFiller_6.cxx L882-980).
        //   Sub-PBs that are too short (no valid range) will be excluded there.

        // �?OCCT-aligned: FilterPavesOnCurves �?cross-curve vertex deduplication
        //   OCCT L796, L2349-2443: for a vertex on multiple curves, keep the best matching curve.
        {
            let a_sin_angle_min: f64 = 0.5;
            let mut vert_curves: std::collections::HashMap<usize, Vec<(usize, bool)>> = std::collections::HashMap::new();
            let all_curve_ids: Vec<usize> = (0..self.ds.intersection_curves.len()).collect();

            for &ci in &all_curve_ids {

                let ic = &self.ds.intersection_curves[ci];

                vert_curves.entry(ic.start_vertex).or_default().push((ci, true));

                if ic.start_vertex != ic.end_vertex {

                    vert_curves.entry(ic.end_vertex).or_default().push((ci, false));

                }

            }

            // 2. Only process vertices appearing on multiple curves

            // OCCT-aligned: PaveBlockDist { ci, n_v, sq_dist, sin_angle }
            struct PaveBlockDist {
                ci: usize,
                n_v: usize,
                sq_dist: f64,
                sin_angle: f64,
            }
            let mut vert_pbs: std::collections::HashMap<usize, Vec<PaveBlockDist>> = std::collections::HashMap::new();

            for (n_v, curves) in &vert_curves {

                if curves.len() < 2 { continue; }

                for &(ci, is_start) in curves {

                    let ic = &self.ds.intersection_curves[ci];

                    let par = if is_start { ic.t_range[0] } else { ic.t_range[1] };

                    let pt = self.ds.vertices[*n_v].point;

                    let curve_pt = ic.curve.point_at(par);

                    let tangent = ic.curve.tangent_at(par);

                    let to_curve = curve_pt - pt;

                    let sq_dist = to_curve.length_squared();

                    let speed_sq = tangent.length_squared();

                    let sin_angle = if sq_dist > 1e-30 && speed_sq > 1e-30 {

                        (to_curve.cross(tangent).length_squared() / (sq_dist * speed_sq)).sqrt()

                    } else { 0.0 };

                    vert_pbs.entry(*n_v).or_default().push(PaveBlockDist { ci, n_v: *n_v, sq_dist, sin_angle });

                }

            }

            // OCCT-aligned: For each vertex, keep the best-matching curve, remove ext paves from others
            //   (aPBD.PB->RemoveExtPave(nV))
            for (n_v, dists) in &vert_pbs {

                let min_dist = dists.iter().map(|d| d.sq_dist).min_by(|a,b| a.partial_cmp(b).unwrap()).unwrap_or(0.0);

                for d in dists {

                    let check_dist = 100.0 * cur_tol.max(min_dist);

                    if d.sq_dist > check_dist && d.sin_angle < a_sin_angle_min {

                        if let Some(pb) = self.ds.intersection_curves[d.ci].change_pave_block1() {
                            pb.remove_ext_pave(d.n_v);
                        }

                    }

                }

            }

        }

        // �?OCCT-aligned: PutBoundPaveOnCurve after FilterPavesOnCurves
        //   OCCT PaveFiller_6.cxx L798-832: PutBoundPaveOnCurve runs AFTER
        //   FilterPavesOnCurves, and creates bound vertices at IC endpoints
        //   that pass IsValidPointForFaces.  rcad: start_vertex/end_vertex are
        //   always assigned by the analytic IC creators, but the IsValidPointForFaces
        //   check still filters spurious IC endpoints.
        for ci in 0..self.ds.intersection_curves.len() {
            let face_idxs = find_face_idxs_for_curve(&self.ds, ci);
            if face_idxs[0] == usize::MAX || face_idxs[1] == usize::MAX { continue; }
            let ic = &self.ds.intersection_curves[ci];
            let tol = self.ff_tol(face_idxs[0], face_idxs[1]);
            // OCCT: aIC.Bounds(aT[0], aT[1], aP[0], aP[1])
            let a_t = ic.t_range;
            let a_p = [ic.curve.point_at(a_t[0]), ic.curve.point_at(a_t[1])];
            let a_pb = &self.ds.intersection_curves[ci].pave_blocks;
            let a_pb1 = a_pb.first().cloned();
            // OCCT: getBoundPaves �?check existing vertex at each end
            let a_sv = ic.start_vertex;
            let a_ev = ic.end_vertex;
            let sv_pos = self.ds.vertices[a_sv].point;
            let ev_pos = self.ds.vertices[a_ev].point;
            let sv_dist = sv_pos.distance_squared(a_p[0]);
            let ev_dist = ev_pos.distance_squared(a_p[1]);
            // OCCT L2324: isClosed check
            let is_closed = a_p[0].distance_squared(a_p[1]) < TOLERANCE_ABS_SQ;
            if is_closed { continue; }
            // OCCT L2340: IsValidPointForFaces for both endpoints, pre-computed to avoid borrow conflict
            let end_tol = tol.max(TOLERANCE_ABS) * 10.0;
            let end0_valid = {
                let s0 = &self.ds.faces[face_idxs[0]].surface;
                let (_, p0) = crate::extrema::closest_point_on_surface(s0, a_p[0]);
                let s1 = &self.ds.faces[face_idxs[1]].surface;
                let (_, p1) = crate::extrema::closest_point_on_surface(s1, a_p[0]);
                p0.distance(a_p[0]) < end_tol && p1.distance(a_p[0]) < end_tol
            };
            let end1_valid = {
                let s0 = &self.ds.faces[face_idxs[0]].surface;
                let (_, p0) = crate::extrema::closest_point_on_surface(s0, a_p[1]);
                let s1 = &self.ds.faces[face_idxs[1]].surface;
                let (_, p1) = crate::extrema::closest_point_on_surface(s1, a_p[1]);
                p0.distance(a_p[1]) < end_tol && p1.distance(a_p[1]) < end_tol
            };
            for j in 0..2 {
                let a_pt = a_p[j];
                let a_tj = if j == 0 { a_t[0] } else { a_t[1] };
                let cur_dist = if j == 0 { sv_dist } else { ev_dist };
                let end_valid = if j == 0 { end0_valid } else { end1_valid };
                // OCCT L2332: if (aBndNV[j] < 0) �?vertex on this end needs creation
                //   rcad: skip if start_vertex matches endpoint position
                if cur_dist < tol * tol { continue; }
                // OCCT L2340: IsValidPointForFaces
                if !end_valid { continue; }
                // OCCT: MakeNewVertex �?append to DS �?AppendExtPave
                let n_vn = self.ds.add_vertex(a_pt);
                let parent_tol = self.ds.faces[face_idxs[0]].geom_tol.max(self.ds.faces[face_idxs[1]].geom_tol).max(self.seam_shift_tol);
                self.ds.vertices[n_vn].geom_tol = parent_tol;
                if let Some(pb) = self.ds.intersection_curves[ci].change_pave_block1() {
                    pb.append_ext_pave1(Pave { vertex_idx: n_vn, param: a_tj });
                }
            }
        }

        // �?OCCT-aligned: MakeSectionEdges from curve PaveBlocks (PaveFiller_6.cxx L882-980).
        //   Creates DSEdges from curve PBs split by ext_paves.
        self.make_section_edges_from_curve_pbs();

        // �?OCCT-aligned: Build edge images from pave blocks (FillImagesEdges)
        self.ds.build_edge_images();

        if std::env::var("RCAD_DEBUG_SPLIT").is_ok() {
            let n_circle = self.ds.intersection_curves.iter().filter(|ic| matches!(ic.curve, Curve3::Circle(_))).count();
            let n_total = self.ds.intersection_curves.len();
            eprintln!("[SPLIT] END_MAKE_BLOCKS total_curves={} circle_curves={}", n_total, n_circle);
            for fi in 0..self.ds.faces.len() {
                let face = &self.ds.faces[fi];
                eprintln!("[SPLIT]   face[{}] curves_sc={} vertices_in={}", fi, face.face_info.curves_sc.len(), face.face_info.vertices_in.len());
            }
        }
        // �?OCCT-aligned: InitPaveBlock1 for all curves (PaveFiller_6.cxx L800).
        //   Curves created by splitting also need initial PaveBlocks.
        for ci in 0..self.ds.intersection_curves.len() {
            self.ds.intersection_curves[ci].init_pave_block1();
        }
    }

    fn check_self_interference(&self) -> Result<(), String> {
        let mut messages: Vec<String> = Vec::new();

        for interference in &self.ds.interferences {
            match interference {
                Interference::EdgeEdge { e1, e2, .. } => {
                    let origin1 = self.ds.edges[*e1].origin;
                    let origin2 = self.ds.edges[*e2].origin;
                    if origin1 == origin2 {
                        messages.push(format!(
                            "  EdgeEdge(e1={}, e2={}) both from {:?}",
                            e1, e2, origin1
                        ));
                    }
                }
                Interference::EdgeFace { edge, face, .. } => {
                    let edge_origin = self.ds.edges[*edge].origin;
                    let face_origin = self.ds.faces[*face].origin;
                    if edge_origin == face_origin {
                        messages.push(format!(
                            "  EdgeFace(edge={}, face={}) both from {:?}",
                            edge, face, edge_origin
                        ));
                    }
                }
                Interference::FaceFace { f1, f2, .. } => {
                    let origin1 = self.ds.faces[*f1].origin;
                    let origin2 = self.ds.faces[*f2].origin;
                    if origin1 == origin2 {
                        messages.push(format!(
                            "  FaceFace(f1={}, f2={}) both from {:?} (non-fatal)",
                            f1, f2, origin1
                        ));
                    }
                }
                _ => {}
            }
        }

        if messages.is_empty() {
            Ok(())
        } else {
            Err(format!("Self-interference detected:\n{}", messages.join("\n")))
        }
    }

    fn make_split_edges(&mut self) {
        self.build_split_edges();
    }

    fn build_split_edges(&mut self) {
        // OCCT L392: UpdateCommonBlocksWithSDVertices �?before creating split edges,
        //   ensure CommonBlocks reference correct (SD-deduplicated) vertex indices.
        self.ds.update_common_blocks_with_sd_vertices();

        // Phase 1: collect PaveBlock data without creating new edges (avoids
        // mutable borrow conflict with self.ds.edges iteration).
        struct BlockData {
            ei: usize,
            sv: usize, ev: usize,
            t_start: f64, t_end: f64,
            curve: Curve3,
            origin: ShapeOrigin,
            geom_tol: f64,
            face_reps: Vec<DSRepOnFace>,
        }
        let mut all_blocks: Vec<BlockData> = Vec::new();
        let n_orig_edges = self.ds.edges.len();

        // �?OCCT-aligned: MakeSplitEdges (PaveFiller_7.cxx) only creates split
        //    edges and sets PaveBlock->Edge() (pb.new_edge).  rcad also initializes
        //    pave_blocks on source edges here so downstream FillImagesEdges can
        //    read pb.new_edge.  my_images / my_origins are NOT populated here �?
        //    that is FillImagesEdges' responsibility (build_edge_images in ds.rs).

        for ei in 0..n_orig_edges {
            let edge = &self.ds.edges[ei];
            if edge.paves.is_empty() {
                // �?OCCT-aligned: no split �?edge stays as-is (no PaveBlock created).
                //   OCCT FillImagesEdges requires HasReference (non-empty pave_blocks)
                //   to create split images.  Empty pave_blocks = un-split edge =
                //   passes through BuildResult unchanged.
                continue;
            }

            let mut all_paves = vec![
                Pave { vertex_idx: edge.start_vertex, param: edge.t_range[0] },
                Pave { vertex_idx: edge.end_vertex, param: edge.t_range[1] },
            ];
            all_paves.extend_from_slice(&edge.paves);
            all_paves.sort_by(|a, b| a.param.partial_cmp(&b.param).unwrap_or(std::cmp::Ordering::Equal));
            all_paves.dedup_by(|a, b| params_equal(a.param, b.param));

            for w in all_paves.windows(2) {
                let pb = PaveBlock::new(ei, w[0], w[1]);
                let t1 = pb.pave1.param;
                let t2 = pb.pave2.param;
                let (t_start, t_end) = if t1 < t2 { (t1, t2) } else { (t2, t1) };
                let split_curve = pb.curve.clone().unwrap_or_else(|| edge.curve.clone());
                all_blocks.push(BlockData {
                    ei,
                    sv: pb.pave1.vertex_idx,
                    ev: pb.pave2.vertex_idx,
                    t_start, t_end,
                    curve: split_curve,
                    origin: edge.origin,
                    geom_tol: edge.geom_tol,
                    face_reps: edge.face_reps.clone(),
                });
            }
        }

        // Phase 2: create new DSEdges for each collected block + set pave_blocks
        // on source edges (MakeSplitEdges).  my_images / my_origins are NOT
        // populated here �?that is FillImagesEdges' job (build_edge_images in ds.rs).
        let mut edge_pbs: std::collections::HashMap<usize, Vec<(usize, usize, f64, f64, usize)>> =
            std::collections::HashMap::new();

        for data in &all_blocks {
            let new_ei = self.ds.edges.len();
            self.ds.edges.push(DSEdge {
                start_vertex: data.sv,
                end_vertex: data.ev,
                curve: data.curve.clone(),
                t_range: [data.t_start, data.t_end],
                origin: data.origin,
                geom_tol: data.geom_tol,
                paves: vec![],
                pave_blocks: vec![],
                face_reps: data.face_reps.clone(),
                is_internal: false,
            });

            // Track for pave_blocks assignment on source edge
            edge_pbs.entry(data.ei).or_default().push((
                data.sv, data.ev, data.t_start, data.t_end, new_ei,
            ));
        }

        // �?OCCT-aligned: Set pave_blocks on source edges that were split,
        //    so Builder::fill_images_edges can read pb.new_edge.
        for (ei, blocks) in &edge_pbs {
            let pbs: Vec<PaveBlock> = blocks.iter().map(|&(sv, ev, t_start, t_end, new_ei)| {
                let mut pb = PaveBlock::new(*ei,
                    Pave { vertex_idx: sv, param: t_start },
                    Pave { vertex_idx: ev, param: t_end },
                );
                pb.new_edge = Some(new_ei);
                pb
            }).collect();
            self.ds.edges[*ei].pave_blocks = pbs;
        }
    }

    // 鈹€鈹€鈹€ Helpers 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

    fn verts_of(&self, origin: ShapeOrigin) -> Vec<usize> {
        self.ds
            .vertices
            .iter()
            .enumerate()
            .filter(|(_, v)| v.origin == Some(origin))
            .map(|(i, _)| i)
            .collect()
    }

    fn edges_of(&self, origin: ShapeOrigin) -> Vec<usize> {
        self.ds
            .edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.origin == origin)
            .map(|(i, _)| i)
            .collect()
    }

    fn faces_of(&self, origin: ShapeOrigin) -> Vec<usize> {
        self.ds
            .faces
            .iter()
            .enumerate()
            .filter(|(_, f)| f.origin == origin)
            .map(|(i, _)| i)
            .collect()
    }
}

// 鈹€鈹€ Phase 2a helpers: vertex �?curve parameter projection 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

fn param_on_line3(pt: DVec3, line: &Line3, tol: f64) -> Option<f64> {
    let dir = line.direction;
    let to_pt = pt - line.origin;
    let t = to_pt.dot(dir);
    let proj = line.origin + t * dir;
    let dist = proj.distance(pt);
    if dist > tol { None } else { Some(t) }
}

fn param_on_circle3(pt: DVec3, circle: &Circle3, tol: f64) -> Option<f64> {
    let r = circle.radius;
    let center = circle.center;
    let normal = circle.normal;
    // �?OCCT-aligned: point must be on the circle's plane (Geom_Circle::Value natural requirement)
    let local = pt - center;
    if local.dot(normal).abs() > tol {
        return None;
    }
    let dist_to_center = local.length();
    if (dist_to_center - r).abs() > tol {
        return None;
    }
    let ref_dir = any_perpendicular(normal);
    let u = ref_dir.normalize();
    let v = normal.cross(u);
    let x = local.dot(u);
    let y = local.dot(v);
    Some(y.atan2(x))
}

fn project_vertex_to_curve(pt: DVec3, curve: &Curve3, tol: f64) -> Option<f64> {
    match curve {
        Curve3::Line(line) => param_on_line3(pt, line, tol),
        Curve3::Circle(circ) => param_on_circle3(pt, circ, tol),
        Curve3::Ellipse(ell) => param_on_ellipse3(pt, ell, tol),
        _ => param_on_curve3_numeric(pt, curve, tol),
    }
}

fn param_on_ellipse3(pt: DVec3, ellipse: &Ellipse3, tol: f64) -> Option<f64> {
    use rcad_kernel::geom::CurveEval;
    let local = pt - ellipse.center;
    if local.dot(ellipse.normal).abs() > tol { return None; }
    let y_ax = ellipse.normal.cross(ellipse.major_dir).normalize();
    let n_sample = 64usize;
    let mut best_t = 0.0f64;
    let mut best_dist = f64::INFINITY;
    for i in 0..n_sample {
        let t = std::f64::consts::TAU * i as f64 / n_sample as f64;
        let p = ellipse.point_at(t);
        let d = p.distance_squared(pt);
        if d < best_dist { best_dist = d; best_t = t; }
    }
    if best_dist.sqrt() > tol * 100.0 { return None; }
    // Newton: f(t) = (C(t)-P)路C'(t)=0
    for _ in 0..6 {
        let (ct, st) = (best_t.cos(), best_t.sin());
        let c = ellipse.center + ellipse.major_radius * ct * ellipse.major_dir
            + ellipse.minor_radius * st * y_ax;
        let d = c - pt;
        let cp = -ellipse.major_radius * st * ellipse.major_dir
            + ellipse.minor_radius * ct * y_ax;
        let f = d.dot(cp);
        let cpp = -ellipse.major_radius * ct * ellipse.major_dir
            - ellipse.minor_radius * st * y_ax;
        let fp = cp.dot(cp) + d.dot(cpp);
        if fp.abs() < 1e-15 { break; }
        best_t -= f / fp;
        if (f/fp).abs() < 1e-12 { break; }
    }
    best_t = best_t.rem_euclid(std::f64::consts::TAU);
    if ellipse.point_at(best_t).distance(pt) > tol { None } else { Some(best_t) }
}

fn param_on_curve3_numeric(pt: DVec3, curve: &Curve3, tol: f64) -> Option<f64> {
    use rcad_kernel::geom::CurveEval;
    let (t0, t1) = match curve {
        Curve3::BSpline(bsp) => { let k = &bsp.knots;
            if k.len() >= 2 { (k[0], k[k.len()-1]) } else { (0.0, 1.0) } }
        _ => (0.0, 1.0),
    };
    if (t1 - t0).abs() < 1e-15 { return None; }
    let n_sample = 128usize;
    let mut best_t = t0; let mut best_dist = f64::INFINITY;
    for i in 0..n_sample {
        let t = t0 + (t1 - t0) * i as f64 / (n_sample - 1) as f64;
        let p = curve.point_at(t);
        let d = p.distance_squared(pt);
        if d < best_dist { best_dist = d; best_t = t; }
    }
    if best_dist.sqrt() > tol * 100.0 { return None; }
    for _ in 0..6 {
        let c = curve.point_at(best_t);
        let cp = curve.tangent_at(best_t);
        let f = (c - pt).dot(cp);
        let eps = 1e-6_f64.max((best_t - t0).abs() * 1e-6);
        let t_eps = (best_t + eps).min(t1);
        let f2 = (curve.point_at(t_eps) - pt).dot(curve.tangent_at(t_eps));
        let fp = (f2 - f) / (t_eps - best_t);
        if fp.abs() < 1e-15 { break; }
        best_t = (best_t - f / fp).clamp(t0, t1);
        if (f/fp).abs() < 1e-12 { break; }
    }
    if curve.point_at(best_t).distance(pt) > tol { None } else { Some(best_t) }
}

// 鈹€鈹€ FindValidRange / ShrunkData helper functions 鈹€鈹€
// OCCT references:
//   IntTools_ShrunkRange::Perform()  IntTools_ShrunkRange.cxx L107-191
//   BRepLib::FindValidRange          BRepLib_1.cxx L173-258
//   findNearestValidPoint            BRepLib_1.cxx L31-148

fn curve_resolution(curve: &Curve3, t: f64, tol: f64) -> f64 {

    use rcad_kernel::geom::CurveEval;

    let speed = curve.tangent_at(t).length();

    if speed < 1e-15 { tol } else { tol / speed }

}

fn find_nearest_valid_point(

    curve: &Curve3, first: f64, last: f64, is_first: bool,

    vert_pt: DVec3, vert_tol: f64, eps: f64,

) -> Option<f64> {

    use rcad_kernel::geom::CurveEval;

    let (start_u, end_u) = if is_first { (first, last) } else { (last, first) };

    let tol_sq = vert_tol * vert_tol;

    // 1. Check if endpoint is inside tolerance sphere

    if curve.point_at(start_u).distance_squared(vert_pt) > tol_sq { return None; }

    // 2. Step until outside tolerance sphere

    let step = curve_resolution(curve, start_u, vert_tol).max(eps);

    let step = if is_first { step } else { -step };

    let (mut u_in, mut u_out) = (start_u, start_u);

    loop {

        u_in = u_out; u_out += step;

        if (is_first && u_out > end_u) || (!is_first && u_out < end_u) {

            if curve.point_at(end_u).distance_squared(vert_pt) <= tol_sq { return None; }

            u_out = end_u; break;

        }

        if curve.point_at(u_out).distance_squared(vert_pt) > tol_sq { break; }

    }

    // 3. Bisection refinement

    while (u_out - u_in).abs() > eps {

        let mid = (u_in + u_out) * 0.5;

        if curve.point_at(mid).distance_squared(vert_pt) > tol_sq {

            u_out = mid;

        } else { u_in = mid; }

    }

    Some(if is_first { u_out } else { u_in })

}

fn find_valid_range(

    curve: &Curve3, t0: f64, t1: f64,

    sv_pt: DVec3, sv_tol: f64, ev_pt: DVec3, ev_tol: f64,

) -> Option<(f64, f64)> {

    use rcad_kernel::geom::CurveEval;

    if (t1 - t0).abs() < 1e-12 { return None; }

    let abs_max = t0.abs().max(t1.abs()).max(1.0);

    let eps = curve_resolution(curve, (t0+t1)*0.5, 1e-7).max(abs_max * 1e-12).max(1e-12);

    // Start point shrunk

    let first = if t0.is_infinite() { t0 } else {

        match find_nearest_valid_point(curve, t0, t1, true, sv_pt, sv_tol, eps) {

            Some(f) => { if (t1 - f).abs() < eps { return None; } f }

            None => { return None; }

        }

    };

    // End point shrunk

    let last = if t1.is_infinite() { t1 } else {

        match find_nearest_valid_point(curve, t0, t1, false, ev_pt, ev_tol, eps) {

            Some(l) => { if (l - t0).abs() < eps { return None; } l }

            None => { return None; }

        }

    };

    if first > last { None } else { Some((first, last)) }

}

// 鈹€鈹€ Seam Edge Shift Struct 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

struct SeamEdgeShift {
    shift_vector: DVec3,
    shift_value: f64,
    shifted_face: u8,
}

// 鈹€鈹€ Free Helper Functions 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

fn apply_shift_to_surface(surface: &Surface3, shift: DVec3) -> Surface3 {
    match *surface {
        Surface3::Plane(p) => Surface3::Plane(Plane {
            origin: p.origin + shift,
            ..p
        }),
        Surface3::Cylinder(c) => Surface3::Cylinder(CylindricalSurface {
            origin: c.origin + shift,
            ..c
        }),
        Surface3::Sphere(s) => Surface3::Sphere(SphericalSurface {
            center: s.center + shift,
            ..s
        }),
        Surface3::Torus(t) => Surface3::Torus(ToroidalSurface {
            center: t.center + shift,
            ..t
        }),
        Surface3::Cone(c) => Surface3::Cone(ConicalSurface {
            apex: c.apex + shift,
            ..c
        }),
        Surface3::BSpline(ref bs) => {
            let mut bs = bs.clone();
            for row in &mut bs.control_points {
                for cp in row {
                    *cp += shift;
                }
            }
            Surface3::BSpline(bs)
        }
        Surface3::Bezier(ref bz) => {
            let mut bz = bz.clone();
            for row in &mut bz.control_points {
                for cp in row {
                    *cp += shift;
                }
            }
            Surface3::Bezier(bz)
        }
        Surface3::LinearExtrusion(ref le) => {
            let mut le = le.clone();
            le.direction = le.direction; // direction unchanged
            // The profile curve's origin is not directly accessible as a field;
            // clone without position modification for now
            Surface3::LinearExtrusion(le)
        }
        ref other => other.clone(),
    }
}

fn translate_curve3(curve: &Curve3, shift: DVec3) -> Curve3 {
    match *curve {
        Curve3::Line(l) => Curve3::Line(Line3 {
            origin: l.origin + shift,
            ..l
        }),
        Curve3::Circle(c) => Curve3::Circle(Circle3 {
            center: c.center + shift,
            ..c
        }),
        Curve3::Ellipse(e) => Curve3::Ellipse(Ellipse3 {
            center: e.center + shift,
            ..e
        }),
        Curve3::BSpline(ref bs) => {
            let mut bs = bs.clone();
            for cp in &mut bs.control_points {
                *cp += shift;
            }
            Curve3::BSpline(bs)
        }
        Curve3::Bezier(ref bz) => {
            let mut bz = bz.clone();
            for cp in &mut bz.control_points {
                *cp += shift;
            }
            Curve3::Bezier(bz)
        }
        Curve3::Hyperbola(h) => Curve3::Hyperbola(Hyperbola3 {
            center: h.center + shift,
            ..h
        }),
        Curve3::Parabola(p) => Curve3::Parabola(Parabola3 {
            vertex: p.vertex + shift,
            ..p
        }),
        Curve3::Offset(ref o) => {
            let mut o = o.clone();
            o.basis = Box::new(translate_curve3(&o.basis, shift));
            Curve3::Offset(o)
        }
        Curve3::CircularHelix(ref h) => {
            let mut h = h.clone();
            h.origin += shift;
            Curve3::CircularHelix(h)
        }
        Curve3::SineWave(ref sw) => {
            let mut sw = sw.clone();
            sw.origin += shift;
            Curve3::SineWave(sw)
        }
    }
}
#[cfg(test)]
mod phase2a_tests {
    use super::*;
    use crate::tolerance::*;
    use rcad_kernel::geom::any_perpendicular;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn test_param_on_line3() {
        let line = Line3 { origin: DVec3::ZERO, direction: DVec3::X };
        let pt = DVec3::new(3.0, 0.0, 0.0);
        let t = param_on_line3(pt, &line, 1e-6).unwrap();
        assert!((t - 3.0).abs() < 1e-6, "expected 3.0, got {}", t);

        // Point not on the line
        let off = DVec3::new(3.0, 1.0, 0.0);
        assert!(param_on_line3(off, &line, 1e-6).is_none());
    }

    #[test]
    fn test_param_on_circle3() {
        let circle = Circle3 { center: DVec3::ZERO, normal: DVec3::Z, radius: 1.0 };
        // Point (1,0,0) �?angle 0
        let pt = DVec3::new(1.0, 0.0, 0.0);
        let t = param_on_circle3(pt, &circle, 1e-6).unwrap();
        assert!(t < 1e-6 || (t - 2.0 * PI).abs() < 1e-6,
            "expected ~0 or 2�? got {}", t);

        // Point (0,1,0) �?angle �?2
        let pt2 = DVec3::new(0.0, 1.0, 0.0);
        let t2 = param_on_circle3(pt2, &circle, 1e-6).unwrap();
        assert!((t2 - FRAC_PI_2).abs() < 1e-6,
            "expected �?2, got {}", t2);

        // Point not on the circle
        let off = DVec3::new(2.0, 0.0, 0.0);
        assert!(param_on_circle3(off, &circle, 1e-6).is_none());
    }
}

// 鈹€鈹€ Phase 2a: MakeBlocks candidate injection helpers 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

fn find_face_idxs_for_curve(ds: &DS, ci: usize) -> [usize; 2] {
    let mut result = [usize::MAX; 2];
    let mut idx = 0;
    for (fi, face) in ds.faces.iter().enumerate() {
        if face.face_info.curves_sc.contains(&ci)
        {
            if idx < 2 {
                result[idx] = fi;
                idx += 1;
            }
        }
    }
    result
}

fn put_bound_pave_on_curve(
    _ds: &DS,
    _curve_idx: usize,
    _face_idxs: &[usize; 2],
    _tol: f64,
) -> Vec<(f64, usize)> {
    // �?OCCT-aligned: PutBoundPaveOnCurve only processes the TWO curve endpoint
    //   positions (aP[0], aP[1]).  In rcad, start_vertex/end_vertex are always set,
    //   so bound vertex injection is a no-op.  OCCT does NOT iterate boundary_verts.
    vec![]
}

fn put_pave_on_curve_full(
    ds: &DS,
    curve_idx: usize,
    face_idxs: &[usize; 2],
) -> Vec<(f64, usize)> {
    let ic = &ds.intersection_curves[curve_idx];
    let [t0, t1] = ic.t_range;
    let a_tol_r3d = ic.geom_tol; // OCCT L2384: max(theNC.Tolerance(), theNC.TangentialTolerance())
    let mut paves: Vec<(f64, usize)> = Vec::new();

    // OCCT-aligned: compute curve bounding box for vertex filtering (L2409: aBoxC.IsOut(aBoxV)).
    let curve_bbox = curve_bounding_box_simple(&ic.curve, a_tol_r3d);

    // OCCT-aligned: GetStickVertices (PaveFiller_6.cxx L2847-2905) collects EF vertex set
    //   per FF pair.  Only EF vertices belonging to this specific pair are added to aMVEF.
    //   rcad: filter EF vertices by checking if the interference's face is in this pair.
    let ef_vertices: std::collections::HashSet<usize> = ds.interferences.iter()
        .filter_map(|inf| {
            if let Interference::EdgeFace { new_vertex, face, .. } = inf {
                // OCCT L2896: aMI.Contains(nS1) && aMI.Contains(nS2)
                //   Both sub-shapes belong to the two faces — the EF vertex involves this pair.
                if *face == face_idxs[0] || *face == face_idxs[1] {
                    Some(*new_vertex)
                } else { None }
            } else { None }
        })
        .collect();

    for &fi in face_idxs.iter().filter(|&&fi| fi != usize::MAX) {
        let face = &ds.faces[fi];

        // OCCT L2386-2392: EF vertices first — only vertices belonging to this FF pair
        for &vi in &face.face_info.vertices_on {
            if !ef_vertices.contains(&vi) { continue; } // OCCT GetStickVertices: skip non-pair EF
            if vi == ic.start_vertex || vi == ic.end_vertex { continue; }
            if paves.iter().any(|&(_, v)| v == vi) { continue; }
            let pt = ds.vertices[vi].point;
            if let Some(t) = project_vertex_to_curve(pt, &ic.curve, a_tol_r3d) {
                if t >= t0 - a_tol_r3d && t <= t1 + a_tol_r3d {
                    paves.push((t, vi));
                }
            }
        }

        // OCCT L2394-2420: ON/IN vertices with BBox + IsNewShape filtering
        for &vi in &face.face_info.vertices_in {
            if vi == ic.start_vertex || vi == ic.end_vertex { continue; }
            // OCCT L2399-2401: skip ON/IN vertices already in EF set
            if ef_vertices.contains(&vi) { continue; }
            if paves.iter().any(|&(_, v)| v == vi) { continue; }

            // OCCT L2404-2412: BBox filtering
            if let Some([c_min, c_max]) = curve_bbox {
                let v_pt = ds.vertices[vi].point;
                let v_tol = ds.vertices[vi].geom_tol.max(a_tol_r3d);
                let v_min = v_pt - DVec3::splat(v_tol);
                let v_max = v_pt + DVec3::splat(v_tol);
                if v_max.x < c_min.x || v_min.x > c_max.x ||
                   v_max.y < c_min.y || v_min.y > c_max.y ||
                   v_max.z < c_min.z || v_min.z > c_max.z {
                    continue;
                }
            }

            // OCCT L2413-2415: IsNewShape filter - skip non-new vertices
            if !ds.is_new_vertex(vi) {
                continue;
            }

            let pt = ds.vertices[vi].point;
            if let Some(t) = project_vertex_to_curve(pt, &ic.curve, a_tol_r3d) {
                if t >= t0 - a_tol_r3d && t <= t1 + a_tol_r3d {
                    paves.push((t, vi));
                }
            }
        }
    }

    // Sort by parameter
    paves.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    // Deduplicate by parameter or vertex idx
    paves.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-12 || a.1 == b.1);

    // PutClosingPaveOnCurve for closed curves
    put_closing_pave_on_curve(&mut paves, matches!(&ic.curve, Curve3::Circle(_)));

    paves
}

fn curve_bounding_box_simple(curve: &Curve3, tol: f64) -> Option<[DVec3; 2]> {
    let bbox = match curve {
        Curve3::Line(_) => {
            // Lines are infinite; use unit-length box centered at origin.
            Some([DVec3::splat(-1.0), DVec3::splat(1.0)])
        }
        Curve3::Circle(c) => {
            let n = c.normal.normalize();
            let extent = DVec3::new(
                c.radius * (1.0 - n.x * n.x).sqrt(),
                c.radius * (1.0 - n.y * n.y).sqrt(),
                c.radius * (1.0 - n.z * n.z).sqrt(),
            );
            Some([c.center - extent, c.center + extent])
        }
        Curve3::Ellipse(e) => {
            let n = e.normal.normalize();
            let max_r = e.major_radius.max(e.minor_radius);
            let extent = DVec3::new(
                max_r * (1.0 - n.x * n.x).sqrt(),
                max_r * (1.0 - n.y * n.y).sqrt(),
                max_r * (1.0 - n.z * n.z).sqrt(),
            );
            Some([e.center - extent, e.center + extent])
        }
        Curve3::BSpline(b) => {
            let mut mn = DVec3::splat(f64::INFINITY);
            let mut mx = DVec3::splat(f64::NEG_INFINITY);
            for &p in &b.control_points {
                mn = mn.min(p);
                mx = mx.max(p);
            }
            if mn.is_finite() { Some([mn, mx]) } else { None }
        }
        Curve3::Bezier(b) => {
            let mut mn = DVec3::splat(f64::INFINITY);
            let mut mx = DVec3::splat(f64::NEG_INFINITY);
            for &p in &b.control_points {
                mn = mn.min(p);
                mx = mx.max(p);
            }
            if mn.is_finite() { Some([mn, mx]) } else { None }
        }
        _ => None,
    };
    // Expand by tolerance (OCCT: aBox.Enlarge(aBoxExpandValue))
    bbox.map(|[mn, mx]| [mn - DVec3::splat(tol), mx + DVec3::splat(tol)])
}

fn filter_paves_on_curves(
    ds: &DS,
    curve_idx: usize,
    paves: &[(f64, usize)],
) -> Vec<(f64, usize)> {
    let ic = &ds.intersection_curves[curve_idx];
    // OCCT-aligned: curve tolerance + fuzzy (SUM matching PutPaveOnCurve L2976 aTolR3D + myFuzzyValue)
    let tol = ic.geom_tol + ds.fuzzy_tol;
    let tol_sq = tol * tol;
    paves.iter().filter(|&&(_, vi)| {
        let pt = ds.vertices[vi].point;
        let dist_sq = match &ic.curve {
            Curve3::Line(line) => {
                let to_pt = pt - line.origin;
                let proj = line.origin + line.direction * to_pt.dot(line.direction);
                proj.distance_squared(pt)
            }
            Curve3::Circle(circ) => {
                let center_dist = pt.distance(circ.center);
                (center_dist - circ.radius).powi(2)
            }
            _ => 0.0,
        };
        dist_sq < tol_sq
    }).copied().collect()
}

fn put_closing_pave_on_curve(
    paves: &mut Vec<(f64, usize)>,
    is_closed: bool,
) {
    if paves.len() < 2 { return; }
    if is_closed {
        let first_t = paves[0].0;
        let last_t = paves[paves.len() - 1].0;
        let span = last_t - first_t;
        // Only replace if the curve spans at least one full period (�?2�?for circles)
        if (span - std::f64::consts::TAU).abs() < 0.1 {
            let first_vi = paves[0].1;
            let last_idx = paves.len() - 1;
            paves[last_idx].1 = first_vi;
        }
    }
}

fn intersect_line_line(
    l1: &Line3,
    r1: [f64; 2],
    l2: &Line3,
    r2: [f64; 2],
    coincidence_tol: f64,
) -> Option<(f64, f64, DVec3)> {
    let tol = coincidence_tol.max(TOLERANCE_ABS);
    let tol_sq = tol * tol;
    let d1 = l1.direction;
    let d2 = l2.direction;
    let w0 = l1.origin - l2.origin;

    let a = d1.dot(d1);
    let b = d1.dot(d2);
    let c = d2.dot(d2);
    let d = d1.dot(w0);
    let e = d2.dot(w0);

    let denom = a * c - b * b;
    if denom.abs() < TOLERANCE_ABS * TOLERANCE_ABS {
        // Parallel lines. Check if they are colinear (on the same line).
        // If colinear, compute the overlap of their ranges and return the midpoint.
        let cross_sq = d1.cross(w0).length_squared();
        let d1_sq = d1.length_squared();
        if cross_sq > tol_sq * d1_sq.max(1.0) {
            return None; // parallel but not colinear �?no intersection
        }
        // Colinear: map l2's range into l1's parameter space.
        // l1: P(t) = l1.origin + t * d1
        // l2: P(s) = l2.origin + s * d2
        // For colinear lines, d2 = +/- d1 (parallel). Map s-parameter to t:
        // l2.origin + s * d2 = l1.origin + t * d1
        // t = (l2.origin - l1.origin).dot(d1) / d1_sq + s * (d2.dot(d1) / d1_sq)
        let sign = if d1.dot(d2) > 0.0 { 1.0 } else { -1.0 };
        let origin_offset = (l2.origin - l1.origin).dot(d1) / d1_sq;
        let t2_lo = origin_offset + r2[0] * sign;
        let t2_hi = origin_offset + r2[1] * sign;
        let overlap_lo = r1[0].max(t2_lo.min(t2_hi));
        let overlap_hi = r1[1].min(t2_lo.max(t2_hi));
        if overlap_hi <= overlap_lo + tol {
            return None; // no overlap
        }
        let t_mid = (overlap_lo + overlap_hi) * 0.5;
        let s_mid = (t_mid - origin_offset) * sign;
        let p = l1.origin + d1 * t_mid;
        return Some((t_mid, s_mid, p));
    }

    let t1 = (b * e - c * d) / denom;
    let t2 = (a * e - b * d) / denom;

    // Check within ranges
    if t1 < r1[0] - tol || t1 > r1[1] + tol || t2 < r2[0] - tol || t2 > r2[1] + tol {
        return None;
    }

    let p1 = l1.origin + d1 * t1;
    let p2 = l2.origin + d2 * t2;

    if (p1 - p2).length_squared() > tol_sq {
        return None; // skew, don't actually intersect
    }

    Some((t1, t2, (p1 + p2) * 0.5))
}

// 鈹€鈹€ Sampling helpers for marching seed-point generation 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

fn sample_plane(plane: &Plane, half_extent: f64, n: usize) -> Vec<DVec3> {
    let u = rcad_kernel::any_perpendicular(plane.normal);
    let v = plane.normal.cross(u);
    let mut pts = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            let su = -half_extent + 2.0 * half_extent * i as f64 / (n - 1).max(1) as f64;
            let sv = -half_extent + 2.0 * half_extent * j as f64 / (n - 1).max(1) as f64;
            pts.push(plane.origin + u * su + v * sv);
        }
    }
    pts
}

fn sample_cone(
    cone: &ConicalSurface,
    h_min: f64,
    h_max: f64,
    n_theta: usize,
    n_h: usize,
) -> Vec<DVec3> {
    let u = rcad_kernel::any_perpendicular(cone.axis);
    let v = cone.axis.cross(u);
    let tan_h = cone.half_angle_rad.tan();
    let mut pts = Vec::with_capacity(n_theta * n_h);
    for ih in 0..n_h {
        let h = h_min + (h_max - h_min) * ih as f64 / (n_h - 1).max(1) as f64;
        let r = h * tan_h;
        for it in 0..n_theta {
            let theta = 2.0 * std::f64::consts::PI * it as f64 / n_theta as f64;
            let p = cone.apex + cone.axis * h + (u * theta.cos() + v * theta.sin()) * r;
            pts.push(p);
        }
    }
    pts
}

fn sample_circle_arc(circle: &Circle3, t_start: f64, t_end: f64, n: usize) -> Vec<DVec3> {
    use rcad_kernel::CurveEval;
    use rcad_kernel::geom::Curve3;
    let curve = Curve3::Circle(*circle);
    (0..n)
        .map(|i| {
            let t = t_start + (t_end - t_start) * i as f64 / (n - 1).max(1) as f64;
            curve.point_at(t)
        })
        .collect()
}

fn circle_param(point: DVec3, circle: &Circle3) -> f64 {
    let u = rcad_kernel::any_perpendicular(circle.normal);
    let v = circle.normal.cross(u);
    let d = point - circle.center;
    let mut theta = d.dot(v).atan2(d.dot(u));
    if theta < 0.0 {
        theta += std::f64::consts::TAU;
    }
    theta
}

fn intersect_line_circle(
    line: &Line3,
    circle: &Circle3,
    tol: f64,
) -> Vec<(f64, f64, DVec3)> {
    let mut results = Vec::new();
    let d = line.direction;
    let o = line.origin;
    let c = circle.center;
    let n = circle.normal;
    let r = circle.radius;
    let r_sq = r * r;

    // Planarity constraint: every point on the circle satisfies (P - c)路n = 0.
    let dn = d.dot(n);
    let w = o - c;
    let wn = w.dot(n);

    if dn.abs() > tol {
        // Line pierces the circle plane at one point.
        let t = -wn / dn;
        let p = o + d * t;
        // �?OCCT-aligned: check distance to circle circumference, not inside-circle.
        // (p - c).length_squared <= r_sq allows points at the circle CENTER (false positive).
        let dist = (p - c).length();
        if (dist - r).abs() <= tol {
            results.push((t, circle_param(p, circle), p));
        }
    } else if wn.abs() <= tol {
        // Line lies in the circle plane �?solve 2D line-circle.
        let t_closest = -w.dot(d);
        let perp_dist_sq = ((o + d * t_closest) - c).length_squared();

        if perp_dist_sq <= r_sq + tol * tol {
            let along = (r_sq - perp_dist_sq).max(0.0).sqrt();
            let t1 = t_closest - along;
            let p1 = o + d * t1;
            results.push((t1, circle_param(p1, circle), p1));

            let t2 = t_closest + along;
            if (t2 - t1).abs() > tol {
                let p2 = o + d * t2;
                results.push((t2, circle_param(p2, circle), p2));
            }
        }
    }

    results
}

fn intersect_coplanar_circles(c1: &Circle3, c2: &Circle3, tol: f64) -> Vec<DVec3> {
    let d_vec = c2.center - c1.center;
    let d = d_vec.length();
    let r1 = c1.radius;
    let r2 = c2.radius;

    // Disjoint or concentric �?no isolated intersection points
    if d > r1 + r2 + tol || d < (r1 - r2).abs() - tol || d < tol {
        return vec![];
    }

    // 2D circle-circle intersection
    // x = projection of intersection point onto the line of centers
    let x = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let y_sq = r1 * r1 - x * x;

    if y_sq < -tol * tol {
        return vec![];
    }
    let y = y_sq.max(0.0).sqrt();

    let dir = d_vec / d;
    let perp = c1.normal.cross(dir).try_normalize().unwrap_or(DVec3::ZERO);

    let mid = c1.center + dir * x;
    if y < tol || perp == DVec3::ZERO {
        vec![mid]
    } else {
        vec![mid + perp * y, mid - perp * y]
    }
}

fn intersect_circle_circle(
    c1: &Circle3,
    c2: &Circle3,
    tol: f64,
) -> Vec<DVec3> {
    let n1 = c1.normal;
    let n2 = c2.normal;
    let cross = n1.cross(n2);
    let cross_len_sq = cross.length_squared();

    // Parallel/coincident planes �?coplanar circle-circle case
    if cross_len_sq < TOLERANCE_ANG * TOLERANCE_ANG {
        let offset = (c2.center - c1.center).dot(n1).abs();
        if offset > tol {
            return vec![];
        }
        return intersect_coplanar_circles(c1, c2, tol);
    }

    // Planes intersect in a line L along the cross-product direction.
    let line_dir = cross / cross_len_sq.sqrt();
    let b = n1.dot(n2);
    let denom = 1.0 - b * b; // sin虏胃 > 0 (not parallel)
    let h1 = c1.center.dot(n1);
    let h2 = c2.center.dot(n2);
    let alpha = (h1 - h2 * b) / denom;
    let beta = (h2 - h1 * b) / denom;
    let base = n1 * alpha + n2 * beta; // a point on line L

    // Intersect sphere of circle1 (center=c1.center, radius=r1) with line L.
    let w = base - c1.center;
    let a = line_dir.dot(line_dir); // = 1 for unit direction
    let b2 = 2.0 * w.dot(line_dir);
    let c = w.dot(w) - c1.radius * c1.radius;
    let disc = b2 * b2 - 4.0 * a * c;

    if disc < -tol * tol {
        return vec![];
    }
    if disc < tol * tol {
        let t = -b2 / (2.0 * a);
        let p = base + line_dir * t;
        return if (p - c2.center).length_squared() <= (c2.radius + tol) * (c2.radius + tol) {
            vec![p]
        } else {
            vec![]
        };
    }

    let sqrt_disc = disc.sqrt();
    let t1 = (-b2 - sqrt_disc) / (2.0 * a);
    let t2 = (-b2 + sqrt_disc) / (2.0 * a);
    let p1 = base + line_dir * t1;
    let p2 = base + line_dir * t2;

    let r2_tol_sq = (c2.radius + tol) * (c2.radius + tol);
    let mut results = Vec::with_capacity(2);
    if (p1 - c2.center).length_squared() <= r2_tol_sq {
        results.push(p1);
    }
    if (p2 - p1).length_squared() > tol * tol
        && (p2 - c2.center).length_squared() <= r2_tol_sq
    {
        results.push(p2);
    }
    results
}

fn in_range(t: f64, range: [f64; 2], tol: f64) -> bool {
    let lo = range[0].min(range[1]) - tol;
    let hi = range[0].max(range[1]) + tol;
    t >= lo && t <= hi
}
fn point_in_sphere_face(pt: DVec3, boundary_verts: &[DVec3], _ds: &DS) -> bool {
    if boundary_verts.is_empty() {
        return false;
    }
    // OCCT-style single-seam sphere: only two pole vertices. An axis-aligned hull of those
    // poles rejects almost every real point on the sphere (e.g. equator vs poles on 卤Y),
    // so plane鈥搒phere tangent handling never records `FaceFace` points and downstream
    // trimming misses imprint geometry (see OCCT `bcommon_simple/A4`).
    if boundary_verts.len() == 2 {
        let a = boundary_verts[0];
        let b = boundary_verts[1];
        let diam = (a - b).length();
        let r = diam * 0.5;
        if r < TOLERANCE_LEN_MIN {
            return false;
        }
        let c = (a + b) * 0.5;
        let radial_err = ((pt - c).length() - r).abs();
        return radial_err < (TOLERANCE_ABS * 500.0).max(TOLERANCE_COORD_SUB * r);
    }
    // Convex hull approximation for faces with a full boundary polygon.
    let cx = boundary_verts.iter().map(|v| v.x).fold(0.0_f64, f64::min)
        ..(boundary_verts.iter().map(|v| v.x).fold(0.0_f64, f64::max) + TOLERANCE_COORD_SUB);
    let cy = boundary_verts.iter().map(|v| v.y).fold(0.0_f64, f64::min)
        ..(boundary_verts.iter().map(|v| v.y).fold(0.0_f64, f64::max) + TOLERANCE_COORD_SUB);
    let cz = boundary_verts.iter().map(|v| v.z).fold(0.0_f64, f64::min)
        ..(boundary_verts.iter().map(|v| v.z).fold(0.0_f64, f64::max) + TOLERANCE_COORD_SUB);
    cx.contains(&pt.x) && cy.contains(&pt.y) && cz.contains(&pt.z)
}

fn sample_surface_generic(surface: &Surface3, n_u: usize, n_v: usize) -> Vec<DVec3> {
    use rcad_kernel::geom::SurfaceEval;
    let [u0, u1, v0, v1] = surface.default_domain();
    let mut pts = Vec::with_capacity(n_u * n_v);
    for iu in 0..n_u {
        for iv in 0..n_v {
            let u = u0 + (u1 - u0) * iu as f64 / (n_u - 1).max(1) as f64;
            let v = v0 + (v1 - v0) * iv as f64 / (n_v - 1).max(1) as f64;
            let p = surface.point_at(u, v);
            if p.is_finite() {
                pts.push(p);
            }
        }
    }
    pts
}

fn intersect_edge_face_numeric(
    curve: &Curve3,
    surface: &Surface3,
    t_range: [f64; 2],
    geom_tol: f64,
) -> Vec<(DVec3, f64)> {
    use rcad_kernel::geom::SurfaceEval;
    use rcad_kernel::projection::closest_point_on_surface;
    use rcad_kernel::CurveEval;
    const N_SAMPLES: usize = 64;
    const MAX_BISECT: usize = 30;

    let eps = geom_tol.max(TOLERANCE_ABS);
    let zero_tol = (eps * TOLERANCE_AREA_REL).max(TOLERANCE_LEN_MIN);

    let [t0, t1] = t_range;
    let mut values = Vec::with_capacity(N_SAMPLES + 1);
    let mut points = Vec::with_capacity(N_SAMPLES + 1);

    for i in 0..=N_SAMPLES {
        let t = t0 + (t1 - t0) * i as f64 / N_SAMPLES as f64;
        let p = curve.point_at(t);
        if !p.is_finite() {
            values.push(f64::NAN);
            points.push(p);
            continue;
        }
        values.push(inttools::marching::surface_implicit(surface, p));
        points.push(p);
    }

    let mut hits = Vec::new();
    for i in 0..N_SAMPLES {
        let va = values[i];
        let vb = values[i + 1];
        if va.is_nan() || vb.is_nan() {
            continue;
        }
        if va * vb > 0.0 {
            continue;
        }
        // Bisection refinement (Stage 1 �?coarse detection)
        let mut ta = t0 + (t1 - t0) * i as f64 / N_SAMPLES as f64;
        let mut tb = t0 + (t1 - t0) * (i + 1) as f64 / N_SAMPLES as f64;
        let mut fa = va;
        let mut fb = vb;
        for _ in 0..MAX_BISECT {
            let tm = (ta + tb) * 0.5;
            let pm = curve.point_at(tm);
            if !pm.is_finite() {
                break;
            }
            let fm = inttools::marching::surface_implicit(surface, pm);
            if fm.abs() < zero_tol {
                hits.push((pm, tm));
                break;
            }
            if (tb - ta).abs() < zero_tol {
                hits.push((pm, tm));
                break;
            }
            if fa * fm < 0.0 {
                tb = tm;
                fb = fm;
            } else {
                ta = tm;
                fa = fm;
            }
        }
        // If bisection didn't converge well, use midpoint
        let tm = (ta + tb) * 0.5;
        let pm = curve.point_at(tm);
        let dedup_dt = (TOLERANCE_MESH_LEGACY).max(eps * 10.0);
        if pm.is_finite() && !hits.iter().any(|(_, t)| (t - tm).abs() < dedup_dt) {
            hits.push((pm, tm));
        }
    }

    // Stage 2 �?Newton refinement: polish each bisection result
    // �?OCCT-aligned: IntCurveSurface_TheExactHInter two-stage approach
    //   coarse sign-change detection �?Newton-Raphson refinement.
    for (point, t) in hits.iter_mut() {
        let initial_t = *t;
        let initial_point = *point;

        // Get initial UV guess via closest-point projection
        let proj = closest_point_on_surface(surface, initial_point, 8);
        let initial_uv = DVec2::new(proj.params.0, proj.params.1);

        if let Some((refined_t, refined_uv)) =
            inttools::curve_surface::newton_refine_curve_surface(
                curve,
                initial_t,
                surface,
                initial_uv,
                20,
                eps,
            )
        {
            // �?OCCT-aligned validation (Stage 3):
            //   1. t within the curve's parametric range
            if refined_t < t_range[0] - eps || refined_t > t_range[1] + eps {
                continue; // Keep bisection result
            }

            //   2. uv within the surface's natural UV domain (if bounded)
            let [u0, u1, v0, v1] = surface.default_domain();
            let u_ok = u0.is_infinite()
                || u1.is_infinite()
                || (refined_uv.x >= u0 - eps && refined_uv.x <= u1 + eps);
            let v_ok = v0.is_infinite()
                || v1.is_infinite()
                || (refined_uv.y >= v0 - eps && refined_uv.y <= v1 + eps);
            if !u_ok || !v_ok {
                continue; // Keep bisection result
            }

            //   3. Distance |C(t) - S(uv)| within tolerance
            let refined_point = curve.point_at(refined_t);
            let surface_point = surface.point_at(refined_uv.x, refined_uv.y);
            if (refined_point - surface_point).length() > eps * 10.0 {
                continue; // Keep bisection result
            }

            // Newton refinement passed all checks �?replace the hit
            *t = refined_t;
            *point = refined_point;
        }
    }

    hits
}

impl<'a> PaveFiller<'a> {

    // ============================================================
    // Edge Overlap Detection
    // ============================================================

} // close impl

#[cfg(test)]
mod tests;
