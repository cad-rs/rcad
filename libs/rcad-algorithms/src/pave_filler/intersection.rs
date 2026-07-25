use std::collections::{BTreeMap, HashSet};
use std::sync::Arc;

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{any_perpendicular, Curve3, Surface3};
use rcad_kernel::CurveEval;

use crate::bopalgo::{fill_map, make_blocks, GlueEnum};
use crate::bopds::ds::{
    DSVertex, InterferenceEE, InterferenceEF, InterferenceVE, InterferenceVF, InterferenceVV,
    ShapeOrigin, DS,
};
use crate::bopds::pave::Pave;
use crate::inttools;
use crate::inttools::context::VfError;
use crate::inttools::edge_edge::compute_curve_aabb;
use crate::inttools::fclass2d::{FClass2d, State};
use crate::pave_filler::helpers::*;
use crate::tolerance::*;
use rcad_kernel::topods::ShapeType;

/// OCCT BRepLib::BoundingVertex (BRepLib.cxx L3013-3123).
/// Given vertex indices, compute the fused vertex center and tolerance.
///
/// For 2 vertices: if one completely encloses the other (distance <= tolerance
/// difference), the larger tolerance vertex is returned. Otherwise the midpoint
/// is shifted toward the larger-tolerance vertex proportional to dR / distance.
///
/// For 3+ vertices: positions are sorted for deterministic floating-point order,
/// then averaged for the center. Tolerance is max(vertex_dist_to_center + vertex_tolerance).
pub fn bounding_vertex(vi: &[usize], ds: &crate::bopds::ds::DS) -> (glam::DVec3, f64) {
    let n = vi.len();
    if n == 0 {
        return (glam::DVec3::ZERO, 0.0);
    }
    if n == 1 {
        return (ds.vertex_point(vi[0]), ds.vertex_tolerance(vi[0]));
    }
    if n == 2 {
        let p0 = ds.vertex_point(vi[0]);
        let p1 = ds.vertex_point(vi[1]);
        let t0 = ds.vertex_tolerance(vi[0]);
        let t1 = ds.vertex_tolerance(vi[1]);
        // m = vertex with larger tolerance, n = vertex with smaller
        let (pm, pn, rm, rn) = if t0 >= t1 {
            (p0, p1, t0, t1)
        } else {
            (p1, p0, t1, t0)
        };
        let dr = rm - rn; // >= 0
        let vd = pn - pm; // vector from larger-tol vertex to smaller-tol vertex
        let d = vd.length();
        if d <= dr || d < f64::EPSILON {
            return (pm, rm);
        }
        // OCCT: aRr = 0.5*(aR[m]+aR[n]+aD); aXYZr = 0.5*(Pm+Pn - VD*(dR/aD))
        let new_tol = 0.5 * (rm + rn + d);
        let new_center = 0.5 * (pm + pn - vd * (dr / d));
        return (new_center, new_tol);
    }
    // n > 2
    let mut points: Vec<glam::DVec3> = vi.iter().map(|&i| ds.vertex_point(i)).collect();
    // OCCT sorts for deterministic float-sum order (issue 0027540)
    points.sort_by(|a, b| {
        a.x.partial_cmp(&b.x)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.y.partial_cmp(&b.y).unwrap_or(std::cmp::Ordering::Equal))
            .then_with(|| a.z.partial_cmp(&b.z).unwrap_or(std::cmp::Ordering::Equal))
    });
    let center = points.iter().sum::<glam::DVec3>() / points.len() as f64;
    // aDmax = max(vertex_dist + vertex_tolerance)
    let tol = vi
        .iter()
        .map(|&i| (ds.vertex_point(i) - center).length() + ds.vertex_tolerance(i))
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap_or(0.0);
    (center, tol)
}

/// OCCT IntTools_CommonPrt::Type() ??VERTEX or EDGE.

/// BOPAlgo_VertexEdge (PaveFiller_2.cxx L40-134).
/// Solver that projects a vertex onto an edge's curve via IntTools_Context::ComputeVE.
/// Used by IntersectVE to compute the exact parameter for PB splitting.
pub(crate) struct VertexEdgeSolver {
    pub(crate) vi: usize,
    pub(crate) ei: usize,
    pub(crate) param: f64,
    pub(crate) tol_vnew: f64,
    flag: i32,
}

impl VertexEdgeSolver {
    pub fn new() -> Self {
        Self {
            vi: 0,
            ei: 0,
            param: -1.0,
            tol_vnew: -1.0,
            flag: -1,
        }
    }
    pub fn set_data(&mut self, vi: usize, ei: usize) {
        self.vi = vi;
        self.ei = ei;
    }
    /// OCCT L104-121: Perform() = myContext->ComputeVE(myV, myE, myT, myTolVNew, myFuzzyValue)
    pub fn perform(&mut self, ctx: &mut super::IntToolsContext, ds: &DS, fuzz: f64) {
        match ctx.compute_ve(ds, self.vi, self.ei, fuzz) {
            Ok(res) => {
                self.flag = 0;
                self.param = res.param;
                self.tol_vnew = res.tolerance;
            }
            Err(_) => {
                self.flag = -3;
            }
        }
    }
    pub fn is_done(&self) -> bool {
        self.flag == 0
    }
    pub fn flag(&self) -> i32 {
        self.flag
    }
}
#[derive(Clone, Copy)]
enum EfHit {
    Vertex { point: DVec3, param: f64 },
    Edge { t1: f64, t2: f64 },
}

/// OCCT L55-93: BOPAlgo_EdgeFace (architecture diff: rcad equivalent).
// OCCT L55-162: BOPAlgo_EdgeFace — target for future alignment
struct EdgeFaceTask {
    myIE: usize,
    myIF: usize,
    myFlag: i32,
    myPB: usize,
    myNewSR: [f64; 2],
    myRange: [f64; 2],
    myBox1: Option<(DVec3, DVec3)>,
    myBox2: Option<(DVec3, DVec3)>,
    myFuzzyValue: f64,
    bExpressCompute: bool,
    myCommonParts: Vec<EfHit>,
    myMinDist: f64,
    myHasErrors: bool,
    myIsDone: bool,
}

impl EdgeFaceTask {
    fn new() -> Self {
        EdgeFaceTask {
            myIE: usize::MAX,
            myIF: usize::MAX,
            myFlag: -1,
            myPB: usize::MAX,
            myNewSR: [0.0, 0.0],
            myRange: [0.0, 0.0],
            myBox1: None,
            myBox2: None,
            myFuzzyValue: 0.0,
            bExpressCompute: false,
            myCommonParts: Vec::new(),
            myMinDist: f64::MAX,
            myHasErrors: false,
            myIsDone: false,
        }
    }
    fn set_indices(&mut self, nE: usize, nF: usize) {
        self.myIE = nE;
        self.myIF = nF;
    }
    fn indices(&self) -> (usize, usize) {
        (self.myIE, self.myIF)
    }
    fn set_pave_block(&mut self, pb: usize) {
        self.myPB = pb;
    }
    fn pave_block(&self) -> usize {
        self.myPB
    }
    fn set_fuzzy_value(&mut self, f: f64) {
        self.myFuzzyValue = f;
    }
    fn set_boxes(&mut self, b1: Option<(DVec3, DVec3)>, b2: Option<(DVec3, DVec3)>) {
        self.myBox1 = b1;
        self.myBox2 = b2;
    }
    fn use_quick_coincidence_check(&mut self, flag: bool) {
        self.bExpressCompute = flag;
    }
    fn set_new_sr(&mut self, sr: [f64; 2]) {
        self.myNewSR = sr;
    }
    fn set_range(&mut self, r: [f64; 2]) {
        self.myRange = r;
    }
    fn is_done(&self) -> bool {
        self.myIsDone
    }
    fn has_errors(&self) -> bool {
        self.myHasErrors
    }
    fn common_parts(&self) -> &[EfHit] {
        &self.myCommonParts
    }
    fn minimal_distance(&self) -> f64 {
        self.myMinDist
    }
    fn perform(&mut self, ctx: &mut crate::inttools::context::Context, ds: &DS) {
        let ef_range = self.myRange;
        let etr = ds.edges[self.myIE].t_range;
        let range = [ef_range[0].max(etr[0]), ef_range[1].min(etr[1])];
        let tol_ef = ds
            .edge_tolerance(self.myIE)
            .max(ds.face_tolerance(self.myIF))
            .max(CONFUSION);
        if range[1] - range[0] <= tol_ef {
            return;
        }
        let hits = compute_ef_hits(ds, self.myIE, self.myIF, &range, tol_ef);
        self.myCommonParts = hits;
        self.myIsDone = true;
    }
}

/// OCCT IntTools_EdgeFace (PaveFiller_5.cxx L340-480): compute edge-face intersection hits.
/// Returns VERTEX-type (point) and EDGE-type (range) common parts.
fn compute_ef_hits(
    ds: &DS,
    edge_idx: usize,
    face_idx: usize,
    ef_range: &[f64; 2],
    ef_tol: f64,
) -> Vec<EfHit> {
    let edge_curve = &ds.edges[edge_idx].curve;
    let face_surface = &ds.faces[face_idx].surface;
    let mut hits: Vec<EfHit> = match (edge_curve, face_surface) {
        (Curve3::Line(line), Surface3::Plane(plane)) => {
            crate::inttools::edge_face::intersect_line_plane_with_tol(
                line, *ef_range, plane, ef_tol,
            )
            .into_iter()
            .map(|h| EfHit::Vertex {
                point: h.point,
                param: h.edge_param,
            })
            .collect()
        }
        (Curve3::Line(line), Surface3::Cylinder(cyl)) => {
            crate::inttools::curve_surface::intersect_line_cylinder_with_tol(
                line, *ef_range, cyl, ef_tol,
            )
            .into_iter()
            .map(|h| EfHit::Vertex {
                point: h.point,
                param: h.curve_param,
            })
            .collect()
        }
        (Curve3::Line(line), Surface3::Sphere(sph)) => {
            crate::inttools::curve_surface::intersect_line_sphere_with_tol(
                line, *ef_range, sph, ef_tol,
            )
            .into_iter()
            .map(|h| EfHit::Vertex {
                point: h.point,
                param: h.curve_param,
            })
            .collect()
        }
        (Curve3::Line(line), Surface3::Cone(cone)) => {
            crate::inttools::curve_surface::intersect_line_cone_with_tol(
                line, *ef_range, cone, ef_tol,
            )
            .into_iter()
            .map(|h| EfHit::Vertex {
                point: h.point,
                param: h.curve_param,
            })
            .collect()
        }
        _ => Vec::new(),
    };
    // Phase 2: edge coincident with face (OCCT MakeType EDGE)
    if hits.is_empty() && ef_range[1] - ef_range[0] > ef_tol {
        let mid_t = (ef_range[0] + ef_range[1]) * 0.5;
        let mid_pt = edge_curve.point_at(mid_t);
        use crate::inttools::bean_face_intersector::BeanFaceIntersector;
        let mut bfi = BeanFaceIntersector::new();
        bfi.init_curve_surface(
            ds.edges[edge_idx].curve.clone(),
            ef_tol,
            ds.faces[face_idx].surface.clone(),
            ef_tol,
        );
        bfi.set_bean_parameters(ef_range[0], ef_range[1]);
        bfi.perform();
        if bfi.is_done() && !bfi.result().is_empty() {
            for r in bfi.result() {
                let t = (r.first() + r.last()) * 0.5;
                let p = ds.edges[edge_idx].curve.point_at(t);
                hits.push(EfHit::Vertex { point: p, param: t });
            }
        }
    }
    hits
}

// OCCT L37-132: BOPAlgo_VertexFace
struct VertexFaceTask {
    myIV: usize,
    myIF: usize,
    myFlag: i32,
    myT1: f64,
    myT2: f64,
    myTolVNew: f64,
    myHasErrors: bool,
}

impl VertexFaceTask {
    fn new() -> Self {
        VertexFaceTask {
            myIV: usize::MAX,
            myIF: usize::MAX,
            myFlag: -1,
            myT1: -1.0,
            myT2: -1.0,
            myTolVNew: -1.0,
            myHasErrors: false,
        }
    }
    fn set_indices(&mut self, nV: usize, nF: usize) {
        self.myIV = nV;
        self.myIF = nF;
    }
    fn indices(&self) -> (usize, usize) {
        (self.myIV, self.myIF)
    }
    fn flag(&self) -> i32 {
        self.myFlag
    }
    fn has_errors(&self) -> bool {
        self.myHasErrors
    }
    fn parameters(&self) -> (f64, f64) {
        (self.myT1, self.myT2)
    }
    fn vertex_new_tolerance(&self) -> f64 {
        self.myTolVNew
    }
    fn perform(&mut self, ctx: &mut crate::inttools::context::Context, ds: &DS, fuzz: f64) {
        match ctx.compute_vf(ds, self.myIV, self.myIF, fuzz) {
            Ok(res) => {
                self.myFlag = 0;
                self.myT1 = res.u;
                self.myT2 = res.v;
                self.myTolVNew = res.tolerance;
            }
            Err(VfError::ProjectionFailed) => {
                self.myFlag = -1;
            }
            Err(VfError::DistanceTooLarge) => {
                self.myFlag = -2;
            }
            Err(VfError::PointOutsideFace) => {
                self.myFlag = -3;
            }
        }
    }
}

impl<'a> super::PaveFiller<'a> {
    /// OCCT PaveFiller L145: BOPDS_Iterator  ?BVH pair enumeration
    pub(crate) fn build_box_tree(&self, is_a: bool, is_edge: bool) -> crate::bvh::BoxTree {
        use crate::bvh::{Aabb, BoxTree};
        let (ds_start, end) = if is_edge {
            if is_a {
                (0, self.ds.a_edge_count)
            } else {
                (self.ds.a_edge_count, self.ds.edges.len())
            }
        } else {
            if is_a {
                (0, self.ds.a_vertex_count)
            } else {
                (self.ds.a_vertex_count, self.ds.vertices.len())
            }
        };
        let n = end - ds_start;
        let mut indices = Vec::with_capacity(n);
        let mut aabbs = Vec::with_capacity(n);

        for local_i in 0..n {
            let ds_i = ds_start + local_i;
            indices.push(ds_i);
            let aabb = if is_edge {
                let e = &self.ds.edges[ds_i];
                let pts = [
                    self.ds.vertices[e.start_vertex].point,
                    self.ds.vertices[e.end_vertex].point,
                ];
                let mut a = Aabb::empty();
                for &p in &pts {
                    a.expand_point(p);
                }
                // Expand for edge tolerance — OCCT: extents + gap
                let gap_tol = e.geom_tol.max(CONFUSION);
                a.min -= DVec3::splat(gap_tol);
                a.max += DVec3::splat(gap_tol);
                a.gap = gap_tol;
                a
            } else {
                let pt = self.ds.vertex_point(ds_i);
                // OCCT BOPDS_ShapeInfo::Box(): extents = point, gap = tolerance
                // OCCT: use box_gap from shape_info (set at Init, NOT updated during VE)
                let si = self.ds.vertex_shape_idx.get(ds_i).copied().unwrap_or(ds_i);
                let gap = if si < self.ds.shape_info.len() {
                    self.ds.shape_info[si].box_gap
                } else {
                    CONFUSION
                };
                Aabb {
                    min: pt,
                    max: pt,
                    gap,
                }
            };
            aabbs.push(aabb);
        }
        BoxTree::build(indices, aabbs)
    }
    ///  BOPDS_Iterator  ?build a single BVH for all elements
    /// of the given shape type (both operands A and B combined), used for
    /// single-pass cross-operand pair traversal.
    pub(crate) fn build_box_tree_combined(&self, is_edge: bool) -> crate::bvh::BoxTree {
        use crate::bvh::{Aabb, BoxTree};
        let n = if is_edge {
            self.ds.edges.len()
        } else {
            self.ds.vertices.len()
        };
        let mut indices = Vec::with_capacity(n);
        let mut aabbs = Vec::with_capacity(n);
        for ds_i in 0..n {
            indices.push(ds_i);
            let aabb = if is_edge {
                let e = &self.ds.edges[ds_i];
                // OCCT BOPDS_ShapeInfo::Box() computes AABB from full curve geometry,
                // not just endpoint vertices. For Line edges, endpoints suffice.
                // For Circle/Ellipse edges, endpoints alone miss the full extent.
                let mut a = Aabb::empty();
                let curve = &e.curve;
                // Sample curve along its range to capture full geometric extent
                let t_range = e.t_range;
                let n_samples = match curve {
                    Curve3::Line(_) => 2,     // Line: endpoints are sufficient
                    Curve3::Circle(_) => 16,  // Circle: sample around the circle
                    Curve3::Ellipse(_) => 16, // Ellipse: sample around the ellipse
                    _ => 8,                   // Other curves: conservative sampling
                };
                for si in 0..n_samples {
                    let t = t_range[0]
                        + (t_range[1] - t_range[0]) * si as f64 / (n_samples - 1).max(1) as f64;
                    let p = curve.point_at(t);
                    a.expand_point(p);
                }
                let gap_tol = e.geom_tol.max(CONFUSION);
                a.min -= DVec3::splat(gap_tol);
                a.max += DVec3::splat(gap_tol);
                a.gap = gap_tol;
                a
            } else {
                let pt = self.ds.vertex_point(ds_i);
                // OCCT BOPDS_ShapeInfo::Box(): extents = point, gap = tolerance
                // OCCT: use box_gap from shape_info (set at Init, NOT updated during VE)
                let si = self.ds.vertex_shape_idx.get(ds_i).copied().unwrap_or(ds_i);
                let gap = if si < self.ds.shape_info.len() {
                    self.ds.shape_info[si].box_gap
                } else {
                    CONFUSION
                };
                Aabb {
                    min: pt,
                    max: pt,
                    gap,
                }
            };
            aabbs.push(aabb);
        }
        BoxTree::build(indices, aabbs)
    }
    // OCCT PaveFiller_2.cxx L141-208: PerformVE
    // Groups cross-operand (nV, nE) pairs by edge, then calls IntersectVE.
    // rcad: receives pre-computed cross-operand pairs from BOPDS_Iterator.
    // OCCT BOPAlgo_PaveFiller_2.cxx L141: PerformVE
    // OCCT BOPAlgo_PaveFiller_2.cxx L141-207: PerformVE
    pub(crate) fn perform_ve(&mut self) {
        // OCCT L143: FillShrunkData(VERTEX, EDGE)
        self.fill_shrunk_data(ShapeType::Vertex, ShapeType::Edge);

        // OCCT L148-152: iSize = myIterator->ExpectedLength()
        self.my_iterator
            .initialize(ShapeType::Vertex, ShapeType::Edge);
        let i_size = self.my_iterator.expected_length();
        if i_size == 0 {
            return;
        }

        // OCCT L155: NCollection_IndexedDataMap<handle<PaveBlock>, NCollection_List<int>> aMVEPairs
        // rcad: HashMap<edge_idx, Vec<vertex_idx> (edge identified by its first PB)
        let a_vc = self.ds.a_vertex_count;
        let a_ec = self.ds.a_edge_count;
        let mut a_mve_pairs: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();

        // OCCT L156-205: for (; myIterator->More(); myIterator->Next())
        for _ in 0..i_size {
            let (n_v, n_e) = self.my_iterator.value();
            self.my_iterator.next();
            // rcad: cross-operand filter (OCCT: BOPDS_Iterator enforces this)
            if (n_v < a_vc) == (n_e < a_ec) {
                continue;
            }

            // OCCT L165-168: aSIE.HasSubShape(nV) — vertex is subshape of edge?
            if self.ds.edge_has_vertex(n_v, n_e) {
                continue;
            }

            // OCCT L171-174: aSIE.HasFlag()
            if self.ds.edge_has_flag(n_e) {
                continue;
            }

            // OCCT L176-179: myDS->HasInterf(nV, nE)
            if self.ds.has_interf_ve(n_v, n_e) {
                continue;
            }

            // OCCT L181-184: myDS->HasInterfShapeSubShapes(nV, nE)
            if self.ds.has_interf_ve_via_faces(n_v, n_e) {
                continue;
            }

            // OCCT L186-190: aLPB empty (const accessor, no lazy init)
            if self.ds.edge_pave_blocks(n_e).is_empty() {
                continue;
            }

            // OCCT L192-197: first PB not splittable
            if !self.ds.edge_pave_blocks(n_e)[0]
                .0
                .read()
                .unwrap()
                .is_splittable
            {
                continue;
            }

            // OCCT L199-204: group vertices by edge (keyed by first PB)
            a_mve_pairs.entry(n_e).or_default().push(n_v);
        }

        // OCCT L207: IntersectVE(aMVEPairs, ...)
        self.intersect_ve(&a_mve_pairs, true);
    }

    // OCCT PaveFiller_2.cxx L212-395: IntersectVE
    fn intersect_ve(
        &mut self,
        the_ve_pairs: &std::collections::HashMap<usize, Vec<usize>>,
        the_add_interfs: bool,
    ) {
        // OCCT L217-221: aNbVE = theVEPairs.Extent()
        let a_nb_ve = the_ve_pairs.len();
        if a_nb_ve == 0 {
            return;
        }

        // OCCT L223-227: aVEs.SetIncrement(aNbVE)
        if the_add_interfs {
            self.ds.interf_ve.reserve(a_nb_ve);
        }

        // OCCT L230: BOPAlgo_VectorOfVertexEdge aVVE
        // rcad: Vec storing (nV, nE) task data
        struct VeTask {
            n_v: usize,
            n_e: usize,
        }
        let mut a_vve: Vec<VeTask> = Vec::new();

        // OCCT L235: aDMVSD ??map (nVSD, nE) -> list of original vertices
        let mut a_dmv_sd: std::collections::HashMap<(usize, usize), Vec<usize>> =
            std::collections::HashMap::new();

        // OCCT L238-291: for (i = 1; i <= aNbVE; ++i)
        for (&n_e, verts) in the_ve_pairs {
            // OCCT L244: nE = aPB->OriginalEdge() ??rcad: n_e is the edge index directly

            // OCCT L247-254: build aMVPB from all PBs of this edge
            let mut a_mv_pb: std::collections::HashSet<usize> = std::collections::HashSet::new();
            for spb in self.ds.edge_pave_blocks(n_e) {
                let pb = spb.0.read().unwrap();
                a_mv_pb.insert(pb.pave1.vertex_idx);
                a_mv_pb.insert(pb.pave2.vertex_idx);
            }

            // OCCT L256-291: iterate vertex list for this PB
            for &n_v in verts {
                // OCCT L262-263: resolve SD vertex
                let n_vsd = self.ds.has_shape_sd(n_v).unwrap_or(n_v);

                // OCCT L265-268: skip if nVSD is a PB endpoint
                if a_mv_pb.contains(&n_vsd) {
                    continue;
                }

                // OCCT L270-277: check if (nVSD, nE) already in aDMVSD
                let a_pair = (n_vsd, n_e);
                if let Some(p_li) = a_dmv_sd.get_mut(&a_pair) {
                    // Already added ??just append the original vertex
                    p_li.push(n_v);
                    continue;
                }

                // OCCT L279-291: new pair ??create solver task
                a_dmv_sd.insert(a_pair, vec![n_v]);
                a_vve.push(VeTask { n_v: n_vsd, n_e });
            }
        }

        // OCCT L294: aNbVE = aVVE.Length()
        let a_nb_ve = a_vve.len();

        // OCCT L302-304: BOPTools_Parallel::Perform(myRunParallel, aVVE, myContext)
        // rcad: sequential execution

        // OCCT L312: NCollection_Map<int> aMEdges
        let mut a_m_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();

        // OCCT L315-387: for (i = 0; i < aNbVE; ++i)
        for task in &a_vve {
            // OCCT L321-329: if flag != 0 ??skip / warn
            let res = match self.context.compute_ve(
                self.ds,
                task.n_v,
                task.n_e,
                self.fuzzy_tolerance,
            ) {
                Ok(res) => {
                    if std::env::var("RCAD_DEBUG_VE").is_ok() {
                        let pt_v = self.ds.vertex_point(task.n_v);
                        let tol_v = self.ds.vertex_tolerance(task.n_v);
                        let tol_e = self.ds.edge_tolerance(task.n_e);
                        eprintln!(
                            "[VE_OK] nV={} nE={} aT={:.12e} dist={:.12e} tolV={:.12e} tolE={:.12e} V=({:.12e},{:.12e},{:.12e})",
                            task.n_v,
                            task.n_e,
                            res.param,
                            res.tolerance - tol_e,
                            tol_v,
                            tol_e,
                            pt_v.x,
                            pt_v.y,
                            pt_v.z
                        );
                    }
                    res
                }
                Err(e) => {
                    if std::env::var("RCAD_DEBUG_VE").is_ok() {
                        eprintln!("[VE_ERR] nV={} nE={} err={:?}", task.n_v, task.n_e, e);
                    }
                    // OCCT L324-328: HasErrors ??AddIntersectionFailedWarning
                    self.add_intersection_failed_warning(task.n_v, task.n_e);
                    continue;
                }
            };

            // OCCT L332-338: extract result
            let a_t = res.param;
            let a_tol_v_new = res.tolerance;
            // OCCT L338: nVx = UpdateVertex(nV, aTolVNew)
            let n_vx = self.update_vertex(task.n_v, a_tol_v_new);

            // OCCT L341-354: Find PB on edge containing aT
            let a_lpb = self.ds.edge_pave_blocks(task.n_e);
            let pb_idx = a_lpb.iter().position(|spb| {
                let pb = spb.0.read().unwrap();
                let (a_t1, a_t2) = pb.range();
                // OCCT L380: if (aT > aT1 && aT < aT2)
                a_t > a_t1 && a_t < a_t2
            });
            let pb_idx = match pb_idx {
                Some(i) => i,
                None => continue,
            };

            // OCCT L360-363: AppendExtPave
            let a_pave = Pave {
                vertex_idx: n_vx,
                param: a_t,
            };
            a_lpb[pb_idx].0.write().unwrap().append_ext_pave(a_pave);
            a_m_edges.insert(task.n_e);

            // OCCT L366-387: create interferences
            if the_add_interfs {
                // OCCT L369: BOPDS_Pair aPair(nV, nE)
                let a_pair = (task.n_v, task.n_e);
                // OCCT L370: aDMVSD.Find(aPair)
                if let Some(a_li) = a_dmv_sd.get(&a_pair) {
                    // OCCT L371-386: for each original vertex
                    for &n_v_old in a_li {
                        // OCCT L376-378: create VE interference
                        let b_new = self.ds.is_new_vertex(n_vx);
                        self.ds.interf_ve.push(InterferenceVE {
                            vertex: n_v_old,
                            edge: task.n_e,
                            param: a_t,
                            index_new: n_vx,
                        });
                        // OCCT L380: myDS->AddInterf(nVOld, nE)
                        self.ds.try_add_interf(n_v_old, task.n_e);
                        // OCCT L382-385: if new shape, SetIndexNew
                        // rcad: index_new is set directly above
                    }
                }
            }
        }

        // OCCT L394: SplitPaveBlocks(aMEdges, theAddInterfs)
        if !a_m_edges.is_empty() {
            self.split_pave_blocks(&a_m_edges, the_add_interfs);
        }
    }
    // OCCT BOPAlgo_PaveFiller_3.cxx L145-590: PerformEE
    //
    // OCCT structure:
    //   L147: FillShrunkData(EDGE, EDGE)
    //   L149-150: Iterator init, iSize check
    //   L157-175: variable declarations (aEEs, aMEdges, allocators)
    //   L181-267: Phase 1 -- collect BOPAlgo_EdgeEdge tasks
    //   L269-278: Phase 2 -- parallel execution (BOPTools_Parallel)
    //   L285-556: Phase 3 -- process CommonPrt (VERTEX/EDGE types)
    //   L558-585: Phase 4 -- PerformCommonBlocks + PerformNewVertices + SplitPaveBlocks
    //
    // rcad architecture differences:
    //   - intersect_ee combines computation + InterfEE creation (no CommonPrt)
    //   - treat_new_vertices() called separately from perform() (not inside)
    //   - No aMVCPB / aMPBLPB coupling (common blocks handled elsewhere)
    //   - Sequential execution (no BOPTools_Parallel)
    // OCCT BOPAlgo_PaveFiller_3.cxx L914-955: GetPBBox
    /// Translates OCCT's BOPAlgo_PaveFiller::GetPBBox.
    /// Returns false if the pave block's range is degenerate (length ≤ TOLERANCE_LEN_MIN).
    /// On return, `the_box` contains the AABB of the PB's curve segment.
    pub(crate) fn get_pb_box(
        &self,
        the_e: usize,
        the_pb: &crate::bopds::pave::PaveBlock,
        the_pb_box: &mut std::collections::HashMap<usize, (glam::DVec3, glam::DVec3)>,
        the_first: &mut f64,
        the_last: &mut f64,
        the_s_first: &mut f64,
        the_s_last: &mut f64,
        the_box: &mut (glam::DVec3, glam::DVec3),
    ) -> bool {
        // OCCT L916: thePB->Range(theFirst, theLast)
        let (r1, r2) = the_pb.range();
        *the_first = r1;
        *the_last = r2;

        // OCCT L917: bValid = theLast - theFirst > Precision::PConfusion()
        let b_valid = *the_last - *the_first > crate::tolerance::TOLERANCE_LEN_MIN;
        if !b_valid {
            return b_valid;
        }

        // OCCT L918-920: if (HasShrunkData) { ShrunkData(...); return bValid; }
        // NOTE: rcad's shrunk_data() does NOT return a Bnd_Box, so we fall through
        // to compute the AABB (either from cache or by evaluating the curve).
        if the_pb.has_shrunk_data() {
            let (ts1, ts2, _spl) = the_pb.shrunk_data();
            *the_s_first = ts1;
            *the_s_last = ts2;
        } else {
            // OCCT L922: theSFirst = theFirst; theSLast = theLast
            *the_s_first = *the_first;
            *the_s_last = *the_last;
        }

        // OCCT L923-927: cache lookup
        let key = the_pb as *const crate::bopds::pave::PaveBlock as usize;
        if let Some(cached) = the_pb_box.get(&key) {
            *the_box = *cached;
        } else {
            // OCCT L928-933: BRepAdaptor_Curve + BndLib_Add3dCurve
            let curve = &self.ds.edges[the_e].curve;
            let a_tol = self.ds.edge_tolerance(the_e) + crate::tolerance::CONFUSION;
            let aabb = compute_curve_aabb(curve, *the_s_first, *the_s_last, a_tol);
            the_pb_box.insert(key, aabb);
            *the_box = aabb;
        }

        b_valid
    }

    // OCCT BOPAlgo_PaveFiller_3.cxx L145: PerformEE
    // OCCT BOPAlgo_PaveFiller_3.cxx L145-585: PerformEE
    pub(crate) fn perform_ee(&mut self) {
        // OCCT L147: FillShrunkData(EDGE, EDGE)
        self.fill_shrunk_data(ShapeType::Edge, ShapeType::Edge);

        // OCCT L149-150: myIterator->Initialize(EDGE, EDGE)
        // iSize = myIterator->ExpectedLength()
        self.my_iterator
            .initialize(ShapeType::Edge, ShapeType::Edge);
        let i_size = self.my_iterator.expected_length();

        // OCCT L152-155: if (!iSize) return
        if i_size == 0 {
            return;
        }

        // OCCT L157-166: variable declarations
        //   bExpressCompute, bIsPBSplittable1/2, i, iX, nE1, nE2, aNbCPrts, k, aNbEdgeEdge
        //   nV11-22, aTS11-22, aT11-22, aType
        //   NCollection_List::Iterator aIt1, aIt2
        //   NCollection_Map<int> aMEdges  -- modified edges set
        // OCCT L169-176: allocator + data maps
        //   aMPBLPB: IndexedDataMap<PB, List<PB>>
        //   aMVCPB:  IndexedDataMap<Shape, CoupleOfPaveBlocks>
        //   aDMPBBox: DataMap<PB, Bnd_Box>

        // rcad EeTask replaces BOPAlgo_EdgeEdge (no TopoDS_Shape, no handle types)
        struct EeTask {
            pb1: crate::bopds::pave::SharedPB,
            pb2: crate::bopds::pave::SharedPB,
            nE1: usize,
            nE2: usize,
            aT11: f64,
            aT12: f64,
            aTS11: f64,
            aTS12: f64,
            aT21: f64,
            aT22: f64,
            aTS21: f64,
            aTS22: f64,
            nV11: usize,
            nV12: usize,
            nV21: usize,
            nV22: usize,
            b_express_compute: bool,
            b_is_pb_splittable1: bool,
            b_is_pb_splittable2: bool,
        }

        // OCCT L167: NCollection_Map<int> aMEdges
        let mut a_m_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();

        // OCCT L169-176: aMPBLPB, aMVCPB, aDMPBBox
        let mut a_mplpb: std::collections::BTreeMap<usize, Vec<usize>> =
            std::collections::BTreeMap::new();
        let mut a_mvcpb: std::collections::HashMap<usize, crate::bopds::pave::CoupleOfPaveBlocks> =
            std::collections::HashMap::new();
        let mut a_dmpbbox: std::collections::HashMap<usize, (glam::DVec3, glam::DVec3)> =
            std::collections::HashMap::new();

        // OCCT L178-179: aEEs.SetIncrement(iSize)
        self.ds.interf_ee.reserve(i_size);

        // OCCT L181: for (; myIterator->More(); myIterator->Next())
        let a_ec = self.ds.a_edge_count;
        let mut a_vee: Vec<EeTask> = Vec::new();

        for _ in 0..i_size {
            let (nE1, nE2) = self.my_iterator.value();
            self.my_iterator.next();
            // rcad: cross-operand filter (OCCT: done by BOPDS_Iterator)
            if (nE1 < a_ec) == (nE2 < a_ec) {
                continue;
            }

            // OCCT L189-192: myDS->ShapeInfo(nE1).HasFlag()
            if self.ds.edge_has_flag(nE1) || self.ds.edge_has_flag(nE2) {
                continue;
            }

            // OCCT L200-204: myDS->ChangePaveBlocks(nE1).IsEmpty()
            let a_lpb1 = self.ds.edge_pave_blocks(nE1);
            if a_lpb1.is_empty() {
                continue;
            }

            // OCCT L206-210: myDS->ChangePaveBlocks(nE2).IsEmpty()
            let a_lpb2 = self.ds.edge_pave_blocks(nE2);
            if a_lpb2.is_empty() {
                continue;
            }

            // OCCT L215-266: PB pair iteration — using get_pb_box for AABB + cache (OCCT L914-955)
            for pb1 in a_lpb1.iter() {
                let mut aT11 = 0.0;
                let mut aT12 = 0.0;
                let mut aTS11 = 0.0;
                let mut aTS12 = 0.0;
                let mut aBB1 = (glam::DVec3::ZERO, glam::DVec3::ZERO);
                let pb1_r = pb1.0.read().unwrap();
                if !self.get_pb_box(
                    nE1,
                    &pb1_r,
                    &mut a_dmpbbox,
                    &mut aT11,
                    &mut aT12,
                    &mut aTS11,
                    &mut aTS12,
                    &mut aBB1,
                ) {
                    drop(pb1_r);
                    continue;
                }
                let b_is_pb_splittable1 = pb1_r.has_shrunk_data() && pb1_r.is_splittable;

                // OCCT L231: aPB1->Indices(nV11, nV12)
                let (nV11, nV12) = pb1_r.indices();
                drop(pb1_r);

                // OCCT L233-265: aIt2.Initialize(aLPB2); for (; aIt2.More(); aIt2.Next())
                for pb2 in a_lpb2.iter() {
                    let mut aT21 = 0.0;
                    let mut aT22 = 0.0;
                    let mut aTS21 = 0.0;
                    let mut aTS22 = 0.0;
                    let mut aBB2 = (glam::DVec3::ZERO, glam::DVec3::ZERO);
                    let pb2_r = pb2.0.read().unwrap();
                    if !self.get_pb_box(
                        nE2,
                        &pb2_r,
                        &mut a_dmpbbox,
                        &mut aT21,
                        &mut aT22,
                        &mut aTS21,
                        &mut aTS22,
                        &mut aBB2,
                    ) {
                        drop(pb2_r);
                        continue;
                    }
                    let b_is_pb_splittable2 = pb2_r.has_shrunk_data() && pb2_r.is_splittable;

                    // OCCT L245: if (aBB1.IsOut(aBB2)) continue;
                    if aBB1.0.x > aBB2.1.x
                        || aBB1.1.x < aBB2.0.x
                        || aBB1.0.y > aBB2.1.y
                        || aBB1.1.y < aBB2.0.y
                        || aBB1.0.z > aBB2.1.z
                        || aBB1.1.z < aBB2.0.z
                    {
                        drop(pb2_r);
                        continue;
                    }

                    // OCCT L250: aPB2->Indices(nV21, nV22)
                    let (nV21, nV22) = pb2_r.indices();

                    // OCCT L252: bExpressCompute = same vertex bounds
                    let b_express_compute =
                        (nV11 == nV21 && nV12 == nV22) || (nV12 == nV21 && nV11 == nV22);

                    drop(pb2_r);

                    a_vee.push(EeTask {
                        pb1: pb1.clone(),
                        pb2: pb2.clone(),
                        nE1,
                        nE2,
                        aT11,
                        aT12,
                        aTS11,
                        aTS12,
                        aT21,
                        aT22,
                        aTS21,
                        aTS22,
                        nV11,
                        nV12,
                        nV21,
                        nV22,
                        b_express_compute,
                        b_is_pb_splittable1,
                        b_is_pb_splittable2,
                    });
                }
            }
        }

        // OCCT L269: aNbEdgeEdge = aVEdgeEdge.Length()
        let a_nb_edge_edge = a_vee.len();

        // OCCT L271-278: SetProgressRange + BOPTools_Parallel::Perform
        // rcad: sequential execution (no BOPTools_Parallel)

        // OCCT L285-556: Process results
        for k in 0..a_nb_edge_edge {
            // OCCT L291: Bnd_Box aBB1, aBB2;
            // OCCT L293: BOPAlgo_EdgeEdge& anEdgeEdge = aVEdgeEdge(k);
            let task = &a_vee[k];
            let nE1 = task.nE1;
            let nE2 = task.nE2;
            // Use PB's FULL range [aT11, aT12] for intersection computation
            let sr1 = [task.aT11.min(task.aT12), task.aT11.max(task.aT12)];
            let sr2 = [task.aT21.min(task.aT22), task.aT21.max(task.aT22)];

            // OCCT L294-301: if (!IsDone() || HasErrors()) continue;
            let edge1 = &self.ds.edges[nE1];
            let edge2 = &self.ds.edges[nE2];
            let tol = self.ee_tol(nE1, nE2);
            let e1_curve = edge1.curve.clone();
            let e2_curve = edge2.curve.clone();
            drop(edge1);
            drop(edge2);
            // OCCT IntTools_EdgeEdge computes CommonParts from the curve geometry.
            let hits: Vec<(f64, f64, DVec3, [f64; 2], [f64; 2])> = match (&e1_curve, &e2_curve) {
                (Curve3::Line(l1), Curve3::Line(l2)) => intersect_line_line(l1, sr1, l2, sr2, tol)
                    .into_iter()
                    .map(|(t1, t2, p)| (t1, t2, p, [t1, t1], [t2, t2]))
                    .collect(),
                _ => {
                    use crate::inttools::edge_edge::EdgeEdgeIntersector;
                    let mut ee = EdgeEdgeIntersector::new();
                    ee.set_edges(nE1, sr1, nE2, sr2, self.ds);
                    ee.set_fuzzy_value(tol);
                    ee.perform();
                    ee.common_parts()
                        .iter()
                        .filter_map(|cp| {
                            Some((
                                cp.vertex_param1,
                                cp.vertex_param2,
                                cp.bounding_point1,
                                cp.range1,
                                *cp.ranges2.first().unwrap_or(&[0.0, 0.0]),
                            ))
                        })
                        .collect()
                }
            };

            // OCCT L303-308: aCPrts = anEdgeEdge.CommonParts(); aNbCPrts = aCPrts.Length();
            // if (!aNbCPrts) continue;
            let a_nb_cprts = hits.len();
            if a_nb_cprts == 0 {
                continue;
            }

            // OCCT L310-322: aPB1->Range(aT11, aT12); HasShrunkData -> aTS11/12
            let aT11 = task.aT11.min(task.aT12);
            let aT12 = task.aT11.max(task.aT12);
            let aTS11 = task.aTS11.min(task.aTS12);
            let aTS12 = task.aTS11.max(task.aTS12);

            // OCCT L324-336: aPB2->Range(aT21, aT22); HasShrunkData -> aTS21/22
            let aT21 = task.aT21.min(task.aT22);
            let aT22 = task.aT21.max(task.aT22);
            let aTS21 = task.aTS21.min(task.aTS22);
            let aTS22 = task.aTS21.max(task.aTS22);

            // OCCT L339: IntTools_Range aR11(aT11, aTS11), aR12(aTS12, aT12), aR21(aT21, aTS21), aR22(aTS22, aT22);
            let a_r11 = [aT11, aTS11];
            let a_r12 = [aTS12, aT12];
            let a_r21 = [aT21, aTS21];
            let a_r22 = [aTS22, aT22];

            // OCCT L341-353: bAnalytical = (Line && Circle)
            let b_analytical = matches!(
                (&e1_curve, &e2_curve),
                (Curve3::Line(_), Curve3::Circle(_)) | (Curve3::Circle(_), Curve3::Line(_))
            );

            // OCCT L355: for (i = 1; i <= aNbCPrts; ++i)
            for i in 0..a_nb_cprts {
                // OCCT L357-361: get CommonPart
                let (a_t1, a_t2, a_pnew, a_cr1_range, a_cr2_range) = hits[i];
                let a_cr1_first = a_cr1_range[0];
                let a_cr1_last = a_cr1_range[1];
                let a_cr2_first = a_cr2_range[0];
                let a_cr2_last = a_cr2_range[1];

                // OCCT L366-368: aType = aCPart.Type()
                // rcad: determine type based on collinearity
                // OCCT L367: switch (aType)

                // Check if this is EDGE-type (collinear overlap)
                let b_is_edge_type = match (&e1_curve, &e2_curve) {
                    (Curve3::Line(l1), Curve3::Line(l2)) => {
                        let cross = l1.direction.cross(l2.direction);
                        cross.length_squared() <= tol * tol
                            && (l2.origin - l1.origin).cross(l1.direction).length() <= tol
                    }
                    _ => false,
                };

                if b_is_edge_type {
                    // OCCT L530-531: if (aNbCPrts > 1) break;
                    if a_nb_cprts > 1 {
                        continue;
                    }
                    // OCCT L529-550: case TopAbs_EDGE:
                    // OCCT L530-531: if (aNbCPrts > 1) break;
                    // (rcad: collinear lines produce one hit, so aNbCPrts == 1 always)
                    // OCCT L535-536: bHasSameBounds = aPB1->HasSameBounds(aPB2)
                    let b_has_same_bounds = task.nV11 == task.nV21 && task.nV12 == task.nV22
                        || task.nV11 == task.nV22 && task.nV12 == task.nV21;
                    if !b_has_same_bounds {
                        // OCCT L537-539: if (!bHasSameBounds) break;
                        continue;
                    }
                    // OCCT L541-549: create InterfEE + FillMap for common blocks
                    self.ds.interf_ee.push(InterferenceEE {
                        e1: nE1,
                        e2: nE2,
                        point: self.ds.edges[nE1].curve.point_at(a_t1),
                        param1: a_t1,
                        param2: a_t2,
                        new_vertex: usize::MAX,
                        range1: [a_t1, a_t1],
                        range2: [a_t2, a_t2],
                    });
                    a_m_edges.insert(nE1);
                    a_m_edges.insert(nE2);
                    // OCCT L549: BOPAlgo_Tools::FillMap(aPB1, aPB2, aMPBLPB, aAllocator)
                    // rcad: fill bidirectional map
                    let pb1_key = Arc::as_ptr(&task.pb1.0) as usize;
                    let pb2_key = Arc::as_ptr(&task.pb2.0) as usize;
                    a_mplpb.entry(pb1_key).or_default().push(pb2_key);
                    a_mplpb.entry(pb2_key).or_default().push(pb1_key);
                    continue;
                }

                // OCCT L369-526: case TopAbs_VERTEX:

                // OCCT L370-373: if (!bIsPBSplittable1 || !bIsPBSplittable2) continue
                if !task.b_is_pb_splittable1 || !task.b_is_pb_splittable2 {
                    continue;
                }

                // OCCT L374-382: VertexParameters, aTol, aCR1, aCR2
                let a_tol_confusion = crate::tolerance::CONFUSION;

                // OCCT L386-394: bIsOnPave[0..3]
                let b0 = crate::boptools::is_on_pave(a_t1, a_r11, a_tol_confusion)
                    || (a_r11[0] - a_t1).abs() <= a_tol_confusion;
                let b1 = crate::boptools::is_on_pave(a_t1, a_r12, a_tol_confusion)
                    || (a_r12[1] - a_t1).abs() <= a_tol_confusion;
                let b2 = crate::boptools::is_on_pave(a_t2, a_r21, a_tol_confusion)
                    || (a_r21[0] - a_t2).abs() <= a_tol_confusion;
                let b3 = crate::boptools::is_on_pave(a_t2, a_r22, a_tol_confusion)
                    || (a_r22[1] - a_t2).abs() <= a_tol_confusion;

                // OCCT L396-397: aPB1->Indices(nV[0], nV[1]); aPB2->Indices(nV[2], nV[3]);
                let n_v = [task.nV11, task.nV12, task.nV21, task.nV22];

                // OCCT L399-403: if both sides on pave -> continue
                if (b0 && b2) || (b0 && b3) || (b1 && b2) || (b1 && b3) {
                    continue;
                }

                // OCCT L405-417: ForceInterfVE for individual on-pave
                // OCCT: for (j = 0; j < 4; ++j) { if (bIsOnPave[j]) { bIsOnPave[j] = ForceInterfVE(nV[j], aPB, aMEdges); ... } }
                let mut b_is_on_pave = [b0, b1, b2, b3];
                let mut is_v_exists = false;
                for j in 0..4 {
                    if b_is_on_pave[j] {
                        let a_pb = if j < 2 {
                            task.pb1.clone()
                        } else {
                            task.pb2.clone()
                        };
                        b_is_on_pave[j] = self.force_interf_ve(n_v[j], &a_pb, &mut a_m_edges);
                        if b_is_on_pave[j] {
                            is_v_exists = true;
                        }
                    }
                }

                // OCCT L419-451: MakeNewVertex + isVExists check
                // rcad: a_pnew is the intersection point from the hit
                if is_v_exists {
                    // OCCT L422-436: check if this is a real intersection or just touching
                    let e1_p = self.ds.edges[nE1].curve.point_at(a_t1);
                    let e2_p = self.ds.edges[nE2].curve.point_at(a_t2);
                    if (e1_p - e2_p).length() > crate::tolerance::INTERSECTION {
                        continue;
                    }
                    // OCCT L440-451: UpdateVertex for each on-pave vertex
                    for j in 0..4 {
                        if b_is_on_pave[j] {
                            let a_v = self.ds.vertex_point(n_v[j]);
                            let a_dist_pp = (a_pnew - a_v).length();
                            self.update_vertex(n_v[j], a_dist_pp);
                        }
                    }
                }

                // OCCT L454: double aTolVnew = BRep_Tool::Tolerance(aVnew);
                // OCCT BOPTools_AlgoTools::MakeNewVertex (AlgoTools_2.cxx L224-250):
                //   aDist = aPnt1.Distance(aPnt2); aMaxTol = max(tol1, tol2) + 0.5 * aDist
                let e1_pt = self.ds.edges[nE1].curve.point_at(a_t1);
                let e2_pt = self.ds.edges[nE2].curve.point_at(a_t2);
                let a_dist = (e1_pt - e2_pt).length();
                let mut a_tol_vnew =
                    self.ds.edge_tolerance(nE1).max(self.ds.edge_tolerance(nE2)) + 0.5 * a_dist;

                // OCCT L455-466: bAnalytical tolerance adjustment
                // OCCT: if (bAnalytical) {
                //   aTolMin = BRepAdaptor_Curve(aE1).GetType() == GeomAbs_Line
                //           ? (aCR1.Last() - aCR1.First()) / 2.
                //           : (aCR2.Last() - aCR2.First()) / 2.;
                if b_analytical {
                    let a_tol_min = match &e1_curve {
                        Curve3::Line(_) => (a_cr1_last - a_cr1_first).abs() * 0.5,
                        _ => (a_cr2_last - a_cr2_first).abs() * 0.5,
                    };
                    if a_tol_min > a_tol_vnew {
                        a_tol_vnew = a_tol_min;
                    }
                }

                // OCCT L468-510: vertex coincidence check (100x tolerance)
                let mut i_found = false;
                let mut a_mv: std::collections::HashSet<usize> = std::collections::HashSet::new();
                a_mv.insert(n_v[0]);
                a_mv.insert(n_v[1]);
                let mut n_vs: [usize; 2] = [usize::MAX; 2];
                let mut j: i32 = -1;
                if a_mv.contains(&n_v[2]) {
                    j += 1;
                    n_vs[j as usize] = n_v[2];
                }
                if a_mv.contains(&n_v[3]) {
                    j += 1;
                    n_vs[j as usize] = n_v[3];
                }
                // OCCT L490-504: check each shared vertex
                for k1 in 0..=j {
                    let n_vx = n_vs[k1 as usize];
                    if n_vx < self.ds.vertices.len() {
                        let a_px = self.ds.vertex_point(n_vx);
                        let a_tol_vx = self.ds.vertex_tolerance(n_vx);
                        let a_d2: f64 = (a_pnew - a_px).length_squared();
                        let a_dt2 = 100.0 * (a_tol_vnew + a_tol_vx) * (a_tol_vnew + a_tol_vx);
                        if a_d2 < a_dt2 {
                            i_found = true;
                            break;
                        }
                    }
                }
                // OCCT L506-509: if (iFound) continue;
                if i_found {
                    continue;
                }

                // OCCT L512-516: BOPDS_InterfEE& aEE = aEEs.Appended();
                // iX = aEEs.Length() - 1; aEE.SetIndices(nE1, nE2); aEE.SetCommonPart(aCPart);
                let new_v = self.ds.add_vertex(a_pnew);
                self.ds.interf_ee.push(InterferenceEE {
                    e1: nE1,
                    e2: nE2,
                    point: a_pnew,
                    param1: a_t1,
                    param2: a_t2,
                    new_vertex: new_v,
                    range1: [a_t1, a_t1],
                    range2: [a_t2, a_t2],
                });
                let i_x = self.ds.interf_ee.len() - 1;

                // OCCT L518: myDS->AddInterf(nE1, nE2);
                a_m_edges.insert(nE1);
                a_m_edges.insert(nE2);

                // OCCT L520-525: aCPB.SetPaveBlocks(aPB1,aPB2); aCPB.SetIndexInterf(iX);
                // aCPB.SetTolerance(aTolVnew); aMVCPB.Add(aVnew, aCPB);
                a_mvcpb.insert(
                    new_v,
                    crate::bopds::pave::CoupleOfPaveBlocks {
                        interf_idx: i_x,
                        vertex_index: new_v,
                        pb1: task.pb1.clone(),
                        pb2: task.pb2.clone(),
                        tolerance: a_tol_vnew,
                    },
                );
            } // for (i=0; i<aNbCPrts; i++)
        } // for (k=0; k<aNbEdgeEdge; k++)

        // OCCT L558-560: PerformCommonBlocks + UpdateVerticesOfCB
        // (rcad: aMPBLPB populated above, passed to perform_common_blocks)
        if !a_mplpb.is_empty() {
            // rcad: existing perform_common_blocks scans all edges.
            crate::bopds::tools::perform_common_blocks(&mut self.ds);
        }
        self.update_vertices_of_cb();

        // OCCT L565: PerformNewVertices
        self.perform_new_vertices(a_mvcpb, true);

        // OCCT L571-585: SplitPaveBlocks
        if !a_m_edges.is_empty() {
            self.split_pave_blocks(&a_m_edges, false);
        }
    }

    /// OCCT BOPAlgo_PaveFiller::ForceInterfVE (PaveFiller_3.cxx L828-910).
    /// Forces an interference between vertex nV and the edge of aPB.
    /// Returns true if the vertex was already on the edge or if a new
    /// VE interference was created.
    pub(crate) fn force_interf_ve(
        &mut self,
        n_v: usize,
        a_pb: &crate::bopds::pave::SharedPB,
        the_m_edges: &mut std::collections::HashSet<usize>,
    ) -> bool {
        let pb = a_pb.0.read().unwrap();
        let n_e = pb.original_edge;
        drop(pb);
        if self.ds.edge_has_vertex(n_v, n_e) {
            return true;
        }
        if self.ds.has_interf_ve(n_v, n_e) {
            return true;
        }
        if self.ds.has_interf_ve_via_faces(n_v, n_e) {
            return true;
        }
        {
            let pb = a_pb.0.read().unwrap();
            if pb.pave1.vertex_idx == n_v || pb.pave2.vertex_idx == n_v {
                return true;
            }
        }
        let n_vx = self.ds.has_shape_sd(n_v).unwrap_or(n_v);
        let res = match self
            .context
            .compute_ve(self.ds, n_vx, n_e, self.fuzzy_tolerance)
        {
            Ok(r) => r,
            Err(_) => return false,
        };
        let a_t = res.param;
        let a_tol_vnew = res.tolerance;
        let a_pave = crate::bopds::pave::Pave {
            vertex_idx: n_vx,
            param: a_t,
        };
        {
            let mut pb = a_pb.0.write().unwrap();
            pb.append_ext_pave(a_pave);
        }
        self.ds
            .interf_ve
            .push(crate::bopds::ds::types::InterferenceVE {
                vertex: n_v,
                edge: n_e,
                param: a_t,
                index_new: n_vx,
            });
        self.ds.try_add_interf(n_v, n_e);
        self.update_vertex(n_v, a_tol_vnew);
        the_m_edges.insert(n_e);
        return true;
    }

    /// OCCT: PaveBlock range extraction (GetPBBox equivalent)
    pub(crate) fn get_pb_boxes(ds: &DS, edge_idx: usize, edge_t_range: [f64; 2]) -> Vec<[f64; 2]> {
        let paves = &ds.edge_paves(edge_idx);
        if paves.is_empty() {
            return vec![edge_t_range];
        }
        let mut params: Vec<f64> = paves
            .iter()
            .map(|p| p.param)
            .filter(|p| p.is_finite())
            .collect();
        params.sort_by(|a, b| a.partial_cmp(b).unwrap());
        params.dedup();
        let tol = ds
            .edge_tolerance(edge_idx)
            .max(crate::tolerance::TOLERANCE_ABS);
        let mut ranges = Vec::new();
        let mut prev = edge_t_range[0];
        for &p in &params {
            if (p - prev).abs() > tol {
                ranges.push([prev, p]);
            }
            prev = p;
        }
        if (edge_t_range[1] - prev).abs() > tol {
            ranges.push([prev, edge_t_range[1]]);
        }
        ranges
    }
    // OCCT BOPAlgo_PaveFiller_4.cxx L139-301: PerformVF
    pub(crate) fn perform_vf(&mut self) {
        // OCCT L141: myIterator->Initialize(VERTEX, FACE)
        // OCCT L142: iSize = myIterator->ExpectedLength()
        self.my_iterator
            .initialize(ShapeType::Vertex, ShapeType::Face);
        let i_size = self.my_iterator.expected_length();
        //
        // OCCT L147-160: myGlue == GlueFull
        if self.glue == GlueEnum::GlueFull {
            for _ in 0..i_size {
                let (nV, nF) = self.my_iterator.value();
                self.my_iterator.next();
                if !self.ds.is_sub_shape(nV, nF) {
                    self.ds.face_info_mut(nF);
                }
            }
            return;
        }
        //
        // OCCT L162: NCollection_DynamicArray<BOPDS_InterfVF>& aVFs = myDS->InterfVF()
        // OCCT L163-170: if (!iSize)
        if i_size == 0 {
            // OCCT L165: iSize = 10
            // OCCT L166: aVFs.SetIncrement(iSize)
            self.ds.interf_vf.reserve(10);
            // OCCT L168: TreatVerticesEE()
            self.treat_vertices_ee();
            return;
        }
        //
        // OCCT L172-174: variable declarations
        // OCCT L174: BOPAlgo_VectorOfVertexFace aVVF
        let mut a_vv_f: Vec<VertexFaceTask> = Vec::new();
        //
        // OCCT L176: aVFs.SetIncrement(iSize)
        self.ds.interf_vf.reserve(i_size);
        //
        // OCCT L178-180: NCollection_DataMap<BOPDS_Pair, NCollection_Map<int>> aMVFPairs
        let mut a_mvf_pairs: std::collections::HashMap<(usize, usize), Vec<usize>> =
            std::collections::HashMap::new();
        //
        // OCCT L181: for (; myIterator->More(); myIterator->Next())
        for _ in 0..i_size {
            // OCCT L183-186: UserBreak
            //
            // OCCT L187: myIterator->Value(nV, nF)
            let (nV, nF) = self.my_iterator.value();
            self.my_iterator.next();
            //
            // OCCT L189-192: IsSubShape
            let flat_f = self.ds.face_shape_idx[nF];
            let flat_v = self.ds.vertex_shape_idx[nV];
            if self.ds.is_sub_shape(flat_v, flat_f) {
                continue;
            }
            //
            // OCCT L194-197: if (myDS->HasInterf(nV, nF)) continue;
            //   (type-specific indices — matches rcad's interference matrix)
            if self.ds.has_interf(nV, nF) {
                continue;
            }
            //
            // OCCT L199: myDS->ChangeFaceInfo(nF)
            self.ds.face_info_mut(nF);
            //
            // OCCT L200-203: HasInterfShapeSubShapes(nV, nF)
            //   (architecture: face sub_shapes are wire indices; V-Wire interference
            //    never registered in OCCT either — effectively always false)
            //
            // OCCT L205-209: SD resolution
            let nVx = self.ds.has_shape_sd(nV).unwrap_or(nV);
            //
            // OCCT L211-220: aMVFPairs dedup (key = nVx, nF)
            let key = (nVx, nF);
            let entry = a_mvf_pairs.entry(key).or_default();
            entry.push(nV);
            if entry.len() > 1 {
                // OCCT L216: continue - already have a task for this SD-face pair
                continue;
            }
            //
            // OCCT L222-230: Create BOPAlgo_VertexFace task
            let mut a_task = VertexFaceTask::new();
            a_task.set_indices(nVx, nF);
            a_vv_f.push(a_task);
        } // for (; myIterator->More(); myIterator->Next()) {
          //
          // OCCT L234-240: SetProgressRange
          // OCCT L242: BOPTools_Parallel::Perform
        for task in &mut a_vv_f {
            task.perform(&mut self.context, self.ds, self.fuzzy_tolerance);
        }
        // OCCT L244-247: UserBreak check (not ported)
        //
        // OCCT L249: for (k = 0; k < aNbVF; ++k)
        for task in &a_vv_f {
            // OCCT L257: iFlag = aVertexFace.Flag();
            if task.flag() != 0 {
                continue;
            }
            //
            // OCCT L268-270: Indices, Parameters, VertexNewTolerance
            let (mut nVx, nF) = task.indices();
            let (a_t1, a_t2) = task.parameters();
            let a_tol_vnew = task.vertex_new_tolerance();
            //
            // OCCT L272-273: aMVFPairs.Find(aVFPair) - get all original vertices
            let key = (nVx, nF);
            let orig_verts: Vec<usize> = match a_mvf_pairs.get(&key) {
                Some(v) => v.clone(),
                None => vec![nVx],
            };
            // OCCT L275: for (; itMV.More(); itMV.Next())
            for &nV in &orig_verts {
                // OCCT L279-281: BOPDS_InterfVF
                let a_vf_idx = self.ds.interf_vf.len();
                self.ds.interf_vf.push(InterferenceVF {
                    vertex: nV,
                    face: nF,
                    u: a_t1,
                    v: a_t2,
                    index_new: None,
                });
                // OCCT L283: myDS->AddInterf(nV, nF);
                self.ds.try_add_interf(nV, nF);
                // OCCT L286: nVx = UpdateVertex(nV, aTolVNew);
                nVx = self.update_vertex(nV, a_tol_vnew);
                // OCCT L289-292: IsNewShape -> SetIndexNew
                if self.ds.is_new_vertex(nVx) {
                    if let Some(a_vf) = self.ds.interf_vf.get_mut(a_vf_idx) {
                        a_vf.index_new = Some(nVx);
                    }
                }
            }
            // OCCT L295-297: FaceInfo VerticesIn
            let a_fi = self.ds.face_info_mut(nF);
            a_fi.vertices_in.insert(nVx);
        }
        //
        // OCCT L300: TreatVerticesEE()
        self.treat_vertices_ee();
    }
    /// OCCT BOPAlgo_PaveFiller_5.cxx L165-592: PerformEF
    pub(crate) fn perform_ef(&mut self, pairs: &[(usize, usize)]) {
        // OCCT L167: FillShrunkData(TopAbs_EDGE, TopAbs_FACE)
        self.fill_shrunk_data(ShapeType::Edge, ShapeType::Face);
        // OCCT L169: myIterator->Initialize(EDGE, FACE)
        // OCCT L171: iSize = myIterator->ExpectedLength()
        let i_size = pairs.len();
        if i_size == 0 {
            return;
        }
        // OCCT L177: int nE, nF
        //
        // OCCT L179-192: myGlue == GlueFull
        if self.glue == GlueEnum::GlueFull {
            for &(nE, nF) in pairs {
                // OCCT L186: if (!myDS->ShapeInfo(nE).HasFlag())
                if !self.ds.edge_has_flag(nE) {
                    // OCCT L188: myDS->ChangeFaceInfo(nF)
                    self.ds.face_info_mut(nF);
                }
            }
            return;
        }
        //
        // OCCT L194-214: variable declarations
        let mut a_mi_efc: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_v_edge_face: Vec<EdgeFaceTask> = Vec::new();
        // OCCT L209-210: NCollection_IndexedDataMap<TopoDS_Shape, CoupleOfPaveBlocks> aMVCPB
        let mut a_mvcpb: std::collections::HashMap<usize, crate::bopds::pave::CoupleOfPaveBlocks> =
            std::collections::HashMap::new();
        // OCCT L216-217: InterfEF
        self.ds.interf_ef.reserve(i_size);
        //
        // OCCT L219-307: for (; myIterator->More(); myIterator->Next())
        for &(nE, nF) in pairs {
            if self.ds.edge_has_flag(nE) {
                continue;
            }
            let a_tol_e = self.ds.edge_tolerance(nE);
            let a_tol_f = self.ds.face_tolerance(nF);
            let a_fi = self.ds.face_info_mut(nF);
            let a_mv_in: std::collections::HashSet<usize> =
                a_fi.vertices_in.iter().copied().collect();
            let a_mv_on: std::collections::HashSet<usize> =
                a_fi.vertices_on.iter().copied().collect();
            let a_mpbf: Vec<usize> = a_fi.pave_blocks_on.iter().copied().collect();
            drop(a_fi);
            let n_pave_blocks = self.ds.edge_pave_blocks(nE).len();
            if n_pave_blocks == 0 {
                continue;
            }
            for pb_local_idx in 0..n_pave_blocks {
                let pb = &self.ds.edge_pave_blocks(nE)[pb_local_idx];
                let pb_ref = pb.0.read().unwrap();
                let pb_key = pb_ref.new_edge.unwrap_or(pb_ref.original_edge);
                if a_mpbf.contains(&pb_key) {
                    continue;
                }
                let (aT1, aT2) = pb_ref.range();
                let (aTS1, aTS2) = match pb_ref.shrunk_range {
                    Some(sr) => (sr[0], sr[1]),
                    None => continue,
                };
                let (nV1, nV2) = pb_ref.indices();
                // OCCT L268-271: AABB overlap check (aBBF.IsOut(aBBE))
                let a_bbf = {
                    let si = self.ds.face_shape_idx[nF];
                    let info = &self.ds.shape_info[si];
                    info.box_min.zip(info.box_max)
                };
                let a_bbe = pb_ref.my_shrunk_box;
                let ((f_min, f_max), (e_min, e_max)) = match (a_bbf, a_bbe) {
                    (Some(fb), Some(eb)) => (fb, eb),
                    _ => continue,
                };
                if f_min.x > e_max.x
                    || f_max.x < e_min.x
                    || f_min.y > e_max.y
                    || f_max.y < e_min.y
                    || f_min.z > e_max.z
                    || f_max.z < e_min.z
                {
                    continue;
                }
                let b_express_compute = (a_mv_in.contains(&nV1) || a_mv_on.contains(&nV1))
                    && (a_mv_in.contains(&nV2) || a_mv_on.contains(&nV2));
                drop(pb_ref);
                let tol_ef = a_tol_e.max(a_tol_f).max(CONFUSION);
                let a_pb_range = [aT1.min(aT2), aT1.max(aT2)];
                let _a_sr_corrected = Self::correct_range_for_face(
                    &self.ds.edges[nE].curve,
                    tol_ef,
                    [aTS1.min(aTS2), aTS1.max(aTS2)],
                );
                let a_pb_corrected =
                    Self::correct_range_for_face(&self.ds.edges[nE].curve, tol_ef, a_pb_range);
                if a_pb_corrected[1] - a_pb_corrected[0] <= tol_ef {
                    continue;
                }
                self.fpbdone.entry(nF).or_default().insert(pb_key);
                let mut a_task = EdgeFaceTask::new();
                a_task.set_indices(nE, nF);
                a_task.set_pave_block(pb_local_idx);
                a_task.set_fuzzy_value(self.fuzzy_tolerance);
                a_task.use_quick_coincidence_check(b_express_compute);
                a_task.set_new_sr([aTS1, aTS2]);
                a_task.set_range(a_pb_range);
                a_v_edge_face.push(a_task);
            }
        }
        //
        // OCCT L309-317: BOPTools_Parallel::Perform
        for task in &mut a_v_edge_face {
            task.perform(&mut self.context, self.ds);
        }
        //
        // OCCT L324-571: Process results
        for k in 0..a_v_edge_face.len() {
            let task = &mut a_v_edge_face[k];
            if !task.is_done() || task.has_errors() {
                continue;
            }
            let (nE, nF) = task.indices();
            let a_tol_e = self.ds.edge_tolerance(nE);
            let a_tol_f = self.ds.face_tolerance(nF);
            let a_cprts = task.common_parts();
            let a_nb_cprts = a_cprts.len();
            if a_nb_cprts == 0 {
                let md = task.minimal_distance();
                if md < f64::MAX && md > a_tol_e + a_tol_f {
                    let t_range = task.myRange;
                    self.distances.entry((nE, nF)).or_default().push(
                        crate::pave_filler::EdgeRangeDistance {
                            first: t_range[0],
                            last: t_range[1],
                            distance: md,
                        },
                    );
                }
                continue;
            }
            let a_pb_local = task.pave_block();
            let nV = {
                let pb = &self.ds.edge_pave_blocks(nE)[a_pb_local];
                let pb_ref = pb.0.read().unwrap();
                let (a, b) = pb_ref.indices();
                [a, b]
            };
            let mut new_sr = task.myNewSR;
            let pb_range = task.myRange;
            // OCCT L373-380: ReduceIntersectionRange for VERTEX type
            if a_nb_cprts > 0 {
                if let EfHit::Vertex { .. } = &a_cprts[0] {
                    let mut ts1 = new_sr[0];
                    let mut ts2 = new_sr[1];
                    self.reduce_ef_intersection_range(nV[0], nV[1], nE, nF, &mut ts1, &mut ts2);
                    new_sr = [ts1, ts2];
                }
            }
            let a_r1 = [pb_range[0].min(new_sr[0]), pb_range[0].max(new_sr[0])];
            let a_r2 = [new_sr[1].min(pb_range[1]), new_sr[1].max(pb_range[1])];
            let a_fi = self.ds.face_info_mut(nF);
            let a_mif_on: std::collections::HashSet<usize> =
                a_fi.vertices_on.iter().copied().collect();
            let a_mif_in: std::collections::HashSet<usize> =
                a_fi.vertices_in.iter().copied().collect();
            drop(a_fi);
            let b_line_plane = matches!(
                (&self.ds.edges[nE].curve, &self.ds.faces[nF].surface),
                (Curve3::Line(_), Surface3::Plane(_))
            );
            for i in 0..a_nb_cprts {
                let a_cpart = &a_cprts[i];
                match a_cpart {
                    EfHit::Vertex { point, param: a_t } => {
                        // OCCT L406-543: TopAbs_VERTEX
                        let a_tol_to_decide = 5e-8;
                        let mut b_is_on_pave = [false, false];
                        b_is_on_pave[0] = (a_t - a_r1[0]).abs() <= a_tol_to_decide
                            || (a_t - a_r1[1]).abs() <= a_tol_to_decide;
                        b_is_on_pave[1] = (a_t - a_r2[0]).abs() <= a_tol_to_decide
                            || (a_t - a_r2[1]).abs() <= a_tol_to_decide;
                        // OCCT L421-439: both on pave or (bLinePlane && one on pave)
                        if (b_is_on_pave[0] && b_is_on_pave[1])
                            || (b_line_plane && (b_is_on_pave[0] || b_is_on_pave[1]))
                        {
                            // OCCT L423-425: CheckFacePaves(nV[0..1], aMIFOn, aMIFIn)
                            let bv0 = a_mif_on.contains(&nV[0]) || a_mif_in.contains(&nV[0]);
                            let bv1 = a_mif_on.contains(&nV[1]) || a_mif_in.contains(&nV[1]);
                            if bv0 && bv1 {
                                // OCCT L427-437: EDGE-type treatment — set CommonPart
                                self.ds.interf_ef.push(InterferenceEF {
                                    edge: nE,
                                    face: nF,
                                    point: *point,
                                    edge_param: *a_t,
                                    new_vertex: 0,
                                });
                                self.ds.try_add_interf(nE, nF);
                                a_mi_efc.insert(nF);
                                continue;
                            }
                            // OCCT L448-455: one vertex NOT on face, just AddInterf
                            self.ds.try_add_interf(nE, nF);
                        }
                        // OCCT L442-444: splittable — if PB not splittable, skip
                        if !self.ds.edge_pave_blocks(nE)[a_pb_local]
                            .0
                            .read()
                            .unwrap()
                            .is_splittable
                        {
                            continue;
                        }
                        // OCCT L447-457: ForceInterfVF for on-pave vertices
                        for j in 0..2 {
                            if b_is_on_pave[j] {
                                let bv = a_mif_on.contains(&nV[j]) || a_mif_in.contains(&nV[j]);
                                if !bv {
                                    b_is_on_pave[j] = self.force_interf_vf_pair(nV[j], nF);
                                }
                            }
                        }
                        // OCCT L459-502: real intersection check for on-pave vertices
                        if b_is_on_pave[0] || b_is_on_pave[1] {
                            for j in 0..2 {
                                if b_is_on_pave[j] {
                                    let dist_pp = (self.ds.vertex_point(nV[j]) - *point).length();
                                    let a_tol = self.ds.vertex_tolerance(nV[j]);
                                    let mut a_max_dist = 1e4 * a_tol;
                                    if a_tol < 0.01 {
                                        a_max_dist = a_max_dist.min(0.1);
                                    }
                                    if dist_pp < a_max_dist {
                                        self.update_vertex(nV[j], dist_pp);
                                        self.verts_to_avoid_extension.insert(nV[j]);
                                    }
                                }
                            }
                            continue;
                        }
                        // OCCT L505-508: CheckFacePaves(aVnew, aMIFOn)
                        {
                            let near_face_vx = a_mif_on.iter().chain(a_mif_in.iter()).any(|&vi| {
                                vi < self.ds.vertices.len()
                                    && (self.ds.vertex_point(vi) - *point).length()
                                        <= a_tol_e.max(a_tol_f).max(CONFUSION * 10.0)
                            });
                            if near_face_vx {
                                continue;
                            }
                        }
                        // OCCT L510-526: aTolVnew + IsPointInFace
                        let mut a_tol_vnew = a_tol_e.max(a_tol_f);
                        if b_line_plane {
                            // OCCT L513-518: increase tolerance for Line/Plane
                            a_tol_vnew = a_tol_vnew.max(CONFUSION * 100.0);
                        }
                        // OCCT L523-526: IsPointInFace 3D
                        if !self.is_point_in_face(*point, nF, a_tol_vnew) {
                            continue;
                        }
                        // OCCT L528-542: Create EF interference + CPB
                        a_mi_efc.insert(nF);
                        let new_v = self.ds.add_vertex_no_dedup(*point);
                        let interf_idx = self.ds.interf_ef.len();
                        let pb_for_cpb = self.ds.edges[nE].pave_blocks[a_pb_local].clone();
                        self.ds.interf_ef.push(InterferenceEF {
                            edge: nE,
                            face: nF,
                            point: *point,
                            edge_param: *a_t,
                            new_vertex: new_v,
                        });
                        self.ds.try_add_interf(nE, nF);
                        a_mvcpb.insert(
                            new_v,
                            crate::bopds::pave::CoupleOfPaveBlocks {
                                interf_idx,
                                vertex_index: new_v,
                                pb1: pb_for_cpb.clone(),
                                pb2: pb_for_cpb,
                                tolerance: a_tol_vnew,
                            },
                        );
                        self.ds.faces[nF].face_info.vertices_on.insert(new_v);
                        if nE < self.ds.edge_paves.len() {
                            self.add_pave_to_edge(
                                nE,
                                Pave {
                                    vertex_idx: new_v,
                                    param: *a_t,
                                },
                            );
                        }
                    }
                    EfHit::Edge { t1, t2 } => {
                        a_mi_efc.insert(nF);
                        let mid_t = (t1 + t2) * 0.5;
                        let mid_pt = self.ds.edges[nE].curve.point_at(mid_t);
                        let bv0 = a_mif_on.contains(&nV[0]) || a_mif_in.contains(&nV[0]);
                        let bv1 = a_mif_on.contains(&nV[1]) || a_mif_in.contains(&nV[1]);
                        self.ds.interf_ef.push(InterferenceEF {
                            edge: nE,
                            face: nF,
                            point: mid_pt,
                            edge_param: mid_t,
                            new_vertex: 0,
                        });
                        self.ds.try_add_interf(nE, nF);
                    }
                }
            }
        }
        //
        // ==================================================================
        // Phase 4: Post-treatment (OCCT L576-592)
        // ==================================================================
        // OCCT L576: BOPAlgo_Tools::PerformCommonBlocks(aMPBLI, ...)
        // OCCT L577: UpdateVerticesOfCB()
        // OCCT L578: PerformNewVertices(aMVCPB, ...)
        // OCCT L585: myDS->UpdateFaceInfoIn(aMIEFC)
        if !a_v_edge_face.is_empty() {
            crate::bopds::tools::perform_common_blocks(&mut self.ds);
        }
        self.update_vertices_of_cb();
        // OCCT L578+L612+L687: PerformNewVertices → TreatNewVertices → IntersectVE
        self.perform_new_vertices(a_mvcpb, false);
        for &fi in &a_mi_efc {
            self.ds.update_face_info_in(fi);
        }
    }
    /// OCCT PaveFiller_5.cxx L685-768: ReduceIntersectionRange
    fn reduce_ef_intersection_range(
        &self,
        the_v1: usize,
        the_v2: usize,
        the_e: usize,
        the_f: usize,
        the_ts1: &mut f64,
        the_ts2: &mut f64,
    ) {
        if !self.ds.is_new_vertex(the_v1) && !self.ds.is_new_vertex(the_v2) {
            return;
        }
        let has_interf_shape_sub_shapes = {
            let face_edges: std::collections::HashSet<usize> =
                self.ds.face_boundary_edges(the_f).iter().copied().collect();
            self.ds.interf_ee.iter().any(|inf| {
                let involves_our_edge = inf.e1 == the_e || inf.e2 == the_e;
                let involves_face_edge =
                    face_edges.contains(&inf.e1) || face_edges.contains(&inf.e2);
                involves_our_edge && involves_face_edge
            })
        };
        if !has_interf_shape_sub_shapes {
            return;
        }
        let a_nb_ees = self.ds.interf_ee.len();
        if a_nb_ees == 0 {
            return;
        }
        let face_edges: std::collections::HashSet<usize> =
            self.ds.face_boundary_edges(the_f).iter().copied().collect();
        for inf in &self.ds.interf_ee {
            let nv = inf.new_vertex;
            if nv != the_v1 && nv != the_v2 {
                continue;
            }
            let involves_our_edge = inf.e1 == the_e || inf.e2 == the_e;
            let involves_face_edge = face_edges.contains(&inf.e1) || face_edges.contains(&inf.e2);
            if !involves_our_edge || !involves_face_edge {
                continue;
            }
            let ee_param = if inf.e1 == the_e {
                inf.param1
            } else {
                inf.param2
            };
            if nv == the_v1 {
                if *the_ts1 < ee_param {
                    *the_ts1 = ee_param;
                }
            } else {
                if *the_ts2 > ee_param {
                    *the_ts2 = ee_param;
                }
            }
        }
    }
    /// OCCT BOPDS_Iterator: face BVH construction
    pub(crate) fn build_box_tree_face(&self, is_a: bool) -> crate::bvh::BoxTree {
        use crate::bvh::{Aabb, BoxTree};
        let (start, end) = if is_a {
            (0, self.ds.a_face_count)
        } else {
            (self.ds.a_face_count, self.ds.faces.len())
        };
        let n = end - start;
        let mut indices = Vec::with_capacity(n);
        let mut aabbs = Vec::with_capacity(n);
        for local_i in 0..n {
            let fi = start + local_i;
            indices.push(fi);
            let f = &self.ds.faces[fi];
            let mut aabb = Aabb::empty();
            // Boundary vertices
            for &vi in &f.boundary_verts {
                if vi < self.ds.vertices.len() {
                    aabb.expand_point(self.ds.vertex_point(vi));
                }
            }
            // OCCT BndLib_AddSurface: expand AABB for curved surfaces.
            // Sphere: full sphere AABB = center  ?radius (face boundary
            // vertices only cover a patch, not the whole sphere volume).
            // Cylinder/Cone: boundary vertices already span the full
            // parametric extent  ?no extra expansion needed.
            if let Surface3::Sphere(s) = &f.surface {
                let r = s.radius.abs();
                aabb.expand_point(s.center + DVec3::splat(r));
                aabb.expand_point(s.center - DVec3::splat(r));
            }
            let tol = f.geom_tol.max(CONFUSION);
            aabb.min -= DVec3::splat(tol);
            aabb.max += DVec3::splat(tol);
            aabbs.push(aabb);
        }
        BoxTree::build(indices, aabbs)
    }
    ///  BOPDS_Iterator  ?combined face BVH (both operands).
    pub(crate) fn build_box_tree_face_all(&self) -> crate::bvh::BoxTree {
        use crate::bvh::{Aabb, BoxTree};
        let n = self.ds.faces.len();
        let mut indices = Vec::with_capacity(n);
        let mut aabbs = Vec::with_capacity(n);
        for fi in 0..n {
            indices.push(fi);
            let f = &self.ds.faces[fi];
            let mut aabb = Aabb::empty();
            for &vi in &f.boundary_verts {
                if vi < self.ds.vertices.len() {
                    aabb.expand_point(self.ds.vertex_point(vi));
                }
            }
            if let Surface3::Sphere(s) = &f.surface {
                let r = s.radius.abs();
                aabb.expand_point(s.center + DVec3::splat(r));
                aabb.expand_point(s.center - DVec3::splat(r));
            }
            let tol = f.geom_tol.max(CONFUSION);
            aabb.min -= DVec3::splat(tol);
            aabb.max += DVec3::splat(tol);
            aabbs.push(aabb);
        }
        BoxTree::build(indices, aabbs)
    }
    /// rcad glue-mode acceleration (no OCCT equivalent)
    pub(crate) fn should_skip_ve_pass(&self) -> bool {
        if !self.use_glue() {
            return false;
        }

        // If all vertices are shared, skip V-E pass
        let shared_verts = &self.ds.shared_topology.shared_vertices;
        if shared_verts.is_empty() {
            return false;
        }

        // Check if all vertices from shape A have matches in shape B
        let a_verts: std::collections::HashSet<usize> = self
            .ds
            .vertices
            .iter()
            .enumerate()
            .filter(|(_, v)| v.origin == Some(ShapeOrigin::ShapeA))
            .map(|(i, _)| i)
            .collect();

        let matched_a: std::collections::HashSet<usize> =
            shared_verts.iter().map(|(a, _)| *a).collect();

        a_verts == matched_a && !a_verts.is_empty()
    }
    /// rcad glue-mode acceleration (no OCCT equivalent)
    pub(crate) fn should_skip_ee_pass(&self) -> bool {
        if !self.use_glue() {
            return false;
        }

        let shared_edges = &self.ds.shared_topology.shared_edges;
        if shared_edges.is_empty() {
            return false;
        }

        // Check if all edges from shape A have matches in shape B
        let a_edges: std::collections::HashSet<usize> = self
            .ds
            .edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.origin == ShapeOrigin::ShapeA)
            .map(|(i, _)| i)
            .collect();

        let matched_a: std::collections::HashSet<usize> =
            shared_edges.iter().map(|(a, _)| *a).collect();

        a_edges == matched_a && !a_edges.is_empty()
    }
    /// rcad glue-mode acceleration (no OCCT equivalent)
    pub(crate) fn should_skip_vf_pass(&self) -> bool {
        if !self.use_glue() {
            return false;
        }

        // If all faces are fully glued, skip V-F pass
        !self.ds.shared_topology.fully_glued_faces.is_empty()
            && self.ds.shared_topology.fully_glued_faces.len()
                == self.ds.a_face_count * (self.ds.faces.len() - self.ds.a_face_count)
    }
    /// rcad glue-mode acceleration (no OCCT equivalent)
    pub(crate) fn should_skip_ef_pass(&self) -> bool {
        if !self.use_glue() {
            return false;
        }

        // If all faces are fully glued, skip E-F pass
        !self.ds.shared_topology.fully_glued_faces.is_empty()
            && self.ds.shared_topology.fully_glued_faces.len()
                == self.ds.a_face_count * (self.ds.faces.len() - self.ds.a_face_count)
    }
    /// rcad glue-mode acceleration (no OCCT equivalent)
    pub(crate) fn should_skip_ff_pass(&self) -> bool {
        if !self.use_glue() {
            return false;
        }

        // If all faces are fully glued, skip F-F pass
        let total_face_pairs = self.ds.a_face_count * (self.ds.faces.len() - self.ds.a_face_count);
        self.ds.shared_topology.fully_glued_faces.len() == total_face_pairs && total_face_pairs > 0
    }
    // OCCT BOPAlgo_PaveFiller_1.cxx L45-132: PerformVV
    pub(crate) fn perform_vv(&mut self) {
        // L47-48: n1, n2, iFlag, aSize; aAllocator
        // L50-51: myIterator->Initialize(VERTEX, VERTEX); aSize = ExpectedLength()
        self.my_iterator
            .initialize(ShapeType::Vertex, ShapeType::Vertex);
        let a_size = self.my_iterator.expected_length();
        // L52: Message_ProgressScope (rcad: sequential, no progress)
        // L53-56: if (!aSize) return
        if a_size == 0 {
            return;
        }
        // L58-59: myDS->InterfVV().SetIncrement(aSize)
        self.ds.interf_vv.reserve(a_size);
        // L62-64: aAllocator, aMILI, aMBlocks
        //   NCollection_IndexedDataMap<int, NCollection_List<int>> aMILI(100, aAllocator);
        //   NCollection_List<NCollection_List<int>> aMBlocks(aAllocator);
        let mut a_mili: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        // L66-98: 1. Map V/LV
        // L68: Message_ProgressScope aPSLoop (rcad: sequential)
        for _ in 0..a_size {
            // L71-74: UserBreak check (rcad: not ported)
            // L75: myIterator->Value(n1, n2)
            let (n1, n2) = self.my_iterator.value();
            self.my_iterator.next();
            //
            // L77-81: if HasInterf(n1, n2) -> FillMap + continue
            let key = if n1 < n2 { (n1, n2) } else { (n2, n1) };
            if self.ds.interf_tb.contains(&key) {
                fill_map(&mut a_mili, n1, n2);
                continue;
            }
            // L84-88: Resolve SD vertices
            //   int n1SD = n1; myDS->HasShapeSD(n1, n1SD);
            let n1sd: usize = self.ds.has_shape_sd(n1).unwrap_or(n1);
            let n2sd: usize = self.ds.has_shape_sd(n2).unwrap_or(n2);
            // L90-93: ComputeVV(aV1, aV2, myFuzzyValue)
            //   BOPTools_AlgoTools::ComputeVV returns: 0 = interfered, 1 = not
            let a_tol =
                self.ds.vertex_tolerance(n1sd) + self.ds.vertex_tolerance(n2sd) + self.tol();
            let a_sq_dist =
                (self.ds.vertex_point(n1sd) - self.ds.vertex_point(n2sd)).length_squared();
            let i_flag = if a_sq_dist <= a_tol * a_tol { 0 } else { 1 };
            // L94-97: if (!iFlag) -> FillMap(n1, n2)
            if i_flag == 0 {
                fill_map(&mut a_mili, n1, n2);
            }
        }
        // L100-101: 2. Make blocks — BOPAlgo_Tools::MakeBlocks(aMILI, aMBlocks, aAllocator)
        let a_m_blocks: Vec<Vec<usize>> = make_blocks(&a_mili);
        // L103-113: 3. Make SD vertices
        //   NCollection_List<NCollection_List<int>>::Iterator aItB(aMBlocks)
        for block in &a_m_blocks {
            // L107-110: UserBreak check (rcad: not ported)
            // L111-112: MakeSDVertices(aLI)
            self.make_sd_vertices_vv(block);
        }
        // L115-127: 4. InitPaveBlocksForVertex for each SD vertex source
        // L117: myDS->ShapesSD()
        let a_dmii: std::collections::HashSet<usize> = self
            .ds
            .shape_sd
            .sd_vertices_iter()
            .map(|&(k, _)| k)
            .collect();
        // L118-119: aItDMII.Initialize(aDMII)
        for &n1_key in &a_dmii {
            // L121-124: UserBreak check (rcad: not ported)
            // L125-126: myDS->InitPaveBlocksForVertex(n1)
            self.ds.init_pave_blocks_for_vertex(n1_key);
        }
        // L129-131: aMBlocks.Clear(); aMILI.Clear() — handled by Rust Drop
    }

    /// OCCT BOPAlgo_PaveFiller::MakeSDVertices (PaveFiller_1.cxx L136-233).
    /// Merges a connected group of vertices into a single SD vertex.
    /// If any member already has an SD partner (nSD), that SD vertex is
    /// updated in-place.  Otherwise a new vertex is appended to the DS.
    /// Every pair in the block gets AddShapeSD + a VV interference record
    /// pointing to the merged vertex.
    pub(super) fn make_sd_vertices_vv(&mut self, block: &[usize]) {
        // L136-138: return early if fewer than 2 vertices
        if block.len() < 2 {
            return;
        }
        // L141-161: 1. Collect vertices + track existing SD partner
        let mut n_sd: Option<usize> = None;
        let mut a_lv: Vec<usize> = Vec::with_capacity(block.len());
        for &n_x in block {
            // L145-158: check if vertex already has an SD partner
            if let Some(n_sd1) = self.ds.has_shape_sd(n_x) {
                if n_sd.is_none() {
                    // L148-153: keep the first SD vertex as the merge target
                    n_sd = Some(n_sd1);
                }
            }
            // L159-160: add vertex to aLV list
            a_lv.push(n_x);
        }
        // L162: MakeVertex(aLV, aVn) — compute centroid + bounding tolerance.
        // OCCT calls BRepLib::BoundingVertex to compute centroid and tolerance
        // large enough to enclose all input vertices.
        let centroid: DVec3 = a_lv
            .iter()
            .map(|&vi| self.ds.vertex_point(vi))
            .fold(DVec3::ZERO, |acc, p| acc + p)
            / a_lv.len() as f64;
        let bounding_tol: f64 = a_lv
            .iter()
            .map(|&vi| {
                (self.ds.vertex_point(vi) - centroid).length() + self.ds.vertex_tolerance(vi)
            })
            .fold(TOLERANCE_ABS, |acc, d| acc.max(d));
        // L163-179: 2. Determine nV — either update existing SD or append new
        let n_v: usize;
        if let Some(n_sd_idx) = n_sd {
            // L166-171: update existing SD vertex in-place (position + tolerance)
            self.ds.vertex_data_mut(n_sd_idx).point = centroid;
            self.ds.vertex_data_mut(n_sd_idx).tolerance = bounding_tol;
            n_v = n_sd_idx;
        } else {
            // L176-179: append new vertex to DS
            n_v = self.ds.vertices.len();
            self.ds.push_vertex(
                DSVertex {
                    point: centroid,
                    origin: None,
                    geom_tol: bounding_tol,
                    is_internal: true,
                    location: 0,
                },
                None,
            );
            // push_vertex now creates ShapeInfo (keeps shapes:shape_info 1:1).
            // Update the SD vertex box/gap with merged values.
            if n_v < self.ds.vertex_shape_idx.len() {
                let si = self.ds.vertex_shape_idx[n_v];
                if let Some(si_mut) = self.ds.shape_info.get_mut(si) {
                    si_mut.box_min = Some(centroid);
                    si_mut.box_max = Some(centroid);
                    si_mut.box_gap = bounding_tol + TOLERANCE_ABS;
                }
            }
        }
        // L181-184: update bounding box for the SD vertex (both nSD and new).
        // Use shapes-index si (not vertex index n_v) to stay 1:1 with shape_info.
        let si = self
            .ds
            .vertex_shape_idx
            .get(n_v)
            .copied()
            .unwrap_or(usize::MAX);
        if let Some(si_mut) = self.ds.shape_info.get_mut(si) {
            si_mut.box_min = Some(centroid);
            si_mut.box_max = Some(centroid);
            si_mut.box_gap = bounding_tol + TOLERANCE_ABS;
        }
        // L186-231: 3. Record SD mappings + VV interferences for every pair
        for i in 0..block.len() {
            let n1 = block[i];
            // L197: AddShapeSD(n1, nV)
            self.ds.add_shape_sd(n1, n_v);
            // L199-218: self-interfering shape warning
            let i_r1 = self.ds.rank(n1);
            // L221-228: VV interference for each pair (n1, n2)
            for j in (i + 1)..block.len() {
                let n2 = block[j];
                // OCCT L199-218: if same rank, add self-interfering shape warning
                if i_r1 == self.ds.rank(n2) {
                    self.my_report
                        .add_alert(crate::bopalgo::Alert::SelfInterferingShape(n1, n2));
                }
                // OCCT L223: if (myDS->AddInterf(n1, n2)) — fence check
                let key = if n1 < n2 { (n1, n2) } else { (n2, n1) };
                if self.ds.interf_tb.insert(key) {
                    // L225-227: aVV.SetIndices(n1, n2); aVV.SetIndexNew(nV)
                    self.ds.interf_vv.push(InterferenceVV {
                        v1: n1,
                        v2: n2,
                        merged_vertex: n_v,
                    });
                }
            }
        }
    }
    /// OCCT BOPAlgo_PaveFiller_3.cxx L145: PerformEE
    /// OCCT BOPAlgo_PaveFiller_3.cxx L145-585: PerformEE
    /// OCCT BOPAlgo_PaveFiller_3.cxx L580-640: CheckEdgeEdge
    pub(crate) fn intersect_ee(
        &mut self,
        e1: usize,
        e2: usize,
        range1: [f64; 2],
        range2: [f64; 2],
        modified: &mut std::collections::HashSet<usize>,
    ) {
        let edge1 = &self.ds.edges[e1];
        let edge2 = &self.ds.edges[e2];
        let tol = self.ee_tol(e1, e2);
        // Capture all edge data before mutable borrow
        let e1_curve = edge1.curve.clone();
        let e2_curve = edge2.curve.clone();
        drop(edge1);
        drop(edge2);
        // ??FillShrunkData computes shrunk ranges for each pave block.
        // If shrunk_range fails (edge too short), skip this pair entirely
        // (=OCCT BOPAlgo_PaveFiller_3: !aPB->IsSplittable() ??continue).
        let sr1 = match crate::inttools::curve_range::shrunk_range(&e1_curve, range1, tol, tol, tol)
        {
            Some(sr) => sr,
            None => return,
        };
        let sr2 = match crate::inttools::curve_range::shrunk_range(&e2_curve, range2, tol, tol, tol)
        {
            Some(sr) => sr,
            None => return,
        };

        // Compute intersections restricted to shrunk sub-ranges.
        let hits: Vec<(f64, f64, DVec3)> = match (&e1_curve, &e2_curve) {
            (Curve3::Line(l1), Curve3::Line(l2)) => intersect_line_line(l1, sr1, l2, sr2, tol)
                .into_iter()
                .map(|(t1, t2, p)| (t1, t2, p))
                .collect(),
            (Curve3::Line(l), Curve3::Circle(c)) => intersect_line_circle(l, c, tol)
                .into_iter()
                .filter(|(t_line, t_circle, _)| {
                    in_range(*t_line, sr1, tol) && in_range(*t_circle, sr2, tol)
                })
                .map(|(t_line, t_circle, p)| (t_line, t_circle, p))
                .collect(),
            (Curve3::Circle(c), Curve3::Line(l)) => intersect_line_circle(l, c, tol)
                .into_iter()
                .filter(|(t_line, t_circle, _)| {
                    in_range(*t_line, sr2, tol) && in_range(*t_circle, sr1, tol)
                })
                .map(|(t_line, t_circle, p)| (t_circle, t_line, p))
                .collect(),
            (Curve3::Circle(c1), Curve3::Circle(c2)) => intersect_circle_circle(c1, c2, tol)
                .into_iter()
                .filter_map(|p| {
                    let t1 = circle_param(p, c1);
                    let t2 = circle_param(p, c2);
                    if in_range(t1, sr1, tol) && in_range(t2, sr2, tol) {
                        Some((t1, t2, p))
                    } else {
                        None
                    }
                })
                .collect(),
            _ => {
                // Fallback: use EdgeEdgeIntersector for non-analytic curve pairs
                use crate::inttools::edge_edge::EdgeEdgeIntersector;
                let mut ee = EdgeEdgeIntersector::new();
                ee.set_edges(e1, sr1, e2, sr2, self.ds);
                ee.set_fuzzy_value(tol);
                ee.perform();
                ee.common_parts()
                    .iter()
                    .map(|cp| {
                        let t1 = cp.vertex_param1;
                        let t2 = cp.vertex_param2;
                        Some((t1, t2, cp.bounding_point1))
                    })
                    .flatten()
                    .collect()
            }
        };

        // OCCT L529-551: EDGE-type common parts — coincident edges create EE interferences
        if hits.is_empty() {
            let b_coincident = match (&e1_curve, &e2_curve) {
                (Curve3::Line(l1), Curve3::Line(l2)) => {
                    let cross = l1.direction.cross(l2.direction);
                    if cross.length() > tol {
                        false
                    } else {
                        (l2.origin - l1.origin).cross(l1.direction).length() <= tol
                    }
                }
                _ => false,
            };
            if b_coincident {
                let mid_t = (range1[0] + range1[1]) * 0.5;
                let mid_pt = e1_curve.point_at(mid_t);
                let new_v = self.ds.add_vertex(mid_pt);
                self.ds.interf_ee.push(InterferenceEE {
                    e1,
                    e2,
                    point: mid_pt,
                    param1: range1[0],
                    param2: range2[0],
                    new_vertex: new_v,
                    range1,
                    range2,
                });
                self.add_pave_to_edge(
                    e1,
                    Pave {
                        vertex_idx: new_v,
                        param: range1[0],
                    },
                );
                self.add_pave_to_edge(
                    e2,
                    Pave {
                        vertex_idx: new_v,
                        param: range2[0],
                    },
                );
                modified.insert(e1);
                modified.insert(e2);
                return;
            }
            return;
        }

        //  Process each intersection result (PaveFiller_3.cxx L682-750).
        // For each valid intersection, create a new vertex and record EE interference.
        for (t1, t2, point) in hits {
            if t1 < sr1[0] || t1 > sr1[1] || t2 < sr2[0] || t2 > sr2[1] {
                continue;
            }
            let new_v = self.ds.add_vertex(point);
            self.ds.interf_ee.push(InterferenceEE {
                e1,
                e2,
                point,
                param1: t1,
                param2: t2,
                new_vertex: new_v,
                range1: [t1, t1],
                range2: [t2, t2],
            });
            self.add_pave_to_edge(
                e1,
                Pave {
                    vertex_idx: new_v,
                    param: t1,
                },
            );
            self.add_pave_to_edge(
                e2,
                Pave {
                    vertex_idx: new_v,
                    param: t2,
                },
            );
            modified.insert(e1);
            modified.insert(e2);
        }
    }
    /// OCCT PaveFiller_3.cxx L580-640: CheckEdgeEdge
    pub(crate) fn check_edge_edge(&mut self, e1: usize, e2: usize) {
        let edge1 = &self.ds.edges[e1];
        let edge2 = &self.ds.edges[e2];
        let tol = self.ee_tol(e1, e2);

        let hits: Vec<(f64, f64, DVec3)> = match (&edge1.curve, &edge2.curve) {
            (Curve3::Line(l1), Curve3::Line(l2)) => {
                intersect_line_line(l1, edge1.t_range, l2, edge2.t_range, tol)
                    .into_iter()
                    .map(|(t1, t2, p)| (t1, t2, p))
                    .collect()
            }
            (Curve3::Line(l), Curve3::Circle(c)) => intersect_line_circle(l, c, tol)
                .into_iter()
                .filter(|(t_line, t_circle, _)| {
                    in_range(*t_line, edge1.t_range, tol) && in_range(*t_circle, edge2.t_range, tol)
                })
                .map(|(t_line, t_circle, p)| (t_line, t_circle, p))
                .collect(),
            (Curve3::Circle(c), Curve3::Line(l)) => intersect_line_circle(l, c, tol)
                .into_iter()
                .filter(|(t_line, t_circle, _)| {
                    in_range(*t_line, edge2.t_range, tol) && in_range(*t_circle, edge1.t_range, tol)
                })
                .map(|(t_line, t_circle, p)| (t_circle, t_line, p))
                .collect(),
            (Curve3::Circle(c1), Curve3::Circle(c2)) => intersect_circle_circle(c1, c2, tol)
                .into_iter()
                .filter_map(|p| {
                    let t1 = circle_param(p, c1);
                    let t2 = circle_param(p, c2);
                    if in_range(t1, edge1.t_range, tol) && in_range(t2, edge2.t_range, tol) {
                        Some((t1, t2, p))
                    } else {
                        None
                    }
                })
                .collect(),
            _ => vec![],
        };

        for (t1, t2, point) in hits {
            let new_v = self.ds.add_vertex(point);
            self.ds.interf_ee.push(InterferenceEE {
                e1,
                e2,
                point,
                param1: t1,
                param2: t2,
                new_vertex: new_v,
                range1: [t1, t1],
                range2: [t2, t2],
            });
            self.add_pave_to_edge(
                e1,
                Pave {
                    vertex_idx: new_v,
                    param: t1,
                },
            );
            self.add_pave_to_edge(
                e2,
                Pave {
                    vertex_idx: new_v,
                    param: t2,
                },
            );
        }
    }
    /// OCCT PaveFiller L575-590 (L692-723 in source): TreatNewVertices
    pub(crate) fn treat_new_vertices(
        &mut self,
        a_mvcpb: &std::collections::HashMap<usize, crate::bopds::pave::CoupleOfPaveBlocks>,
    ) -> Vec<usize> {
        let a_nb_v = a_mvcpb.len();
        if a_nb_v == 0 {
            return vec![];
        }

        // = =  Phase 1: Collect new vertices from a_mvcpb (OCCT L696-702) = = = = = = = = =
        #[derive(Clone, Copy)]
        struct NewVertInfo {
            idx: usize,
            pos: DVec3,
            tol: f64,
        }
        let mut new_verts: Vec<NewVertInfo> = Vec::with_capacity(a_nb_v);
        for (&vi, cpb) in a_mvcpb {
            let pos = self.ds.vertex_point(vi);
            let tol = cpb.tolerance.max(self.ds.fuzzy_tol);
            new_verts.push(NewVertInfo { idx: vi, pos, tol });
        }
        if new_verts.len() < 2 {
            return new_verts.iter().map(|v| v.idx).collect();
        }

        // = =  Phase 2: IntersectVertices — group nearby vertices into chains = = =
        // OCCT L704-708: BOPAlgo_Tools::IntersectVertices(aVerts, myFuzzyValue, aChains)
        // rcad: BVH-based grouping (equivalent to IntersectVertices)
        let gap = self.ds.fuzzy_tol / 2.0;
        use crate::bvh::{Aabb, BoxTree};
        let nv = new_verts.len();
        let mut bvh_indices: Vec<usize> = Vec::with_capacity(nv);
        let mut bvh_aabbs: Vec<Aabb> = Vec::with_capacity(nv);
        for i in 0..nv {
            let v = &new_verts[i];
            bvh_indices.push(i);
            let half = v.tol + gap;
            bvh_aabbs.push(Aabb {
                min: v.pos - DVec3::splat(half),
                max: v.pos + DVec3::splat(half),
                gap: 0.0,
            });
        }
        let bvh = BoxTree::build(bvh_indices, bvh_aabbs);
        let pairs = BoxTree::candidate_pairs(&bvh, &bvh);
        let mut a_mili: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for &(ia, ib) in &pairs {
            let tol_a = new_verts[ia].tol;
            let tol_b = new_verts[ib].tol;
            let merge_tol = tol_a + tol_b + self.ds.fuzzy_tol;
            if (new_verts[ia].pos - new_verts[ib].pos).length() <= merge_tol {
                fill_map(&mut a_mili, ia, ib);
            }
        }
        let a_blocks = make_blocks(&a_mili);
        let mut groups: Vec<Vec<usize>> = a_blocks
            .iter()
            .map(|block| block.iter().map(|&i| new_verts[i].idx).collect())
            .collect();
        {
            let mut taken: HashSet<usize> = HashSet::new();
            for group in &groups {
                for &vi in group {
                    taken.insert(vi);
                }
            }
            for i in 0..nv {
                if !taken.contains(&new_verts[i].idx) {
                    groups.push(vec![new_verts[i].idx]);
                }
            }
        }

        // = =  Phase 3: MakeVertex for each chain (BOPTools_AlgoTools::MakeVertex) = = =
        // OCCT L714-721: for each chain, create a new TopoDS_Vertex via MakeVertex
        // rcad: groups with >1 member get a new fused vertex; singletons reuse the old index.
        let mut survivors: Vec<usize> = Vec::with_capacity(groups.len());
        for members in &groups {
            if members.len() < 2 {
                survivors.push(members[0]);
                continue;
            }

            // OCCT: BRepLib::BoundingVertex computes center + tolerance.
            let (center, tol) = bounding_vertex(members, &self.ds);

            // OCCT BRep_Builder::MakeVertex: create new vertex at center.
            let new_vi = self.ds.add_vertex_no_dedup(center);
            self.ds.vertex_data_mut(new_vi).tolerance = tol;
            // OCCT: myIncreasedSS.Add(nV)  ?mark tolerance as increased.
            self.my_increased_ss.insert(new_vi);

            // Update interferences and paves to point to the new vertex
            for &old_vi in members {
                if old_vi == new_vi {
                    continue;
                }
                for edge in &mut self.ds.edges {
                    for pave in &mut edge.paves {
                        if pave.vertex_idx == old_vi {
                            pave.vertex_idx = new_vi;
                        }
                    }
                }
                for inf in &mut self.ds.interf_ee {
                    if inf.new_vertex == old_vi {
                        inf.new_vertex = new_vi;
                    }
                }
                for inf in &mut self.ds.interf_ef {
                    if inf.new_vertex == old_vi {
                        inf.new_vertex = new_vi;
                    }
                }
                for face in &mut self.ds.faces {
                    if face.face_info.vertices_on.remove(&old_vi) {
                        face.face_info.vertices_on.insert(new_vi);
                    }
                    if face.face_info.vertices_in.remove(&old_vi) {
                        face.face_info.vertices_in.insert(new_vi);
                    }
                }
            }
            survivors.push(new_vi);
        }
        survivors
    }
    // OCCT BOPAlgo_PaveFiller_3.cxx L594-687: PerformNewVertices
    pub(crate) fn perform_new_vertices(
        &mut self,
        a_mvcpb: std::collections::HashMap<usize, crate::bopds::pave::CoupleOfPaveBlocks>,
        b_is_ee_intersection: bool,
    ) {
        let a_nb_v = a_mvcpb.len();
        if a_nb_v == 0 {
            return;
        }

        // Step 1: Fuse the new vertices via TreatNewVertices (OCCT L609-612)
        let _survivors = self.treat_new_vertices(&a_mvcpb);

        // Steps 2-4: Build edge→vertices map and call IntersectVE (OCCT L655-687)
        // OCCT builds aMPBLI (PB → vertices) then calls IntersectVE(aMPBLI, ..., false).
        // rcad: build edge→vertices map and call self.intersect_ve().
        //
        // For EF intersections where the vertex coincides with a PB endpoint,
        // OCCT's TreatNewVertices creates a new vertex index so IntersectVE
        // doesn't skip it. In rcad, if the vertex wasn't fused (singleton group),
        // it keeps its old index → IntersectVE skips it. Handle this case manually.
        let mut edge_verts: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        let mut ef_endpoint_cases: Vec<(usize, usize)> = Vec::new(); // (edge, vertex)

        for (_v_key, cpb) in &a_mvcpb {
            let (n_v, n_e) = if b_is_ee_intersection {
                let inf = &self.ds.interf_ee[cpb.interf_idx];
                (inf.new_vertex, inf.e1)
            } else {
                let inf = &self.ds.interf_ef[cpb.interf_idx];
                if inf.new_vertex == usize::MAX {
                    continue;
                }
                (inf.new_vertex, inf.edge)
            };
            if n_v == usize::MAX || n_e >= self.ds.edges.len() {
                continue;
            }

            // Check if vertex is a PB endpoint
            let is_endpoint = {
                let pbs = self.ds.edge_pave_blocks(n_e);
                if pbs.is_empty() {
                    continue;
                }
                let mut ep_set: HashSet<usize> = HashSet::new();
                for spb in pbs {
                    let pb = spb.0.read().unwrap();
                    ep_set.insert(pb.pave1.vertex_idx);
                    ep_set.insert(pb.pave2.vertex_idx);
                }
                ep_set.contains(&n_v)
            };

            if is_endpoint && !b_is_ee_intersection {
                // EF vertex at PB endpoint: handle manually
                ef_endpoint_cases.push((n_e, n_v));
            } else if !is_endpoint {
                edge_verts.entry(n_e).or_default().push(n_v);
            }
            // EE endpoint case: skip (IntersectVE skips them too)
        }

        // Step 4a: Split PBs at new vertices via IntersectVE (OCCT L687)
        if !edge_verts.is_empty() {
            self.intersect_ve(&edge_verts, false);
        }

        // Step 4b: Handle EF endpoint cases (rcad-specific, preserves existing behavior)
        if !ef_endpoint_cases.is_empty() {
            use crate::bopds::pave::{Pave, PaveBlock, SharedPB};
            for (n_e, n_v) in ef_endpoint_cases {
                if self.ds.edges[n_e].pave_blocks.len() != 1 {
                    continue;
                }
                let mut a_ve = VertexEdgeSolver::new();
                a_ve.set_data(n_v, n_e);
                a_ve.perform(&mut self.context, &self.ds, self.ds.fuzzy_tol);
                if !a_ve.is_done() {
                    continue;
                }
                let a_t = a_ve.param;
                let (pb_sv, pb_ev, pb_t1, pb_t2, orig_e) = {
                    let pb = self.ds.edges[n_e].pave_blocks[0].0.read().unwrap();
                    (
                        pb.pave1.vertex_idx,
                        pb.pave2.vertex_idx,
                        pb.pave1.param,
                        pb.pave2.param,
                        pb.original_edge,
                    )
                };
                let a_t_split = a_t.clamp(pb_t1 + 1e-12, pb_t2 - 1e-12);
                let pv1 = Pave {
                    vertex_idx: pb_sv,
                    param: pb_t1,
                };
                let pv_new = Pave {
                    vertex_idx: n_v,
                    param: a_t_split,
                };
                let pv2 = Pave {
                    vertex_idx: pb_ev,
                    param: pb_t2,
                };
                let mut pb1 = PaveBlock::new(orig_e, pv1, pv_new);
                let mut pb2 = PaveBlock::new(orig_e, pv_new, pv2);
                pb1.set_shrunk_data(pb_t1, a_t_split, true);
                pb2.set_shrunk_data(a_t_split, pb_t2, true);
                self.ds.edges[n_e].pave_blocks.clear();
                self.ds.edges[n_e].pave_blocks.push(SharedPB::new(pb1));
                self.ds.edges[n_e].pave_blocks.push(SharedPB::new(pb2));
            }
        }
    }
    // OCCT BOPAlgo_PaveFiller.cxx L383-448: RepeatIntersection
    pub(crate) fn repeat_intersection(&mut self) {
        // OCCT L385-386: NCollection_Map<int> anExtraInterfMap;
        let mut an_extra_interf_map: HashSet<usize> = HashSet::new();
        // OCCT L387: const int aNbS = myDS->NbSourceShapes();
        let a_nb_s = self.ds.nb_source_shapes;
        // OCCT L388: Message_ProgressScope aPS(theRange, "Repeat intersection", 3);
        // OCCT L389-414: for (int i = 0; i < aNbS; ++i)
        for i in 0..a_nb_s {
            // OCCT L391-395: if ShapeType != VERTEX, continue
            if self.ds.shape_type_of(i) != ShapeType::Vertex {
                continue;
            }
            // OCCT L397-401: if (myIncreasedSS.Contains(i)) { anExtraInterfMap.Add(i); continue; }
            if self.my_increased_ss.contains(&i) {
                an_extra_interf_map.insert(i);
                continue;
            }
            // OCCT L404-408: int nVSD; if (!myDS->HasShapeSD(i, nVSD)) { continue; }
            if let Some(n_vsd) = self.ds.has_shape_sd(i) {
                // OCCT L410-413: if (myIncreasedSS.Contains(nVSD)) { anExtraInterfMap.Add(i); }
                if self.my_increased_ss.contains(&n_vsd) {
                    an_extra_interf_map.insert(i);
                }
            } // else: OCCT L405-407 continue (handled by if-let None)
        }
        // OCCT L416-419: if (anExtraInterfMap.IsEmpty()) return;
        if an_extra_interf_map.is_empty() {
            return;
        }

        // OCCT L422: myIterator->IntersectExt(anExtraInterfMap);
        self.my_iterator.intersect_ext(&an_extra_interf_map);

        // OCCT L426-430: PerformVV(aPS.Next());
        self.perform_vv();
        if self.my_report.has_errors() {
            return;
        }
        // OCCT L431: UpdatePaveBlocksWithSDVertices();
        self.ds.update_pave_blocks_with_sd_vertices();

        // OCCT L433-438: PerformVE(aPS.Next());
        self.perform_ve();
        if self.my_report.has_errors() {
            return;
        }
        // OCCT L438: UpdatePaveBlocksWithSDVertices();
        self.ds.update_pave_blocks_with_sd_vertices();

        // OCCT L440-444: PerformVF(aPS.Next());
        self.perform_vf();
        if self.my_report.has_errors() {
            return;
        }

        // OCCT L446-447: UpdatePaveBlocksWithSDVertices(); UpdateInterfsWithSDVertices();
        self.ds.update_pave_blocks_with_sd_vertices();
        self.update_interfs_with_sd_vertices();
    }
    /// OCCT PaveFiller_4.cxx: PerformVF
    /// Shrinks the range by the face tolerance converted to parametric space.
    pub(crate) fn correct_range_for_face(
        edge_curve: &Curve3,
        etf: f64,
        range: [f64; 2],
    ) -> [f64; 2] {
        const DT: f64 = 1e-12;
        let a_tf = range[0];
        let a_tl = range[1];
        let mut a_new_first = a_tf;
        let mut a_new_last = a_tl;
        // OCCT L387-433: for (i = 0; i < 2; ++i)
        for i in 0..2 {
            let t = if i == 0 { a_tf } else { a_tl };
            // OCCT L389: aRes = aTolF; then convert to parametric space
            let a_res = match edge_curve {
                // OCCT L416-417: analytic ??aBC.Resolution(aRes)
                Curve3::Line(l) => {
                    let dir_len = l.direction.length();
                    if dir_len > 1e-12 {
                        etf / dir_len
                    } else {
                        etf
                    }
                }
                Curve3::Circle(c) => etf / c.radius.max(TOLERANCE_ABS),
                Curve3::Ellipse(e) => etf / e.major_radius.max(TOLERANCE_ABS),
                // OCCT L391-413: BSpline/Bezier ??aRes / |derivative|
                _ => {
                    let dt = 1e-7;
                    let p1 = edge_curve.point_at(t - dt);
                    let p2 = edge_curve.point_at(t + dt);
                    let dm = (p2 - p1).length() / (2.0 * dt);
                    if dm > 1e-12 {
                        etf / dm
                    } else {
                        etf
                    }
                }
            };
            // OCCT L420-427: shrink endpoint
            if i == 0 {
                a_new_first = a_tf + a_res;
            } else {
                a_new_last = a_tl - a_res;
            }
            // OCCT L429-432: if too small, restore original
            if (a_new_last - a_new_first) < DT {
                return range;
            }
        }
        [a_new_first, a_new_last]
    }
    /// OCCT IntTools_FClass2d: point-in-face check.
    /// Uses fclass2d for all surface types to properly check
    /// against the face's trimming wires (not just parametric bounds).
    pub(crate) fn is_point_in_face(&mut self, point: DVec3, face_idx: usize, tol: f64) -> bool {
        let face = &self.ds.faces[face_idx];
        match &face.surface {
            Surface3::Plane(plane) => {
                let verts = self.ds.face_boundary_points(face_idx);
                inttools::edge_face::point_in_planar_face_with_tol(point, plane, &verts, tol)
            }
            _ => {
                // OCCT IntTools_FClass2d: project to UV, check within face wires.
                self.context
                    .is_point_in_face_3d(&self.ds, face_idx, point, tol)
            }
        }
    }
    /// OCCT PaveFiller_5.cxx L340-480: IntersectEdgeFace
    pub(crate) fn intersect_ef(&mut self, edge_idx: usize, face_idx: usize, pb_range: &[f64; 2]) {
        let edge_curve = self.ds.edges[edge_idx].curve.clone();
        let edge_t_range = self.ds.edge_range(edge_idx);

        // Use PaveBlock range to constrain intersection interval (OCCT L262: SetRange(aPBRange))
        let ef_range = [
            pb_range[0].max(edge_t_range[0]),
            pb_range[1].min(edge_t_range[1]),
        ];
        let etf = self.ef_tol(edge_idx, face_idx);
        if ef_range[1] - ef_range[0] <= etf {
            return;
        }
        let ef_range = Self::correct_range_for_face(&edge_curve, etf, ef_range);
        if ef_range[1] - ef_range[0] <= etf {
            return;
        }
        let face_surface = self.ds.faces[face_idx].surface.clone();

        // Dispatch based on curve type  ?surface type
        let hits: Vec<(DVec3, f64)> = match (&edge_curve, &face_surface) {
            (Curve3::Line(line), Surface3::Plane(plane)) => {
                inttools::edge_face::intersect_line_plane_with_tol(line, ef_range, plane, etf)
                    .into_iter()
                    .map(|h| (h.point, h.edge_param))
                    .collect()
            }
            (Curve3::Line(line), Surface3::Cylinder(cyl)) => {
                inttools::curve_surface::intersect_line_cylinder_with_tol(line, ef_range, cyl, etf)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            (Curve3::Line(line), Surface3::Sphere(sph)) => {
                inttools::curve_surface::intersect_line_sphere_with_tol(line, ef_range, sph, etf)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            (Curve3::Line(line), Surface3::Cone(cone)) => {
                inttools::curve_surface::intersect_line_cone_with_tol(line, ef_range, cone, etf)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            (Curve3::Circle(circle), Surface3::Plane(plane)) => {
                // Use edge start vertex as reference direction for  ?0
                let sv = self.ds.edge_start_vertex_ds(edge_idx);
                let ref_dir = (self.ds.vertex_point(sv) - circle.center).normalize();
                inttools::curve_surface::intersect_circle_plane_with_ref(
                    circle,
                    ef_range,
                    plane,
                    etf,
                    Some(ref_dir),
                )
                .into_iter()
                .map(|h| (h.point, h.curve_param))
                .collect()
            }
            (Curve3::Circle(circle), Surface3::Cylinder(cyl)) => {
                inttools::curve_surface::intersect_circle_cylinder_with_tol(
                    circle, ef_range, cyl, etf,
                )
                .into_iter()
                .map(|h| (h.point, h.curve_param))
                .collect()
            }
            (Curve3::Circle(circle), Surface3::Sphere(sph)) => {
                inttools::curve_surface::intersect_circle_sphere_with_tol(
                    circle, ef_range, sph, etf,
                )
                .into_iter()
                .map(|h| (h.point, h.curve_param))
                .collect()
            }
            (Curve3::Circle(circle), Surface3::Cone(cone)) => {
                inttools::curve_surface::intersect_circle_cone_with_tol(circle, ef_range, cone, etf)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            (Curve3::Ellipse(ellipse), Surface3::Plane(plane)) => {
                //  IntAna_IntConicQuad Ellipse  ?Plane
                inttools::ellipse_intersection::intersect_ellipse_plane_with_tol(
                    ellipse, ef_range, plane, etf,
                )
                .into_iter()
                .map(|h| (h.point, h.ellipse_param))
                .collect()
            }
            (Curve3::Ellipse(ellipse), Surface3::Cylinder(cyl)) => {
                //  ?Partially aligned: numeric fallback, same as OCCT for rare cases
                inttools::ellipse_intersection::intersect_ellipse_cylinder_with_tol(
                    ellipse, ef_range, cyl, etf,
                )
                .into_iter()
                .map(|h| (h.point, h.ellipse_param))
                .collect()
            }
            (Curve3::Ellipse(ellipse), Surface3::Sphere(sph)) => {
                inttools::ellipse_intersection::intersect_ellipse_sphere_with_tol(
                    ellipse, ef_range, sph, etf,
                )
                .into_iter()
                .map(|h| (h.point, h.ellipse_param))
                .collect()
            }
            (Curve3::Ellipse(ellipse), Surface3::Cone(cone)) => {
                inttools::ellipse_intersection::intersect_ellipse_cone_with_tol(
                    ellipse, ef_range, cone, etf,
                )
                .into_iter()
                .map(|h| (h.point, h.ellipse_param))
                .collect()
            }
            (Curve3::Parabola(parabola), Surface3::Plane(plane)) => {
                //  IntAna_IntConicQuad Parabola  ?Plane
                inttools::parabola_intersection::intersect_parabola_plane_with_tol(
                    parabola, ef_range, plane, etf,
                )
                .into_iter()
                .map(|h| (h.point, h.parabola_param))
                .collect()
            }
            (Curve3::Parabola(parabola), Surface3::Cylinder(cyl)) => {
                //  ?Partially aligned: numeric fallback
                inttools::parabola_intersection::intersect_parabola_cylinder_with_tol(
                    parabola, ef_range, cyl, etf,
                )
                .into_iter()
                .map(|h| (h.point, h.parabola_param))
                .collect()
            }
            (Curve3::Parabola(parabola), Surface3::Sphere(sph)) => {
                inttools::parabola_intersection::intersect_parabola_sphere_with_tol(
                    parabola, ef_range, sph, etf,
                )
                .into_iter()
                .map(|h| (h.point, h.parabola_param))
                .collect()
            }
            (Curve3::Parabola(parabola), Surface3::Cone(cone)) => {
                inttools::parabola_intersection::intersect_parabola_cone_with_tol(
                    parabola, ef_range, cone, etf,
                )
                .into_iter()
                .map(|h| (h.point, h.parabola_param))
                .collect()
            }
            (Curve3::Hyperbola(hyperbola), Surface3::Plane(plane)) => {
                //  IntAna_IntConicQuad Hyperbola  ?Plane
                inttools::hyperbola_intersection::intersect_hyperbola_plane_with_tol(
                    hyperbola, ef_range, plane, etf,
                )
                .into_iter()
                .map(|h| (h.point, h.hyperbola_param))
                .collect()
            }
            (Curve3::Hyperbola(hyperbola), Surface3::Cylinder(cyl)) => {
                //  ?Partially aligned: numeric fallback
                inttools::hyperbola_intersection::intersect_hyperbola_cylinder_with_tol(
                    hyperbola, ef_range, cyl, etf,
                )
                .into_iter()
                .map(|h| (h.point, h.hyperbola_param))
                .collect()
            }
            (Curve3::Hyperbola(hyperbola), Surface3::Sphere(sph)) => {
                inttools::hyperbola_intersection::intersect_hyperbola_sphere_with_tol(
                    hyperbola, ef_range, sph, etf,
                )
                .into_iter()
                .map(|h| (h.point, h.hyperbola_param))
                .collect()
            }
            (Curve3::Hyperbola(hyperbola), Surface3::Cone(cone)) => {
                inttools::hyperbola_intersection::intersect_hyperbola_cone_with_tol(
                    hyperbola, ef_range, cone, etf,
                )
                .into_iter()
                .map(|h| (h.point, h.hyperbola_param))
                .collect()
            }
            _ => {
                // Use BeanFaceIntersector for generic curve-surface pairs.
                use crate::inttools::bean_face_intersector::BeanFaceIntersector;
                use rcad_kernel::geom::SurfaceEval;
                let mut bfi = BeanFaceIntersector::from_curve_surface(
                    edge_curve.clone(),
                    face_surface.clone(),
                );
                bfi.set_bean_parameters(ef_range[0], ef_range[1]);
                let [u_min, u_max, v_min, v_max] = face_surface.default_domain();
                bfi.set_surface_parameters(u_min, u_max, v_min, v_max);
                bfi.init_curve_surface(edge_curve.clone(), etf, face_surface.clone(), etf);
                bfi.set_bean_parameters(ef_range[0], ef_range[1]);
                bfi.perform();
                let mut pts = Vec::new();
                if bfi.is_done() {
                    for r in bfi.result() {
                        let t = (r.first() + r.last()) * 0.5;
                        let p = edge_curve.point_at(t);
                        pts.push((p, t));
                    }
                }
                pts
            }
        };

        for (point, edge_param) in hits {
            //  IsPointInFace check for ALL surface types (PaveFiller_5.cxx L523)
            let in_face = self.is_point_in_face(point, face_idx, etf);
            if !in_face {
                let near_face_vert = self
                    .ds
                    .face_boundary_points(face_idx)
                    .iter()
                    .any(|&vp| (vp - point).length() <= etf * 2.0);
                if !near_face_vert {
                    continue;
                }
            }
            // OCCT IntTools_EdgeFace creates a new vertex for each hit, even when
            // the hit coincides with an existing edge endpoint.  SD vertex merging
            // handles near-coincident vertices later (MakeSDVerticesFF in PostTreat).
            // rcad: do NOT skip endpoint-coincident hits  ?they are needed for
            // PutPaveOnCurve to split intersection curve pave blocks.
            let new_v = self.ds.add_vertex(point);
            // Register vertices_on for the new vertex if it's near the edge boundary
            let sv = self.ds.edge_start_vertex_ds(edge_idx);
            let ev = self.ds.edge_end_vertex_ds(edge_idx);
            let tol = etf
                .max(self.ds.vertex_tolerance(sv))
                .max(self.ds.vertex_tolerance(ev));
            if (point - self.ds.vertex_point(sv)).length() <= tol
                || (point - self.ds.vertex_point(ev)).length() <= tol
            {
                self.ds.faces[face_idx].face_info.vertices_on.insert(new_v);
            }
            //  Create EF interference for EVERY hit, even at edge endpoints.
            // OCCT IntTools_EdgeFace creates a new vertex for each hit (no dedup).
            // rcad: remove the vertices_on skip check  ?always push interference.
            self.ds.interf_ef.push(InterferenceEF {
                edge: edge_idx,
                face: face_idx,
                point,
                edge_param,
                new_vertex: new_v,
            });
            // Only mark vertices_on if actually inserted (avoid duplicate insert msg)
            if !self.ds.faces[face_idx]
                .face_info
                .vertices_on
                .contains(&new_v)
            {
                self.ds.faces[face_idx].face_info.vertices_on.insert(new_v);
            }
            self.add_pave_to_edge(
                edge_idx,
                Pave {
                    vertex_idx: new_v,
                    param: edge_param,
                },
            );
        }
    }
    /// OCCT HasInterf: skip already-processed pairs
    pub fn skip_redundant_interferences(&self) -> std::collections::HashSet<(usize, usize, u8)> {
        let mut skip_set = std::collections::HashSet::new();

        if !self.use_glue() {
            return skip_set;
        }

        // Skip V-V for shared vertices
        for &(va, vb) in &self.ds.shared_topology.shared_vertices {
            skip_set.insert((va, vb, 0)); // 0 = V-V
        }

        // Skip E-E for shared edges
        for &(ea, eb) in &self.ds.shared_topology.shared_edges {
            skip_set.insert((ea, eb, 2)); // 2 = E-E
        }

        // Skip F-F for fully glued faces
        for &(fa, fb) in &self.ds.shared_topology.fully_glued_faces {
            skip_set.insert((fa, fb, 5)); // 5 = F-F
        }

        skip_set
    }

    /// OCCT myDS->HasInterf: check existing EE interference
    pub(crate) fn has_ee_interf(&self, e1: usize, e2: usize) -> bool {
        self.ds
            .interf_ee
            .iter()
            .any(|inf| (inf.e1 == e1 && inf.e2 == e2) || (inf.e1 == e2 && inf.e2 == e1))
    }
}
