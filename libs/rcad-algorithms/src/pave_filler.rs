use glam::DVec3;
use rcad_kernel::geom::*;

use crate::bopds::ds::*;
use crate::bopds::pave::*;
use crate::bvh::Bvh;
use crate::inttools;
use crate::tolerance::*;

/// Minimum total face count before BVH acceleration is used.
/// Below this threshold, brute-force O(n²) is faster due to BVH build overhead.
const BVH_THRESHOLD: usize = 20;

/// PaveFiller executes the six intersection passes (OCCT: BOPAlgo_PaveFiller).
pub struct PaveFiller<'a> {
    pub ds: &'a mut DS,
    bvh_a: Option<&'a Bvh>,
    bvh_b: Option<&'a Bvh>,
    use_glue: bool,
    glue_tolerance: f64,
}

impl<'a> PaveFiller<'a> {
    pub fn new(ds: &'a mut DS) -> Self {
        Self {
            ds,
            bvh_a: None,
            bvh_b: None,
            use_glue: false,
            glue_tolerance: TOLERANCE_ABS,
        }
    }

    /// Create a PaveFiller with optional BVH acceleration for face-face intersection.
    ///
    /// `bvh_a` and `bvh_b` must be built from the same BReps that were used to
    /// construct the DS. Face indices in the BVHs map directly to DS face indices
    /// (A faces come first, then B faces).
    pub fn with_bvh(ds: &'a mut DS, bvh_a: &'a Bvh, bvh_b: &'a Bvh) -> Self {
        let total_faces = ds.faces.len();
        let use_bvh = total_faces >= BVH_THRESHOLD;
        Self {
            ds,
            bvh_a: if use_bvh { Some(bvh_a) } else { None },
            bvh_b: if use_bvh { Some(bvh_b) } else { None },
            use_glue: false,
            glue_tolerance: TOLERANCE_ABS,
        }
    }

    /// Configure shared-face glue detection for the face-face pass.
    pub fn configure_glue(&mut self, enable: bool, tolerance: f64) {
        self.use_glue = enable;
        self.glue_tolerance = tolerance.max(TOLERANCE_ABS);
    }

    /// Configure glue with adaptive tolerance based on input geometry.
    ///
    /// This function analyzes the input shapes and computes an appropriate
    /// glue tolerance based on geometry characteristics such as:
    /// - Minimum feature size
    /// - Face area distribution
    /// - Edge length distribution
    ///
    /// # Arguments
    /// * `enable` - Whether to enable glue detection.
    /// * `base_tolerance` - Base tolerance to start with.
    /// * `adaptive` - Whether to use adaptive tolerance adjustment.
    ///
    /// # Returns
    /// The computed adaptive glue tolerance.
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

    /// Compute adaptive glue tolerance based on geometry characteristics.
    fn compute_adaptive_glue_tolerance(&self, base_tolerance: f64) -> f64 {
        let mut min_feature_size = f64::INFINITY;
        let mut min_edge_length = f64::INFINITY;
        let mut min_face_area = f64::INFINITY;

        // Analyze edge lengths
        for edge in &self.ds.edges {
            let p1 = edge.curve.point_at(edge.t_range[0]);
            let p2 = edge.curve.point_at(edge.t_range[1]);
            let length = (p2 - p1).length();
            if length > 1e-10 {
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
                if diag > 1e-10 {
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

    /// Effective tolerance for coincidence tests in all passes.
    ///
    /// Returns the DS `fuzzy_tol` (already clamped to ≥ `TOLERANCE_ABS`).
    #[inline]
    fn tol(&self) -> f64 {
        self.ds.fuzzy_tol
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
        let snap_tol = (local_scale * 4.0).max(1e-4);

        if (start - chain[0]).length() <= snap_tol {
            chain[0] = start;
        }
        if (end - chain[chain.len() - 1]).length() <= snap_tol {
            let last = chain.len() - 1;
            chain[last] = end;
        }
    }

    /// Execute all intersection passes.
    pub fn perform(&mut self) {
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

        if !skip_ve {
            self.perform_ve();
        }

        if !skip_ee {
            self.perform_ee();
        }

        if !skip_vf {
            self.perform_vf();
        }

        if !skip_ef {
            self.perform_ef();
        }

        if !skip_ff {
            self.perform_ff();
        }

        self.build_split_edges();
    }

    /// Determine if Vertex-Edge pass can be skipped.
    ///
    /// Returns true when all shared vertices are connected to shared edges,
    /// meaning no additional V-E intersections are needed.
    fn should_skip_ve_pass(&self) -> bool {
        if !self.use_glue {
            return false;
        }

        // If all vertices are shared, skip V-E pass
        let shared_verts = &self.ds.shared_topology.shared_vertices;
        if shared_verts.is_empty() {
            return false;
        }

        // Check if all vertices from shape A have matches in shape B
        let a_verts: std::collections::HashSet<usize> = self.ds.vertices
            .iter()
            .enumerate()
            .filter(|(_, v)| v.origin == Some(ShapeOrigin::ShapeA))
            .map(|(i, _)| i)
            .collect();

        let matched_a: std::collections::HashSet<usize> = shared_verts
            .iter()
            .map(|(a, _)| *a)
            .collect();

        a_verts == matched_a && !a_verts.is_empty()
    }

    /// Determine if Edge-Edge pass can be skipped.
    ///
    /// Returns true when all shared edges are detected, meaning no additional
    /// E-E intersections are needed.
    fn should_skip_ee_pass(&self) -> bool {
        if !self.use_glue {
            return false;
        }

        let shared_edges = &self.ds.shared_topology.shared_edges;
        if shared_edges.is_empty() {
            return false;
        }

        // Check if all edges from shape A have matches in shape B
        let a_edges: std::collections::HashSet<usize> = self.ds.edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.origin == ShapeOrigin::ShapeA)
            .map(|(i, _)| i)
            .collect();

        let matched_a: std::collections::HashSet<usize> = shared_edges
            .iter()
            .map(|(a, _)| *a)
            .collect();

        a_edges == matched_a && !a_edges.is_empty()
    }

    /// Determine if Vertex-Face pass can be skipped.
    fn should_skip_vf_pass(&self) -> bool {
        if !self.use_glue {
            return false;
        }

        // If all faces are fully glued, skip V-F pass
        self.ds.shared_topology.fully_glued_faces.len() > 0
            && self.ds.shared_topology.fully_glued_faces.len()
                == self.ds.a_face_count * (self.ds.faces.len() - self.ds.a_face_count)
    }

    /// Determine if Edge-Face pass can be skipped.
    fn should_skip_ef_pass(&self) -> bool {
        if !self.use_glue {
            return false;
        }

        // If all faces are fully glued, skip E-F pass
        self.ds.shared_topology.fully_glued_faces.len() > 0
            && self.ds.shared_topology.fully_glued_faces.len()
                == self.ds.a_face_count * (self.ds.faces.len() - self.ds.a_face_count)
    }

    /// Determine if Face-Face pass can be skipped.
    ///
    /// Returns true when all faces have been detected as fully glued,
    /// meaning no F-F intersections are needed.
    fn should_skip_ff_pass(&self) -> bool {
        if !self.use_glue {
            return false;
        }

        // If all faces are fully glued, skip F-F pass
        let total_face_pairs = self.ds.a_face_count * (self.ds.faces.len() - self.ds.a_face_count);
        self.ds.shared_topology.fully_glued_faces.len() == total_face_pairs && total_face_pairs > 0
    }

    /// Skip redundant interferences based on pre-detected shared topology.
    ///
    /// This function identifies interference computations that can be skipped
    /// because the involved sub-shapes are already known to share topology.
    ///
    /// # Returns
    /// A set of (subshape_a, subshape_b, interference_type) pairs that can be skipped.
    pub fn skip_redundant_interferences(&self) -> std::collections::HashSet<(usize, usize, u8)> {
        let mut skip_set = std::collections::HashSet::new();

        if !self.use_glue {
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

    // ─── Pass 1: Vertex-Vertex ─────────────────────────────────────────

    fn perform_vv(&mut self) {
        // Use pre-detected shared vertices if glue is enabled
        if self.use_glue && !self.ds.shared_topology.shared_vertices.is_empty() {
            for &(vi_a, vi_b) in &self.ds.shared_topology.shared_vertices {
                self.ds.interferences.push(Interference::VertexVertex {
                    v1: vi_a,
                    v2: vi_b,
                    merged_vertex: vi_a,
                });
            }
            return;
        }

        // Fallback: brute-force search
        let a_verts: Vec<usize> = self
            .ds
            .vertices
            .iter()
            .enumerate()
            .filter(|(_, v)| v.origin == Some(ShapeOrigin::ShapeA))
            .map(|(i, _)| i)
            .collect();
        let b_verts: Vec<usize> = self
            .ds
            .vertices
            .iter()
            .enumerate()
            .filter(|(_, v)| v.origin == Some(ShapeOrigin::ShapeB))
            .map(|(i, _)| i)
            .collect();

        for &ai in &a_verts {
            for &bi in &b_verts {
                let tol = self.tol();
                let dist = (self.ds.vertices[ai].point - self.ds.vertices[bi].point).length();
                if dist <= tol {
                    self.ds.interferences.push(Interference::VertexVertex {
                        v1: ai,
                        v2: bi,
                        merged_vertex: ai,
                    });
                }
            }
        }
    }

    // ─── Pass 2: Vertex-Edge ───────────────────────────────────────────

    fn perform_ve(&mut self) {
        let a_verts: Vec<usize> = self.verts_of(ShapeOrigin::ShapeA);
        let b_edges: Vec<usize> = self.edges_of(ShapeOrigin::ShapeB);

        for &vi in &a_verts {
            for &ei in &b_edges {
                self.check_vertex_edge(vi, ei);
            }
        }

        let b_verts: Vec<usize> = self.verts_of(ShapeOrigin::ShapeB);
        let a_edges: Vec<usize> = self.edges_of(ShapeOrigin::ShapeA);

        for &vi in &b_verts {
            for &ei in &a_edges {
                self.check_vertex_edge(vi, ei);
            }
        }
    }

    fn check_vertex_edge(&mut self, vi: usize, ei: usize) {
        let point = self.ds.vertices[vi].point;
        let edge = &self.ds.edges[ei];
        match &edge.curve {
            Curve3::Line(line) => {
                if let Some(t) = inttools::vertex_ops::vertex_on_line(point, line, edge.t_range) {
                    self.ds.interferences.push(Interference::VertexEdge {
                        vertex: vi,
                        edge: ei,
                        param: t,
                    });
                    self.ds.edges[ei].paves.push(Pave {
                        vertex_idx: vi,
                        param: t,
                    });
                }
            }
            Curve3::Circle(circle) => {
                // Check if point lies on the circle arc
                let v = point - circle.center;
                let dist = v.length();
                if (dist - circle.radius).abs() < TOLERANCE_ABS {
                    let on_plane = v.dot(circle.normal).abs() < TOLERANCE_ABS;
                    if on_plane {
                        // Compute angular parameter
                        let u = if circle.normal.x.abs() < 0.9 {
                            circle.normal.cross(DVec3::X).normalize()
                        } else {
                            circle.normal.cross(DVec3::Y).normalize()
                        };
                        let w = circle.normal.cross(u);
                        let theta = w.dot(v).atan2(u.dot(v));
                        let t_range = edge.t_range;
                        if theta >= t_range[0] - TOLERANCE_ABS
                            && theta <= t_range[1] + TOLERANCE_ABS
                        {
                            self.ds.interferences.push(Interference::VertexEdge {
                                vertex: vi,
                                edge: ei,
                                param: theta,
                            });
                            self.ds.edges[ei].paves.push(Pave {
                                vertex_idx: vi,
                                param: theta,
                            });
                        }
                    }
                }
            }
            _ => {}
        }
    }

    // ─── Pass 3: Edge-Edge ─────────────────────────────────────────────

    fn perform_ee(&mut self) {
        let a_edges: Vec<usize> = self.edges_of(ShapeOrigin::ShapeA);
        let b_edges: Vec<usize> = self.edges_of(ShapeOrigin::ShapeB);

        // Build a set of shared edge pairs for fast lookup when glue is enabled
        let shared_edge_set: std::collections::HashSet<(usize, usize)> = if self.use_glue {
            self.ds
                .shared_topology
                .shared_edges
                .iter()
                .map(|(e1, e2)| (*e1, *e2))
                .collect()
        } else {
            std::collections::HashSet::new()
        };

        for &ae in &a_edges {
            for &be in &b_edges {
                // Skip shared edges when glue is enabled
                if self.use_glue && shared_edge_set.contains(&(ae, be)) {
                    // Add interference for the shared edge but skip geometric intersection
                    self.ds.interferences.push(Interference::EdgeEdge {
                        e1: ae,
                        e2: be,
                        point: self.ds.vertices[self.ds.edges[ae].start_vertex].point,
                        param1: self.ds.edges[ae].t_range[0],
                        param2: self.ds.edges[be].t_range[0],
                        new_vertex: self.ds.edges[ae].start_vertex,
                    });
                    continue;
                }
                self.check_edge_edge(ae, be);
            }
        }
    }

    fn check_edge_edge(&mut self, e1: usize, e2: usize) {
        let edge1 = &self.ds.edges[e1];
        let edge2 = &self.ds.edges[e2];

        if let (Curve3::Line(l1), Curve3::Line(l2)) = (&edge1.curve, &edge2.curve)
            && let Some((t1, t2, point)) = intersect_line_line(l1, edge1.t_range, l2, edge2.t_range)
        {
            let new_v = self.ds.add_vertex(point);
            self.ds.interferences.push(Interference::EdgeEdge {
                e1,
                e2,
                point,
                param1: t1,
                param2: t2,
                new_vertex: new_v,
            });
            self.ds.edges[e1].paves.push(Pave {
                vertex_idx: new_v,
                param: t1,
            });
            self.ds.edges[e2].paves.push(Pave {
                vertex_idx: new_v,
                param: t2,
            });
        }
    }

    // ─── Pass 4: Vertex-Face ───────────────────────────────────────────

    fn perform_vf(&mut self) {
        let a_verts = self.verts_of(ShapeOrigin::ShapeA);
        let b_faces = self.faces_of(ShapeOrigin::ShapeB);

        for &vi in &a_verts {
            for &fi in &b_faces {
                self.check_vertex_face(vi, fi);
            }
        }

        let b_verts = self.verts_of(ShapeOrigin::ShapeB);
        let a_faces = self.faces_of(ShapeOrigin::ShapeA);

        for &vi in &b_verts {
            for &fi in &a_faces {
                self.check_vertex_face(vi, fi);
            }
        }
    }

    fn check_vertex_face(&mut self, vi: usize, fi: usize) {
        let point = self.ds.vertices[vi].point;
        let face = &self.ds.faces[fi];

        if let Surface3::Plane(plane) = &face.surface
            && inttools::vertex_ops::vertex_on_plane(point, plane)
        {
            let face_verts = self.ds.face_boundary_points(fi);
            if inttools::edge_face::point_in_planar_face(point, plane, &face_verts) {
                self.ds.interferences.push(Interference::VertexFace {
                    vertex: vi,
                    face: fi,
                });
                self.ds.faces[fi].face_info.vertices_on.insert(vi);
            }
        } else {
            // For curved surfaces, use closest-point projection to check if
            // the vertex lies on the surface within tolerance.
            let surface = face.surface.clone();
            if !matches!(surface, Surface3::Plane(_)) {
                let proj =
                    rcad_kernel::projection::closest_point_on_surface(&surface, point, 16);
                if proj.distance < self.tol() {
                    self.ds.interferences.push(Interference::VertexFace {
                        vertex: vi,
                        face: fi,
                    });
                }
            }
        }
    }

    // ─── Pass 5: Edge-Face ─────────────────────────────────────────────

    fn perform_ef(&mut self) {
        let a_edges = self.edges_of(ShapeOrigin::ShapeA);
        let b_faces = self.faces_of(ShapeOrigin::ShapeB);

        for &ei in &a_edges {
            for &fi in &b_faces {
                self.intersect_edge_face(ei, fi);
            }
        }

        let b_edges = self.edges_of(ShapeOrigin::ShapeB);
        let a_faces = self.faces_of(ShapeOrigin::ShapeA);

        for &ei in &b_edges {
            for &fi in &a_faces {
                self.intersect_edge_face(ei, fi);
            }
        }
    }

    fn intersect_edge_face(&mut self, edge_idx: usize, face_idx: usize) {
        let edge_curve = self.ds.edges[edge_idx].curve.clone();
        let edge_t_range = self.ds.edges[edge_idx].t_range;
        let face_surface = self.ds.faces[face_idx].surface.clone();

        // Dispatch based on curve type × surface type
        let hits: Vec<(DVec3, f64)> = match (&edge_curve, &face_surface) {
            (Curve3::Line(line), Surface3::Plane(plane)) => {
                inttools::edge_face::intersect_line_plane(line, edge_t_range, plane)
                    .into_iter()
                    .map(|h| (h.point, h.edge_param))
                    .collect()
            }
            (Curve3::Line(line), Surface3::Cylinder(cyl)) => {
                inttools::curve_surface::intersect_line_cylinder(line, edge_t_range, cyl)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            (Curve3::Line(line), Surface3::Sphere(sph)) => {
                inttools::curve_surface::intersect_line_sphere(line, edge_t_range, sph)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            (Curve3::Line(line), Surface3::Cone(cone)) => {
                inttools::curve_surface::intersect_line_cone(line, edge_t_range, cone)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            (Curve3::Circle(circle), Surface3::Plane(plane)) => {
                inttools::curve_surface::intersect_circle_plane(circle, edge_t_range, plane)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            (Curve3::Circle(circle), Surface3::Cylinder(cyl)) => {
                inttools::curve_surface::intersect_circle_cylinder(circle, edge_t_range, cyl)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            (Curve3::Circle(circle), Surface3::Sphere(sph)) => {
                inttools::curve_surface::intersect_circle_sphere(circle, edge_t_range, sph)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            (Curve3::Circle(circle), Surface3::Cone(cone)) => {
                inttools::curve_surface::intersect_circle_cone(circle, edge_t_range, cone)
                    .into_iter()
                    .map(|h| (h.point, h.curve_param))
                    .collect()
            }
            _ => {
                // Numeric fallback: sample the curve, find sign changes of the
                // surface implicit function. Works for any Curve3 × Surface3 pair.
                intersect_edge_face_numeric(&edge_curve, &face_surface, edge_t_range)
            }
        };

        for (point, edge_param) in hits {
            // Verify hit is within face boundary (for planar faces)
            let in_face = match &face_surface {
                Surface3::Plane(plane) => {
                    let face_verts = self.ds.face_boundary_points(face_idx);
                    inttools::edge_face::point_in_planar_face(point, plane, &face_verts)
                }
                _ => true,
            };

            if !in_face {
                continue;
            }

            // Skip if point is an edge endpoint
            let sv = self.ds.edges[edge_idx].start_vertex;
            let ev = self.ds.edges[edge_idx].end_vertex;
            let tol = self.tol();
            if (point - self.ds.vertices[sv].point).length() <= tol
                || (point - self.ds.vertices[ev].point).length() <= tol
            {
                continue;
            }

            let new_v = self.ds.add_vertex(point);
            self.ds.interferences.push(Interference::EdgeFace {
                edge: edge_idx,
                face: face_idx,
                point,
                edge_param,
                new_vertex: new_v,
            });
            self.ds.faces[face_idx].face_info.vertices_on.insert(new_v);
            self.ds.edges[edge_idx].paves.push(Pave {
                vertex_idx: new_v,
                param: edge_param,
            });
        }
    }

    // ─── Pass 6: Face-Face ─────────────────────────────────────────────

    fn perform_ff(&mut self) {
        let a_faces = self.faces_of(ShapeOrigin::ShapeA);
        let b_faces = self.faces_of(ShapeOrigin::ShapeB);

        if a_faces.is_empty() || b_faces.is_empty() {
            return;
        }

        if let (Some(bvh_a), Some(bvh_b)) = (self.bvh_a, self.bvh_b) {
            // Build reverse maps: BRep face index → position in a_faces/b_faces
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
            for (fa_brep, fb_brep) in candidates {
                if let (Some(&ai), Some(&bi)) = (a_rev.get(fa_brep), b_rev.get(fb_brep)) {
                    if ai != usize::MAX && bi != usize::MAX {
                        let af = a_faces[ai];
                        let bf = b_faces[bi];
                        if self.should_skip_glued_face_pair(af, bf) {
                            continue;
                        }
                        self.intersect_face_face(af, bf);
                    }
                }
            }
        } else {
            // Brute-force: all A-face × B-face pairs
            for &af in &a_faces {
                for &bf in &b_faces {
                    if self.should_skip_glued_face_pair(af, bf) {
                        continue;
                    }
                    self.intersect_face_face(af, bf);
                }
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

    fn surfaces_glue_compatible(&self, s1: &Surface3, s2: &Surface3) -> bool {
        let tol = self.glue_tolerance;
        let axis_parallel = |a: DVec3, b: DVec3| {
            let la = a.length();
            let lb = b.length();
            if la <= TOLERANCE_ABS || lb <= TOLERANCE_ABS {
                return false;
            }
            (a / la).dot(b / lb).abs() >= 0.999
        };
        match (s1, s2) {
            (Surface3::Plane(p1), Surface3::Plane(p2)) => {
                if !axis_parallel(p1.normal, p2.normal) {
                    return false;
                }
                let n = p1.normal.normalize_or_zero();
                (p2.origin - p1.origin).dot(n).abs() <= tol * 2.0
            }
            (Surface3::Sphere(s1), Surface3::Sphere(s2)) => {
                (s1.center - s2.center).length() <= tol * 2.0
                    && (s1.radius - s2.radius).abs() <= tol
            }
            (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => {
                if !axis_parallel(c1.axis, c2.axis) {
                    return false;
                }
                let a = c1.axis.normalize_or_zero();
                (c2.origin - c1.origin).cross(a).length() <= tol * 2.0
                    && (c1.radius - c2.radius).abs() <= tol
            }
            (Surface3::Cone(c1), Surface3::Cone(c2)) => {
                axis_parallel(c1.axis, c2.axis)
                    && (c1.apex - c2.apex).length() <= tol * 2.0
                    && (c1.radius - c2.radius).abs() <= tol
                    && (c1.half_angle_rad - c2.half_angle_rad).abs() <= tol
            }
            (Surface3::Torus(t1), Surface3::Torus(t2)) => {
                axis_parallel(t1.axis, t2.axis)
                    && (t1.center - t2.center).length() <= tol * 2.0
                    && (t1.major_radius - t2.major_radius).abs() <= tol
                    && (t1.minor_radius - t2.minor_radius).abs() <= tol
            }
            _ => false,
        }
    }

    fn boundaries_fully_overlap(&self, f1: usize, f2: usize) -> bool {
        let pts1 = self.ds.face_boundary_points(f1);
        let pts2 = self.ds.face_boundary_points(f2);
        if pts1.len() < 3 || pts2.len() < 3 || pts1.len() != pts2.len() {
            return false;
        }
        let tol = self.glue_tolerance;
        let mut used = vec![false; pts2.len()];
        for p1 in &pts1 {
            let mut found = false;
            for (j, p2) in pts2.iter().enumerate() {
                if used[j] {
                    continue;
                }
                if (*p1 - *p2).length() <= tol {
                    used[j] = true;
                    found = true;
                    break;
                }
            }
            if !found {
                return false;
            }
        }
        true
    }

    /// Detect partially shared edges between two faces (for enhanced glue detection).
    /// Returns a list of (edge_idx_in_f1, edge_idx_in_f2) pairs for shared edges.
    fn detect_shared_edges_between_faces(&self, f1: usize, f2: usize) -> Vec<(usize, usize)> {
        let tol = self.glue_tolerance;
        let tol_sq = tol * tol;
        let mut shared_edges = Vec::new();

        let edges1: Vec<usize> = self.ds.faces[f1].boundary_edges.iter().copied().collect();
        let edges2: Vec<usize> = self.ds.faces[f2].boundary_edges.iter().copied().collect();

        for &e1 in &edges1 {
            let edge1 = match self.ds.edges.get(e1) {
                Some(e) => e,
                None => continue,
            };

            // Get endpoints of edge1
            let pts1: Vec<DVec3> = (0..=8)
                .map(|k| {
                    let t = edge1.t_range[0] + (edge1.t_range[1] - edge1.t_range[0]) * k as f64 / 8.0;
                    edge1.curve.point_at(t)
                })
                .collect();

            for &e2 in &edges2 {
                let edge2 = match self.ds.edges.get(e2) {
                    Some(e) => e,
                    None => continue,
                };

                // Check if edges are geometrically compatible (same curve direction or reverse)
                let curve_compatible = self.edges_curve_compatible(e1, e2, tol);
                if !curve_compatible {
                    continue;
                }

                // Sample edge2 and check proximity
                let pts2: Vec<DVec3> = (0..=8)
                    .map(|k| {
                        let t = edge2.t_range[0] + (edge2.t_range[1] - edge2.t_range[0]) * k as f64 / 8.0;
                        edge2.curve.point_at(t)
                    })
                    .collect();

                // Check if sampled points are close
                let mut all_close = true;
                for p1 in &pts1 {
                    let min_dist = pts2.iter().map(|p2| (*p1 - *p2).length_squared()).fold(f64::INFINITY, f64::min);
                    if min_dist > tol_sq {
                        all_close = false;
                        break;
                    }
                }

                if all_close {
                    shared_edges.push((e1, e2));
                    break; // Each edge in f1 matches at most one in f2
                }
            }
        }

        shared_edges
    }

    /// Check if two edges have compatible curves (same geometry, possibly reversed direction).
    fn edges_curve_compatible(&self, e1: usize, e2: usize, tol: f64) -> bool {
        let edge1 = match self.ds.edges.get(e1) {
            Some(e) => e,
            None => return false,
        };
        let edge2 = match self.ds.edges.get(e2) {
            Some(e) => e,
            None => return false,
        };

        match (&edge1.curve, &edge2.curve) {
            (Curve3::Line(l1), Curve3::Line(l2)) => {
                // Check if lines are parallel (or anti-parallel)
                let d1 = l1.direction.normalize_or_zero();
                let d2 = l2.direction.normalize_or_zero();
                if d1.dot(d2).abs() < 0.999 {
                    return false;
                }
                // Check if origins are on the same line
                let v = l2.origin - l1.origin;
                let perp = v - d1 * v.dot(d1);
                perp.length() <= tol
            }
            (Curve3::Circle(c1), Curve3::Circle(c2)) => {
                // Check if circles are the same
                (c1.center - c2.center).length() <= tol
                    && c1.normal.dot(c2.normal).abs() >= 0.999
                    && (c1.radius - c2.radius).abs() <= tol
            }
            (Curve3::Ellipse(e1), Curve3::Ellipse(e2)) => {
                // Simplified ellipse compatibility check
                (e1.center - e2.center).length() <= tol
                    && e1.normal.dot(e2.normal).abs() >= 0.999
                    && (e1.major_radius - e2.major_radius).abs() <= tol
                    && (e1.minor_radius - e2.minor_radius).abs() <= tol
            }
            // For other curve types, return false (conservative)
            _ => false,
        }
    }

    /// Check if two faces have partial glue (share some edges but not full boundary).
    fn has_partial_glue(&self, f1: usize, f2: usize) -> bool {
        if !self.use_glue {
            return false;
        }

        let face1 = &self.ds.faces[f1];
        let face2 = &self.ds.faces[f2];

        // Faces must come from different original shapes
        if face1.origin == face2.origin {
            return false;
        }

        // Surfaces must be glue-compatible
        if !self.surfaces_glue_compatible(&face1.surface, &face2.surface) {
            return false;
        }

        // Check for shared edges
        let shared = self.detect_shared_edges_between_faces(f1, f2);
        !shared.is_empty()
    }

    fn intersect_face_face(&mut self, f1: usize, f2: usize) {
        let s1 = self.ds.faces[f1].surface.clone();
        let s2 = self.ds.faces[f2].surface.clone();

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
            // ── Torus × * ─────────────────────────────────────────────────
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
            _ => {
                // General case: numerical marching
                self.intersect_ff_by_marching(f1, f2);
            }
        }
    }

    fn intersect_plane_plane_faces(&mut self, f1: usize, f2: usize, p1: &Plane, p2: &Plane) {
        use inttools::pcurve_derive::line_pcurve_on_plane;

        match inttools::plane_plane::intersect_plane_plane(p1, p2) {
            inttools::plane_plane::PlanePlaneResult::Parallel => {}
            inttools::plane_plane::PlanePlaneResult::Coincident => {
                // Coplanar — handled via coplanar analysis
                self.handle_coplanar_faces(f1, f2, p1);
            }
            inttools::plane_plane::PlanePlaneResult::Line(line) => {
                let verts1 = self.ds.face_boundary_points(f1);
                let verts2 = self.ds.face_boundary_points(f2);

                let range1 = inttools::edge_face::clip_line_to_convex_polygon(&line, p1, &verts1);
                let range2 = inttools::edge_face::clip_line_to_convex_polygon(&line, p2, &verts2);

                if let (Some((t1_min, t1_max)), Some((t2_min, t2_max))) = (range1, range2) {
                    let t_min = t1_min.max(t2_min);
                    let t_max = t1_max.min(t2_max);
                    if t_max - t_min < TOLERANCE_ABS {
                        return;
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
                    });

                    self.ds.interferences.push(Interference::FaceFace {
                        f1,
                        f2,
                        curves: vec![curve_idx],
                        points: vec![],
                    });

                    self.ds.faces[f1].face_info.curves_in.insert(curve_idx);
                    self.ds.faces[f2].face_info.curves_in.insert(curve_idx);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f1].face_info.vertices_in.insert(v_end);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_start);
                    self.ds.faces[f2].face_info.vertices_in.insert(v_end);
                }
            }
        }
    }

    fn handle_coplanar_faces(&mut self, f1: usize, f2: usize, plane: &Plane) {
        let verts1 = self.ds.face_boundary_points(f1);
        let verts2 = self.ds.face_boundary_points(f2);

        let result = inttools::coplanar::analyze_coplanar_faces(&verts1, &verts2, plane);

        if !result.overlap.is_empty() {
            // Record as a FaceFace interference with no curves (coplanar overlap)
            self.ds.interferences.push(Interference::FaceFace {
                f1,
                f2,
                curves: vec![],
                points: vec![],
            });
        }
    }

    // ── Plane × Sphere analytic face-face intersection ─────────────────────────

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

        match intersect_plane_sphere(plane, sphere) {
            PlaneSphereResult::NoIntersection => {}
            PlaneSphereResult::TangentPoint(pt) => {
                let verts1 = self.ds.face_boundary_points(f1);
                let verts2 = self.ds.face_boundary_points(f2);
                if inttools::edge_face::point_in_planar_face(pt, plane, &verts1)
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
                // Sample the circle and clip to both face boundaries
                let pts = sample_circle_arc(&circle, 0.0, std::f64::consts::TAU, 32);
                if pts.len() < 2 {
                    return;
                }

                let pcurve_plane = circle_pcurve_on_plane(&circle, plane);
                // `circle_pcurve_on_sphere` is only analytically correct when the
                // sphere axis is parallel to the cutting plane normal (i.e. the
                // intersection circle is a latitude line in the sphere's UV domain).
                // When the axis is not aligned with the plane normal, we fall back to
                // projection-based sampling — exactly as `intersect_sphere_sphere_faces`
                // already does — to obtain the correct parameter-space curve.
                let axis_dot_normal = sphere
                    .axis
                    .normalize()
                    .dot(plane.normal.normalize())
                    .abs();
                let pcurve_sphere = if (axis_dot_normal - 1.0).abs() < 1e-6 {
                    // Axis is parallel to plane normal → latitude line is exact.
                    circle_pcurve_on_sphere(&circle, sphere)
                } else {
                    // Axis is NOT aligned → use projection fallback.
                    fallback_pcurve_by_projection(
                        &Curve3::Circle(circle),
                        &[0.0, std::f64::consts::TAU],
                        &Surface3::Sphere(*sphere),
                    )
                };
                let (pcurve_on_a, pcurve_on_b) = if plane_is_f1 {
                    (Some(pcurve_plane), Some(pcurve_sphere))
                } else {
                    (Some(pcurve_sphere), Some(pcurve_plane))
                };

                let v_start = self.ds.add_vertex(pts[0]);
                let v_end = self.ds.add_vertex(pts[pts.len() - 1]);

                let curve_idx = self.ds.intersection_curves.len();
                self.ds.intersection_curves.push(IntersectionCurve {
                    curve: Curve3::Circle(circle),
                    polyline: vec![],
                    start_vertex: v_start,
                    end_vertex: v_end,
                    t_range: [0.0, std::f64::consts::TAU],
                    pcurve_on_a,
                    pcurve_on_b,
                });

                self.ds.faces[f1].face_info.curves_in.insert(curve_idx);
                self.ds.faces[f2].face_info.curves_in.insert(curve_idx);
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

    // ── Sphere × Sphere analytic face-face intersection ───────────────────────

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
        if d < 1e-14 || d >= sph1.radius + sph2.radius || d <= (sph1.radius - sph2.radius).abs() {
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
        });

        self.ds.faces[f1].face_info.curves_in.insert(curve_idx);
        self.ds.faces[f2].face_info.curves_in.insert(curve_idx);
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

    // ── Sphere × Cylinder analytic face-face intersection ─────────────────────

    fn intersect_sphere_cylinder_faces(
        &mut self,
        f1: usize,
        f2: usize,
        sphere: &SphericalSurface,
        cyl: &CylindricalSurface,
    ) {
        use inttools::pcurve_derive::{circle_pcurve_on_cylinder, circle_pcurve_on_sphere};
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
                });
                ds.faces[f1].face_info.curves_in.insert(curve_idx);
                ds.faces[f2].face_info.curves_in.insert(curve_idx);
                ds.faces[f1].face_info.vertices_in.insert(v_start);
                ds.faces[f1].face_info.vertices_in.insert(v_end);
                ds.faces[f2].face_info.vertices_in.insert(v_start);
                ds.faces[f2].face_info.vertices_in.insert(v_end);
                curve_idx
            };

        // Closure to compute pcurves for one intersection circle.
        // The intersection circle is always a latitude line on the sphere
        // (φ = acos((h − h_c) / R)), so `circle_pcurve_on_sphere` is exact
        // here regardless of whether the sphere and cylinder axes are parallel.
        let make_circle_pcurves = |circle: &Circle3| -> (Option<Curve2d>, Option<Curve2d>) {
            let pcurve_sph = circle_pcurve_on_sphere(circle, sphere);
            let pcurve_cyl = circle_pcurve_on_cylinder(circle, cyl);
            make_pcurves(pcurve_sph, pcurve_cyl)
        };

        match intersect_sphere_cylinder(sphere, cyl) {
            SphereCylinderResult::NoIntersection => return,
            SphereCylinderResult::General => {
                // Fall back to numeric marching for the quartic case.
                self.intersect_ff_by_marching(f1, f2);
                return;
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
        }
    }

    // ── Cylinder × Cylinder analytic face-face intersection ──────────────────

    fn intersect_cylinder_cylinder_faces(
        &mut self,
        f1: usize,
        f2: usize,
        cyl1: &CylindricalSurface,
        cyl2: &CylindricalSurface,
    ) {
        use inttools::cylinder_cylinder::{CylinderCylinderResult, intersect_cylinder_cylinder};
        use inttools::pcurve_derive::{
            fallback_pcurve_by_projection, line_pcurve_on_cylinder,
        };
        use std::f64::consts::TAU;

        // Determine which face is cyl1 (for pcurve_on_a/b ordering)
        let cyl1_is_f1 = {
            if let Surface3::Cylinder(c) = &self.ds.faces[f1].surface {
                (c.origin - cyl1.origin).length_squared() < 1e-10
                    && (c.axis - cyl1.axis).length_squared() < 1e-10
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
                });
                ds.faces[f1].face_info.curves_in.insert(ci);
                ds.faces[f2].face_info.curves_in.insert(ci);
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
            });
            ds.faces[f1].face_info.curves_in.insert(ci);
            ds.faces[f2].face_info.curves_in.insert(ci);
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
            });
            ds.faces[f1].face_info.curves_in.insert(ci);
            ds.faces[f2].face_info.curves_in.insert(ci);
            ds.faces[f1].face_info.vertices_in.insert(v_start);
            ds.faces[f1].face_info.vertices_in.insert(v_end);
            ds.faces[f2].face_info.vertices_in.insert(v_start);
            ds.faces[f2].face_info.vertices_in.insert(v_end);
            ci
        };

        let extent = 20.0_f64;
        let mut curve_indices = Vec::new();

        match intersect_cylinder_cylinder(cyl1, cyl2) {
            CylinderCylinderResult::NoIntersection | CylinderCylinderResult::Coaxial => return,

            CylinderCylinderResult::General => {
                // Fall back to numeric marching for skew/oblique cases.
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
                // PCurves for the cylinder surfaces use projection fallback since
                // these circles are not latitude or generator lines.
                let pca1 = fallback_pcurve_by_projection(
                    &Curve3::Circle(c1),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl1),
                );
                let pcb1 = fallback_pcurve_by_projection(
                    &Curve3::Circle(c1),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl2),
                );
                let (pca1, pcb1) = make_pcurves(pca1, pcb1);
                let ci1 = add_circle(self.ds, &c1, pca1, pcb1, f1, f2);

                let pca2 = fallback_pcurve_by_projection(
                    &Curve3::Circle(c2),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl1),
                );
                let pcb2 = fallback_pcurve_by_projection(
                    &Curve3::Circle(c2),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl2),
                );
                let (pca2, pcb2) = make_pcurves(pca2, pcb2);
                let ci2 = add_circle(self.ds, &c2, pca2, pcb2, f1, f2);

                curve_indices.push(ci1);
                curve_indices.push(ci2);
            }

            CylinderCylinderResult::TwoEllipses(e1, e2) => {
                // Perpendicular Steinmetz unequal-radii.
                // Each ellipse lies in a plane; use ellipse_pcurve_on_plane for the
                // plane PCurve and fallback for the cylinder PCurve.
                let plane1 = rcad_kernel::geom::Plane {
                    origin: e1.center,
                    normal: e1.normal,
                };
                let plane2 = rcad_kernel::geom::Plane {
                    origin: e2.center,
                    normal: e2.normal,
                };

                let pca1 = fallback_pcurve_by_projection(
                    &Curve3::Ellipse(e1),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl1),
                );
                let pcb1 = fallback_pcurve_by_projection(
                    &Curve3::Ellipse(e1),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl2),
                );
                let (pca1, pcb1) = make_pcurves(pca1, pcb1);
                let ci1 = add_ellipse(self.ds, &e1, pca1, pcb1, f1, f2);

                let pca2 = fallback_pcurve_by_projection(
                    &Curve3::Ellipse(e2),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl1),
                );
                let pcb2 = fallback_pcurve_by_projection(
                    &Curve3::Ellipse(e2),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl2),
                );
                let (pca2, pcb2) = make_pcurves(pca2, pcb2);
                let ci2 = add_ellipse(self.ds, &e2, pca2, pcb2, f1, f2);

                curve_indices.push(ci1);
                curve_indices.push(ci2);

                let _ = (plane1, plane2); // suppress warnings
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

    // ── Plane × Cylinder analytic face-face intersection ──────────────────────

    fn intersect_plane_cylinder_faces(
        &mut self,
        f1: usize,
        f2: usize,
        plane: &Plane,
        cyl: &CylindricalSurface,
    ) {
        use inttools::pcurve_derive::{
            circle_pcurve_on_cylinder, circle_pcurve_on_plane, ellipse_pcurve_on_plane,
            fallback_pcurve_by_projection, line_pcurve_on_cylinder, line_pcurve_on_plane,
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
            });
            ds.faces[f1].face_info.curves_in.insert(curve_idx);
            ds.faces[f2].face_info.curves_in.insert(curve_idx);
            ds.faces[f1].face_info.vertices_in.insert(v_start);
            ds.faces[f1].face_info.vertices_in.insert(v_end);
            ds.faces[f2].face_info.vertices_in.insert(v_start);
            ds.faces[f2].face_info.vertices_in.insert(v_end);
            curve_idx
        };

        let mut curve_indices = Vec::new();

        match result {
            PlaneCylinderResult::NoIntersection => return,
            PlaneCylinderResult::TangentLine(_) => return, // zero-area intersection
            PlaneCylinderResult::TwoLines(l1, l2) => {
                // Clip each line to the face bounding-box extent
                let extent = 20.0_f64;
                let (pca1, pcb1) = make_pcurves(
                    line_pcurve_on_plane(&l1, plane),
                    line_pcurve_on_cylinder(&l1, cyl),
                );
                let ci1 = add_curve(
                    self.ds,
                    Curve3::Line(l1),
                    [-extent, extent],
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
                    [-extent, extent],
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
            }
            PlaneCylinderResult::Ellipse(ellipse) => {
                let pca_plane = ellipse_pcurve_on_plane(&ellipse, plane);
                let pcb_cyl = fallback_pcurve_by_projection(
                    &Curve3::Ellipse(ellipse),
                    &[0.0, TAU],
                    &Surface3::Cylinder(*cyl),
                );
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

    // ── Plane × Cone analytic face-face intersection ──────────────────────────

    fn intersect_plane_cone_faces(
        &mut self,
        f1: usize,
        f2: usize,
        plane: &Plane,
        cone: &ConicalSurface,
    ) {
        use inttools::pcurve_derive::{
            circle_pcurve_on_plane, ellipse_pcurve_on_plane, fallback_pcurve_by_projection,
            line_pcurve_on_plane,
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
            });
            ds.faces[f1].face_info.curves_in.insert(ci);
            ds.faces[f2].face_info.curves_in.insert(ci);
            ds.faces[f1].face_info.vertices_in.insert(v_start);
            ds.faces[f1].face_info.vertices_in.insert(v_end);
            ds.faces[f2].face_info.vertices_in.insert(v_start);
            ds.faces[f2].face_info.vertices_in.insert(v_end);
            ci
        };

        let mut curve_indices = Vec::new();

        match intersect_plane_cone(plane, cone) {
            PlaneConicalResult::NoIntersection | PlaneConicalResult::Point(_) => return,

            PlaneConicalResult::SingleLine(line) => {
                let extent = 20.0_f64;
                let pca_plane = line_pcurve_on_plane(&line, plane);
                let pcb_cone = fallback_pcurve_by_projection(
                    &Curve3::Line(line),
                    &[-extent, extent],
                    &Surface3::Cone(*cone),
                );
                let (pca, pcb) = make_pcurves(pca_plane, pcb_cone);
                let ci = add_curve(self.ds, Curve3::Line(line), [-extent, extent], pca, pcb, f1, f2);
                curve_indices.push(ci);
            }

            PlaneConicalResult::TwoLines(l1, l2) => {
                let extent = 20.0_f64;
                let pca1 = line_pcurve_on_plane(&l1, plane);
                let pcb1 = fallback_pcurve_by_projection(
                    &Curve3::Line(l1),
                    &[-extent, extent],
                    &Surface3::Cone(*cone),
                );
                let (pca1, pcb1) = make_pcurves(pca1, pcb1);
                let ci1 = add_curve(self.ds, Curve3::Line(l1), [-extent, extent], pca1, pcb1, f1, f2);

                let pca2 = line_pcurve_on_plane(&l2, plane);
                let pcb2 = fallback_pcurve_by_projection(
                    &Curve3::Line(l2),
                    &[-extent, extent],
                    &Surface3::Cone(*cone),
                );
                let (pca2, pcb2) = make_pcurves(pca2, pcb2);
                let ci2 = add_curve(self.ds, Curve3::Line(l2), [-extent, extent], pca2, pcb2, f1, f2);
                curve_indices.push(ci1);
                curve_indices.push(ci2);
            }

            PlaneConicalResult::Circle(circle) => {
                let pca_plane = circle_pcurve_on_plane(&circle, plane);
                let pcb_cone = fallback_pcurve_by_projection(
                    &Curve3::Circle(circle),
                    &[0.0, TAU],
                    &Surface3::Cone(*cone),
                );
                let (pca, pcb) = make_pcurves(pca_plane, pcb_cone);
                let ci = add_curve(self.ds, Curve3::Circle(circle), [0.0, TAU], pca, pcb, f1, f2);
                curve_indices.push(ci);
            }

            PlaneConicalResult::Ellipse(ellipse) => {
                let pca_plane = ellipse_pcurve_on_plane(&ellipse, plane);
                let pcb_cone = fallback_pcurve_by_projection(
                    &Curve3::Ellipse(ellipse),
                    &[0.0, TAU],
                    &Surface3::Cone(*cone),
                );
                let (pca, pcb) = make_pcurves(pca_plane, pcb_cone);
                let ci = add_curve(self.ds, Curve3::Ellipse(ellipse), [0.0, TAU], pca, pcb, f1, f2);
                curve_indices.push(ci);
            }

            PlaneConicalResult::Parabola(parabola) => {
                // Parabola: use fallback projection on both surfaces.
                // Parameterise the parabola over a finite t range centred at the vertex.
                let t_range = [-20.0_f64, 20.0_f64];
                let pca_plane = fallback_pcurve_by_projection(
                    &Curve3::Parabola(parabola),
                    &t_range,
                    &Surface3::Plane(*plane),
                );
                let pcb_cone = fallback_pcurve_by_projection(
                    &Curve3::Parabola(parabola),
                    &t_range,
                    &Surface3::Cone(*cone),
                );
                let (pca, pcb) = make_pcurves(pca_plane, pcb_cone);
                let ci =
                    add_curve(self.ds, Curve3::Parabola(parabola), t_range, pca, pcb, f1, f2);
                curve_indices.push(ci);
            }

            PlaneConicalResult::Hyperbola(hyperbola) => {
                // Hyperbola: use fallback projection on both surfaces.
                let t_range = [-20.0_f64, 20.0_f64];
                let pca_plane = fallback_pcurve_by_projection(
                    &Curve3::Hyperbola(hyperbola),
                    &t_range,
                    &Surface3::Plane(*plane),
                );
                let pcb_cone = fallback_pcurve_by_projection(
                    &Curve3::Hyperbola(hyperbola),
                    &t_range,
                    &Surface3::Cone(*cone),
                );
                let (pca, pcb) = make_pcurves(pca_plane, pcb_cone);
                let ci =
                    add_curve(self.ds, Curve3::Hyperbola(hyperbola), t_range, pca, pcb, f1, f2);
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

    // ── Cylinder × Cone analytic face-face intersection ───────────────────────

    fn intersect_cylinder_cone_faces(
        &mut self,
        f1: usize,
        f2: usize,
        cyl: &CylindricalSurface,
        cone: &ConicalSurface,
    ) {
        use inttools::cylinder_cone::{CylinderConeResult, intersect_cylinder_cone};
        use inttools::pcurve_derive::{circle_pcurve_on_cylinder, fallback_pcurve_by_projection};
        use std::f64::consts::TAU;

        // Determine which face carries the cylinder (for pcurve_on_a/b ordering).
        let cyl_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Cylinder(_));

        let make_pcurves = |pca: Curve2d, pcb: Curve2d| -> (Option<Curve2d>, Option<Curve2d>) {
            if cyl_is_f1 { (Some(pca), Some(pcb)) } else { (Some(pcb), Some(pca)) }
        };

        match intersect_cylinder_cone(cyl, cone) {
            CylinderConeResult::NoIntersection => return,

            CylinderConeResult::General => {
                self.intersect_ff_by_marching(f1, f2);
                return;
            }

            CylinderConeResult::CoaxialCircle(circle) => {
                let pca_cyl = circle_pcurve_on_cylinder(&circle, cyl);
                let pcb_cone = fallback_pcurve_by_projection(
                    &Curve3::Circle(circle),
                    &[0.0, TAU],
                    &Surface3::Cone(*cone),
                );
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
                });
                self.ds.faces[f1].face_info.curves_in.insert(ci);
                self.ds.faces[f2].face_info.curves_in.insert(ci);
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

    // ── Cone × Cone analytic face-face intersection ────────────────────────────

    fn intersect_cone_cone_faces(
        &mut self,
        f1: usize,
        f2: usize,
        cone1: &ConicalSurface,
        cone2: &ConicalSurface,
    ) {
        use inttools::cone_cone::{ConeConeResult, intersect_cone_cone};
        use inttools::pcurve_derive::fallback_pcurve_by_projection;
        use std::f64::consts::TAU;

        // Determine which face is cone1 (for pcurve_on_a/b ordering).
        let cone1_is_f1 = {
            if let Surface3::Cone(c) = &self.ds.faces[f1].surface {
                (c.apex - cone1.apex).length_squared() < 1e-10
                    && (c.axis - cone1.axis).length_squared() < 1e-10
            } else {
                false
            }
        };

        let make_pcurves = |pca: Curve2d, pcb: Curve2d| -> (Option<Curve2d>, Option<Curve2d>) {
            if cone1_is_f1 { (Some(pca), Some(pcb)) } else { (Some(pcb), Some(pca)) }
        };

        match intersect_cone_cone(cone1, cone2) {
            ConeConeResult::NoIntersection | ConeConeResult::Coaxial => return,

            ConeConeResult::CoaxialPoint(_pt) => {
                // Single shared apex — a point contact, not a curve.
                return;
            }

            ConeConeResult::General => {
                self.intersect_ff_by_marching(f1, f2);
                return;
            }

            ConeConeResult::CoaxialCircle(circle) => {
                let pca_cone1 = fallback_pcurve_by_projection(
                    &Curve3::Circle(circle),
                    &[0.0, TAU],
                    &Surface3::Cone(*cone1),
                );
                let pcb_cone2 = fallback_pcurve_by_projection(
                    &Curve3::Circle(circle),
                    &[0.0, TAU],
                    &Surface3::Cone(*cone2),
                );
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
                });
                self.ds.faces[f1].face_info.curves_in.insert(ci);
                self.ds.faces[f2].face_info.curves_in.insert(ci);
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

    // ── Torus intersection helpers ─────────────────────────────────────────────

    /// Generic helper: call `intersect_surfaces` and wire all results into the DS.
    /// `torus_is_f1` controls pcurve ordering.
    fn register_torus_intersection(
        &mut self,
        f1: usize,
        f2: usize,
        s1: &Surface3,
        s2: &Surface3,
        torus_is_f1: bool,
    ) {
        use inttools::intss::{intersect_surfaces, SurfaceCurve};
        use inttools::pcurve_derive::polyline_pcurve_by_projection;

        let result = intersect_surfaces(s1, s2);
        if result.is_empty() {
            return;
        }

        for sir in &result.curves {
            match &sir.curve_3d {
                SurfaceCurve::Circle(circle) => {
                    let pts = sample_circle_arc(circle, 0.0, std::f64::consts::TAU, 32);
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
                        curve: Curve3::Circle(*circle),
                        polyline: vec![],
                        start_vertex: v_start,
                        end_vertex: v_end,
                        t_range: [0.0, std::f64::consts::TAU],
                        pcurve_on_a: pca,
                        pcurve_on_b: pcb,
                    });

                    self.ds.faces[f1].face_info.curves_in.insert(curve_idx);
                    self.ds.faces[f2].face_info.curves_in.insert(curve_idx);
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
                        t_range: [0.0, arc_len.max(1e-10)],
                        pcurve_on_a: pca,
                        pcurve_on_b: pcb,
                    });

                    self.ds.faces[f1].face_info.curves_in.insert(curve_idx);
                    self.ds.faces[f2].face_info.curves_in.insert(curve_idx);
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
                    });

                    self.ds.faces[f1].face_info.curves_in.insert(curve_idx);
                    self.ds.faces[f2].face_info.curves_in.insert(curve_idx);
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
                    });

                    self.ds.faces[f1].face_info.curves_in.insert(curve_idx);
                    self.ds.faces[f2].face_info.curves_in.insert(curve_idx);
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

    /// For curved×curved face pairs, use numeric_intss_with_density (sign-change
    /// edge marching) which returns ordered polylines without the closure/drift
    /// issues of the gradient marcher.
    fn intersect_ff_by_numeric_intss(
        &mut self,
        f1: usize,
        f2: usize,
        s1: &Surface3,
        s2: &Surface3,
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

        let result = numeric_intss_with_domains(s1, s2, 64, dom1, dom2);
        if result.is_empty() {
            return;
        }

        let mut curve_indices = Vec::new();
        for sir in &result.curves {
            let mut chain = match &sir.curve_3d {
                crate::inttools::intss::SurfaceCurve::Polyline(pts) => pts.clone(),
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
                curve: Curve3::Line(Line3 {
                    origin: chain[0],
                    direction: if dir.length_squared() > 0.5 {
                        dir
                    } else {
                        DVec3::X
                    },
                }),
                polyline: chain,
                start_vertex: v_start,
                end_vertex: v_end,
                t_range: [0.0, arc_len.max(1e-10)],
                pcurve_on_a: pcurve_a,
                pcurve_on_b: pcurve_b,
            });

            self.ds.faces[f1].face_info.curves_in.insert(curve_idx);
            self.ds.faces[f2].face_info.curves_in.insert(curve_idx);
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
        use inttools::pcurve_derive::polyline_pcurve_by_projection;

        let s1 = self.ds.faces[f1].surface.clone();
        let s2 = self.ds.faces[f2].surface.clone();

        // For curved-curved surface pairs (e.g. Cylinder × Cylinder), use the
        // sign-change grid marching from numeric_intss_with_density, which produces
        // ordered polylines directly without the closure/drift issues of the gradient
        // marcher.
        let both_curved = !matches!(&s1, Surface3::Plane(_)) && !matches!(&s2, Surface3::Plane(_));
        if both_curved {
            self.intersect_ff_by_numeric_intss(f1, f2, &s1, &s2);
            return;
        }

        // Use adaptive sampling density based on surface geometry
        let base_density = 16usize;
        let sampling1 = adaptive_sampling_density(&s1, base_density);
        let sampling2 = adaptive_sampling_density(&s2, base_density);
        // Use the higher density to ensure we don't miss narrow intersections
        let n_u = sampling1.n_u.max(sampling2.n_u);
        let n_v = sampling1.n_v.max(sampling2.n_v);

        let samples = self.generate_surface_samples_grid(&s1, n_u, n_v);
        // Use multi-scale seed detection for improved robustness
        // Scales: coarse (8x8), medium (16x16), fine (32x32)
        let base_step = self.estimate_step_size(&s1, &s2);
        let seeds = inttools::marching::find_seed_points_multiscale(
            &s1,
            &s2,
            |nu, nv| self.generate_surface_samples_grid(&s1, nu, nv),
            &[8, 16, 32],
            base_step * 2.0, // dedup tolerance
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
        let step_size = base_step.min(char_len * 0.5).max(1e-6);

        // Configure marching with convergence monitoring
        let marching_config = MarchingConfig {
            step_size,
            min_step_size: step_size * 0.01,
            max_steps: 500,
            max_oscillations: 3,
            step_reduction_factor: 0.5,
            multiscale_seeds: true,
        };

        let mut curve_indices = Vec::new();
        // Track all points already covered by marched curves, to deduplicate
        // seeds that trace the same intersection curve.
        let mut covered_points: Vec<DVec3> = Vec::new();
        let dedup_tol = step_size * 3.0;

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
            let pcurve_a = polyline_pcurve_by_projection(&curve.points, &s1);
            let pcurve_b = polyline_pcurve_by_projection(&curve.points, &s2);
            self.ds.intersection_curves.push(IntersectionCurve {
                curve: Curve3::Line(Line3 {
                    origin: curve.points[0],
                    direction: if dir.length_squared() > 0.5 {
                        dir
                    } else {
                        DVec3::X
                    },
                }),
                polyline: curve.points.clone(),
                start_vertex: v_start,
                end_vertex: v_end,
                t_range: [0.0, arc_len.max(1e-10)],
                pcurve_on_a: pcurve_a,
                pcurve_on_b: pcurve_b,
            });

            self.ds.faces[f1].face_info.curves_in.insert(curve_idx);
            self.ds.faces[f2].face_info.curves_in.insert(curve_idx);
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

    fn generate_surface_samples(&self, surface: &Surface3, n1: usize, n2: usize) -> Vec<DVec3> {
        match surface {
            Surface3::Cylinder(cyl) => {
                inttools::marching::sample_cylinder(cyl, [-5.0, 5.0], n1, n2)
            }
            Surface3::Sphere(sph) => inttools::marching::sample_sphere(sph, n1, n2),
            Surface3::Torus(torus) => inttools::marching::sample_torus(torus, n1, n2),
            Surface3::Plane(plane) => sample_plane(plane, 5.0, n1),
            Surface3::Cone(cone) => sample_cone(cone, 0.01, 5.0, n1, n2),
            // Generic fallback: sample via surface.default_domain() UV grid.
            // Works for BSpline, Bezier, Offset, Revolution, Trimmed, LinearExtrusion.
            _ => sample_surface_generic(surface, n1, n2),
        }
    }

    /// Like `generate_surface_samples` but returns a structured `n_u × n_v` grid
    /// (row-major) so callers can use grid-aware adjacency for seed detection.
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
                // Rebuild in (n_u azimuth) × (n_v height) order.
                let height_range = [-5.0_f64, 5.0_f64];
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
            | Surface3::Trimmed(_)
            | Surface3::Gordon(_) => 1.0,
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
            | Surface3::Trimmed(_)
            | Surface3::Gordon(_) => 1.0,
        };
        size1.min(size2) * 0.1
    }

    // ─── Edge splitting ────────────────────────────────────────────────

    fn build_split_edges(&mut self) {
        for ei in 0..self.ds.edges.len() {
            let edge = &self.ds.edges[ei];
            if edge.paves.is_empty() {
                // No splits — single pave block spanning entire edge
                let pb = PaveBlock::new(
                    ei,
                    Pave {
                        vertex_idx: edge.start_vertex,
                        param: edge.t_range[0],
                    },
                    Pave {
                        vertex_idx: edge.end_vertex,
                        param: edge.t_range[1],
                    },
                );
                self.ds.edges[ei].pave_blocks = vec![pb];
                continue;
            }

            // Collect all paves including endpoints, sort by parameter
            let mut all_paves = vec![
                Pave {
                    vertex_idx: edge.start_vertex,
                    param: edge.t_range[0],
                },
                Pave {
                    vertex_idx: edge.end_vertex,
                    param: edge.t_range[1],
                },
            ];
            all_paves.extend_from_slice(&edge.paves);
            all_paves.sort_by(|a, b| a.param.partial_cmp(&b.param).unwrap_or(std::cmp::Ordering::Equal));

            // Deduplicate paves at the same parameter
            all_paves.dedup_by(|a, b| params_equal(a.param, b.param));

            // Create pave blocks between consecutive paves
            let mut blocks = Vec::new();
            for w in all_paves.windows(2) {
                blocks.push(PaveBlock::new(ei, w[0], w[1]));
            }
            self.ds.edges[ei].pave_blocks = blocks;
        }
    }

    // ─── Helpers ───────────────────────────────────────────────────────

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

/// Intersect two bounded line segments in 3D. Returns (t1, t2, point) if they
/// cross within tolerance.
fn intersect_line_line(
    l1: &Line3,
    r1: [f64; 2],
    l2: &Line3,
    r2: [f64; 2],
) -> Option<(f64, f64, DVec3)> {
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
        return None; // parallel
    }

    let t1 = (b * e - c * d) / denom;
    let t2 = (a * e - b * d) / denom;

    // Check within ranges
    if t1 < r1[0] - TOLERANCE_ABS
        || t1 > r1[1] + TOLERANCE_ABS
        || t2 < r2[0] - TOLERANCE_ABS
        || t2 > r2[1] + TOLERANCE_ABS
    {
        return None;
    }

    let p1 = l1.origin + d1 * t1;
    let p2 = l2.origin + d2 * t2;

    if !points_coincide(p1, p2) {
        return None; // skew, don't actually intersect
    }

    Some((t1, t2, (p1 + p2) * 0.5))
}

// ── Sampling helpers for marching seed-point generation ──────────────────────

/// Sample a flat plane (infinite) over a 2D square of side `half_extent*2`
/// centred at `plane.origin`.
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

/// Sample a cone surface between heights `h_min` and `h_max` along its axis.
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

/// Sample `n` points on a circular arc from `t_start` to `t_end`.
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

/// Check whether a point lies within the boundary of a sphere-face, defined by
/// the sphere face boundary vertices (used for tangent-point containment check).
fn point_in_sphere_face(pt: DVec3, boundary_verts: &[DVec3], _ds: &DS) -> bool {
    // Simple bounding-box check: the point should be within the convex hull
    // of the boundary vertices on the sphere surface (rough approximation).
    if boundary_verts.is_empty() {
        return true;
    }
    let cx = boundary_verts.iter().map(|v| v.x).fold(0.0_f64, f64::min)
        ..(boundary_verts.iter().map(|v| v.x).fold(0.0_f64, f64::max) + 1e-9);
    let cy = boundary_verts.iter().map(|v| v.y).fold(0.0_f64, f64::min)
        ..(boundary_verts.iter().map(|v| v.y).fold(0.0_f64, f64::max) + 1e-9);
    let cz = boundary_verts.iter().map(|v| v.z).fold(0.0_f64, f64::min)
        ..(boundary_verts.iter().map(|v| v.z).fold(0.0_f64, f64::max) + 1e-9);
    cx.contains(&pt.x) && cy.contains(&pt.y) && cz.contains(&pt.z)
}

/// Generic UV-grid sampling for any surface type via `SurfaceEval::default_domain()`.
/// Works for BSpline, Bezier, Offset, Revolution, Trimmed, LinearExtrusion.
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

/// Numeric edge-face intersection: sample the curve, find sign changes of the
/// surface implicit function, then refine via bisection.
///
/// Used as fallback for unsupported curve×surface combinations (Ellipse,
/// Hyperbola, Parabola, BSpline, Bezier, OffsetCurve × any surface).
fn intersect_edge_face_numeric(
    curve: &Curve3,
    surface: &Surface3,
    t_range: [f64; 2],
) -> Vec<(DVec3, f64)> {
    use rcad_kernel::CurveEval;
    const N_SAMPLES: usize = 64;
    const MAX_BISECT: usize = 30;

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
        // Bisection refinement
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
            if fm.abs() < 1e-12 {
                hits.push((pm, tm));
                break;
            }
            if (tb - ta).abs() < 1e-12 {
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
        if pm.is_finite() && !hits.iter().any(|(_, t)| (t - tm).abs() < 1e-6) {
            hits.push((pm, tm));
        }
    }
    hits
}

/// Result of partial face overlap analysis.
#[derive(Debug, Clone)]
pub struct PartialOverlapInfo {
    /// Face index in shape A.
    pub face_a: usize,
    /// Face index in shape B.
    pub face_b: usize,
    /// Estimated overlap ratio (0.0 to 1.0).
    pub overlap_ratio: f64,
    /// Overlap type.
    pub overlap_type: PartialOverlapType,
}

/// Type of partial overlap between faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PartialOverlapType {
    /// Faces are coplanar with partial boundary overlap.
    CoplanarBoundary,
    /// Faces share an edge partially.
    /// TODO: Implement edge overlap detection
    EdgeOverlap,
    /// One face is contained within another.
    /// TODO: Implement containment detection
    Contained,
}

impl<'a> PaveFiller<'a> {
    /// Detect partial overlaps between faces for Glue mode.
    ///
    /// This method identifies face pairs where the boundaries partially overlap,
    /// as opposed to `should_skip_glued_face_pair` which only detects complete overlaps.
    ///
    /// # Returns
    /// A vector of `PartialOverlapInfo` describing the detected partial overlaps.
    pub fn detect_partial_glue_overlaps(&self) -> Vec<PartialOverlapInfo> {
        let mut overlaps = Vec::new();
        let tol = self.tol();

        // Iterate over all face pairs from different shapes
        for f1_idx in 0..self.ds.a_face_count {
            for f2_idx in self.ds.a_face_count..self.ds.faces.len() {
                if let Some(overlap) = self.check_partial_overlap(f1_idx, f2_idx, tol) {
                    overlaps.push(overlap);
                }
            }
        }

        overlaps
    }

    fn check_partial_overlap(
        &self,
        f1_idx: usize,
        f2_idx: usize,
        tol: f64,
    ) -> Option<PartialOverlapInfo> {
        // First check if surfaces are compatible for overlap
        let face1 = &self.ds.faces[f1_idx];
        let face2 = &self.ds.faces[f2_idx];

        // Skip if same origin
        if face1.origin == face2.origin {
            return None;
        }

        // Check surface compatibility
        if !self.surfaces_glue_compatible(&face1.surface, &face2.surface) {
            return None;
        }

        // Get boundary points for both faces
        let pts1 = self.sampled_face_boundary_points(f1_idx, 8);
        let pts2 = self.sampled_face_boundary_points(f2_idx, 8);

        if pts1.is_empty() || pts2.is_empty() {
            return None;
        }

        // Compute overlap ratio by counting points near the other face's boundary
        let overlap_ratio = self.compute_boundary_overlap_ratio(&pts1, &pts2, tol);

        // Partial overlap: some but not complete
        if overlap_ratio > 0.1 && overlap_ratio < 0.99 {
            return Some(PartialOverlapInfo {
                face_a: f1_idx,
                face_b: f2_idx,
                overlap_ratio,
                overlap_type: PartialOverlapType::CoplanarBoundary,
            });
        }

        None
    }

    fn compute_boundary_overlap_ratio(&self, pts1: &[DVec3], pts2: &[DVec3], tol: f64) -> f64 {
        let proximity_tol = tol * 100.0; // More lenient for overlap detection

        // Count points from pts1 that are near pts2
        let in_2 = pts1
            .iter()
            .filter(|p| pts2.iter().any(|b| (*b - **p).length() < proximity_tol))
            .count();

        // Count points from pts2 that are near pts1
        let in_1 = pts2
            .iter()
            .filter(|p| pts1.iter().any(|b| (*b - **p).length() < proximity_tol))
            .count();

        let total = pts1.len() + pts2.len();
        if total == 0 {
            return 0.0;
        }

        (in_2 + in_1) as f64 / total as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::{BRep, PrimitiveSolid};
    use crate::bopds::ds::DS;

    #[test]
    fn glue_detects_partial_face_overlap() {
        // Two boxes that partially overlap on one face
        let box1 = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });
        let box2 = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 2.0,
            depth: 2.0,
        });

        // Translate box2 so it partially overlaps box1's face
        let mut box2_moved = box2.clone();
        for v in &mut box2_moved.vertices {
            v.point.x += 1.5; // Partial overlap
        }

        let mut ds = DS::new(&box1, &box2_moved);
        let filler = PaveFiller::new(&mut ds);

        // Should detect partial overlap on faces
        let overlaps = filler.detect_partial_glue_overlaps();
        assert!(
            !overlaps.is_empty(),
            "Should detect partial face overlaps"
        );

        // Verify the detected overlap makes sense
        for overlap in &overlaps {
            // Overlap ratio should be in partial range
            assert!(
                overlap.overlap_ratio > 0.0 && overlap.overlap_ratio < 1.0,
                "Overlap ratio should be partial, got {}",
                overlap.overlap_ratio
            );
            // Type should be CoplanarBoundary for box-box overlap
            assert_eq!(overlap.overlap_type, PartialOverlapType::CoplanarBoundary);
        }
    }
}
