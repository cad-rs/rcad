use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval};
use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_CLAMP_MIN};

/// �?OCCT-aligned: IntTools_CommonPrt (IntTools_CommonPrt.hxx L32-128).
/// Describes a common part between two edges: either a VERTEX-type (point)
/// or EDGE-type (overlapping segment) intersection.
#[derive(Debug, Clone)]
pub struct CommonPrt {
    /// Type of common part: `false` = VERTEX (point), `true` = EDGE (segment).
    pub is_edge_type: bool,
    /// Range on first edge `[t1, t2]`.
    pub range1: [f64; 2],
    /// Ranges on second edge �?sequence of `[t1, t2]` supporting 1-to-N mapping.
    /// OCCT: `NCollection_Sequence<IntTools_Range> myRanges2`.
    pub ranges2: Vec<[f64; 2]>,
    /// Parameter of the first vertex on edge1 (for VERTEX-type).
    pub vertex_param1: f64,
    /// Parameter of the second vertex on edge2 (for VERTEX-type).
    pub vertex_param2: f64,
    /// Bounding points of the common part (3D).
    pub bounding_point1: DVec3,
    pub bounding_point2: DVec3,
    /// All-null flag: true when both ranges are degenerate/null.
    pub all_null_flag: bool,
}

impl CommonPrt {
    /// Create a VERTEX-type common part (point intersection).
    pub fn new_vertex(t1: f64, t2: f64, p: DVec3) -> Self {
        Self {
            is_edge_type: false,
            range1: [t1, t1],
            ranges2: vec![[t2, t2]],
            vertex_param1: t1,
            vertex_param2: t2,
            bounding_point1: p,
            bounding_point2: p,
            all_null_flag: false,
        }
    }

    /// Create an EDGE-type common part (overlapping segment).
    pub fn new_edge(range1: [f64; 2], ranges2: Vec<[f64; 2]>, p1: DVec3, p2: DVec3) -> Self {
        let vp2 = ranges2.first().map_or(0.0, |r| r[0]);
        Self {
            is_edge_type: true,
            range1,
            ranges2,
            vertex_param1: range1[0],
            vertex_param2: vp2,
            bounding_point1: p1,
            bounding_point2: p2,
            all_null_flag: false,
        }
    }
}

/// �?OCCT-aligned: IntTools_EdgeEdge::TypeToInteger (cxx L1456-1482).
/// Maps curve type to integer priority for edge swapping (lower = simpler).
pub fn curve_type_to_integer(curve: &Curve3) -> i32 {
    match curve {
        Curve3::Line(_) => 0,
        Curve3::Hyperbola(_) | Curve3::Parabola(_) => 1,
        Curve3::Circle(_) | Curve3::Ellipse(_) => 2,
        Curve3::BSpline(_) | Curve3::Bezier(_) => 3,
        _ => 4,
    }
}

/// �?OCCT-aligned: IntTools_EdgeEdge::PointBoxDistance (cxx L1423-1452).
/// Computes min distance from a point to an axis-aligned bounding box.
pub fn point_box_distance(p: DVec3, box_min: DVec3, box_max: DVec3) -> f64 {
    let mut dist = 0.0;
    for i in 0..3 {
        let c = [p.x, p.y, p.z][i];
        let bmin = [box_min.x, box_min.y, box_min.z][i];
        let bmax = [box_max.x, box_max.y, box_max.z][i];
        if c < bmin {
            let d = bmin - c;
            dist += d * d;
        } else if c > bmax {
            let d = c - bmax;
            dist += d * d;
        }
    }
    dist.sqrt()
}

/// �?OCCT-aligned: IntTools_EdgeEdge::SplitRangeOnSegments (cxx L1366-1406).
/// Splits range [aT1, aT2] into segments based on resolution. Returns number of segments.
pub fn split_range_on_segments(t1: f64, t2: f64, resolution: f64, nb_seg: i32) -> (i32, Vec<[f64; 2]>) {
    let diff = t2 - t1;
    if diff < resolution || nb_seg == 1 {
        return (1, vec![[t1, t2]]);
    }
    let mut a_nb_segments = nb_seg;
    let mut a_dt = diff / a_nb_segments as f64;
    if a_dt < resolution {
        let seg = (diff / resolution) as i32;
        a_nb_segments = seg + 1;
        a_dt = diff / a_nb_segments as f64;
    }
    let mut segments = Vec::new();
    let mut t1x = t1;
    for _ in 1..a_nb_segments {
        let t2x = t1x + a_dt;
        segments.push([t1x, t2x]);
        t1x = t2x;
    }
    segments.push([t1x, t2]);
    (a_nb_segments, segments)
}

/// �?OCCT-aligned: IntTools_EdgeEdge::Resolution (cxx L1561-1607).
/// Computes curve resolution (parameter step for a given 3D tolerance).
/// For lines: returns theR3D directly.
/// For circles: 2*asin(res_coeff * theR3D).
/// For BSpline/Bezier: delegates to curve resolution method.
pub fn curve_resolution_edge(curve: &Curve3, res_coeff: f64, r3d: f64) -> f64 {
    match curve {
        Curve3::Line(_) => r3d,
        Curve3::Circle(c) => {
            let dt = res_coeff * r3d;
            if dt <= 1.0 { 2.0 * dt.asin() } else { std::f64::consts::TAU }
        }
        Curve3::Ellipse(e) => {
            let dt = res_coeff * r3d;
            if dt <= 1.0 { 2.0 * dt.asin() } else { std::f64::consts::TAU }
        }
        Curve3::BSpline(_) | Curve3::Bezier(_) => {
            // OCCT: theCurve->Resolution(theR3D, aRes)
            // rcad: approximate using curve resolution
            crate::inttools::curve_range::curve_resolution(curve, 0.0, r3d)
        }
        _ => res_coeff * r3d,
    }
}

/// �?OCCT-aligned: IntTools_EdgeEdge::ResolutionCoeff (cxx L1486-1557).
/// Computes the resolution coefficient for a curve type.
/// For circles: 1/(2*radius). For ellipses: 1/major_radius.
pub fn resolution_coeff(curve: &Curve3, t_range: [f64; 2]) -> f64 {
    match curve {
        Curve3::Circle(c) => 1.0 / (2.0 * c.radius.max(TOLERANCE_LEN_SQ_DIV_SAFE)),
        Curve3::Ellipse(e) => 1.0 / e.major_radius.max(TOLERANCE_LEN_SQ_DIV_SAFE),
        _ => {
            // OCCT: sample 30 points, find min dt/dist ratio
            let nb_p = 30;
            let t1 = t_range[0];
            let t2 = t_range[1];
            let dt = (t2 - t1) / nb_p as f64;
            let mut t = t1;
            let mut p1 = curve.point_at(t1);
            let mut k_min = 10.0;
            for _ in 1..=nb_p {
                t += dt;
                let p2 = curve.point_at(t);
                let dist = (p1 - p2).length();
                if dist > TOLERANCE_LEN_SQ_DIV_SAFE {
                    let k = dt / dist;
                    if k < k_min { k_min = k; }
                }
                p1 = p2;
            }
            k_min
        }
    }
}

/// �?OCCT-aligned: IntTools_EdgeEdge::CurveDeflection (cxx L1611-1638).
/// Computes total angular deflection of a curve over its range by sampling.
pub fn curve_deflection(curve: &Curve3, t_range: [f64; 2]) -> f64 {
    let nb_p = 10;
    let t1 = t_range[0];
    let t2 = t_range[1];
    let dt = (t2 - t1) / nb_p as f64;
    let mut t = t1;
    let mut v1 = curve.tangent_at(t1);
    let mut defl = 0.0;
    for _ in 1..=nb_p {
        t += dt;
        let v2 = curve.tangent_at(t);
        let len1 = v1.length_squared();
        let len2 = v2.length_squared();
        if len1 > TOLERANCE_LEN_SQ_DIV_SAFE && len2 > TOLERANCE_LEN_SQ_DIV_SAFE {
            let d1 = v1 / len1.sqrt();
            let d2 = v2 / len2.sqrt();
            defl += d1.dot(d2).acos();
        }
        v1 = v2;
    }
    defl
}

/// �?OCCT-aligned: IsClosed (IntTools_EdgeEdge.cxx L1642-1659).
/// Checks if the curve segment between aT1 and aT2 is closed.
pub fn is_curve_segment_closed(curve: &Curve3, t1: f64, t2: f64, tol: f64, res: f64) -> bool {
    if (t1 - t2).abs() < res { return false; }
    let p1 = curve.point_at(t1);
    let p2 = curve.point_at(t2);
    (p1 - p2).length() < tol
}

/// �?OCCT-aligned: ComputeLineLine common part detection (cxx L902-1056).
/// Determines if two line segments intersect, returning parameters if they do.
pub fn intersect_line_line_3d(
    l1_origin: DVec3, l1_dir: DVec3, t1_range: [f64; 2],
    l2_origin: DVec3, l2_dir: DVec3, t2_range: [f64; 2],
    tol: f64,
) -> Option<([f64; 2], [f64; 2], bool)> {
    let d1 = l1_dir.normalize();
    let d2 = l2_dir.normalize();
    let angle = d1.dot(d2).acos();
    let ang_tol = 1e-12;
    let is_coincide = angle < ang_tol || (std::f64::consts::PI - angle).abs() < ang_tol;

    if is_coincide {
        // OCCT L916-919: check distance between lines
        let dist = (l2_origin - l1_origin).cross(d1).length();
        if dist > tol { return None; }
        // Project both endpoints onto line1
        let t21 = (l2_origin - l1_origin).dot(d1);
        let t22 = t21 + (l2_origin + d2 - l1_origin).dot(d1);
        let (mut t21, mut t22) = if t21 < t22 { (t21, t22) } else { (t22, t21) };
        let [t11, t12] = t1_range;
        if (t21 > t12 && t22 > t12) || (t21 < t11 && t22 < t11) { return None; }
        t21 = t21.max(t11);
        t22 = t22.min(t12);
        let range1 = [t11, t12];
        let range2 = [t21, t22];
        return Some((range1, range2, true));
    }

    // Non-coincident lines: find intersection point
    let cross = d1.cross(d2);
    let cross_len2 = cross.length_squared();
    if cross_len2 < TOLERANCE_LEN_SQ_DIV_SAFE { return None; }
    let o1o2 = l2_origin - l1_origin;
    let dist_ll = o1o2.dot(cross / cross_len2.sqrt()).abs();
    if dist_ll > tol { return None; }

    // Find parameters of closest approach
    let a = d1.dot(d2);
    let b = d1.dot(o1o2);
    let c = d2.dot(o1o2);
    let denom = 1.0 - a * a;
    if denom.abs() < TOLERANCE_CLAMP_MIN { return None; }
    let t1 = (b - a * c) / denom;
    let t2 = (a * b - c) / denom;

    let [t11, t12] = t1_range;
    let [t21, t22] = t2_range;
    if t1 < t11 || t1 > t12 || t2 < t21 || t2 > t22 { return None; }

    let p1 = l1_origin + d1 * t1;
    let p2 = l2_origin + d2 * t2;
    if (p1 - p2).length_squared() > tol * tol { return None; }

    // OCCT L1047-1055: compute intersection range with ComputeIntRange
    let a_tol_1 = tol;
    let a_tol_2 = tol;
    let a_dt1 = crate::boptools::compute_int_range(a_tol_1, a_tol_2, angle);
    let a_dt2 = crate::boptools::compute_int_range(a_tol_2, a_tol_1, angle);

    let range1 = [t1 - a_dt1, t1 + a_dt1];
    let range2 = [t2 - a_dt2, t2 + a_dt2];
    Some((range1, range2, false))
}

/// �?OCCT-aligned: IntTools_EdgeEdge full class.
///
/// Edge/Edge intersection engine based on bounding box refinement.
/// Algorithm:
///   1. CheckData + Prepare (swap edges, set tolerances/resolution)
///   2. ComputeLineLine if both lines
///   3. FindSolutions: recursive range refinement using bounding boxes
///      �?FindParameters (golden-section search for range within box)
///   4. MergeSolutions: merge overlapping ranges
///   5. AddSolution: create CommonPart for each solution
///
/// rcad: uses DS edge indices instead of TopoDS_Edge handles.
pub struct EdgeEdgeIntersector {
    /// First edge DS index
    edge1: usize,
    /// Second edge DS index
    edge2: usize,
    /// Curve of first edge (3D)
    curve1: Curve3,
    /// Curve of second edge (3D)
    curve2: Curve3,
    /// Tolerance for edge1 (incl. fuzzy/2)
    tol1: f64,
    /// Tolerance for edge2 (incl. fuzzy/2)
    tol2: f64,
    /// Combined tolerance
    tol: f64,
    /// Fuzzy value
    fuzzy_value: f64,
    /// Resolution for edge1
    res1: f64,
    /// Resolution for edge2
    res2: f64,
    /// Resolution coefficient for edge1
    res_coeff1: f64,
    /// Resolution coefficient for edge2
    res_coeff2: f64,
    /// Parameter tolerance for edge1
    p_tol1: f64,
    /// Parameter tolerance for edge2
    p_tol2: f64,
    /// Range on edge1
    range1: [f64; 2],
    /// Range on edge2
    range2: [f64; 2],
    /// True if edges were swapped (simpler curve first)
    swapped: bool,
    /// Error status (0=ok)
    error_status: i32,
    /// Result common parts as `CommonPrt` entries.
    common_parts: Vec<CommonPrt>,
    /// Quick coincidence check flag
    quick_coincidence_check: bool,
    /// DS reference for edge data
    ds_edges: Vec<Curve3>,
    ds_t_ranges: Vec<[f64; 2]>,
}

impl EdgeEdgeIntersector {
    /// Create empty intersector.
    pub fn new() -> Self {
        Self {
            edge1: usize::MAX,
            edge2: usize::MAX,
            curve1: Curve3::Line(rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X }),
            curve2: Curve3::Line(rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X }),
            tol1: 0.0, tol2: 0.0, tol: 0.0,
            fuzzy_value: 0.0,
            res1: 0.0, res2: 0.0,
            res_coeff1: 0.0, res_coeff2: 0.0,
            p_tol1: 0.0, p_tol2: 0.0,
            range1: [0.0, 0.0], range2: [0.0, 0.0],
            swapped: false,
            error_status: 0,
            common_parts: Vec::new(),
            quick_coincidence_check: false,
            ds_edges: Vec::new(),
            ds_t_ranges: Vec::new(),
        }
    }

    /// Set edges from DS with ranges. Returns self for chaining.
    pub fn set_edges(
        &mut self,
        ei1: usize,
        range1: [f64; 2],
        ei2: usize,
        range2: [f64; 2],
        ds: &crate::bopds::ds::DS,
    ) -> &mut Self {
        self.edge1 = ei1;
        self.edge2 = ei2;
        self.curve1 = ds.edges[ei1].curve.clone();
        self.curve2 = ds.edges[ei2].curve.clone();
        self.range1 = range1;
        self.range2 = range2;
        self.ds_edges = ds.edges.iter().map(|e| e.curve.clone()).collect();
        self.ds_t_ranges = ds.edges.iter().map(|e| e.t_range).collect();
        self
    }

    /// Set fuzzy value.
    pub fn set_fuzzy_value(&mut self, fuzz: f64) { self.fuzzy_value = fuzz; }

    /// Enable/disable quick coincidence check.
    pub fn use_quick_coincidence_check(&mut self, b: bool) { self.quick_coincidence_check = b; }

    /// Perform the intersection.
    pub fn perform(&mut self) {
        self.error_status = 0;
        self.common_parts.clear();
        // 1. CheckData
        if self.edge1 == usize::MAX || self.edge2 == usize::MAX {
            self.error_status = 1;
            return;
        }
        // 2. Prepare
        self.prepare();
        if self.error_status != 0 { return; }

        // 3.1. Line/Line case
        if matches!(self.curve1, Curve3::Line(_)) && matches!(self.curve2, Curve3::Line(_)) {
            self.compute_line_line();
            return;
        }

        // 3.2. Quick coincidence check
        if self.quick_coincidence_check {
            if self.is_coincident() {
                self.add_solution(self.range1[0], self.range1[1],
                    self.range2[0], self.range2[1], true);
                return;
            }
        }

        // 3.3. Find solutions via range refinement
        let mut ranges1: Vec<[f64; 2]> = Vec::new();
        let mut ranges2: Vec<[f64; 2]> = Vec::new();
        let mut b_split2 = false;
        self.find_solutions(&mut ranges1, &mut ranges2, &mut b_split2);

        // 4. Merge solutions
        self.merge_solutions(&ranges1, &ranges2, b_split2);
    }

    /// Returns true if common parts found.
    pub fn is_done(&self) -> bool { !self.common_parts.is_empty() }

    /// Returns common parts.
    pub fn common_parts(&self) -> &[CommonPrt] { &self.common_parts }

    // ===== Private methods =====

    fn prepare(&mut self) {
        // OCCT L94-147: set curve adaptors, swap edges by type priority
        let ct1 = curve_type_to_integer(&self.curve1);
        let ct2 = curve_type_to_integer(&self.curve2);

        // OCCT L99-107: set default ranges if needed
        if self.range1[1] <= self.range1[0] {
            self.range1 = [0.0, 1.0]; // fallback
        }
        if self.range2[1] <= self.range2[0] {
            self.range2 = [0.0, 1.0];
        }

        // OCCT L115-130: adjust type priority by curve deflection
        let mut i_ct1 = ct1;
        let mut i_ct2 = ct2;
        if i_ct1 == i_ct2 && i_ct1 != 0 {
            let c2 = curve_deflection(&self.curve2, self.range2);
            let c1 = if c2 > TOLERANCE_ABS { curve_deflection(&self.curve1, self.range1) } else { 1.0 };
            if c1 < c2 { i_ct1 -= 1; }
        }

        // OCCT L132-147: swap so simpler curve is first
        if i_ct1 < i_ct2 {
            std::mem::swap(&mut self.edge1, &mut self.edge2);
            std::mem::swap(&mut self.curve1, &mut self.curve2);
            std::mem::swap(&mut self.range1, &mut self.range2);
            self.swapped = true;
        }

        // OCCT L149-152: set tolerances
        let a_tol_add = self.fuzzy_value / 2.0;
        let edge_tol1 = TOLERANCE_ABS; // approximate geom_tol
        let edge_tol2 = TOLERANCE_ABS;
        self.tol1 = edge_tol1 + a_tol_add;
        self.tol2 = edge_tol2 + a_tol_add;
        self.tol = self.tol1 + self.tol2;

        // OCCT L154-180: compute resolution coefficients and resolutions
        if i_ct1 != 0 || i_ct2 != 0 {
            self.res_coeff1 = resolution_coeff(&self.curve1, self.range1);
            self.res_coeff2 = resolution_coeff(&self.curve2, self.range2);
            self.res1 = curve_resolution_edge(&self.curve1, self.res_coeff1, self.tol1);
            self.res2 = curve_resolution_edge(&self.curve2, self.res_coeff2, self.tol2);

            // Parameter tolerances
            self.p_tol1 = 5e-13;
            let tm1 = self.range1[0].abs().max(self.range1[1].abs());
            if tm1 > 999.0 { self.p_tol1 = 5e-16 * tm1; }
            self.p_tol2 = 5e-13;
            let tm2 = self.range2[0].abs().max(self.range2[1].abs());
            if tm2 > 999.0 { self.p_tol2 = 5e-16 * tm2; }
        }
    }

    /// OCCT L247-286: IsCoincident �?sample 23 points, check projection distance.
    fn is_coincident(&self) -> bool {
        let a_tresh = 0.5;
        let a_nb_seg = 23usize;
        let [t11, t12] = self.range1;
        let [t21, t22] = self.range2;
        let dt = (t12 - t11) / a_nb_seg as f64;

        let mut i_cnt = 0;
        for i in 0..=a_nb_seg {
            let t1 = t11 + i as f64 * dt;
            let p1 = self.curve1.point_at(t1);
            // Project onto curve2
            let proj = rcad_kernel::closest_point_on_curve(&self.curve2, p1, 16);
            let d = proj.distance;
            if d < self.tol { i_cnt += 1; }
        }
        let coeff = i_cnt as f64 / (a_nb_seg + 1) as f64;
        coeff > a_tresh
    }

    /// OCCT L902-1056: ComputeLineLine
    fn compute_line_line(&mut self) {
        let Curve3::Line(l1) = &self.curve1 else { return; };
        let Curve3::Line(l2) = &self.curve2 else { return; };
        let d1 = l1.direction.normalize();
        let d2 = l2.direction.normalize();
        let angle = d1.dot(d2).acos();
        let ang_tol = 1e-12;
        let is_coincide = angle < ang_tol || (std::f64::consts::PI - angle).abs() < ang_tol;

        let [t11, t12] = self.range1;
        let [t21, t22] = self.range2;
        let a_tol = self.tol * self.tol;

        if is_coincide {
            // OCCT L916-919: check distance between lines
            let dist = (l2.origin - l1.origin).cross(d1).length();
            if dist > self.tol { return; }
            // OCCT L958-993: project endpoints
            let p11 = l1.origin + d1 * t11;
            let p12 = l1.origin + d1 * t12;
            let t2a = (p11 - l2.origin).dot(d2);
            let t2b = (p12 - l2.origin).dot(d2);
            let (mut t21p, mut t22p) = if t2a < t2b { (t2a, t2b) } else { (t2b, t2a) };
            if (t21p > t22 && t22p > t22) || (t21p < t21 && t22p < t21) { return; }
            t21p = t21p.max(t21).min(t22);
            t22p = t22p.max(t21).min(t22);
            self.add_solution(t11, t12, t21p, t22p, true);
            return;
        }

        // OCCT L922-947: check coincidence at endpoints
        let p11 = l1.origin + d1 * t11;
        let p12 = l1.origin + d1 * t12;
        let o2 = l2.origin;
        let v1 = (o2 - p11).cross(d2);
        let v2 = (o2 - p12).cross(d2);
        let d_sq1 = v1.length_squared();
        let d_sq2 = v2.length_squared();
        let is_coincide2 = d_sq1 <= a_tol && d_sq2 <= a_tol;
        if !is_coincide2 && v1.dot(v2) > 0.0 { return; }

        // OCCT L996-1055: find exact intersection
        let cross = d1.cross(d2);
        let cross_len2 = cross.length_squared();
        if cross_len2 < TOLERANCE_LEN_SQ_DIV_SAFE { return; }
        let denom = 1.0 - (d1.dot(d2)).powi(2);
        if denom.abs() < TOLERANCE_CLAMP_MIN { return; }
        let o1o2 = l2.origin - l1.origin;
        let t2 = (d1.dot(o1o2) * d1.dot(d2) - o1o2.dot(d2)) / denom;
        if t2 < t21 || t2 > t22 { return; }
        let p2 = l2.origin + d2 * t2;
        let t1 = (p2 - l1.origin).dot(d1);
        if t1 < t11 || t1 > t12 { return; }
        let p1 = l1.origin + d1 * t1;
        let dist = (p1 - p2).length_squared();
        if dist > a_tol { return; }

        let a_dt1 = crate::boptools::compute_int_range(self.tol1, self.tol2, angle);
        let a_dt2 = crate::boptools::compute_int_range(self.tol2, self.tol1, angle);
        self.add_solution(t1 - a_dt1, t1 + a_dt1, t2 - a_dt2, t2 + a_dt2, false);
    }

    /// OCCT L290-349: FindSolutions (top-level dispatch)
    fn find_solutions(
        &mut self,
        ranges1: &mut Vec<[f64; 2]>,
        ranges2: &mut Vec<[f64; 2]>,
        b_split2: &mut bool,
    ) {
        let [t11, t12] = self.range1;
        let [t21, t22] = self.range2;
        *b_split2 = false;

        let b_is_closed2 = is_curve_segment_closed(&self.curve2, t21, t22, self.tol2, self.res2);
        if b_is_closed2 {
            // Build box for curve1, check if curve2 start point is inside
            let box1_min = self.curve_aabb(&self.curve1, t11, t12, self.tol1);
            let p2_start = self.curve2.point_at(t21);
            if !point_in_aabb(p2_start, box1_min.0, box1_min.1) {
                // Not closed for this intersection
            }
        }

        // OCCT L312-317: direct call if not closed
        let box1 = self.curve_aabb(&self.curve1, t11, t12, self.tol1);
        let box2 = self.curve_aabb(&self.curve2, t21, t22, self.tol2);
        self.find_solutions_rec(t11, t12, box1, t21, t22, box2, ranges1, ranges2);
    }

    /// OCCT L353-549: FindSolutions (recursive)
    fn find_solutions_rec(
        &self,
        t11: f64, t12: f64,
        box1: (DVec3, DVec3),
        t21: f64, t22: f64,
        box2: (DVec3, DVec3),
        ranges1: &mut Vec<[f64; 2]>,
        ranges2: &mut Vec<[f64; 2]>,
    ) {
        let mut a_t11 = t11; let mut a_t12 = t12;
        let mut a_t21 = t21; let mut a_t22 = t22;
        let mut a_b1 = box1;
        let mut a_b2 = box2;
        let mut b_stop = false;
        let mut b_out = false;
        let mut thin = false;

        let mut iter = 0i32;
        while !b_stop && iter < 50 {
            iter += 1;
            let tb11 = a_t11; let tb12 = a_t12;
            let tb21 = a_t21; let tb22 = a_t22;

            // 1. Check box overlap
            if !aabb_overlap(a_b1, a_b2) { b_out = true; break; }

            // 2. Check thin status
            let is_thin = (a_t12 - a_t11) < self.res1
                || (a_b1.1.x - a_b1.0.x < self.tol && a_b1.1.y - a_b1.0.y < self.tol
                    && a_b1.1.z - a_b1.0.z < self.tol);
            thin = is_thin;

            // 3. Find parameters of curve2 within box1
            let found = self.find_parameters(&self.curve2, tb21, tb22, self.tol2, self.res2,
                self.p_tol2, self.res_coeff2, a_b1, &mut a_t21, &mut a_t22);
            if !found || thin { b_out = !found; break; }

            // 4. Rebuild box2 and check
            a_b2 = self.curve_aabb(&self.curve2, a_t21, a_t22, self.tol2);
            if !aabb_overlap(a_b1, a_b2) { b_out = true; break; }

            let is_thin2 = (a_t22 - a_t21) < self.res2
                || (a_b2.1.x - a_b2.0.x < self.tol && a_b2.1.y - a_b2.0.y < self.tol
                    && a_b2.1.z - a_b2.0.z < self.tol);
            thin = is_thin2;

            // 5. Find parameters of curve1 within box2
            let found2 = self.find_parameters(&self.curve1, tb11, tb12, self.tol1, self.res1,
                self.p_tol1, self.res_coeff1, a_b2, &mut a_t11, &mut a_t12);
            if !found2 || thin { b_out = !found2; break; }

            // 6. Check convergence
            let small_step1 = (tb12 - tb11) / 250.0;
            let small_step2 = (tb22 - tb21) / 250.0;
            let ss1 = small_step1.max(self.res1);
            let ss2 = small_step2.max(self.res2);
            if ((a_t11 - tb11).abs() < ss1 && (tb12 - a_t12).abs() < ss1
                && (a_t21 - tb21).abs() < ss2 && (tb22 - a_t22).abs() < ss2) {
                b_stop = true;
            } else {
                a_b1 = self.curve_aabb(&self.curve1, a_t11, a_t12, self.tol1);
            }
        }

        if b_out { return; }

        // OCCT L468-476: CheckCoincidence on refined ranges
        let mut i_com = 0;
        if !thin {
            i_com = self.check_coincidence(a_t11, a_t12, a_t21, a_t22, self.tol, self.res1);
            if i_com == 0 { thin = true; }
        }

        if thin {
            // OCCT L478-520: check intermediate point
            if i_com != 0 {
                let at = (a_t11 + a_t12) * 0.5;
                let p1 = self.curve1.point_at(at);
                let proj = rcad_kernel::closest_point_on_curve(&self.curve2, p1, 16);
                let b_sol = proj.distance <= self.tol;
                if !b_sol { return; }
            }
            ranges1.push([a_t11, a_t12]);
            ranges2.push([a_t21, a_t22]);
            return;
        }

        // OCCT L522-549: split and recurse
        if !self.is_intersection(a_t11, a_t12, a_t21, a_t22) { return; }

        let (nb1, segs1) = split_range_on_segments(a_t11, a_t12, self.res1, 3);
        let a_b2 = self.curve_aabb(&self.curve2, a_t21, a_t22, self.tol2);
        let a_b1_sq_extent = aabb_sq_extent(a_b1);

        for i in 0..nb1 as usize {
            let [r1, r2] = segs1[i];
            let a_b1_seg = self.curve_aabb(&self.curve1, r1, r2, self.tol1);
            if aabb_overlap(a_b1_seg, a_b2) && (nb1 == 1 || aabb_sq_extent(a_b1_seg) < a_b1_sq_extent) {
                self.find_solutions_rec(r1, r2, a_b1_seg, a_t21, a_t22, a_b2,
                    ranges1, ranges2);
            }
        }
    }

    /// OCCT L553-671: FindParameters �?golden-section search for range within a bounding box.
    fn find_parameters(
        &self,
        curve: &Curve3,
        t1: f64, t2: f64,
        tol_val: f64, res: f64, p_tol: f64, res_coeff: f64,
        cbox: (DVec3, DVec3),
        tb1: &mut f64, tb2: &mut f64,
    ) -> bool {
        const CF: f64 = 0.6180339887498948482045868343656;
        let a_cfx = cbox;
        let a_tol_box = tol_val;
        let max_dt = (t2 - t1) * 0.01;

        for side in 0..2i32 {
            let mut a_tb = if side == 0 { t1 } else { t2 };
            let a_t_end = if side == 0 { t2 } else { t1 };
            let a_c: f64 = if side == 0 { 1.0 } else { -1.0 };
            let mut a_dt = res;
            let mut a_dist_p = 0.0;
            let mut b_ret = false;
            let mut k = 1.0;

            while a_c * (a_t_end - a_tb) >= 0.0 {
                let p = curve.point_at(a_tb);
                let a_dist = point_box_distance(p, a_cfx.0, a_cfx.1);
                if a_dist > a_tol_box {
                    if a_dist_p > 0.0 {
                        if (a_dist_p - a_dist).abs() / a_dist_p < 0.1 {
                            let a_dt_new = curve_resolution_edge(curve, res_coeff, k * a_dist);
                            if a_dt_new < max_dt { k *= 2.0; a_dt = a_dt_new; }
                            else { k = 1.0; a_dt = curve_resolution_edge(curve, res_coeff, a_dist); }
                        } else {
                            k = 1.0;
                            a_dt = curve_resolution_edge(curve, res_coeff, a_dist);
                        }
                    }
                    a_tb += a_c * a_dt;
                } else {
                    b_ret = true;
                    break;
                }
                a_dist_p = a_dist;
            }
            if !b_ret {
                if side == 0 { return false; }
                else { b_ret = true; a_tb = *tb1; a_dt = t2 - *tb1; }
            }
            let a_t = if side == 0 { t1 } else { t2 };
            if (a_tb - a_t).abs() > TOLERANCE_LEN_SQ_DIV_SAFE {
                // Refine with golden section
                let mut t_in = a_tb;
                let mut t_out = a_tb - a_c * a_dt;
                let mut diff = t_in - t_out;
                while diff.abs() > p_tol {
                    let t_mid = t_out + diff * CF;
                    let p = curve.point_at(t_mid);
                    if point_box_distance(p, a_cfx.0, a_cfx.1) > a_tol_box {
                        t_out = t_mid;
                    } else {
                        t_in = t_mid;
                    }
                    diff = t_in - t_out;
                }
                a_tb = t_in;
            }
            if side == 0 { *tb1 = a_tb; } else { *tb2 = a_tb; }
        }
        true
    }

    /// OCCT L675-776: MergeSolutions
    fn merge_solutions(
        &mut self,
        ranges1: &[[f64; 2]],
        ranges2: &[[f64; 2]],
        b_split2: bool,
    ) {
        let a_nb_cp = ranges1.len();
        if a_nb_cp == 0 { return; }
        let [t11, t12] = self.range1;
        let [t21, t22] = self.range2;
        let a_res1 = curve_resolution_edge(&self.curve1, self.res_coeff1, self.tol);
        let a_res2 = curve_resolution_edge(&self.curve2, self.res_coeff2, self.tol);
        let d_tr1 = 20.0 * a_res1;
        let d_tr2 = 20.0 * a_res2;

        let mut taken = vec![false; a_nb_cp];
        for i in 0..a_nb_cp {
            if taken[i] { continue; }
            let mut ti1 = ranges1[i][0]; let mut ti2 = ranges1[i][1];
            let mut tj1 = ranges2[i][0]; let mut tj2 = ranges2[i][1];
            taken[i] = true;
            for j in (i + 1)..a_nb_cp {
                if taken[j] { continue; }
                let t_j1 = ranges1[j][0]; let t_j2 = ranges1[j][1];
                let t_jj1 = ranges2[j][0]; let t_jj2 = ranges2[j][1];
                let b_cond = (ti2 - t_j1).abs() < d_tr1
                    || (t_j1 > ti1 && t_j1 < ti2)
                    || (ti1 > t_j1 && ti1 < t_j2)
                    || (b_split2 && (t_j2 - ti1).abs() < d_tr1);
                if b_cond && b_split2 {
                    // rcad: skip the complex second condition (approximated)
                }
                if b_cond {
                    ti1 = ti1.min(t_j1);
                    ti2 = ti2.max(t_j2);
                    tj1 = tj1.min(t_jj1);
                    tj2 = tj2.max(t_jj2);
                    taken[j] = true;
                } else if !b_split2 {
                    break;
                }
            }
            // OCCT L758-766: check if whole range
            let whole1 = (t11 - ti1).abs() < self.res1 && (t12 - ti2).abs() < self.res1;
            let whole2 = (t21 - tj1).abs() < self.res2 && (t22 - tj2).abs() < self.res2;
            let is_edge = whole1 || whole2;
            if is_edge {
                self.common_parts.clear();
                self.add_solution(ti1, ti2, tj1, tj2, true);
                return;
            }
            self.add_solution(ti1, ti2, tj1, tj2, false);
        }
    }

    /// OCCT L780-822: AddSolution
    fn add_solution(&mut self, t11: f64, t12: f64, t21: f64, t22: f64, is_edge: bool) {
        let p1 = self.curve1.point_at(t11);
        let p2 = self.curve2.point_at(t22);
        let cp = if is_edge {
            CommonPrt::new_edge(
                if self.swapped { [t21, t22] } else { [t11, t12] },
                vec![if self.swapped { [t11, t12] } else { [t21, t22] }],
                p1, p2,
            )
        } else {
            let t1 = if self.swapped { t21 } else { t11 };
            let t2 = if self.swapped { t11 } else { t21 };
            CommonPrt::new_vertex(t1, t2, p1)
        };
        self.common_parts.push(cp);
    }

    /// OCCT L1150-1206: CheckCoincidence �?checks whether refined ranges
    /// from `find_solutions_rec` are truly coincident, using golden-section
    /// search to find the max distance between two curves on the interval.
    fn check_coincidence(&self, t11: f64, t12: f64, t21: f64, t22: f64,
        criteria: f64, curve_res1: f64) -> i32 {
        // Step 1: quick rejection �?project 10 sample points from edge1 to edge2
        let a_nb = 10usize;
        let dt1 = (t12 - t11) / a_nb as f64;
        let mut t = t11;
        for _ in 0..a_nb {
            t += dt1;
            if t > t12 { break; }
            let p = self.curve1.point_at(t);
            let r = rcad_kernel::closest_point_on_curve(&self.curve2, p, 16);
            if r.distance > criteria * 100.0 { return 2; }
        }

        // Step 2: golden-section search for max distance on [t11, t12]
        // OCCT uses FindDistPC which finds the point on edge2 closest to each
        // sample on edge1, then does golden-section to find the MAX of these
        // distances over the range �?a one-dimensional min-max optimization.
        const CF: f64 = 0.6180339887498948482045868343656;
        let mut a = t11;
        let mut b = t12;
        let fa = |t: f64| -> f64 {
            let p = self.curve1.point_at(t);
            let r = rcad_kernel::closest_point_on_curve(&self.curve2, p, 16);
            -r.distance // negate because we're minimizing -distance = maximizing distance
        };
        // Golden-section: find the point that MAXIMIZES distance (minimizes -distance)
        let mut x1 = b - CF * (b - a);
        let mut x2 = a + CF * (b - a);
        let mut f1 = fa(x1);
        let mut f2 = fa(x2);
        let eps = curve_res1 * 0.1;
        while (b - a).abs() > eps {
            if f1 < f2 {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = b - CF * (b - a);
                f1 = fa(x1);
            } else {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = a + CF * (b - a);
                f2 = fa(x2);
            }
        }
        let max_dist = -f1; // re-negate to get actual distance

        if max_dist <= criteria { 0 } else { 1 }
    }

    /// OCCT L826-898: FindBestSolution
    #[allow(unused)]
    fn find_best_solution(&self, t11: f64, t12: f64, t21: f64, t22: f64) -> (f64, f64) {
        let a_nb = 10usize;
        let mut a_d_min = f64::MAX;
        let mut a_t1 = 0.0;
        let mut a_t2 = 0.0;
        let dt = (t12 - t11) / a_nb as f64;
        let mut t = t11;
        for _ in 0..=a_nb {
            let p = self.curve1.point_at(t);
            let proj = rcad_kernel::closest_point_on_curve(&self.curve2, p, 16);
            if proj.distance < a_d_min {
                a_d_min = proj.distance;
                a_t1 = t;
                a_t2 = proj.param;
            }
            t += dt;
        }
        (a_t1, a_t2)
    }

    /// OCCT L1060-1146: IsIntersection
    fn is_intersection(&self, t11: f64, t12: f64, t21: f64, t22: f64) -> bool {
        let p11 = self.curve1.point_at(t11);
        let p12 = self.curve1.point_at(t12);
        let p21 = self.curve2.point_at(t21);
        let p22 = self.curve2.point_at(t22);
        let criteria = 100.0 * self.tol;
        let criteria2 = criteria * criteria;
        let d11_21 = (p11 - p21).length_squared();
        let d11_22 = (p11 - p22).length_squared();
        let d12_21 = (p12 - p21).length_squared();
        let d12_22 = (p12 - p22).length_squared();
        let b11_21 = d11_21 < criteria2;
        let b11_22 = d11_22 < criteria2;
        let b12_21 = d12_21 < criteria2;
        let b12_22 = d12_22 < criteria2;
        if (b11_21 && b12_22) || (b11_22 && b12_21) {
            // Check tangent angles
            let v11 = self.curve1.tangent_at(t11);
            let v12 = self.curve1.tangent_at(t12);
            let v21 = self.curve2.tangent_at(t21);
            let v22 = self.curve2.tangent_at(t22);
            if v11.length_squared() > TOLERANCE_CLAMP_MIN && v12.length_squared() > TOLERANCE_CLAMP_MIN
                && v21.length_squared() > TOLERANCE_CLAMP_MIN && v22.length_squared() > TOLERANCE_CLAMP_MIN {
                let a1 = if b11_21 && b12_22 {
                    (v11.normalize().dot(v21.normalize())).acos()
                } else {
                    (v11.normalize().dot(v22.normalize())).acos()
                };
                let a2 = if b11_21 && b12_22 {
                    (v12.normalize().dot(v22.normalize())).acos()
                } else {
                    (v12.normalize().dot(v21.normalize())).acos()
                };
                let ang_criteria = 5e-3;
                if a1 < ang_criteria || (std::f64::consts::PI - a1).abs() < ang_criteria
                    || a2 < ang_criteria || (std::f64::consts::PI - a2).abs() < ang_criteria {
                    // Check min distance
                    let tm = (t11 + t12) * 0.5;
                    let pm = self.curve1.point_at(tm);
                    let r = rcad_kernel::closest_point_on_curve(&self.curve2, pm, 16);
                    return r.distance <= self.tol;
                }
            }
        }
        true
    }

    // ---- Helper: axis-aligned bounding box of a curve segment ----
    fn curve_aabb(&self, curve: &Curve3, t1: f64, t2: f64, tol: f64) -> (DVec3, DVec3) {
        let n = 23usize;
        let dt = (t2 - t1) / n as f64;
        let mut min_p = DVec3::splat(f64::MAX);
        let mut max_p = DVec3::splat(f64::NEG_INFINITY);
        for i in 0..=n {
            let t = t1 + i as f64 * dt;
            let p = curve.point_at(t);
            min_p = min_p.min(p);
            max_p = max_p.max(p);
        }
        (min_p - DVec3::splat(tol), max_p + DVec3::splat(tol))
    }
}

impl Default for EdgeEdgeIntersector {
    fn default() -> Self { Self::new() }
}

/// Point-in-axis-aligned-bounding-box check.
fn point_in_aabb(p: DVec3, box_min: DVec3, box_max: DVec3) -> bool {
    p.x >= box_min.x && p.x <= box_max.x
        && p.y >= box_min.y && p.y <= box_max.y
        && p.z >= box_min.z && p.z <= box_max.z
}

/// Check if two AABBs overlap.
fn aabb_overlap(a: (DVec3, DVec3), b: (DVec3, DVec3)) -> bool {
    a.0.x <= b.1.x && a.1.x >= b.0.x
        && a.0.y <= b.1.y && a.1.y >= b.0.y
        && a.0.z <= b.1.z && a.1.z >= b.0.z
}

/// Square extent of an AABB (max side length squared).
fn aabb_sq_extent(b: (DVec3, DVec3)) -> f64 {
    let d = b.1 - b.0;
    d.length_squared()
}
