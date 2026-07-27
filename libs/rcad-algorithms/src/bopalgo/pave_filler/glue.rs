use crate::bopalgo::pave_filler::helpers::*;
use crate::bopds::ds::NearTangentType;
use crate::bopds::ds::{DSEdge, ShapeOrigin};
use crate::tolerance::*;
use glam::DVec3;
use rcad_kernel::geom::*;

impl<'a> super::PaveFiller<'a> {
    /// OCCT BOPTools_AlgoTools: check surface compatibility for glue
    pub(crate) fn surfaces_glue_compatible(&self, s1: &Surface3, s2: &Surface3) -> bool {
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

    /// OCCT BOPTools_AlgoTools: full boundary overlap check
    pub(crate) fn boundaries_fully_overlap(&self, f1: usize, f2: usize) -> bool {
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

    /// OCCT: find shared edges between faces
    pub(crate) fn detect_shared_edges_between_faces(
        &self,
        f1: usize,
        f2: usize,
    ) -> Vec<(usize, usize)> {
        let tol = self.glue_tolerance;
        let mut shared_edges = Vec::new();

        let edges1: Vec<usize> = self.ds.face_boundary_edges(f1).to_vec();
        let edges2: Vec<usize> = self.ds.face_boundary_edges(f2).to_vec();

        for &e1 in &edges1 {
            for &e2 in &edges2 {
                // Use the new edge overlap detection
                if let Some(overlap) = self.detect_edge_overlap(e1, e2, tol) {
                    // Only consider edges that have at least partial overlap
                    if overlap.overlap_type != EdgeOverlapType::None
                        && overlap.overlap_ratio_a > 0.01
                        && overlap.max_distance < tol * 10.0
                    {
                        shared_edges.push((e1, e2));
                        break; // Each edge in f1 matches at most one in f2
                    }
                }
            }
        }

        shared_edges
    }

    /// OCCT BOPTools_AlgoTools: curve compatibility
    pub(crate) fn edges_curve_compatible(&self, e1: usize, e2: usize, tol: f64) -> bool {
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

    /// OCCT: check partial glue overlap
    pub(crate) fn has_partial_glue(&self, f1: usize, f2: usize) -> bool {
        if !self.use_glue() {
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

    /// OCCT: detect partial overlaps
    pub(crate) fn detect_partial_glue_overlaps(&self) -> Vec<PartialOverlapInfo> {
        let mut overlaps = Vec::new();

        // Iterate over all face pairs from different shapes
        let a_fcount = self.ds.a_face_count;
        let mut pit = crate::bopds::ds::PairIterator::prepare_ab(a_fcount, self.ds.face_count());
        while pit.more() {
            let pk = pit.value();
            let f1_idx = pk.i1;
            let f2_idx = pk.i2;
            let tol = self.ff_tol(f1_idx, f2_idx);
            if let Some(overlap) = self.check_partial_overlap(f1_idx, f2_idx, tol) {
                overlaps.push(overlap);
            }
            pit.next();
        }

        overlaps
    }

    /// OCCT: check partial overlap between faces
    pub(crate) fn check_partial_overlap(
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

        // Check for edge overlap between faces
        let shared_edges = self.detect_shared_edges_between_faces(f1_idx, f2_idx);
        let has_edge_overlap = !shared_edges.is_empty();

        // Check for edge containment
        let mut has_containment = false;
        for &(e1, e2) in &shared_edges {
            if let Some(containment) = self.detect_edge_containment(e1, e2, tol)
                && containment.is_exact
            {
                has_containment = true;
                break;
            }
        }

        // Determine overlap type
        let overlap_type = if has_containment {
            PartialOverlapType::Contained
        } else if has_edge_overlap {
            PartialOverlapType::EdgeOverlap
        } else {
            PartialOverlapType::CoplanarBoundary
        };

        // Partial overlap: some but not complete
        if overlap_ratio > 0.1 && overlap_ratio < 0.99 {
            return Some(PartialOverlapInfo {
                face_a: f1_idx,
                face_b: f2_idx,
                overlap_ratio,
                overlap_type,
            });
        }

        None
    }

    /// OCCT: compute boundary overlap ratio
    pub(crate) fn compute_boundary_overlap_ratio(
        &self,
        pts1: &[DVec3],
        pts2: &[DVec3],
        tol: f64,
    ) -> f64 {
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

    /// OCCT PaveFiller_11: detect edge overlaps
    /// OCCT PaveFiller_11: detect edge overlap between two edges
    pub(crate) fn detect_edge_overlaps(&self) -> Vec<EdgeOverlapResult> {
        let mut overlaps = Vec::new();

        // Iterate over all edge pairs from different shapes
        let a_ecount = self.ds.a_edge_count;
        let mut eit = crate::bopds::ds::PairIterator::prepare_ab(a_ecount, self.ds.edge_count());
        while eit.more() {
            let pk = eit.value();
            let e1_idx = pk.i1;
            let e2_idx = pk.i2;
            let tol = self.ee_tol(e1_idx, e2_idx);
            if let Some(overlap) = self.detect_edge_overlap(e1_idx, e2_idx, tol)
                && overlap.overlap_type != EdgeOverlapType::None
            {
                overlaps.push(overlap);
            }
            eit.next();
        }

        overlaps
    }

    /// OCCT PaveFiller_11: detect edge overlap between two edges
    pub(crate) fn detect_edge_overlap(
        &self,
        e1_idx: usize,
        e2_idx: usize,
        tol: f64,
    ) -> Option<EdgeOverlapResult> {
        let edge1 = self.ds.edges.get(e1_idx)?;
        let edge2 = self.ds.edges.get(e2_idx)?;

        // First check if the curves are compatible (same supporting curve)
        let curve_match = self.curves_are_collinear(&edge1.curve, &edge2.curve, tol);
        if !curve_match {
            return Some(EdgeOverlapResult {
                edge_a: e1_idx,
                edge_b: e2_idx,
                overlap_type: EdgeOverlapType::None,
                overlap_ratio_a: 0.0,
                overlap_ratio_b: 0.0,
                param_range_a: None,
                param_range_b: None,
                max_distance: f64::INFINITY,
            });
        }

        // Compute parameter range overlap in a common parameter space
        let param_overlap = self.compute_param_overlap_for_edges(edge1, edge2, tol);

        // Sample points to compute max distance in overlap region
        let max_distance = if param_overlap.overlap_range.is_some() {
            self.compute_max_edge_distance_in_range(edge1, edge2, &param_overlap, tol)
        } else {
            f64::INFINITY
        };

        let overlap_type = match param_overlap.overlap_type {
            ParamOverlapType::None => EdgeOverlapType::None,
            ParamOverlapType::Partial => EdgeOverlapType::Partial,
            ParamOverlapType::AContainsB => EdgeOverlapType::BContainedInA,
            ParamOverlapType::BContainsA => EdgeOverlapType::AContainedInB,
            ParamOverlapType::Exact => EdgeOverlapType::Full,
        };

        Some(EdgeOverlapResult {
            edge_a: e1_idx,
            edge_b: e2_idx,
            overlap_type,
            overlap_ratio_a: param_overlap.ratio_a,
            overlap_ratio_b: param_overlap.ratio_b,
            param_range_a: param_overlap.overlap_range,
            param_range_b: param_overlap.overlap_range,
            max_distance,
        })
    }

    /// OCCT BOPTools_AlgoTools: curve collinearity check
    pub(crate) fn curves_are_collinear(&self, c1: &Curve3, c2: &Curve3, tol: f64) -> bool {
        match (c1, c2) {
            (Curve3::Line(l1), Curve3::Line(l2)) => self.lines_are_collinear(l1, l2, tol),
            (Curve3::Circle(c1), Curve3::Circle(c2)) => self.circles_are_collinear(c1, c2, tol),
            (Curve3::Ellipse(e1), Curve3::Ellipse(e2)) => self.ellipses_are_collinear(e1, e2, tol),
            (Curve3::BSpline(b1), Curve3::BSpline(b2)) => self.bsplines_are_collinear(b1, b2, tol),
            (Curve3::Bezier(b1), Curve3::Bezier(b2)) => self.beziers_are_collinear(b1, b2, tol),
            // Mixed types could potentially represent the same curve
            // For simplicity, we return false for mixed types
            _ => false,
        }
    }

    /// OCCT BOPTools_AlgoTools: line collinearity
    pub(crate) fn lines_are_collinear(&self, l1: &Line3, l2: &Line3, tol: f64) -> bool {
        let d1 = l1.direction.normalize_or_zero();
        let d2 = l2.direction.normalize_or_zero();

        // Check if directions are parallel (or anti-parallel)
        let dot = d1.dot(d2);
        if dot.abs() < 0.999999 {
            return false;
        }

        // Check if origins are on the same line
        // l2.origin should lie on l1's line
        let v = l2.origin - l1.origin;
        let perp = v - d1 * v.dot(d1);
        perp.length() <= tol * 2.0
    }

    /// OCCT BOPTools_AlgoTools: circle collinearity
    pub(crate) fn circles_are_collinear(&self, c1: &Circle3, c2: &Circle3, tol: f64) -> bool {
        // Centers must be the same
        let center_dist = (c1.center - c2.center).length();
        if center_dist > tol {
            return false;
        }

        // Normals must be parallel (or anti-parallel)
        let normal_dot = c1
            .normal
            .normalize_or_zero()
            .dot(c2.normal.normalize_or_zero());
        if normal_dot.abs() < 0.999999 {
            return false;
        }

        // Radii must be equal
        (c1.radius - c2.radius).abs() <= tol
    }

    /// OCCT BOPTools_AlgoTools: ellipse collinearity
    pub(crate) fn ellipses_are_collinear(&self, e1: &Ellipse3, e2: &Ellipse3, tol: f64) -> bool {
        // Centers must be the same
        let center_dist = (e1.center - e2.center).length();
        if center_dist > tol {
            return false;
        }

        // Normals must be parallel
        let normal_dot = e1
            .normal
            .normalize_or_zero()
            .dot(e2.normal.normalize_or_zero());
        if normal_dot.abs() < 0.999999 {
            return false;
        }

        // Major directions must be parallel (or anti-parallel if normal is flipped)
        let major_dot = e1
            .major_dir
            .normalize_or_zero()
            .dot(e2.major_dir.normalize_or_zero());
        if major_dot.abs() < 0.999999 {
            return false;
        }

        // Radii must be equal
        (e1.major_radius - e2.major_radius).abs() <= tol
            && (e1.minor_radius - e2.minor_radius).abs() <= tol
    }

    /// OCCT BOPTools_AlgoTools: BSpline collinearity
    pub(crate) fn bsplines_are_collinear(
        &self,
        b1: &BSplineCurve3,
        b2: &BSplineCurve3,
        tol: f64,
    ) -> bool {
        // Degrees must match
        if b1.degree != b2.degree {
            return false;
        }

        // Knot vectors should have similar structure
        if b1.knots.len() != b2.knots.len() {
            return false;
        }

        // Control points should match (allowing for reparameterization)
        if b1.control_points.len() != b2.control_points.len() {
            return false;
        }

        // Compare control points with tolerance
        for (p1, p2) in b1.control_points.iter().zip(b2.control_points.iter()) {
            if (*p1 - *p2).length() > tol {
                return false;
            }
        }

        // Compare weights if rational
        for (w1, w2) in b1.weights.iter().zip(b2.weights.iter()) {
            if (w1 - w2).abs() > tol {
                return false;
            }
        }

        true
    }

    /// OCCT BOPTools_AlgoTools: Bezier collinearity
    pub(crate) fn beziers_are_collinear(
        &self,
        b1: &BezierCurve3,
        b2: &BezierCurve3,
        tol: f64,
    ) -> bool {
        // Control point counts must match
        if b1.control_points.len() != b2.control_points.len() {
            return false;
        }

        // Compare control points
        for (p1, p2) in b1.control_points.iter().zip(b2.control_points.iter()) {
            if (*p1 - *p2).length() > tol {
                return false;
            }
        }

        // Compare weights
        for (w1, w2) in b1.weights.iter().zip(b2.weights.iter()) {
            if (w1 - w2).abs() > tol {
                return false;
            }
        }

        true
    }

    /// OCCT: parameter overlap computation
    pub(crate) fn compute_param_overlap_for_edges(
        &self,
        edge1: &DSEdge,
        edge2: &DSEdge,
        tol: f64,
    ) -> ParamOverlap {
        // For collinear edges, we need to map both parameter ranges to a common space
        // The approach depends on the curve type

        match (&edge1.curve, &edge2.curve) {
            (Curve3::Line(l1), Curve3::Line(l2)) => {
                self.compute_line_param_overlap(l1, edge1.t_range, l2, edge2.t_range, tol)
            }
            (Curve3::Circle(c1), Curve3::Circle(c2)) => {
                self.compute_circle_param_overlap(c1, edge1.t_range, c2, edge2.t_range, tol)
            }
            (Curve3::Ellipse(e1), Curve3::Ellipse(e2)) => {
                self.compute_ellipse_param_overlap(e1, edge1.t_range, e2, edge2.t_range, tol)
            }
            (Curve3::BSpline(b1), Curve3::BSpline(b2)) => {
                self.compute_bspline_param_overlap(b1, edge1.t_range, b2, edge2.t_range, tol)
            }
            (Curve3::Bezier(b1), Curve3::Bezier(b2)) => {
                self.compute_bezier_param_overlap(b1, edge1.t_range, b2, edge2.t_range, tol)
            }
            _ => ParamOverlap {
                overlap_type: ParamOverlapType::None,
                overlap_range: None,
                ratio_a: 0.0,
                ratio_b: 0.0,
            },
        }
    }

    /// OCCT BOPTools_AlgoTools: line param overlap
    pub(crate) fn compute_line_param_overlap(
        &self,
        l1: &Line3,
        range1: [f64; 2],
        l2: &Line3,
        range2: [f64; 2],
        tol: f64,
    ) -> ParamOverlap {
        let d1 = l1.direction.normalize_or_zero();
        let d2 = l2.direction.normalize_or_zero();

        // Determine if directions are same or opposite
        let dot = d1.dot(d2);
        let same_direction = dot >= 0.0;

        // Project l2's origin onto l1's parameter space
        // l1: P(t) = l1.origin + t * d1
        // For point p on l2 at parameter s: p = l2.origin + s * d2
        // We need to find t such that: l1.origin + t * d1 = l2.origin + s * d2
        // t = (l2.origin - l1.origin) . d1 + s * (d2 . d1)
        // Since d2 . d1 =  ? (same or opposite direction), we have:
        // t = offset + s * sign

        let offset = (l2.origin - l1.origin).dot(d1);
        let sign = if same_direction { 1.0 } else { -1.0 };

        // Convert range2 to l1's parameter space
        let range2_on_1 = if same_direction {
            [offset + range2[0] * sign, offset + range2[1] * sign]
        } else {
            // Reverse the range when direction is opposite
            [offset + range2[1] * sign, offset + range2[0] * sign]
        };

        // Now compute overlap between range1 and range2_on_1
        self.compute_interval_overlap(range1, range2_on_1, tol)
    }

    /// OCCT: circle param overlap
    pub(crate) fn compute_circle_param_overlap(
        &self,
        c1: &Circle3,
        range1: [f64; 2],
        c2: &Circle3,
        range2: [f64; 2],
        tol: f64,
    ) -> ParamOverlap {
        // For circles, parameters are angles [0, 2
        // Since we already verified circles are the same, we just compare angle ranges
        // But we need to handle periodicity

        let period = 2.0 * std::f64::consts::PI;

        // Check if circles have the same orientation
        let normal_dot = c1
            .normal
            .normalize_or_zero()
            .dot(c2.normal.normalize_or_zero());
        let same_orientation = normal_dot >= 0.0;

        // Normalize ranges to [0, 2
        let r1 = self.normalize_angle_range(range1, period);
        let r2 = self.normalize_angle_range(range2, period);

        // Handle periodic overlap
        if same_orientation {
            self.compute_periodic_interval_overlap(r1, r2, period, tol)
        } else {
            // Flip the range for opposite orientation
            let r2_flipped = [period - r2[1], period - r2[0]];
            self.compute_periodic_interval_overlap(r1, r2_flipped, period, tol)
        }
    }

    /// OCCT: ellipse param overlap
    pub(crate) fn compute_ellipse_param_overlap(
        &self,
        e1: &Ellipse3,
        range1: [f64; 2],
        e2: &Ellipse3,
        range2: [f64; 2],
        tol: f64,
    ) -> ParamOverlap {
        let period = 2.0 * std::f64::consts::PI;

        // Check if ellipses have the same orientation
        let normal_dot = e1
            .normal
            .normalize_or_zero()
            .dot(e2.normal.normalize_or_zero());
        let same_orientation = normal_dot >= 0.0;

        let r1 = self.normalize_angle_range(range1, period);
        let r2 = self.normalize_angle_range(range2, period);

        if same_orientation {
            self.compute_periodic_interval_overlap(r1, r2, period, tol)
        } else {
            let r2_flipped = [period - r2[1], period - r2[0]];
            self.compute_periodic_interval_overlap(r1, r2_flipped, period, tol)
        }
    }

    /// OCCT: BSpline param overlap
    pub(crate) fn compute_bspline_param_overlap(
        &self,
        _b1: &BSplineCurve3,
        range1: [f64; 2],
        _b2: &BSplineCurve3,
        range2: [f64; 2],
        tol: f64,
    ) -> ParamOverlap {
        // For BSplines that have been verified as collinear,
        // we assume the same parameterization and compare ranges directly
        self.compute_interval_overlap(range1, range2, tol)
    }

    /// OCCT: Bezier param overlap
    pub(crate) fn compute_bezier_param_overlap(
        &self,
        _b1: &BezierCurve3,
        range1: [f64; 2],
        _b2: &BezierCurve3,
        range2: [f64; 2],
        tol: f64,
    ) -> ParamOverlap {
        // Bezier curves have domain [0, 1]
        self.compute_interval_overlap(range1, range2, tol)
    }

    /// OCCT BOPTools_AlgoTools: interval overlap
    pub(crate) fn compute_interval_overlap(
        &self,
        a: [f64; 2],
        b: [f64; 2],
        tol: f64,
    ) -> ParamOverlap {
        let a_len = (a[1] - a[0]).abs();
        let b_len = (b[1] - b[0]).abs();

        if a_len < tol || b_len < tol {
            // Degenerate interval
            return ParamOverlap {
                overlap_type: ParamOverlapType::None,
                overlap_range: None,
                ratio_a: 0.0,
                ratio_b: 0.0,
            };
        }

        // Compute overlap range
        let overlap_start = a[0].max(b[0]);
        let overlap_end = a[1].min(b[1]);

        if overlap_start >= overlap_end - tol {
            // No overlap
            return ParamOverlap {
                overlap_type: ParamOverlapType::None,
                overlap_range: None,
                ratio_a: 0.0,
                ratio_b: 0.0,
            };
        }

        let overlap_len = overlap_end - overlap_start;
        let ratio_a = overlap_len / a_len;
        let ratio_b = overlap_len / b_len;

        // Determine overlap type
        let overlap_type = if ratio_a >= 0.999999 && ratio_b >= 0.999999 {
            ParamOverlapType::Exact
        } else if ratio_a >= 0.999999 {
            ParamOverlapType::BContainsA
        } else if ratio_b >= 0.999999 {
            ParamOverlapType::AContainsB
        } else {
            ParamOverlapType::Partial
        };

        ParamOverlap {
            overlap_type,
            overlap_range: Some([overlap_start, overlap_end]),
            ratio_a,
            ratio_b,
        }
    }

    /// OCCT: periodic interval overlap
    pub(crate) fn compute_periodic_interval_overlap(
        &self,
        a: [f64; 2],
        b: [f64; 2],
        period: f64,
        tol: f64,
    ) -> ParamOverlap {
        // Handle wraparound for interval a
        let a_wraps = a[1] > a[0] + period / 2.0 || a[1] < a[0];
        let b_wraps = b[1] > b[0] + period / 2.0 || b[1] < b[0];

        // Simple case: neither wraps
        if !a_wraps && !b_wraps {
            return self.compute_interval_overlap(a, b, tol);
        }

        // For wrapping intervals, we need to handle periodicity
        // Unwrap both intervals to a continuous representation
        let a_unwrapped = if a_wraps {
            vec![[a[0], period], [0.0, a[1]]]
        } else {
            vec![a]
        };

        let b_unwrapped = if b_wraps {
            vec![[b[0], period], [0.0, b[1]]]
        } else {
            vec![b]
        };

        // Compute overlap for each combination
        let mut total_overlap_len = 0.0;
        let mut overlap_ranges = Vec::new();

        for a_seg in &a_unwrapped {
            for b_seg in &b_unwrapped {
                let overlap = self.compute_interval_overlap(*a_seg, *b_seg, tol);
                if let Some(range) = overlap.overlap_range {
                    total_overlap_len += range[1] - range[0];
                    overlap_ranges.push(range);
                }
            }
        }

        let a_len = a_unwrapped.iter().map(|s| s[1] - s[0]).sum::<f64>();
        let b_len = b_unwrapped.iter().map(|s| s[1] - s[0]).sum::<f64>();

        if total_overlap_len < tol {
            return ParamOverlap {
                overlap_type: ParamOverlapType::None,
                overlap_range: None,
                ratio_a: 0.0,
                ratio_b: 0.0,
            };
        }

        let ratio_a = total_overlap_len / a_len;
        let ratio_b = total_overlap_len / b_len;

        let overlap_type = if ratio_a >= 0.999999 && ratio_b >= 0.999999 {
            ParamOverlapType::Exact
        } else if ratio_a >= 0.999999 {
            ParamOverlapType::BContainsA
        } else if ratio_b >= 0.999999 {
            ParamOverlapType::AContainsB
        } else {
            ParamOverlapType::Partial
        };

        // Return the first overlap range (simplified for periodic case)
        ParamOverlap {
            overlap_type,
            overlap_range: overlap_ranges.first().copied(),
            ratio_a,
            ratio_b,
        }
    }

    /// OCCT: normalize angle to [0, period)
    pub(crate) fn normalize_angle_range(&self, range: [f64; 2], period: f64) -> [f64; 2] {
        let mut r1 = range[0] % period;
        let mut r2 = range[1] % period;

        if r1 < 0.0 {
            r1 += period;
        }
        if r2 < 0.0 {
            r2 += period;
        }

        [r1, r2]
    }

    /// OCCT: max edge distance in range
    pub(crate) fn compute_max_edge_distance_in_range(
        &self,
        edge1: &DSEdge,
        edge2: &DSEdge,
        param_overlap: &ParamOverlap,
        _tol: f64,
    ) -> f64 {
        let overlap_range = match param_overlap.overlap_range {
            Some(r) => r,
            None => return f64::INFINITY,
        };

        // Sample points in the overlap region
        let num_samples = 10;
        let mut max_dist = 0.0_f64;

        for i in 0..=num_samples {
            let t = overlap_range[0]
                + (overlap_range[1] - overlap_range[0]) * i as f64 / num_samples as f64;

            let p1 = edge1.curve.point_at(t);

            // Find corresponding point on edge2
            // For now, use simple distance check
            let t2_start = edge2.t_range[0];
            let t2_end = edge2.t_range[1];

            // Sample edge2 and find closest point
            let mut min_dist = f64::INFINITY;
            for j in 0..=num_samples {
                let t2 = t2_start + (t2_end - t2_start) * j as f64 / num_samples as f64;
                let p2 = edge2.curve.point_at(t2);
                let dist = (p1 - p2).length();
                min_dist = min_dist.min(dist);
            }

            max_dist = max_dist.max(min_dist);
        }

        max_dist
    }

    /// OCCT: detect edge containment
    pub(crate) fn detect_edge_containment(
        &self,
        e1_idx: usize,
        e2_idx: usize,
        tol: f64,
    ) -> Option<EdgeContainmentResult> {
        let overlap = self.detect_edge_overlap(e1_idx, e2_idx, tol)?;

        match overlap.overlap_type {
            EdgeOverlapType::AContainedInB => Some(EdgeContainmentResult {
                contained_edge: e1_idx,
                containing_edge: e2_idx,
                containment_ratio: overlap.overlap_ratio_a,
                is_exact: overlap.overlap_ratio_a >= 0.999999,
            }),
            EdgeOverlapType::BContainedInA => Some(EdgeContainmentResult {
                contained_edge: e2_idx,
                containing_edge: e1_idx,
                containment_ratio: overlap.overlap_ratio_b,
                is_exact: overlap.overlap_ratio_b >= 0.999999,
            }),
            _ => None,
        }
    }

    /// OCCT: detect all edge containments
    pub(crate) fn detect_all_edge_containments(&self) -> Vec<EdgeContainmentResult> {
        let mut containments = Vec::new();

        let a_ecount = self.ds.a_edge_count;
        let mut eit = crate::bopds::ds::PairIterator::prepare_ab(a_ecount, self.ds.edge_count());
        while eit.more() {
            let pk = eit.value();
            let e1_idx = pk.i1;
            let e2_idx = pk.i2;
            let tol = self.ee_tol(e1_idx, e2_idx);
            if let Some(containment) = self.detect_edge_containment(e1_idx, e2_idx, tol) {
                containments.push(containment);
            }
            eit.next();
        }

        containments
    }

    /// OCCT: handle near-tangent faces
    pub(crate) fn handle_near_tangent_faces(&self) -> Vec<NearTangentFaceInfo> {
        let mut tangent_faces = Vec::new();

        // Iterate over all face pairs from different shapes
        let a_fcount = self.ds.a_face_count;
        let mut fit = crate::bopds::ds::PairIterator::prepare_ab(a_fcount, self.ds.face_count());
        while fit.more() {
            let pk = fit.value();
            let f1_idx = pk.i1;
            let f2_idx = pk.i2;
            let tangent_threshold = self.ff_tol(f1_idx, f2_idx) * 100.0;
            if let Some(info) = self.check_near_tangent_faces(f1_idx, f2_idx, tangent_threshold) {
                tangent_faces.push(info);
            }
            fit.next();
        }

        tangent_faces
    }

    /// OCCT: check near-tangent faces
    pub(crate) fn check_near_tangent_faces(
        &self,
        f1_idx: usize,
        f2_idx: usize,
        tangent_threshold: f64,
    ) -> Option<NearTangentFaceInfo> {
        let face1 = &self.ds.faces[f1_idx];
        let face2 = &self.ds.faces[f2_idx];

        // Skip if same origin
        if face1.origin == face2.origin {
            return None;
        }

        // Check for near-tangency based on surface types
        match (&face1.surface, &face2.surface) {
            (Surface3::Plane(p1), Surface3::Plane(p2)) => {
                self.check_plane_plane_tangent(f1_idx, f2_idx, p1, p2, tangent_threshold)
            }
            (Surface3::Plane(pl), Surface3::Cylinder(cyl))
            | (Surface3::Cylinder(cyl), Surface3::Plane(pl)) => {
                self.check_plane_cylinder_tangent(f1_idx, f2_idx, pl, cyl, tangent_threshold)
            }
            (Surface3::Plane(pl), Surface3::Sphere(sph))
            | (Surface3::Sphere(sph), Surface3::Plane(pl)) => {
                self.check_plane_sphere_tangent(f1_idx, f2_idx, pl, sph, tangent_threshold)
            }
            (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => {
                self.check_cylinder_cylinder_tangent(f1_idx, f2_idx, c1, c2, tangent_threshold)
            }
            _ => None, // General case not implemented
        }
    }

    /// OCCT: plane-plane tangent check
    pub(crate) fn check_plane_plane_tangent(
        &self,
        f1_idx: usize,
        f2_idx: usize,
        p1: &Plane,
        p2: &Plane,
        tangent_threshold: f64,
    ) -> Option<NearTangentFaceInfo> {
        // Check if normals are nearly parallel (or anti-parallel)
        let n1 = p1.normal.normalize_or_zero();
        let n2 = p2.normal.normalize_or_zero();
        let dot = n1.dot(n2).abs();

        if dot < 0.9999 {
            return None; // Not nearly parallel
        }

        // Compute distance between planes
        let distance = (p2.origin - p1.origin).dot(n1).abs();

        if distance > tangent_threshold {
            return None; // Too far apart
        }

        // Check if faces overlap in XY projection
        let pts1 = self.ds.face_boundary_points(f1_idx);
        let pts2 = self.ds.face_boundary_points(f2_idx);

        if !self.faces_boundaries_overlap(&pts1, &pts2, tangent_threshold) {
            return None;
        }

        Some(NearTangentFaceInfo {
            face_a: f1_idx,
            face_b: f2_idx,
            distance,
            tangent_type: NearTangentType::PlaneParallel,
            should_merge: distance < tangent_threshold * 0.1,
        })
    }

    /// OCCT: plane-cylinder tangent check
    pub(crate) fn check_plane_cylinder_tangent(
        &self,
        f1_idx: usize,
        f2_idx: usize,
        plane: &Plane,
        cyl: &CylindricalSurface,
        tangent_threshold: f64,
    ) -> Option<NearTangentFaceInfo> {
        // A plane is tangent to a cylinder if:
        // 1. Plane normal is perpendicular to cylinder axis
        // 2. Distance from cylinder axis to plane equals radius

        let axis = cyl.axis.normalize_or_zero();
        let normal = plane.normal.normalize_or_zero();

        // Check perpendicularity
        let axis_normal_dot = axis.dot(normal).abs();
        if axis_normal_dot > 0.01 {
            return None; // Not perpendicular
        }

        // Compute distance from cylinder axis to plane
        let axis_point = cyl.origin;
        let dist_to_plane = (axis_point - plane.origin).dot(normal).abs();
        let radius_dist = (dist_to_plane - cyl.radius).abs();

        if radius_dist > tangent_threshold {
            return None; // Not tangent
        }

        Some(NearTangentFaceInfo {
            face_a: f1_idx,
            face_b: f2_idx,
            distance: radius_dist,
            tangent_type: NearTangentType::CylinderPlane,
            should_merge: radius_dist < tangent_threshold * 0.1,
        })
    }

    /// OCCT: plane-sphere tangent check
    pub(crate) fn check_plane_sphere_tangent(
        &self,
        f1_idx: usize,
        f2_idx: usize,
        plane: &Plane,
        sph: &SphericalSurface,
        tangent_threshold: f64,
    ) -> Option<NearTangentFaceInfo> {
        // A plane is tangent to a sphere if distance from center to plane equals radius
        let normal = plane.normal.normalize_or_zero();
        let dist_to_plane = (sph.center - plane.origin).dot(normal).abs();
        let radius_dist = (dist_to_plane - sph.radius).abs();

        if radius_dist > tangent_threshold {
            return None; // Not tangent
        }

        // Check if tangent point is within face boundaries
        let tangent_point = sph.center - normal * sph.radius * dist_to_plane.signum();
        let pts1 = self.ds.face_boundary_points(f1_idx);
        let pts2 = self.ds.face_boundary_points(f2_idx);

        // Simple bounding box check for tangent point
        if !self.point_near_boundary(&tangent_point, &pts1, tangent_threshold * 10.0)
            && !self.point_near_boundary(&tangent_point, &pts2, tangent_threshold * 10.0)
        {
            return None;
        }

        Some(NearTangentFaceInfo {
            face_a: f1_idx,
            face_b: f2_idx,
            distance: radius_dist,
            tangent_type: NearTangentType::SpherePlane,
            should_merge: radius_dist < tangent_threshold * 0.1,
        })
    }

    /// OCCT: cylinder-cylinder tangent check
    pub(crate) fn check_cylinder_cylinder_tangent(
        &self,
        f1_idx: usize,
        f2_idx: usize,
        c1: &CylindricalSurface,
        c2: &CylindricalSurface,
        tangent_threshold: f64,
    ) -> Option<NearTangentFaceInfo> {
        // Check if cylinders have parallel axes
        let a1 = c1.axis.normalize_or_zero();
        let a2 = c2.axis.normalize_or_zero();

        if a1.dot(a2).abs() < 0.999 {
            return None; // Axes not parallel
        }

        // Compute distance between axes
        let v = c2.origin - c1.origin;
        let perp = v - a1 * v.dot(a1);
        let axis_distance = perp.length();

        // Check if tangent (distance equals sum or difference of radii)
        let dist_to_sum = (axis_distance - (c1.radius + c2.radius)).abs();
        let dist_to_diff = (axis_distance - (c1.radius - c2.radius).abs()).abs();
        let min_dist = dist_to_sum.min(dist_to_diff);

        if min_dist > tangent_threshold {
            return None; // Not tangent
        }

        Some(NearTangentFaceInfo {
            face_a: f1_idx,
            face_b: f2_idx,
            distance: min_dist,
            tangent_type: NearTangentType::CylinderCylinder,
            should_merge: min_dist < tangent_threshold * 0.1,
        })
    }

    /// OCCT: face boundary overlap test
    pub(crate) fn faces_boundaries_overlap(
        &self,
        pts1: &[DVec3],
        pts2: &[DVec3],
        tol: f64,
    ) -> bool {
        if pts1.is_empty() || pts2.is_empty() {
            return false;
        }

        // Simple bounding box overlap check
        let mut min1 = DVec3::splat(f64::INFINITY);
        let mut max1 = DVec3::splat(f64::NEG_INFINITY);
        let mut min2 = DVec3::splat(f64::INFINITY);
        let mut max2 = DVec3::splat(f64::NEG_INFINITY);

        for p in pts1 {
            min1 = min1.min(*p);
            max1 = max1.max(*p);
        }
        for p in pts2 {
            min2 = min2.min(*p);
            max2 = max2.max(*p);
        }

        // Check if bounding boxes overlap in all dimensions
        for i in 0..3 {
            if max1[i] + tol < min2[i] || max2[i] + tol < min1[i] {
                return false;
            }
        }

        true
    }

    /// OCCT: point near boundary test
    pub(crate) fn point_near_boundary(&self, point: &DVec3, boundary: &[DVec3], tol: f64) -> bool {
        // Check bounding box first
        let mut min_pt = DVec3::splat(f64::INFINITY);
        let mut max_pt = DVec3::splat(f64::NEG_INFINITY);
        for p in boundary {
            min_pt = min_pt.min(*p);
            max_pt = max_pt.max(*p);
        }

        for i in 0..3 {
            if point[i] < min_pt[i] - tol || point[i] > max_pt[i] + tol {
                return false;
            }
        }

        true
    }

    /// OCCT: handle near-coincident faces
    pub(crate) fn handle_near_coincident_faces(&self) -> Vec<NearCoincidentFaceInfo> {
        let mut coincident_faces = Vec::new();

        let a_fcount = self.ds.a_face_count;
        let mut fit = crate::bopds::ds::PairIterator::prepare_ab(a_fcount, self.ds.face_count());
        while fit.more() {
            let pk = fit.value();
            let f1_idx = pk.i1;
            let f2_idx = pk.i2;
            let coincident_threshold = self.ff_tol(f1_idx, f2_idx) * 10.0;
            if let Some(info) =
                self.check_near_coincident_faces(f1_idx, f2_idx, coincident_threshold)
            {
                coincident_faces.push(info);
            }
            fit.next();
        }

        coincident_faces
    }

    /// OCCT: check near-coincident faces
    pub(crate) fn check_near_coincident_faces(
        &self,
        f1_idx: usize,
        f2_idx: usize,
        coincident_threshold: f64,
    ) -> Option<NearCoincidentFaceInfo> {
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

        // Get boundary points
        let pts1 = self.ds.face_boundary_points(f1_idx);
        let pts2 = self.ds.face_boundary_points(f2_idx);

        // Sample interior points
        let interior1 = self.sample_face_interior(f1_idx, 4);
        let interior2 = self.sample_face_interior(f2_idx, 4);

        // Check distances
        let mut max_distance = 0.0_f64;
        let mut overlap_count = 0;
        let total_points = interior1.len() + interior2.len();

        if total_points == 0 {
            return None;
        }

        // Check interior points of face1 against face2 surface
        for p in &interior1 {
            let dist = self.point_to_surface_distance(*p, &face2.surface);
            if dist < coincident_threshold {
                overlap_count += 1;
            }
            max_distance = max_distance.max(dist);
        }

        // Check interior points of face2 against face1 surface
        for p in &interior2 {
            let dist = self.point_to_surface_distance(*p, &face1.surface);
            if dist < coincident_threshold {
                overlap_count += 1;
            }
            max_distance = max_distance.max(dist);
        }

        // If most points are within threshold, consider faces coincident
        let overlap_ratio = overlap_count as f64 / total_points as f64;
        if overlap_ratio < 0.5 {
            return None;
        }

        // Compute approximate overlap area
        let overlap_area = self.compute_approximate_overlap_area(&pts1, &pts2);

        Some(NearCoincidentFaceInfo {
            face_a: f1_idx,
            face_b: f2_idx,
            max_distance,
            overlap_area,
            should_merge: max_distance < coincident_threshold * 0.1,
        })
    }

    /// OCCT: sample face interior points
    pub(crate) fn sample_face_interior(
        &self,
        face_idx: usize,
        samples_per_dim: usize,
    ) -> Vec<DVec3> {
        let _face = &self.ds.faces[face_idx];
        let boundary = self.ds.face_boundary_points(face_idx);

        if boundary.len() < 3 {
            return Vec::new();
        }

        // Compute centroid
        let centroid: DVec3 = boundary.iter().sum::<DVec3>() / boundary.len() as f64;

        // Sample points along lines from centroid to boundary midpoints
        let mut interior_points = Vec::new();

        for i in 0..boundary.len() {
            let p1 = boundary[i];
            let p2 = boundary[(i + 1) % boundary.len()];
            let mid = (p1 + p2) * 0.5;

            for j in 1..=samples_per_dim {
                let t = j as f64 / (samples_per_dim + 1) as f64;
                let sample = centroid + (mid - centroid) * t;
                interior_points.push(sample);
            }
        }

        interior_points
    }

    /// OCCT: distance from point to surface
    pub(crate) fn point_to_surface_distance(&self, point: DVec3, surface: &Surface3) -> f64 {
        match surface {
            Surface3::Plane(p) => {
                let normal = p.normal.normalize_or_zero();
                (point - p.origin).dot(normal).abs()
            }
            Surface3::Sphere(s) => {
                let dist_to_center = (point - s.center).length();
                (dist_to_center - s.radius).abs()
            }
            Surface3::Cylinder(c) => {
                let axis = c.axis.normalize_or_zero();
                let v = point - c.origin;
                let axial = v.dot(axis);
                let radial = v - axis * axial;
                (radial.length() - c.radius).abs()
            }
            Surface3::Cone(cone) => {
                // Simplified: distance to cone surface
                let axis = cone.axis_dir();
                let v = point - cone.apex;
                let axial = v.dot(axis);
                let radial = (v - axis * axial).length();
                let expected_radius = axial * cone.half_angle_rad.tan();
                (radial - expected_radius).abs()
            }
            Surface3::Torus(t) => {
                // Simplified: distance to torus surface
                let axis = t.axis.normalize_or_zero();
                let v = point - t.center;
                let axial = v.dot(axis);
                let in_plane = v - axis * axial;
                let in_plane_dist = in_plane.length();
                let tube_center_dist = (in_plane_dist - t.major_radius).abs();
                let tube_dist = (tube_center_dist * tube_center_dist + axial * axial).sqrt();
                (tube_dist - t.minor_radius).abs()
            }
            _ => {
                // For other surfaces, use projection
                let proj = rcad_kernel::projection::closest_point_on_surface(surface, point, 16);
                proj.distance
            }
        }
    }

    /// OCCT: approximate overlap area
    pub(crate) fn compute_approximate_overlap_area(&self, pts1: &[DVec3], pts2: &[DVec3]) -> f64 {
        // Compute area of each face
        let area1 = self.compute_polygon_area(pts1);
        let area2 = self.compute_polygon_area(pts2);

        // Return the smaller area as an approximation of overlap
        area1.min(area2)
    }

    /// OCCT: polygon area using Newell's method
    pub(crate) fn compute_polygon_area(&self, pts: &[DVec3]) -> f64 {
        if pts.len() < 3 {
            return 0.0;
        }

        // Find best-fit plane and compute 2D area
        let centroid: DVec3 = pts.iter().sum::<DVec3>() / pts.len() as f64;

        // Use Newell's method to find normal
        let mut normal = DVec3::ZERO;
        for i in 0..pts.len() {
            let p1 = pts[i];
            let p2 = pts[(i + 1) % pts.len()];
            normal.x += (p1.y - p2.y) * (p1.z + p2.z);
            normal.y += (p1.z - p2.z) * (p1.x + p2.x);
            normal.z += (p1.x - p2.x) * (p1.y + p2.y);
        }
        let normal = normal.normalize_or_zero();

        // Project to 2D and compute area
        let (u_dir, v_dir) = if normal.x.abs() > 0.9 {
            (DVec3::Y, DVec3::Z)
        } else {
            (DVec3::X, DVec3::Y)
        };

        let mut area = 0.0;
        for i in 0..pts.len() {
            let p1 = pts[i] - centroid;
            let p2 = pts[(i + 1) % pts.len()] - centroid;
            let u1 = p1.dot(u_dir);
            let v1 = p1.dot(v_dir);
            let u2 = p2.dot(u_dir);
            let v2 = p2.dot(v_dir);
            area += u1 * v2 - u2 * v1;
        }

        area.abs() * 0.5
    }

    /// OCCT: handle micro gaps between edges
    pub(crate) fn handle_micro_gaps(&self) -> Vec<MicroGapInfo> {
        let mut gaps = Vec::new();

        // Check edge-to-edge gaps
        let a_edges: Vec<usize> = self
            .ds
            .edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.origin == ShapeOrigin::ShapeA)
            .map(|(i, _)| i)
            .collect();

        let b_edges: Vec<usize> = self
            .ds
            .edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.origin == ShapeOrigin::ShapeB)
            .map(|(i, _)| i)
            .collect();

        for &ea in &a_edges {
            for &eb in &b_edges {
                let ee = self.ee_tol(ea, eb);
                let gap_threshold = ee * 1000.0;
                if let Some(gap) = self.check_micro_gap(ea, eb, gap_threshold, ee) {
                    gaps.push(gap);
                }
            }
        }

        gaps
    }

    /// OCCT: check micro gap between two edges
    pub(crate) fn check_micro_gap(
        &self,
        e1: usize,
        e2: usize,
        gap_threshold: f64,
        coincident_tol: f64,
    ) -> Option<MicroGapInfo> {
        let _edge1 = &self.ds.edges[e1];
        let _edge2 = &self.ds.edges[e2];

        // Sample points along both edges
        let pts1 = self.sample_edge_points(e1, 8);
        let pts2 = self.sample_edge_points(e2, 8);

        if pts1.is_empty() || pts2.is_empty() {
            return None;
        }

        // Find minimum distance between edges
        let mut min_gap = f64::INFINITY;
        for p1 in &pts1 {
            for p2 in &pts2 {
                let dist = (*p1 - *p2).length();
                min_gap = min_gap.min(dist);
            }
        }

        // Check if it's a micro-gap (within threshold but not coincident)
        if min_gap <= coincident_tol {
            return None; // Already coincident
        }
        if min_gap > gap_threshold {
            return None; // Too large for micro-gap handling
        }

        // Check if edges are approximately parallel
        let parallel = self.edges_approximately_parallel(e1, e2, 0.1);

        Some(MicroGapInfo {
            edge_a: e1,
            edge_b: e2,
            gap_distance: min_gap,
            can_bridge: min_gap < gap_threshold && parallel,
        })
    }

    /// OCCT: sample points along edge
    pub(crate) fn sample_edge_points(&self, edge_idx: usize, n_samples: usize) -> Vec<DVec3> {
        let edge = &self.ds.edges[edge_idx];
        let [t0, t1] = edge.t_range;

        (0..n_samples)
            .map(|i| {
                let t = t0 + (t1 - t0) * i as f64 / (n_samples - 1).max(1) as f64;
                edge.curve.point_at(t)
            })
            .filter(|p| p.is_finite())
            .collect()
    }

    /// OCCT: check if edges are approximately parallel
    pub(crate) fn edges_approximately_parallel(
        &self,
        e1: usize,
        e2: usize,
        angle_tol: f64,
    ) -> bool {
        let edge1 = &self.ds.edges[e1];
        let edge2 = &self.ds.edges[e2];

        // Get edge directions
        let dir1 = match &edge1.curve {
            Curve3::Line(l) => l.direction.normalize_or_zero(),
            Curve3::Circle(_) | Curve3::Ellipse(_) => {
                // For curved edges, check tangent at midpoint
                let t = (edge1.t_range[0] + edge1.t_range[1]) * 0.5;
                let tangent = edge1.curve.tangent_at(t);
                tangent.normalize_or_zero()
            }
            _ => return false,
        };

        let dir2 = match &edge2.curve {
            Curve3::Line(l) => l.direction.normalize_or_zero(),
            Curve3::Circle(_) | Curve3::Ellipse(_) => {
                let t = (edge2.t_range[0] + edge2.t_range[1]) * 0.5;
                let tangent = edge2.curve.tangent_at(t);
                tangent.normalize_or_zero()
            }
            _ => return false,
        };

        // Check parallelism
        let cross = dir1.cross(dir2);
        let sin_angle = cross.length();

        sin_angle < angle_tol
    }

    /// OCCT: handle coincident edges
    pub(crate) fn handle_coincident_edges(&self) -> Vec<CoincidentEdgeInfo> {
        let mut coincident_edges = Vec::new();

        let a_edges: Vec<usize> = self
            .ds
            .edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.origin == ShapeOrigin::ShapeA)
            .map(|(i, _)| i)
            .collect();

        let b_edges: Vec<usize> = self
            .ds
            .edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.origin == ShapeOrigin::ShapeB)
            .map(|(i, _)| i)
            .collect();

        for &ea in &a_edges {
            for &eb in &b_edges {
                let coincident_threshold = self.ee_tol(ea, eb) * 10.0;
                if let Some(info) = self.check_coincident_edges(ea, eb, coincident_threshold) {
                    coincident_edges.push(info);
                }
            }
        }

        coincident_edges
    }

    /// OCCT: check coincident edges between shapes
    pub(crate) fn check_coincident_edges(
        &self,
        e1: usize,
        e2: usize,
        coincident_threshold: f64,
    ) -> Option<CoincidentEdgeInfo> {
        let edge1 = &self.ds.edges[e1];
        let edge2 = &self.ds.edges[e2];

        // Skip if same origin
        if edge1.origin == edge2.origin {
            return None;
        }

        // Check if curves are compatible
        if !self.edges_curve_compatible(e1, e2, coincident_threshold) {
            return None;
        }

        // Sample points and check distances
        let pts1 = self.sample_edge_points(e1, 16);
        let pts2 = self.sample_edge_points(e2, 16);

        if pts1.is_empty() || pts2.is_empty() {
            return None;
        }

        // Compute maximum distance and overlap ratio
        let mut max_distance = 0.0_f64;
        let mut close_count = 0;

        for p1 in &pts1 {
            let min_dist = pts2
                .iter()
                .map(|p2| (*p1 - *p2).length())
                .fold(f64::INFINITY, f64::min);
            max_distance = max_distance.max(min_dist);
            if min_dist < coincident_threshold {
                close_count += 1;
            }
        }

        if max_distance > coincident_threshold {
            return None;
        }

        let overlap_ratio = close_count as f64 / pts1.len() as f64;

        Some(CoincidentEdgeInfo {
            edge_a: e1,
            edge_b: e2,
            max_distance,
            overlap_ratio,
            should_merge: max_distance < coincident_threshold * 0.1 && overlap_ratio > 0.9,
        })
    }
}
