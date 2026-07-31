// OCCT IntTools_EdgeEdge 1:1 translation.
// OCCT ref: IntTools_EdgeEdge.cxx / .hxx / .lxx
//   Perform       L185-243    Prepare      L246-330
//   FindSolutions L290-349    FindSolutions(rec) L353-546
//   FindParameters L553-671   MergeSolutions L675-776
//   AddSolution   L780-822    FindBestSolution L826-898
//   IsIntersection L1060-1146 CheckCoincidence L1150-1206
//   FindDistPC    L1210-1297  DistPC L1301-1362
//   SplitRangeOnSegments L1366-1406  BndBuildBox L1410-1419
//   PointBoxDistance L1423-1452      TypeToInteger L1456-1482
//   ResolutionCoeff L1486-1557       Resolution L1561-1607
//   CurveDeflection L1611-1638       IsClosed L1642-1659
use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval};
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::curve_bounding_box_range;

/// OCCT IntTools_CommonPrt — a common part of two edges.
#[derive(Debug, Clone)]
pub struct CommonPrt {
    pub is_edge: bool,
    pub range1: [f64; 2],
    pub ranges2: Vec<[f64; 2]>,
    pub vertex_param1: f64,
    pub vertex_param2: f64,
    pub bounding_point1: DVec3,
    pub bounding_point2: DVec3,
}

pub struct EdgeEdgeIntersector {
    // OCCT members (myCurve1/2, myRange1/2, myCommonParts, myDone,
    //              myFuzzyValue, myQuickCoincidenceCheck)
    curve1: Curve3, curve2: Curve3,
    range1: [f64; 2], range2: [f64; 2],
    edge1_tol: f64, edge2_tol: f64,
    common_parts: Vec<CommonPrt>, done: bool,
    fuzzy_value: f64,
    quick_coincidence_check: bool,
    // Prepared members (set in Prepare)
    my_tol1: f64, my_tol2: f64, my_tol: f64,
    my_res1: f64, my_res2: f64,
    my_res_coeff1: f64, my_res_coeff2: f64,
    my_p_tol1: f64, my_p_tol2: f64,
    my_swap: bool,
}

// OCCT TypeToInteger (IntTools_EdgeEdge.cxx L1456-1482)
fn type_to_integer(the_c_type: &Curve3) -> i32 {
    match the_c_type {
        Curve3::Line(_) => 0,
        Curve3::Circle(_) | Curve3::Ellipse(_) => 2,
        Curve3::BSpline(_) | Curve3::Bezier(_) => 3,
        _ => 4,
    }
}

impl EdgeEdgeIntersector {
    pub fn new() -> Self {
        EdgeEdgeIntersector {
            curve1: Curve3::Line(rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X }),
            curve2: Curve3::Line(rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X }),
            range1: [0.0, 1.0], range2: [0.0, 1.0],
            edge1_tol: 1e-7, edge2_tol: 1e-7,
            common_parts: Vec::new(), done: false,
            fuzzy_value: 1e-7,
            quick_coincidence_check: false,
            my_tol1: 0.0, my_tol2: 0.0, my_tol: 0.0,
            my_res1: 0.0, my_res2: 0.0,
            my_res_coeff1: 0.0, my_res_coeff2: 0.0,
            my_p_tol1: 0.0, my_p_tol2: 0.0,
            my_swap: false,
        }
    }
    pub fn use_quick_coincidence_check(&mut self, b: bool) { self.quick_coincidence_check = b; }
    pub fn set_edges(&mut self, ei1: usize, r1: [f64; 2], ei2: usize, r2: [f64; 2], ds: &crate::bop::ds::DS) -> &mut Self {
        if let Some(c) = ds.edge_curve(ei1) { self.curve1 = c.clone(); }
        if let Some(c) = ds.edge_curve(ei2) { self.curve2 = c.clone(); }
        self.edge1_tol = ds.edge_tolerance(ei1);
        self.edge2_tol = ds.edge_tolerance(ei2);
        self.range1 = r1; self.range2 = r2; self
    }
    pub fn set_fuzzy_value(&mut self, f: f64) { self.fuzzy_value = f.max(1e-7); }

    /// OCCT IntTools_EdgeEdge::Perform (L185-243).
    pub fn perform(&mut self) {
        // 1. Check data (rcad: edges/curves are set via set_edges; skip CheckData)
        // 2. Prepare
        self.prepare();
        // 3.1 Check Line/Line case
        if type_to_integer(&self.curve1) == 0 && type_to_integer(&self.curve2) == 0 {
            self.compute_line_line();
            return;
        }
        // 3.2 Quick coincidence check
        if self.quick_coincidence_check {
            if self.is_coincident() {
                let a_t11 = self.range1[0]; let a_t12 = self.range1[1];
                let a_t21 = self.range2[0]; let a_t22 = self.range2[1];
                self.add_solution(a_t11, a_t12, a_t21, a_t22, TopAbs::EDGE);
                return;
            }
        }
        // 3.3 Line + analytical fast rejection (OCCT L217-233, DistShapeShape).
        // rcad: box distance fast rejection for line + conic.
        if (type_to_integer(&self.curve1) == 0 || type_to_integer(&self.curve2) == 0)
            && type_to_integer(&self.curve1) <= 2 && type_to_integer(&self.curve2) <= 2
        {
            let b1 = bnd_build_box(&self.curve1, self.range1[0], self.range1[1], self.my_tol1);
            let b2 = bnd_build_box(&self.curve2, self.range2[0], self.range2[1], self.my_tol2);
            if b1.is_out_box(&b2) {
                return;
            }
        }
        // 4. FindSolutions + MergeSolutions
        let mut a_ranges1 = Vec::new();
        let mut a_ranges2 = Vec::new();
        let b_split2 = self.find_solutions(&mut a_ranges1, &mut a_ranges2);
        self.merge_solutions(&a_ranges1, &a_ranges2, b_split2);
        self.done = true;
    }

    /// OCCT IntTools_EdgeEdge::Prepare (L246-330).
    fn prepare(&mut self) {
        let a_ct1 = type_to_integer(&self.curve1);
        let a_ct2 = type_to_integer(&self.curve2);
        let mut i_ct1 = a_ct1;
        let mut i_ct2 = a_ct2;
        // OCCT L274-284: same type (non-line) — deflection-based ordering
        if i_ct1 == i_ct2 && i_ct1 != 0 {
            let a_c2 = curve_deflection(&self.curve2, self.range2);
            let a_c1 = if a_c2 > 1e-7 { curve_deflection(&self.curve1, self.range1) } else { 1.0 };
            if a_c1 < a_c2 { i_ct1 -= 1; }
        }
        // OCCT L286-299: swap so curve1 has the higher type
        if i_ct1 < i_ct2 {
            std::mem::swap(&mut self.curve1, &mut self.curve2);
            std::mem::swap(&mut self.range1, &mut self.range2);
            std::mem::swap(&mut self.edge1_tol, &mut self.edge2_tol);
            self.my_swap = true;
        }
        // OCCT L301-309: tolerances
        let a_tol_add = self.fuzzy_value / 2.0;
        self.my_tol1 = self.edge1_tol + a_tol_add;
        self.my_tol2 = self.edge2_tol + a_tol_add;
        self.my_tol = self.my_tol1 + self.my_tol2;
        // OCCT L311-322: resolutions
        if i_ct1 != 0 || i_ct2 != 0 {
            self.my_res_coeff1 = resolution_coeff(&self.curve1, self.range1);
            self.my_res_coeff2 = resolution_coeff(&self.curve2, self.range2);
            self.my_res1 = resolution(&self.curve1, self.my_res_coeff1, self.my_tol1);
            self.my_res2 = resolution(&self.curve2, self.my_res_coeff2, self.my_tol2);
            self.my_p_tol1 = 5e-13;
            let a_tm = self.range1[0].abs().max(self.range1[1].abs());
            if a_tm > 999.0 { self.my_p_tol1 = 5e-16 * a_tm; }
            self.my_p_tol2 = 5e-13;
            let a_tm = self.range2[0].abs().max(self.range2[1].abs());
            if a_tm > 999.0 { self.my_p_tol2 = 5e-16 * a_tm; }
        }
    }

    /// OCCT IntTools_EdgeEdge::FindSolutions (L290-349).
    fn find_solutions(&self, the_ranges1: &mut Vec<[f64; 2]>, the_ranges2: &mut Vec<[f64; 2]>) -> bool {
        let (a_t11, a_t12) = (self.range1[0], self.range1[1]);
        let (a_t21, a_t22) = (self.range2[0], self.range2[1]);
        let mut b_is_closed2 = is_closed(&self.curve2, a_t21, a_t22, self.my_tol2, self.my_res2);
        if b_is_closed2 {
            let a_b1 = bnd_build_box(&self.curve1, a_t11, a_t12, self.my_tol1);
            let a_p = self.curve2.point_at(a_t21);
            b_is_closed2 = !point_in_box(&a_b1, a_p);
        }
        if !b_is_closed2 {
            let a_b1 = bnd_build_box(&self.curve1, a_t11, a_t12, self.my_tol1);
            let a_b2 = bnd_build_box(&self.curve2, a_t21, a_t22, self.my_tol2);
            self.find_solutions_rec([a_t11, a_t12], &a_b1, [a_t21, a_t22], &a_b2, the_ranges1, the_ranges2, 0);
            return false;
        }
        // OCCT L322-326: coincident closed curve — whole ranges are a common part
        if check_coincidence(&self.curve1, &self.curve2, a_t11, a_t12, a_t21, a_t22, self.my_tol, self.my_res1) == 0 {
            the_ranges1.push([a_t11, a_t12]);
            the_ranges2.push([a_t21, a_t22]);
            return false;
        }
        // OCCT L328-349: split both ranges and recurse
        let a_nb1 = if is_closed(&self.curve1, a_t11, a_t12, self.my_tol1, self.my_res1) { 2 } else { 1 };
        let a_nb2 = 2;
        let a_segments1 = split_range_on_segments(a_t11, a_t12, self.my_res1, a_nb1);
        let a_segments2 = split_range_on_segments(a_t21, a_t22, self.my_res2, a_nb2);
        let a_nb1 = a_segments1.len();
        let a_nb2 = a_segments2.len();
        let a_b2 = bnd_build_box(&self.curve2, a_t21, a_t22, self.my_tol2);
        for i in 0..a_nb1 {
            let a_r1 = a_segments1[i];
            let a_b1 = bnd_build_box(&self.curve1, a_r1[0], a_r1[1], self.my_tol1);
            for j in 0..a_nb2 {
                let a_r2 = a_segments2[j];
                self.find_solutions_rec(a_r1, &a_b1, a_r2, &a_b2, the_ranges1, the_ranges2, 0);
            }
        }
        a_nb2 > 1
    }

    /// OCCT IntTools_EdgeEdge::FindSolutions (L353-546) — recursive.
    fn find_solutions_rec(
        &self,
        the_r1: [f64; 2], the_box1: &BndBox,
        the_r2: [f64; 2], the_box2: &BndBox,
        the_ranges1: &mut Vec<[f64; 2]>, the_ranges2: &mut Vec<[f64; 2]>,
        _depth: usize,
    ) {
        let (mut a_t11, mut a_t12) = (the_r1[0], the_r1[1]);
        let (mut a_t21, mut a_t22) = (the_r2[0], the_r2[1]);
        let mut a_b1 = the_box1.clone();
        let mut a_b2 = the_box2.clone();
        let mut b_out = false;
        let mut b_stop = false;
        let mut b_thin = false;
        let mut i_com = 1i32;
        // OCCT L377-460: do-while loop
        loop {
            let (a_tb11, a_tb12, a_tb21, a_tb22) = (a_t11, a_t12, a_t21, a_t22);
            // OCCT L385-389
            if a_b1.is_out_box(&a_b2) { b_out = true; break; }
            // OCCT L391-392
            b_thin = (a_t12 - a_t11) < self.my_res1 || box_is_thin(&a_b1, self.my_tol);
            // OCCT L394-407
            if !find_parameters(&self.curve2, a_tb21, a_tb22, self.my_tol2, self.my_res2, self.my_p_tol2, self.my_res_coeff2, &a_b1, &mut a_t21, &mut a_t22) {
                b_out = true; break;
            }
            if b_out || b_thin { break; }
            // OCCT L411-419
            a_b2 = bnd_build_box(&self.curve2, a_t21, a_t22, self.my_tol2);
            if a_b1.is_out_box(&a_b2) { b_out = true; break; }
            b_thin = (a_t22 - a_t21) < self.my_res2 || box_is_thin(&a_b2, self.my_tol);
            // OCCT L421-435
            if !find_parameters(&self.curve1, a_tb11, a_tb12, self.my_tol1, self.my_res1, self.my_p_tol1, self.my_res_coeff1, &a_b2, &mut a_t11, &mut a_t12) {
                b_out = true; break;
            }
            if b_out || b_thin { break; }
            // OCCT L438-458: stop if ranges did not shrink
            let mut a_small_step1 = (a_tb12 - a_tb11) / 250.0;
            let mut a_small_step2 = (a_tb22 - a_tb21) / 250.0;
            if a_small_step1 < self.my_res1 { a_small_step1 = self.my_res1; }
            if a_small_step2 < self.my_res2 { a_small_step2 = self.my_res2; }
            if (a_t11 - a_tb11) < a_small_step1 && (a_tb12 - a_t12) < a_small_step1
                && (a_t21 - a_tb21) < a_small_step2 && (a_tb22 - a_t22) < a_small_step2 {
                b_stop = true;
            } else {
                a_b1 = bnd_build_box(&self.curve1, a_t11, a_t12, self.my_tol1);
            }
            if b_stop { break; }
        }
        if b_out { return; }
        // OCCT L468-476
        if !b_thin {
            i_com = check_coincidence(&self.curve1, &self.curve2, a_t11, a_t12, a_t21, a_t22, self.my_tol, self.my_res1);
            if i_com == 0 { b_thin = true; }
        }
        if b_thin {
            // OCCT L480-513: verify intermediate point
            if i_com != 0 {
                let a_t1 = (a_t11 + a_t12) * 0.5;
                let a_p1 = self.curve1.point_at(a_t1);
                if let Some((a_d, _)) = project_on_range(&self.curve2, a_p1, a_t21, a_t22) {
                    if a_d > self.my_tol { return; }
                } else {
                    let a_t2 = (a_t21 + a_t22) * 0.5;
                    let a_p2 = self.curve2.point_at(a_t2);
                    if a_p1.distance(a_p2) > self.my_tol { return; }
                }
            }
            the_ranges1.push([a_t11, a_t12]);
            the_ranges2.push([a_t21, a_t22]);
            return;
        }
        // OCCT L522-546: split ranges on segments and recurse
        if !is_intersection(&self.curve1, &self.curve2, a_t11, a_t12, a_t21, a_t22, self.my_tol, self.my_res1, self.my_res2) {
            return;
        }
        let a_b1_full = bnd_build_box(&self.curve1, a_t11, a_t12, self.my_tol1);
        let a_b1_sq_extent = box_square_extent(&a_b1_full);
        let a_r2 = [a_t21, a_t22];
        let a_b2 = bnd_build_box(&self.curve2, a_t21, a_t22, self.my_tol2);
        let a_segments1 = split_range_on_segments(a_t11, a_t12, self.my_res1, 3);
        let a_nb1 = a_segments1.len();
        for i in 0..a_nb1 {
            let a_r1 = a_segments1[i];
            let a_b1 = bnd_build_box(&self.curve1, a_r1[0], a_r1[1], self.my_tol1);
            if !a_b1.is_out_box(&a_b2) && (a_nb1 == 1 || box_square_extent(&a_b1) < a_b1_sq_extent) {
                self.find_solutions_rec(a_r1, &a_b1, a_r2, &a_b2, the_ranges1, the_ranges2, _depth + 1);
            }
        }
    }

    /// OCCT IntTools_EdgeEdge::MergeSolutions (L675-776).
    fn merge_solutions(&mut self, the_ranges1: &Vec<[f64; 2]>, the_ranges2: &Vec<[f64; 2]>, b_split2: bool) {
        let a_nb_cp = the_ranges1.len();
        if a_nb_cp == 0 { return; }
        let a_res1 = resolution(&self.curve1, self.my_res_coeff1, self.my_tol);
        let a_res2 = resolution(&self.curve2, self.my_res_coeff2, self.my_tol);
        let (a_t11, a_t12) = (self.range1[0], self.range1[1]);
        let (a_t21, a_t22) = (self.range2[0], self.range2[1]);
        let d_tr1 = 20.0 * a_res1;
        let d_tr2 = 20.0 * a_res2;
        let mut a_type_is_edge = false;
        let mut a_mi: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut i = 1usize;
        while i <= a_nb_cp {
            if a_mi.contains(&i) { i += 1; continue; }
            let (mut a_ti11, mut a_ti12) = (the_ranges1[i - 1][0], the_ranges1[i - 1][1]);
            let (mut a_ti21, mut a_ti22) = (the_ranges2[i - 1][0], the_ranges2[i - 1][1]);
            a_mi.insert(i);
            let mut j = i + 1;
            while j <= a_nb_cp {
                if a_mi.contains(&j) { j += 1; continue; }
                let (a_tj11, a_tj12) = (the_ranges1[j - 1][0], the_ranges1[j - 1][1]);
                let (a_tj21, a_tj22) = (the_ranges2[j - 1][0], the_ranges2[j - 1][1]);
                let mut b_cond = (a_ti12 - a_tj11).abs() < d_tr1
                    || (a_tj11 > a_ti11 && a_tj11 < a_ti12)
                    || (a_ti11 > a_tj11 && a_ti11 < a_tj12)
                    || (b_split2 && (a_tj12 - a_ti11).abs() < d_tr1);
                if b_cond && b_split2 {
                    b_cond = ((a_ti22.max(a_tj22) - a_ti21.min(a_tj21))
                              - ((a_ti22 - a_ti21) + (a_tj22 - a_tj21))).abs() < d_tr2
                        || (a_tj21 > a_ti21 && a_tj21 < a_ti22)
                        || (a_ti21 > a_tj21 && a_ti21 < a_tj22);
                }
                if b_cond {
                    a_ti11 = a_ti11.min(a_tj11);
                    a_ti12 = a_ti12.max(a_tj12);
                    a_ti21 = a_ti21.min(a_tj21);
                    a_ti22 = a_ti22.max(a_tj22);
                    a_mi.insert(j);
                } else if !b_split2 {
                    i = j;
                    break;
                }
                j += 1;
            }
            // OCCT L758-763: EDGE type if the merged range spans the whole edge
            if ((a_t11 - a_ti11).abs() < self.my_res1 && (a_t12 - a_ti12).abs() < self.my_res1)
                || ((a_t21 - a_ti21).abs() < self.my_res2 && (a_t22 - a_ti22).abs() < self.my_res2) {
                a_type_is_edge = true;
                self.common_parts.clear();
            }
            self.add_solution(a_ti11, a_ti12, a_ti21, a_ti22, if a_type_is_edge { TopAbs::EDGE } else { TopAbs::VERTEX });
            if a_type_is_edge { break; }
            if b_split2 { i += 1; }
        }
    }

    /// OCCT IntTools_EdgeEdge::AddSolution (L780-822).
    fn add_solution(&mut self, a_t11: f64, a_t12: f64, a_t21: f64, a_t22: f64, the_type: TopAbs) {
        let is_edge = the_type == TopAbs::EDGE;
        let mut cp = CommonPrt {
            is_edge,
            range1: [a_t11, a_t12],
            ranges2: vec![[a_t21, a_t22]],
            vertex_param1: a_t11,
            vertex_param2: a_t21,
            bounding_point1: self.curve1.point_at((a_t11 + a_t12) * 0.5),
            bounding_point2: self.curve2.point_at((a_t21 + a_t22) * 0.5),
        };
        if !is_edge {
            let (a_t1, a_t2) = find_best_solution(
                &self.curve1, &self.curve2, self.my_res_coeff1, self.my_tol, self.my_p_tol1,
                a_t11, a_t12, a_t21, a_t22);
            cp.vertex_param1 = a_t1;
            cp.vertex_param2 = a_t2;
            cp.range1 = [a_t1, a_t1];
            cp.ranges2 = vec![[a_t2, a_t2]];
            cp.bounding_point1 = self.curve1.point_at(a_t1);
            cp.bounding_point2 = self.curve2.point_at(a_t2);
        }
        if self.my_swap {
            std::mem::swap(&mut cp.vertex_param1, &mut cp.vertex_param2);
            std::mem::swap(&mut cp.bounding_point1, &mut cp.bounding_point2);
        }
        self.common_parts.push(cp);
    }

    /// OCCT IntTools_EdgeEdge::IsCoincident (L247-285): 24pt sampling.
    fn is_coincident(&self) -> bool {
        let a_nb_seg = 23usize;
        let t11 = self.range1[0]; let t12 = self.range1[1];
        let t21 = self.range2[0]; let t22 = self.range2[1];
        let dt = (t12 - t11) / a_nb_seg as f64;
        let mut i_cnt = 0;
        for i in 0..=a_nb_seg {
            let a_t1 = t11 + i as f64 * dt;
            let p1 = self.curve1.point_at(a_t1);
            let mut best_d = f64::MAX;
            for j in 0..=a_nb_seg {
                let a_t2 = t21 + (t22 - t21) * j as f64 / a_nb_seg as f64;
                let d = (p1 - self.curve2.point_at(a_t2)).length();
                if d < best_d { best_d = d; }
            }
            if best_d < self.fuzzy_value { i_cnt += 1; }
        }
        i_cnt as f64 / (a_nb_seg + 1) as f64 > 0.5
    }

    /// OCCT IntTools_EdgeEdge::ComputeLineLine (L902+).
    fn compute_line_line(&mut self) {
        let (l1, l2) = match (&self.curve1, &self.curve2) {
            (Curve3::Line(l1), Curve3::Line(l2)) => (l1, l2),
            _ => return,
        };
        let o1 = l1.origin; let d1 = l1.direction;
        let o2 = l2.origin; let d2 = l2.direction;
        let r = o2 - o1;
        let d1xd2 = d1.cross(d2);
        let denom = d1xd2.length_squared();
        if denom < 1e-30 {
            let dist = r.cross(d1).length() / d1.length().max(1e-30);
            if dist <= self.fuzzy_value {
                let s1 = self.range2[0]; let s2 = self.range2[1];
                let d1_sq = d1.length_squared().max(1e-30);
                let proj_s1 = (r.dot(d1) + s1 * d2.dot(d1)) / d1_sq;
                let proj_s2 = (r.dot(d1) + s2 * d2.dot(d1)) / d1_sq;
                let ov_start = self.range1[0].max(proj_s1.min(proj_s2));
                let ov_end = self.range1[1].min(proj_s1.max(proj_s2));
                if ov_end > ov_start + 1e-12 {
                    self.common_parts.push(CommonPrt {
                        is_edge: true, range1: [ov_start, ov_end],
                        ranges2: vec![[s1, s2]],
                        vertex_param1: ov_start, vertex_param2: s1,
                        bounding_point1: self.curve1.point_at(ov_start),
                        bounding_point2: self.curve2.point_at(s1),
                    });
                }
            }
            return;
        }
        let t = r.cross(d2).dot(d1xd2) / denom;
        let s = r.cross(d1).dot(d1xd2) / denom;
        if t >= self.range1[0] - 1e-12 && t <= self.range1[1] + 1e-12
            && s >= self.range2[0] - 1e-12 && s <= self.range2[1] + 1e-12 {
            let p = o1 + d1 * t;
            self.common_parts.push(CommonPrt {
                is_edge: false, range1: [t, t], ranges2: vec![[s, s]],
                vertex_param1: t, vertex_param2: s,
                bounding_point1: p, bounding_point2: p,
            });
        }
    }

    pub fn is_done(&self) -> bool { self.done }
    pub fn common_parts(&self) -> &[CommonPrt] { &self.common_parts }
}

/// OCCT TopAbs_ShapeEnum (VERTEX/EDGE) — rcad enum for common-part types.
#[derive(Clone, Copy, PartialEq)]
enum TopAbs {
    VERTEX,
    EDGE,
}

// === OCCT IntTools_EdgeEdge static helpers ===

/// OCCT BndBuildBox (L1410-1419): BndLib_Add3dCurve box for a curve range + tol.
fn bnd_build_box(curve: &Curve3, a_t1: f64, a_t2: f64, the_tol: f64) -> BndBox {
    if let Some([mn, mx]) = curve_bounding_box_range(curve, a_t1, a_t2, the_tol) {
        let mut b = BndBox::from_corners(mn.x, mn.y, mn.z, mx.x, mx.y, mx.z);
        b.set_gap(b.get_gap() + the_tol);
        b
    } else {
        BndBox::new()
    }
}

/// OCCT PointBoxDistance (L1423-1452).
fn point_box_distance(a_b: &BndBox, a_p: DVec3) -> f64 {
    if let Some((a_b_min_0, a_b_min_1, a_b_min_2, a_b_max_0, a_b_max_1, a_b_max_2)) = a_b.get() {
        let a_pc = [a_p.x, a_p.y, a_p.z];
        let a_bmin = [a_b_min_0, a_b_min_1, a_b_min_2];
        let a_bmax = [a_b_max_0, a_b_max_1, a_b_max_2];
        let mut a_dist = 0.0;
        for i in 0..3 {
            let a_r1 = a_bmin[i] - a_pc[i];
            if a_r1 > 0.0 {
                a_dist += a_r1 * a_r1;
                continue;
            }
            let a_r2 = a_pc[i] - a_bmax[i];
            if a_r2 > 0.0 { a_dist += a_r2 * a_r2; }
        }
        a_dist.sqrt()
    } else {
        f64::MAX
    }
}

/// OCCT Bnd_Box::IsIn (point inside the box).
fn point_in_box(a_b: &BndBox, a_p: DVec3) -> bool {
    point_box_distance(a_b, a_p) <= 1e-12
}

/// OCCT Bnd_Box::IsXThin/IsYThin/IsZThin — raw corner extents within tolerance.
fn box_is_thin(a_b: &BndBox, the_tol: f64) -> bool {
    if a_b.is_void() { return false; }
    let mn = a_b.raw_min();
    let mx = a_b.raw_max();
    (mx.x - mn.x) <= the_tol && (mx.y - mn.y) <= the_tol && (mx.z - mn.z) <= the_tol
}

/// OCCT Bnd_Box::SquareExtent — diagonal length squared.
fn box_square_extent(a_b: &BndBox) -> f64 {
    if let Some((x1, y1, z1, x2, y2, z2)) = a_b.get() {
        (x2 - x1).powi(2) + (y2 - y1).powi(2) + (z2 - z1).powi(2)
    } else {
        f64::MAX
    }
}

/// OCCT SplitRangeOnSegments (L1366-1406).
fn split_range_on_segments(a_t1: f64, a_t2: f64, the_resolution: f64, the_nb_seg: i32) -> Vec<[f64; 2]> {
    let a_diff = a_t2 - a_t1;
    let mut the_segments = Vec::new();
    if a_diff < the_resolution || the_nb_seg == 1 {
        the_segments.push([a_t1, a_t2]);
        return the_segments;
    }
    let mut a_nb_segments = the_nb_seg as usize;
    let mut a_dt = a_diff / a_nb_segments as f64;
    if a_dt < the_resolution {
        let a_seg = a_diff / the_resolution;
        a_nb_segments = (a_seg as usize) + 1;
        a_dt = a_diff / a_nb_segments as f64;
    }
    let mut a_t1x = a_t1;
    for _ in 1..a_nb_segments {
        let a_t2x = a_t1x + a_dt;
        the_segments.push([a_t1x, a_t2x]);
        a_t1x = a_t2x;
    }
    the_segments.push([a_t1x, a_t2]);
    the_segments
}

/// OCCT ResolutionCoeff (L1486-1557).
fn resolution_coeff(the_bac: &Curve3, the_range: [f64; 2]) -> f64 {
    let a_curve = the_bac;
    match a_curve {
        Curve3::Circle(c) => 1.0 / (2.0 * c.radius),
        Curve3::Ellipse(e) => 1.0 / e.major_radius,
        _ => {
            // OCCT L1522-1550: sample-based for Hyperbola/Parabola/Other
            let a_nb_p = 30usize;
            let a_dt = (the_range[1] - the_range[0]) / a_nb_p as f64;
            let mut a_t = the_range[0];
            let mut a_p1 = the_bac.point_at(the_range[0]);
            let mut k_min = 10.0f64;
            for _ in 1..=a_nb_p {
                a_t += a_dt;
                let a_p2 = the_bac.point_at(a_t);
                let a_dist = a_p1.distance(a_p2);
                let k = a_dt / a_dist.max(1e-12);
                if k < k_min { k_min = k; }
                a_p1 = a_p2;
            }
            k_min
        }
    }
}

/// OCCT Resolution (L1561-1607).
fn resolution(the_curve: &Curve3, the_res_coeff: f64, the_r3d: f64) -> f64 {
    match the_curve {
        Curve3::Line(_) => the_r3d,
        Curve3::Circle(_) => {
            let a_dt = the_res_coeff * the_r3d;
            if a_dt <= 1.0 { 2.0 * a_dt.asin() } else { 2.0 * std::f64::consts::PI }
        }
        _ => the_res_coeff * the_r3d,
    }
}

/// OCCT CurveDeflection (L1611-1638).
fn curve_deflection(the_bac: &Curve3, the_range: [f64; 2]) -> f64 {
    let a_nb_p = 10usize;
    let a_dt = (the_range[1] - the_range[0]) / a_nb_p as f64;
    let mut a_t = the_range[0];
    let mut a_defl = 0.0f64;
    let mut a_v1 = the_bac.tangent_at(the_range[0]);
    for _ in 1..=a_nb_p {
        a_t += a_dt;
        let a_v2 = the_bac.tangent_at(a_t);
        if a_v1.length_squared() > 1e-24 && a_v2.length_squared() > 1e-24 {
            let a_d1 = a_v1.normalize();
            let a_d2 = a_v2.normalize();
            a_defl += a_d1.dot(a_d2).clamp(-1.0, 1.0).acos();
        }
        a_v1 = a_v2;
    }
    a_defl
}

/// OCCT IsClosed (L1642-1659).
fn is_closed(the_curve: &Curve3, a_t1: f64, a_t2: f64, the_tol: f64, the_res: f64) -> bool {
    if (a_t1 - a_t2).abs() < the_res { return false; }
    let a_p1 = the_curve.point_at(a_t1);
    let a_p2 = the_curve.point_at(a_t2);
    a_p1.distance(a_p2) <= the_tol
}

/// OCCT IntTools_EdgeEdge::FindParameters (L553-671).
#[allow(clippy::too_many_arguments)]
fn find_parameters(
    the_bac: &Curve3,
    a_t1: f64, a_t2: f64,
    the_tol: f64, the_res: f64, the_p_tol: f64, the_res_coeff: f64,
    the_c_box: &BndBox,
    a_tb1: &mut f64, a_tb2: &mut f64,
) -> bool {
    let a_cf = 0.6180339887498948482045868343656; // OCCT: =0.5*(1.+sqrt(5.))/2.
    let mut a_cbx = the_c_box.clone();
    a_cbx.set_gap(a_cbx.get_gap() + the_tol);
    let a_curve = the_bac;
    let a_max_dt = (a_t2 - a_t1) * 0.01;
    let mut tb1 = *a_tb1;
    let mut tb2 = *a_tb2;
    let mut b_ret = false;
    for i in 0..2 {
        let mut a_tb = if i == 0 { a_t1 } else { a_t2 };
        let a_t = if i == 0 { a_t2 } else { tb1 };
        let a_c: f64 = if i == 0 { 1.0 } else { -1.0 };
        let mut a_dt = the_res;
        let mut a_dist_p = 0.0f64;
        let mut b_ret = false;
        let mut k = 1.0f64;
        // OCCT L590-622: looking for the point on the edge which is in the box
        while a_c * (a_t - a_tb) >= 0.0 {
            let a_p = a_curve.point_at(a_tb);
            let a_dist = point_box_distance(&a_cbx, a_p);
            if a_dist > the_tol {
                if a_dist_p > 0.0 {
                    let mut to_grow = false;
                    if (a_dist_p - a_dist).abs() / a_dist_p < 0.1 {
                        a_dt = resolution(a_curve, the_res_coeff, k * a_dist);
                        if a_dt < a_max_dt {
                            to_grow = true;
                            k *= 2.0;
                        }
                    }
                    if !to_grow {
                        k = 1.0;
                        a_dt = resolution(a_curve, the_res_coeff, a_dist);
                    }
                }
                a_tb += a_c * a_dt;
            } else {
                b_ret = true;
                break;
            }
            a_dist_p = a_dist;
        }
        // OCCT L624-637
        if !b_ret {
            if i == 0 {
                // edge is out of the box
                return false;
            } else {
                a_tb = tb1;
                a_dt = a_t2 - tb1;
            }
        }
        // OCCT L639-660: golden-section bisection
        if a_tb != a_t {
            let mut a_t_in = a_tb;
            let mut a_t_out = a_tb - a_c * a_dt;
            let mut a_diff = a_t_in - a_t_out;
            while a_diff.abs() > the_p_tol {
                a_tb = a_t_out + a_diff * a_cf;
                let a_p = a_curve.point_at(a_tb);
                if a_cbx.is_out_point(a_p) {
                    a_t_out = a_tb;
                } else {
                    a_t_in = a_tb;
                }
                a_diff = a_t_in - a_t_out;
            }
        }
        if i == 0 { tb1 = a_tb; } else { tb2 = a_tb; }
    }
    if tb2 < tb1 { return false; }
    *a_tb1 = tb1;
    *a_tb2 = tb2;
    true
}

/// OCCT GeomAPI_ProjectPointOnCurve — project a point onto a curve within [t1, t2].
fn project_on_range(curve: &Curve3, point: DVec3, t1: f64, t2: f64) -> Option<(f64, f64)> {
    if let Curve3::Line(l) = curve {
        let dir_sq = l.direction.dot(l.direction);
        if dir_sq > 1e-24 {
            let t = ((point - l.origin).dot(l.direction) / dir_sq).clamp(t1, t2);
            let p = l.origin + l.direction * t;
            return Some(((p - point).length(), t));
        }
    }
    if let Curve3::Circle(c) = curve {
        let r = c.radius;
        let center = c.center;
        let d = point - center;
        let radial = d - c.normal * d.dot(c.normal);
        let radial_len = radial.length();
        if radial_len > 1e-30 {
            let proj = center + radial * (r / radial_len);
            let p_rel = proj - center;
            let mut ang = p_rel.dot(c.x_dir).atan2(p_rel.dot(c.y_dir));
            if ang < 0.0 { ang += std::f64::consts::TAU; }
            let ta = t1.rem_euclid(std::f64::consts::TAU);
            let tb = if t2 > t1 { t2 } else { t2 + std::f64::consts::TAU };
            let mut ang_wrap = ang;
            if ang_wrap < ta { ang_wrap += std::f64::consts::TAU; }
            if ang_wrap >= ta && ang_wrap <= tb {
                return Some((proj.distance(point), ang_wrap));
            }
            let d1 = (curve.point_at(t1) - point).length();
            let d2 = (curve.point_at(t2) - point).length();
            if d1 <= d2 { return Some((d1, t1)); } else { return Some((d2, t2)); }
        }
    }
    // General curve: sample + local refine within [t1, t2].
    let n = 64usize;
    let dt = (t2 - t1) / n as f64;
    let mut best_t = t1;
    let mut best_d = f64::MAX;
    for i in 0..=n {
        let t = t1 + dt * i as f64;
        let d = (curve.point_at(t) - point).length();
        if d < best_d { best_d = d; best_t = t; }
    }
    let mut w = dt * 2.0;
    for _ in 0..12 {
        let mut improved = false;
        for i in 0..=8 {
            let t = (best_t - w + 2.0 * w * i as f64 / 8.0).clamp(t1, t2);
            let d = (curve.point_at(t) - point).length();
            if d < best_d { best_d = d; best_t = t; improved = true; }
        }
        if !improved { break; }
        w *= 0.5;
    }
    Some((best_d, best_t))
}

/// OCCT DistPC (L1332-1362).
fn dist_pc(
    the_c1: &Curve3, a_t1: f64, the_criteria: f64,
    the_proj_pc: &mut Projector, a_d: &mut f64, a_t2: &mut f64, i_c: i32,
) -> i32 {
    let a_p1 = the_c1.point_at(a_t1);
    match the_proj_pc.perform(a_p1) {
        None => 1,
        Some((d, t2)) => {
            *a_d = d;
            *a_t2 = t2;
            if i_c as f64 * (d - the_criteria) > 0.0 { 2 } else { 0 }
        }
    }
}

/// OCCT DistPC (L1301-1328) with max tracking.
#[allow(clippy::too_many_arguments)]
fn dist_pc_max(
    the_c1: &Curve3, a_t1: f64, the_criteria: f64,
    the_proj_pc: &mut Projector,
    a_d: &mut f64, a_t2: &mut f64,
    a_dmax: &mut f64, a_t1max: &mut f64, a_t2max: &mut f64, i_c: i32,
) -> i32 {
    let i_err = dist_pc(the_c1, a_t1, the_criteria, the_proj_pc, a_d, a_t2, i_c);
    if i_err == 1 { return i_err; }
    if i_c as f64 * (*a_d - *a_dmax) > 0.0 {
        *a_dmax = *a_d;
        *a_t1max = a_t1;
        *a_t2max = *a_t2;
    }
    i_err
}

/// OCCT FindDistPC (L1210-1297): golden-section search for the extrema distance.
#[allow(clippy::too_many_arguments)]
fn find_dist_pc(
    a_t1a: f64, a_t1b: f64,
    the_c1: &Curve3, the_criteria: f64, the_eps: f64,
    the_proj_pc: &mut Projector,
    a_dmax: &mut f64, a_t1max: &mut f64, a_t2max: &mut f64,
    b_max_dist: bool,
) -> i32 {
    let i_c = if b_max_dist { 1 } else { -1 };
    let a_gs = 0.6180339887498948482045868343656;
    let mut a_a = a_t1a;
    let mut a_b = a_t1b;
    let (mut a_yp, mut a_t2p) = (0.0f64, 0.0f64);
    let (mut a_yl, mut a_t2l) = (0.0f64, 0.0f64);
    // check bounds
    let mut i_err = dist_pc_max(the_c1, a_a, the_criteria, the_proj_pc, &mut a_yp, &mut a_t2p, a_dmax, a_t1max, a_t2max, i_c);
    if i_err == 2 { return i_err; }
    i_err = dist_pc_max(the_c1, a_b, the_criteria, the_proj_pc, &mut a_yl, &mut a_t2l, a_dmax, a_t1max, a_t2max, i_c);
    if i_err == 2 { return i_err; }
    let mut a_xp = a_a + (a_b - a_a) * a_gs;
    let mut a_xl = a_b - (a_b - a_a) * a_gs;
    i_err = dist_pc_max(the_c1, a_xp, the_criteria, the_proj_pc, &mut a_yp, &mut a_t2p, a_dmax, a_t1max, a_t2max, i_c);
    if i_err != 0 { return i_err; }
    i_err = dist_pc_max(the_c1, a_xl, the_criteria, the_proj_pc, &mut a_yl, &mut a_t2l, a_dmax, a_t1max, a_t2max, i_c);
    if i_err != 0 { return i_err; }
    let an_eps = the_eps.max(a_a.abs().max(a_b.abs()) * 1e-14);
    loop {
        if i_c as f64 * (a_yp - a_yl) > 0.0 {
            a_a = a_xl;
            a_xl = a_xp;
            a_yl = a_yp;
            a_xp = a_a + (a_b - a_a) * a_gs;
            i_err = dist_pc_max(the_c1, a_xp, the_criteria, the_proj_pc, &mut a_yp, &mut a_t2p, a_dmax, a_t1max, a_t2max, i_c);
        } else {
            a_b = a_xp;
            a_xp = a_xl;
            a_yp = a_yl;
            a_xl = a_b - (a_b - a_a) * a_gs;
            i_err = dist_pc_max(the_c1, a_xl, the_criteria, the_proj_pc, &mut a_yl, &mut a_t2l, a_dmax, a_t1max, a_t2max, i_c);
        }
        if i_err != 0 {
            if i_err == 2 && !b_max_dist {
                let a_xp = (a_a + a_b) * 0.5;
                let _ = dist_pc_max(the_c1, a_xp, the_criteria, the_proj_pc, &mut a_yp, &mut a_t2p, a_dmax, a_t1max, a_t2max, i_c);
            }
            return i_err;
        }
        if (a_b - a_a) < an_eps { break; }
    }
    i_err
}

/// OCCT IntTools_EdgeEdge::CheckCoincidence (L1150-1206).
fn check_coincidence(
    the_c1: &Curve3, the_c2: &Curve3,
    a_t11: f64, a_t12: f64, a_t21: f64, a_t22: f64,
    the_criteria: f64, the_curve_res1: f64,
) -> i32 {
    let mut the_proj_pc = Projector::new(the_c2, a_t21, a_t22);
    let mut a_dmax = -1.0f64;
    let (mut a_t1max, mut a_t2max) = (0.0f64, 0.0f64);
    // 1. Express evaluation
    let a_nb = 10;
    let a_ranges = split_range_on_segments(a_t11, a_t12, the_curve_res1, a_nb);
    let a_nb1 = a_ranges.len();
    let mut i_err = 0i32;
    for i in 1..a_nb1 {
        let (a_t1a, a_t1b) = (a_ranges[i][0], a_ranges[i][1]);
        let mut a_d = 0.0;
        let mut a_t2 = 0.0;
        i_err = dist_pc_max(the_c1, a_t1b, the_criteria, &mut the_proj_pc, &mut a_d, &mut a_t2, &mut a_dmax, &mut a_t1max, &mut a_t2max, 1);
        if i_err != 0 { return i_err; }
    }
    // if the ranges are fewer than requested, no deep evaluation needed
    if a_nb1 < a_nb as usize { return i_err; }
    // 2. Deep evaluation
    for i in 2..a_nb1 {
        let (a_t1a, a_t1b) = (a_ranges[i][0], a_ranges[i][1]);
        i_err = find_dist_pc(a_t1a, a_t1b, the_c1, the_criteria, the_curve_res1, &mut the_proj_pc, &mut a_dmax, &mut a_t1max, &mut a_t2max, false);
        if i_err != 0 { return i_err; }
    }
    i_err
}

/// OCCT IntTools_EdgeEdge::IsIntersection (L1060-1146).
fn is_intersection(
    the_c1: &Curve3, the_c2: &Curve3,
    a_t11: f64, a_t12: f64, a_t21: f64, a_t22: f64,
    my_tol: f64, my_res1: f64, my_res2: f64,
) -> bool {
    let mut a_coef = 1e5f64;
    if (a_t12 - a_t11) > a_coef * my_res1 && (a_t22 - a_t21) > a_coef * my_res2 {
        a_coef = 5000.0;
    } else {
        let a_tr_min = ((a_t12 - a_t11) / my_res1).min((a_t22 - a_t21) / my_res2);
        a_coef = a_tr_min / 100.0;
        if a_coef < 1.0 { a_coef = 1.0; }
    }
    let a_criteria = a_coef * my_tol;
    let a_criteria = a_criteria * a_criteria;
    let a_p11 = the_c1.point_at(a_t11);
    let a_p12 = the_c1.point_at(a_t12);
    let a_p21 = the_c2.point_at(a_t21);
    let a_p22 = the_c2.point_at(a_t22);
    let a_d11_21 = a_p11.distance_squared(a_p21);
    let a_d11_22 = a_p11.distance_squared(a_p22);
    let a_d12_21 = a_p12.distance_squared(a_p21);
    let a_d12_22 = a_p12.distance_squared(a_p22);
    let b_small_11_21 = a_d11_21 < a_criteria;
    let b_small_11_22 = a_d11_22 < a_criteria;
    let b_small_12_21 = a_d12_21 < a_criteria;
    let b_small_12_22 = a_d12_22 < a_criteria;
    if (b_small_11_21 && b_small_12_22) || (b_small_11_22 && b_small_12_21) {
        if a_coef == 1.0 { return true; }
        let an_angle_criteria = 5e-3f64;
        let (mut an_angle1, mut an_angle2) = (0.0f64, 0.0f64);
        let a_v11 = the_c1.tangent_at(a_t11);
        let a_v12 = the_c1.tangent_at(a_t12);
        let a_v21 = the_c2.tangent_at(a_t21);
        let a_v22 = the_c2.tangent_at(a_t22);
        let vv = |v: DVec3| v.length_squared() > 1e-24;
        if vv(a_v11) && vv(a_v12) && vv(a_v21) && vv(a_v22) {
            if b_small_11_21 && b_small_12_22 {
                an_angle1 = a_v11.normalize().dot(a_v21.normalize()).clamp(-1.0, 1.0).acos();
                an_angle2 = a_v12.normalize().dot(a_v22.normalize()).clamp(-1.0, 1.0).acos();
            } else {
                an_angle1 = a_v11.normalize().dot(a_v22.normalize()).clamp(-1.0, 1.0).acos();
                an_angle2 = a_v12.normalize().dot(a_v21.normalize()).clamp(-1.0, 1.0).acos();
            }
        }
        if (an_angle1 < an_angle_criteria || (std::f64::consts::PI - an_angle1) < an_angle_criteria)
            || (an_angle2 < an_angle_criteria || (std::f64::consts::PI - an_angle2) < an_angle_criteria) {
            let mut the_proj_pc = Projector::new(the_c2, a_t21, a_t22);
            let mut a_d = f64::MAX;
            let (mut a_t1_min, mut a_t2_min) = (0.0f64, 0.0f64);
            let i_err = find_dist_pc(a_t11, a_t12, the_c1, my_tol, my_res1, &mut the_proj_pc, &mut a_d, &mut a_t1_min, &mut a_t2_min, false);
            return i_err == 2;
        }
    }
    true
}

/// OCCT IntTools_EdgeEdge::FindBestSolution (L826-898).
fn find_best_solution(
    the_c1: &Curve3, the_c2: &Curve3,
    the_res_coeff1: f64, my_tol: f64, the_p_tol1: f64,
    a_t11: f64, a_t12: f64, a_t21: f64, a_t22: f64,
) -> (f64, f64) {
    let mut a_d_min = f64::MAX;
    let a_sol_criteria = 5e-16f64;
    let a_touch_criteria = 5e-13f64;
    let mut b_touch = false;
    let mut b_touch_confirm = false;
    let a_res1 = resolution(the_c1, the_res_coeff1, my_tol);
    let a_nb_s = split_range_on_segments(a_t11, a_t12, 3.0 * a_res1, 10);
    let mut the_proj_pc = Projector::new(the_c2, a_t21, a_t22);
    let mut a_t1 = a_t11;
    let mut a_t2 = a_t21;
    let (mut a_t11_touch, mut a_t12_touch) = (a_t11, a_t12);
    let (mut a_t21_touch, mut a_t22_touch) = (a_t21, a_t22);
    let mut is_sol_found = false;
    for r in &a_nb_s {
        let (a_t1a, a_t1b) = (r[0], r[1]);
        let mut a_d = my_tol;
        let (mut a_t1_min, mut a_t2_min) = (0.0f64, 0.0f64);
        let i_err = find_dist_pc(a_t1a, a_t1b, the_c1, a_sol_criteria, the_p_tol1, &mut the_proj_pc, &mut a_d, &mut a_t1_min, &mut a_t2_min, false);
        if i_err != 1 {
            if a_d < a_d_min {
                a_t1 = a_t1_min;
                a_t2 = a_t2_min;
                a_d_min = a_d;
                is_sol_found = true;
            }
            if a_d < a_touch_criteria {
                if b_touch {
                    a_t12_touch = a_t1_min;
                    a_t22_touch = a_t2_min;
                    b_touch_confirm = true;
                } else {
                    a_t11_touch = a_t1_min;
                    a_t21_touch = a_t2_min;
                    b_touch = true;
                }
            }
        }
    }
    if !is_sol_found || b_touch_confirm {
        a_t1 = (a_t11_touch + a_t12_touch) * 0.5;
        let mut a_d = 0.0;
        let mut a_t2_out = 0.0;
        let i_err = dist_pc(the_c1, a_t1, a_sol_criteria, &mut the_proj_pc, &mut a_d, &mut a_t2_out, -1);
        if i_err == 1 {
            a_t2 = (a_t21_touch + a_t22_touch) * 0.5;
        } else {
            a_t2 = a_t2_out;
        }
    }
    (a_t1, a_t2)
}

/// OCCT GeomAPI_ProjectPointOnCurve — projects points onto a curve within a range.
struct Projector<'a> {
    curve: &'a Curve3,
    t1: f64,
    t2: f64,
}

impl<'a> Projector<'a> {
    fn new(curve: &'a Curve3, t1: f64, t2: f64) -> Self {
        Projector { curve, t1, t2 }
    }
    fn perform(&mut self, point: DVec3) -> Option<(f64, f64)> {
        project_on_range(self.curve, point, self.t1, self.t2)
    }
}
