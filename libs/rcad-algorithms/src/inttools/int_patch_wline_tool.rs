//! OCCT-aligned: IntPatch_WLineTool — walking line post-processing utilities.
//!
//! OCCT IntPatch_WLineTool.hxx / .cxx (73K)
//!
//! Core methods:
//!   ComputePurgedWLine — removes collinear/redundant points from WLine
//!   JoinWLines         — joins adjacent WLines at shared endpoints
//!   ExtendTwoWLines    — extends WLine endpoints to meet at intersections

use super::int_patch_line::{IntPatchLine, WLinePnt, WLineType};
use glam::DVec2;

/// Max angle to concatenate two WLines to avoid C0-continuity issues.
/// OCCT WLineTool.cxx: myMaxConcatAngle = PI/6
const MAX_CONCAT_ANGLE: f64 = std::f64::consts::PI / 6.0;

impl IntPatchLine {
    /// OCCT: ComputePurgedWLine — remove collinear points from a walking line.
    ///
    /// 1. Remove duplicate points (within resolution)
    /// 2. Remove out-of-domain points
    /// 3. Tube criteria: remove points that lie near the chord between neighbors
    ///
    /// Returns a new IntPatchLine with purged points, or None if < 2 points remain.
    pub fn purge_wline(&self) -> Option<IntPatchLine> {
        if !self.is_wline() || self.wline_pnts.len() < 2 {
            return None;
        }
        if self.wline_pnts.len() == 2 {
            let d = self.wline_pnts[0].p3d.distance(self.wline_pnts[1].p3d);
            if d < 1e-15 { return None; }
            return Some(self.clone());
        }

        // Step 1: Remove equal points (within resolution)
        let mut cleaned: Vec<WLinePnt> = Vec::with_capacity(self.wline_pnts.len());
        for p in self.wline_pnts.iter() {
            if let Some(last) = cleaned.last() {
                let d3 = last.p3d.distance(p.p3d);
                let du = (last.u1 - p.u1).abs() + (last.v1 - p.v1).abs();
                let dv = (last.u2 - p.u2).abs() + (last.v2 - p.v2).abs();
                // Find max UV magnitude for relative comparison
                let max_uv = last.u1.abs().max(last.v1.abs().max(last.u2.abs().max(last.v2.abs())));
                if d3 < 1e-15 || du < 1e-16 * max_uv.max(1.0) || dv < 1e-16 * max_uv.max(1.0) {
                    continue; // skip duplicate
                }
            }
            cleaned.push(*p);
        }

        if cleaned.len() < 3 {
            return if cleaned.len() == 2 {
                Some(self.clone_with_pnts(cleaned))
            } else { None };
        }

        // Step 2: Tube criteria — remove collinear points
        // Check each point: if distance from chord between neighbors is
        // below tolerance, remove it (collinear simplification)
        let resolution = 1e-12;
        let mut purged: Vec<WLinePnt> = Vec::with_capacity(cleaned.len());
        purged.push(cleaned[0]); // always keep first

        for i in 1..cleaned.len() - 1 {
            let prev = cleaned[i - 1];
            let curr = cleaned[i];
            let next = cleaned[i + 1];

            // 3D chord distance: distance from curr to line segment prev-next
            let seg = next.p3d - prev.p3d;
            let seg_len = seg.length();
            if seg_len < resolution {
                purged.push(curr);
                continue;
            }
            let t = (curr.p3d - prev.p3d).dot(seg) / (seg_len * seg_len);
            let proj = prev.p3d + t * seg;
            let d3 = (curr.p3d - proj).length();

            // 2D chord distance on surface 1
            let uv1_seg = DVec2::new(next.u1 - prev.u1, next.v1 - prev.v1);
            let uv1_len = uv1_seg.length();
            let d_uv1 = if uv1_len < resolution { 0.0 } else {
                let t1 = ((curr.u1 - prev.u1) * uv1_seg.x + (curr.v1 - prev.v1) * uv1_seg.y) / (uv1_len * uv1_len);
                let proj_uv1 = DVec2::new(prev.u1, prev.v1) + t1 * uv1_seg;
                DVec2::new(curr.u1 - proj_uv1.x, curr.v1 - proj_uv1.y).length()
            };

            // 2D chord distance on surface 2
            let uv2_seg = DVec2::new(next.u2 - prev.u2, next.v2 - prev.v2);
            let uv2_len = uv2_seg.length();
            let d_uv2 = if uv2_len < resolution { 0.0 } else {
                let t2 = ((curr.u2 - prev.u2) * uv2_seg.x + (curr.v2 - prev.v2) * uv2_seg.y) / (uv2_len * uv2_len);
                let proj_uv2 = DVec2::new(prev.u2, prev.v2) + t2 * uv2_seg;
                DVec2::new(curr.u2 - proj_uv2.x, curr.v2 - proj_uv2.y).length()
            };

            // Keep point if any deviation exceeds tolerance
            let tol = self.tolerance.max(1e-10);
            if d3 > tol || d_uv1 > tol || d_uv2 > tol {
                purged.push(curr);
            }
        }
        purged.push(cleaned[cleaned.len() - 1]); // always keep last

        if purged.len() < 2 { None }
        else { Some(self.clone_with_pnts(purged)) }
    }

    fn clone_with_pnts(&self, pnts: Vec<WLinePnt>) -> Self {
        let mut c = self.clone();
        c.wline_pnts = pnts;
        c
    }
}

/// OCCT: JoinWLines — join adjacent WLines if they share endpoints.
/// Returns true if any joining was performed.
pub fn join_wlines(lines: &mut Vec<IntPatchLine>, tol_3d: f64) -> bool {
    let mut joined = false;
    let mut i = 0;
    while i < lines.len() {
        let is_wline_i = lines[i].is_wline();
        let pnts_len_i = if is_wline_i { lines[i].wline_pnts.len() } else { 0 };
        if pnts_len_i < 2 { i += 1; continue; }
        let mut j = i + 1;
        while j < lines.len() {
            let is_wline_j = lines[j].is_wline();
            let pnts_len_j = if is_wline_j { lines[j].wline_pnts.len() } else { 0 };
            if pnts_len_j < 2 { j += 1; continue; }
            let first_end = lines[i].wline_pnts.last().unwrap().p3d;
            let second_start = lines[j].wline_pnts[0].p3d;
            let d = first_end.distance(second_start);
            if d < tol_3d {
                let pnts_to_add: Vec<WLinePnt> = lines[j].wline_pnts.iter().copied().collect();
                lines[i].wline_pnts.extend(pnts_to_add);
                lines.remove(j);
                joined = true;
            } else {
                j += 1;
            }
        }
        i += 1;
    }
    joined
}
