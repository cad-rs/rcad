//! Structured healing pipeline for B-Rep analysis and repair.
//!
//! This module provides an analyze -> repair -> recheck workflow similar in
//! spirit to OCCT ShapeAnalysis/ShapeFix orchestration.

use rcad_kernel::BRep;

use crate::brep_check::{CheckIssue, CheckResult, check};
use crate::brep_repair::{RepairReport, repair};
use crate::tolerance::TOLERANCE_ABS;

/// Healing execution mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HealingMode {
    /// Only analyze; no repair pass will run.
    AnalyzeOnly,
    /// Analyze and run repair passes.
    AnalyzeAndRepair,
}

/// Options controlling healing execution.
#[derive(Debug, Clone, Copy)]
pub struct HealingOptions {
    /// Repair tolerance used by [`repair`].
    pub tolerance: f64,
    /// Maximum number of repair passes.
    pub max_passes: usize,
    /// Execution mode for the pipeline.
    pub mode: HealingMode,
}

impl Default for HealingOptions {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_ABS,
            max_passes: 2,
            mode: HealingMode::AnalyzeAndRepair,
        }
    }
}

/// Structured issue counters for checker output.
#[derive(Debug, Clone, Default)]
pub struct HealingIssueStats {
    pub open_wire: usize,
    pub zero_normal: usize,
    pub degenerate_face: usize,
    pub invalid_edge_index: usize,
    pub invalid_vertex_index: usize,
    pub non_manifold_edge: usize,
    pub self_intersecting_wire: usize,
    pub geometric_self_intersection: usize,
}

impl HealingIssueStats {
    pub fn total(&self) -> usize {
        self.open_wire
            + self.zero_normal
            + self.degenerate_face
            + self.invalid_edge_index
            + self.invalid_vertex_index
            + self.non_manifold_edge
            + self.self_intersecting_wire
            + self.geometric_self_intersection
    }

    pub fn from_check_result(result: &CheckResult) -> Self {
        let mut s = Self::default();
        for issue in &result.issues {
            match issue {
                CheckIssue::OpenWire { .. } => s.open_wire += 1,
                CheckIssue::ZeroNormal { .. } => s.zero_normal += 1,
                CheckIssue::DegenerateFace { .. } => s.degenerate_face += 1,
                CheckIssue::InvalidEdgeIndex { .. } => s.invalid_edge_index += 1,
                CheckIssue::InvalidVertexIndex { .. } => s.invalid_vertex_index += 1,
                CheckIssue::NonManifoldEdge { .. } => s.non_manifold_edge += 1,
                CheckIssue::SelfIntersectingWire { .. } => s.self_intersecting_wire += 1,
                CheckIssue::GeometricSelfIntersection { .. } => s.geometric_self_intersection += 1,
            }
        }
        s
    }
}

/// Summary report for analyze/heal workflow.
#[derive(Debug, Clone)]
pub struct HealingReport {
    /// Issues found before any repair.
    pub initial: CheckResult,
    /// Issues after the final pass.
    pub final_result: CheckResult,
    /// Per-pass repair reports.
    pub passes: Vec<RepairReport>,
    /// Structured issue counters before healing.
    pub initial_stats: HealingIssueStats,
    /// Structured issue counters after healing.
    pub final_stats: HealingIssueStats,
}

impl HealingReport {
    pub fn initial_issue_count(&self) -> usize {
        self.initial.issues.len()
    }

    pub fn final_issue_count(&self) -> usize {
        self.final_result.issues.len()
    }

    pub fn fixed_issue_count(&self) -> usize {
        self.initial_issue_count().saturating_sub(self.final_issue_count())
    }

    pub fn is_improved(&self) -> bool {
        self.final_issue_count() < self.initial_issue_count()
    }

    pub fn is_clean(&self) -> bool {
        self.final_result.is_valid()
    }

    pub fn has_issue_kind(&self, pred: impl Fn(&CheckIssue) -> bool) -> bool {
        self.final_result.issues.iter().any(pred)
    }
}

/// Analyze and heal a BRep using the provided options.
pub fn analyze_and_heal(brep: &BRep, options: HealingOptions) -> (BRep, HealingReport) {
    let initial = check(brep);
    let initial_stats = HealingIssueStats::from_check_result(&initial);

    if matches!(options.mode, HealingMode::AnalyzeOnly) {
        return (
            brep.clone(),
            HealingReport {
                initial: initial.clone(),
                final_result: initial,
                passes: Vec::new(),
                initial_stats: initial_stats.clone(),
                final_stats: initial_stats,
            },
        );
    }

    if initial.is_valid() {
        return (
            brep.clone(),
            HealingReport {
                initial: initial.clone(),
                final_result: initial,
                passes: Vec::new(),
                initial_stats: initial_stats.clone(),
                final_stats: initial_stats,
            },
        );
    }

    let mut current = brep.clone();
    let mut passes = Vec::new();
    let pass_count = options.max_passes.max(1);

    for _ in 0..pass_count {
        let (next, rep) = repair(&current, options.tolerance);
        current = next;
        let no_changes = rep.vertices_merged == 0
            && rep.degenerate_faces_removed == 0
            && rep.normals_recomputed == 0
            && rep.wires_fixed == 0;
        passes.push(rep);

        let chk = check(&current);
        if chk.is_valid() || no_changes {
            let final_stats = HealingIssueStats::from_check_result(&chk);
            return (
                current,
                HealingReport {
                    initial,
                    final_result: chk,
                    passes,
                    initial_stats,
                    final_stats,
                },
            );
        }
    }

    let final_result = check(&current);
    let final_stats = HealingIssueStats::from_check_result(&final_result);
    (
        current,
        HealingReport {
            initial,
            final_result,
            passes,
            initial_stats,
            final_stats,
        },
    )
}

/// Convenience wrapper using default options.
pub fn heal(brep: &BRep) -> (BRep, HealingReport) {
    analyze_and_heal(brep, HealingOptions::default())
}

#[cfg(test)]
mod tests {
    use glam::DVec3;
    use rcad_kernel::PrimitiveSolid;

    use super::*;
    use crate::geom_populate;

    fn unit_box() -> BRep {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        geom_populate::populate_box_geom(&mut brep);
        brep
    }

    #[test]
    fn heal_valid_box_is_noop() {
        let b = unit_box();
        let (out, report) = heal(&b);
        assert!(report.initial.is_valid());
        assert!(report.final_result.is_valid());
        assert!(report.passes.is_empty());
        assert_eq!(out.vertices.len(), b.vertices.len());
        assert_eq!(out.edges.len(), b.edges.len());
    }

    #[test]
    fn heal_zero_normal_face_gets_fixed() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let (out, report) = heal(&b);
        assert!(report.initial_issue_count() >= 1);
        assert!(report.is_improved() || report.is_clean());
        assert!(report.initial_stats.zero_normal >= 1);
        assert_eq!(report.initial_stats.total(), report.initial_issue_count());
        assert_eq!(report.final_stats.total(), report.final_issue_count());
        assert!(!out.solids[0].shells[0].faces[0].normal.abs_diff_eq(DVec3::ZERO, 0.0));
    }

    #[test]
    fn analyze_only_preserves_input_and_reports_issues() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let (out, report) = analyze_and_heal(
            &b,
            HealingOptions {
                mode: HealingMode::AnalyzeOnly,
                ..HealingOptions::default()
            },
        );

        assert!(report.initial_issue_count() >= 1);
        assert_eq!(report.initial_issue_count(), report.final_issue_count());
        assert!(report.passes.is_empty());
        assert_eq!(out.solids[0].shells[0].faces[0].normal, DVec3::ZERO);
    }
}
