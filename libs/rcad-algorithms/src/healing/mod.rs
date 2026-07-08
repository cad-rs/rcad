//! Structured healing pipeline for B-Rep analysis and repair.
//!
//! This module provides an analyze -> repair -> recheck workflow similar in
//! spirit to OCCT ShapeAnalysis/ShapeFix orchestration.


use rcad_kernel::BRep;

use glam::DVec3;
use rcad_kernel::topods;
use rcad_kernel::{BSplineSurface, BezierSurface};

use crate::brep_check::{
    CheckIssue, CheckResult, brep_check_analyze, diagnose_same_parameter, diagnose_same_range,
};
use crate::brep_repair::{
    MakeConnectedReport, RepairReport, fix_same_parameter_with_scan,
    fix_same_range_with_scan, make_connected_iterative_with_growth_cap, repair,
};
use crate::tolerance::{
    TOLERANCE_ABS, TOLERANCE_ADAPTIVE_MAX, TOLERANCE_COORD_SUB, TOLERANCE_LEN_MIN,
    TOLERANCE_LINEAR_ULTRA_STRICT, TOLERANCE_MESH_LEGACY, TOLERANCE_RETRY_LADDER_COARSE,
    TOLERANCE_RETRY_LADDER_MID, TOLERANCE_VOL_CUBE_FACTOR,
};

/// Healing execution mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HealingMode {
    /// Only analyze; no repair pass will run.
    AnalyzeOnly,
    /// Analyze and run repair passes.
    AnalyzeAndRepair,
}

/// Policy controlling whether a make-connected prepass is executed before
/// regular repair passes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MakeConnectedPrepassMode {
    /// Never run a prepass.
    Disabled,
    /// Run only when initial checker issues indicate connectivity stress.
    IssueDriven,
    /// Always run before repair passes.
    Always,
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
    /// Control whether to run make-connected before normal repair passes.
    pub make_connected_prepass_mode: MakeConnectedPrepassMode,
    /// Run SameRange/SameParameter scan+fix pass as a prepass.
    pub run_parametric_consistency_prepass: bool,
    /// Re-run parametric consistency pass after each repair iteration when
    /// remaining issues still indicate parametric inconsistency.
    pub run_parametric_consistency_iterative: bool,
    /// When a repair pass makes no progress while issues remain, run a
    /// MakeConnected-style connectivity rebuild pass.
    pub run_make_connected_on_stall: bool,
    /// Base tolerance used by make-connected fallback passes.
    pub make_connected_tolerance: f64,
    /// Maximum number of iterative make-connected passes.
    pub make_connected_max_passes: usize,
    /// Per-pass tolerance growth factor for make-connected fallback.
    pub make_connected_tolerance_growth: f64,
    /// Upper cap for make-connected tolerance growth.
    pub make_connected_tolerance_cap: f64,
}

impl Default for HealingOptions {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_ABS,
            max_passes: 2,
            mode: HealingMode::AnalyzeAndRepair,
            make_connected_prepass_mode: MakeConnectedPrepassMode::Disabled,
            run_parametric_consistency_prepass: true,
            run_parametric_consistency_iterative: true,
            run_make_connected_on_stall: false,
            make_connected_tolerance: TOLERANCE_ABS,
            make_connected_max_passes: 3,
            make_connected_tolerance_growth: 1.5,
            make_connected_tolerance_cap: TOLERANCE_ABS * 1000.0,
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
                // Handle all other variants - they don't have specific counters yet
                _ => {}
            }
        }
        s
    }
}

/// Comprehensive diagnosis report combining all analysis types.
///
/// Analogous to running all ShapeAnalysis tools in OCCT:
/// - ShapeAnalysis_Surface (UV consistency)
/// - ShapeAnalysis_Wire (wire quality)
/// - ShapeAnalysis_ShapeTolerance (tolerance consistency)
/// - BRepCheck_Analyzer (topology validity)
#[derive(Debug, Clone)]
pub struct ComprehensiveDiagnosis {
    /// Basic topology check result.
    pub topology: CheckResult,
    /// Surface UV consistency analysis.
    pub surface_uv: crate::brep_check::SurfaceAnalysisReport,
    /// Wire quality analysis.
    pub wire_quality: crate::brep_check::WireQualityReport,
    /// SameParameter diagnosis.
    pub same_parameter: crate::brep_check::SameParameterDiagnosis,
    /// SameRange diagnosis.
    pub same_range: crate::brep_check::SameRangeDiagnosis,
}

impl Default for ComprehensiveDiagnosis {
    fn default() -> Self {
        Self {
            topology: CheckResult { issues: Vec::new() },
            surface_uv: Default::default(),
            wire_quality: Default::default(),
            same_parameter: Default::default(),
            same_range: Default::default(),
        }
    }
}

impl ComprehensiveDiagnosis {
    /// Returns true if all diagnoses are clean (no issues found).
    pub fn is_clean(&self) -> bool {
        self.topology.is_valid()
            && self.surface_uv.is_clean()
            && self.wire_quality.is_clean()
            && self.same_parameter.is_clean()
            && self.same_range.is_clean()
    }

    /// Returns a summary string of all issues found.
    pub fn summary(&self) -> String {
        if self.is_clean() {
            return "All diagnoses clean: topology valid, no UV issues, all wires closed, parametric consistency OK".to_string();
        }

        let mut parts = Vec::new();

        if !self.topology.is_valid() {
            parts.push(format!("topology: {} issues", self.topology.issues.len()));
        }
        if !self.surface_uv.is_clean() {
            parts.push(format!("UV bounds: {} violations", self.surface_uv.total_issues));
        }
        if !self.wire_quality.is_clean() {
            parts.push(format!(
                "wires: {} open, {} self-intersecting",
                self.wire_quality.open_wires,
                self.wire_quality.self_intersecting_wires
            ));
        }
        if !self.same_parameter.is_clean() {
            parts.push(format!("SameParameter: {} suspect", self.same_parameter.suspect_edges.len()));
        }
        if !self.same_range.is_clean() {
            parts.push(format!("SameRange: {} suspect", self.same_range.suspect_edges.len()));
        }

        parts.join("; ")
    }

    /// Returns total count of all issues found.
    pub fn total_issues(&self) -> usize {
        self.topology.issues.len()
            + self.surface_uv.total_issues
            + (if self.wire_quality.is_clean() { 0 } else { 1 })
            + self.same_parameter.suspect_edges.len()
            + self.same_range.suspect_edges.len()
    }
}

/// Run all available diagnoses on a BRep.
///
/// This is a convenience function that runs all analysis tools and returns
/// a combined report.
///
/// # Arguments
/// * `brep` - The BRep to analyze.
/// * `tolerance` - The tolerance to use for geometric comparisons.
///
/// # Returns
/// A `ComprehensiveDiagnosis` containing all analysis results.
pub fn diagnose_all(brep: &BRep, tolerance: f64) -> ComprehensiveDiagnosis {
    use crate::brep_check::{
        analyze_surface_uv_consistency, analyze_wire_quality,
        diagnose_same_parameter, diagnose_same_range,
    };

    ComprehensiveDiagnosis {
        topology: brep_check_analyze(brep),
        surface_uv: analyze_surface_uv_consistency(brep, tolerance),
        wire_quality: analyze_wire_quality(brep, tolerance),
        same_parameter: diagnose_same_parameter(brep, tolerance),
        same_range: diagnose_same_range(brep, tolerance),
    }
}

/// Summary report for analyze/heal workflow.
/// Result of a healing operation.
#[derive(Debug, Clone, Default)]
pub struct HealingReport {
    /// Issues found before any repair.
    pub initial: CheckResult,
    /// Issues after the final pass.
    pub final_result: CheckResult,
    /// Per-pass repair reports.
    pub passes: Vec<RepairReport>,
    /// Parametric consistency passes (SameRange/SameParameter scan+fix).
    pub parametric_passes: Vec<ParametricConsistencyReport>,
    /// MakeConnected fallback reports executed when repair stalls.
    pub make_connected_passes: Vec<MakeConnectedReport>,
    /// Structured issue counters before healing.
    pub initial_stats: HealingIssueStats,
    /// Structured issue counters after healing.
    pub final_stats: HealingIssueStats,
    /// Stage-by-stage issue counts for analyze/repair pipeline.
    pub stages: Vec<HealingStageReport>,
}

/// Stage marker for healing pipeline reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HealingStage {
    InitialCheck,
    PreprocessPass,
    GeometryRepairPass,
    TopologyRepairPass,
    PreMakeConnected,
    OperatorChainStep,
    ParametricConsistencyPass,
    RepairPass,
    WireGapRepairPass,
    UvBoundsRepairPass,
    MakeConnectedPass,
    FinalizePass,
    FinalCheck,
}

/// ShapeProcess-like healing operators that can be composed into a custom
/// batch pipeline.
#[derive(Debug, Clone, PartialEq)]
pub enum HealingOperator {
    /// MakeConnected-style connectivity rebuild pass.
    MakeConnected,
    /// SameRange/SameParameter consistency pass.
    ParametricConsistency,
    /// General repair pass (`repair`).
    Repair,
    /// Wire gap repair pass.
    WireGapRepair,
    /// UV bounds repair pass.
    UvBoundsRepair,
    /// Stop pipeline execution if the current shape is checker-clean.
    StopIfClean,
    /// Remove faces with area below a threshold.
    FixSmallAreaFaces,
    /// Fix sliver (thin elongated) faces by merging or removal.
    FixSliverFaces,
    /// Repair non-manifold topology by splitting multi-face edges.
    FixNonManifold,
    /// Propagate tolerances through the shape hierarchy.
    PropagateTolerances,
    /// Merge faces that share the same underlying surface geometry.
    UnifySameDomain,
    /// Remove internal faces (faces inside the solid volume).
    RemoveInternalFaces,
    /// Split cylindrical faces at angle thresholds.
    /// Useful for meshing constraints where element size limits angular extent.
    SplitAngle(SplitAngleOperator),
    /// Split edges at continuity breaks (C0/C1/C2 discontinuities).
    SplitContinuity(SplitContinuityOperator),
    /// Convert analytic geometry to BSpline representation.
    ConvertToBSpline(ConvertToBSplineOperator),
    /// Convert BSpline surfaces to Bezier patches by splitting at knot lines.
    SurfaceToBezier(SurfaceToBezierOperator),
    /// Apply uniform or non-uniform scaling transformation.
    ScaleShape(ScaleShapeOperator),
    /// Convert indirect faces to direct (fix face orientation issues).
    DirectFaces(DirectFacesOperator),
    /// Fix SameParameter issues on edges.
    SameParameter(SameParameterOperator),
    /// Remove internal faces after boolean operations.
    RemoveInternalFacesOp(RemoveInternalFacesOperator),
    /// Comprehensive geometry healing combining multiple operations.
    HealGeometry(HealGeometryOperator),
}

/// Operator that splits faces at specified angle thresholds.
///
/// This is particularly useful for:
/// - Cylindrical faces: split into sectors for meshing constraints
/// - Torus faces: split at major/minor angle limits
/// - Conical faces: split into angular sectors
///
/// Analogous to OCCT `ShapeUpgrade_ShapeDivideAngle`.
#[derive(Debug, Clone, PartialEq)]
pub struct SplitAngleOperator {
    /// Maximum angular span in radians for any resulting face sector.
    pub max_angle: f64,
    /// Whether to split cylindrical faces.
    pub split_cylinders: bool,
    /// Whether to split torus faces.
    pub split_tori: bool,
    /// Whether to split conical faces.
    pub split_cones: bool,
    /// Whether to split spherical faces.
    pub split_spheres: bool,
    /// Starting angle offset in radians (for alignment with specific directions).
    pub start_angle: f64,
}

impl Default for SplitAngleOperator {
    fn default() -> Self {
        Self {
            max_angle: std::f64::consts::PI / 2.0, // 90 degrees default
            split_cylinders: true,
            split_tori: true,
            split_cones: true,
            split_spheres: true,
            start_angle: 0.0,
        }
    }
}

/// Operator that splits edges at continuity breaks.
///
/// Detects C0 (position), C1 (tangent), and C2 (curvature) discontinuities
/// and splits edges at those points. This is essential for downstream
/// operations that require specific continuity levels.
///
/// Analogous to OCCT `ShapeUpgrade_ShapeDivideContinuity`.
#[derive(Debug, Clone, PartialEq)]
pub struct SplitContinuityOperator {
    /// Minimum continuity level required (C0, C1, or C2).
    pub min_continuity: ContinuityLevel,
    /// Tolerance for detecting discontinuities.
    pub tolerance: f64,
    /// Whether to check curve continuity.
    pub check_curves: bool,
    /// Whether to check surface continuity at edges.
    pub check_surfaces: bool,
    /// Maximum number of split points per edge.
    pub max_splits_per_edge: usize,
}

/// Continuity level for geometric analysis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default)]
pub enum ContinuityLevel {
    /// C0: position continuous only.
    C0,
    /// C1: tangent continuous (default).
    #[default]
    C1,
    /// C2: curvature continuous.
    C2,
}

impl Default for SplitContinuityOperator {
    fn default() -> Self {
        Self {
            min_continuity: ContinuityLevel::C1,
            tolerance: TOLERANCE_MESH_LEGACY,
            check_curves: true,
            check_surfaces: true,
            max_splits_per_edge: 100,
        }
    }
}

/// Operator that converts analytic geometry to BSpline representation.
///
/// Converts planes, cylinders, spheres, cones, tori, and other analytic
/// surfaces to NURBS (BSpline) representation. This is useful for:
/// - Exporting to formats that only support NURBS
/// - Applying NURBS-specific operations
/// - Ensuring uniform representation for downstream algorithms
///
/// Analogous to OCCT `ShapeUpgrade_ShapeConvertToBSpline`.
#[derive(Debug, Clone, PartialEq)]
pub struct ConvertToBSplineOperator {
    /// Maximum degree for resulting BSpline geometry.
    pub max_degree: usize,
    /// Whether to convert curves to BSpline.
    pub convert_curves: bool,
    /// Whether to convert surfaces to BSpline.
    pub convert_surfaces: bool,
    /// Whether to convert planes (usually kept as analytic).
    pub convert_planes: bool,
    /// Whether to convert elementary surfaces (cylinders, spheres, cones, tori).
    pub convert_elementary: bool,
    /// Number of samples for approximating transcendental surfaces.
    pub approximation_samples: usize,
}

impl Default for ConvertToBSplineOperator {
    fn default() -> Self {
        Self {
            max_degree: 3,
            convert_curves: true,
            convert_surfaces: true,
            convert_planes: false,
            convert_elementary: true,
            approximation_samples: 20,
        }
    }
}

/// Operator that converts BSpline surfaces to Bezier patches.
///
/// Splits BSpline surfaces at all interior knot lines, converting each
/// span into a separate Bezier patch. This is useful for:
/// - Export to formats requiring Bezier patches
/// - Isogeometric analysis workflows
/// - Simplified surface representation
///
/// Analogous to OCCT `ShapeUpgrade_ShapeConvertToBezier`.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceToBezierOperator {
    /// Whether to convert surfaces.
    pub convert_surfaces: bool,
    /// Whether to convert 2D curves (PCurves).
    pub convert_pcurves: bool,
    /// Whether to convert 3D curves.
    pub convert_curves: bool,
    /// Maximum degree for resulting Bezier patches.
    pub max_degree: usize,
}

impl Default for SurfaceToBezierOperator {
    fn default() -> Self {
        Self {
            convert_surfaces: true,
            convert_pcurves: true,
            convert_curves: true,
            max_degree: 25, // High degree allowed for exact conversion
        }
    }
}

/// Operator that applies scaling transformation to a shape.
///
/// Supports both uniform scaling (same factor in all directions) and
/// non-uniform scaling (different factors for X, Y, Z axes).
///
/// Tolerances are scaled appropriately to maintain geometric validity.
///
/// Analogous to OCCT `BRepBuilderAPI_GTransform` for scaling.
#[derive(Debug, Clone, PartialEq)]
pub struct ScaleShapeOperator {
    /// Scale factor for X direction.
    pub scale_x: f64,
    /// Scale factor for Y direction.
    pub scale_y: f64,
    /// Scale factor for Z direction.
    pub scale_z: f64,
    /// Origin point for scaling (default is origin).
    pub origin: Option<glam::DVec3>,
    /// Whether to scale tolerances.
    pub scale_tolerances: bool,
    /// Whether to preserve vertex tolerances on degenerate scaling.
    pub preserve_degenerate_tolerances: bool,
}

impl Default for ScaleShapeOperator {
    fn default() -> Self {
        Self {
            scale_x: 1.0,
            scale_y: 1.0,
            scale_z: 1.0,
            origin: None,
            scale_tolerances: true,
            preserve_degenerate_tolerances: true,
        }
    }
}

impl ScaleShapeOperator {
    /// Create a uniform scaling operator.
    pub fn uniform(scale: f64) -> Self {
        Self {
            scale_x: scale,
            scale_y: scale,
            scale_z: scale,
            ..Default::default()
        }
    }

    /// Create a non-uniform scaling operator.
    pub fn non_uniform(scale_x: f64, scale_y: f64, scale_z: f64) -> Self {
        Self {
            scale_x,
            scale_y,
            scale_z,
            ..Default::default()
        }
    }

    /// Returns true if the scaling is uniform (same factor in all directions).
    pub fn is_uniform(&self) -> bool {
        (self.scale_x - self.scale_y).abs() < TOLERANCE_LEN_MIN
            && (self.scale_y - self.scale_z).abs() < TOLERANCE_LEN_MIN
    }
}

/// Operator that converts indirect faces to direct.
///
/// An indirect face is one where the face orientation does not match
/// the natural surface orientation. This operator ensures all faces
/// are "direct" by correcting orientation flags and surface references.
///
/// This is analogous to OCCT `ShapeFix_Face::FixOrientation` combined
/// with surface orientation adjustments.
#[derive(Debug, Clone, PartialEq)]
pub struct DirectFacesOperator {
    /// Tolerance for geometric comparisons.
    pub tolerance: f64,
    /// Whether to update surface references when fixing orientation.
    pub update_surface_references: bool,
    /// Whether to recompute face normals after orientation fix.
    pub recompute_normals: bool,
    /// Whether to also fix wire orientation on the face.
    pub fix_wire_orientation: bool,
}

impl Default for DirectFacesOperator {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_ABS,
            update_surface_references: true,
            recompute_normals: true,
            fix_wire_orientation: true,
        }
    }
}

impl DirectFacesOperator {
    /// Create a new DirectFacesOperator with specified tolerance.
    pub fn new(tolerance: f64) -> Self {
        Self {
            tolerance,
            ..Default::default()
        }
    }
}

/// Operator that fixes SameParameter issues on edges.
///
/// SameParameter ensures that the 3D curve and 2D PCurves of an edge
/// are parameterized consistently. When violated, the edge's geometry
/// may not match the face's surface geometry at the same parameter value.
///
/// This operator uses the existing `fix_same_parameter_with_scan` function
/// but adds configurable tolerance and additional options.
///
/// Analogous to OCCT `BRepLib::SameParameter` and `ShapeFix_Edge::FixSameParameter`.
#[derive(Debug, Clone, PartialEq)]
pub struct SameParameterOperator {
    /// Tolerance for SameParameter diagnosis and repair.
    pub tolerance: f64,
    /// Maximum number of sampling points for curve comparison.
    pub max_samples: usize,
    /// Whether to enforce SameParameter even on already-flagged edges.
    pub enforce: bool,
    /// Whether to also update PCurve ranges to match 3D curve range.
    pub update_pcurve_ranges: bool,
}

impl Default for SameParameterOperator {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_ABS,
            max_samples: 23,
            enforce: false,
            update_pcurve_ranges: true,
        }
    }
}

impl SameParameterOperator {
    /// Create a new SameParameterOperator with specified tolerance.
    pub fn new(tolerance: f64) -> Self {
        Self {
            tolerance,
            ..Default::default()
        }
    }

    /// Create an enforcing SameParameterOperator that repairs all edges.
    pub fn enforced(tolerance: f64) -> Self {
        Self {
            tolerance,
            enforce: true,
            ..Default::default()
        }
    }
}

/// Operator that removes internal faces after boolean operations.
///
/// Internal faces are partition faces that are completely inside the
/// resulting solid volume after a boolean operation. These faces do not
/// contribute to the outer boundary and should be removed for a clean result.
///
/// This operator detects internal faces by analyzing material sides and
/// connectivity, then removes them while maintaining valid topology.
///
/// Analogous to OCCT `ShapeFix_Shape::FixRemoveInternalFaces` and related
/// post-boolean cleanup operations.
#[derive(Debug, Clone, PartialEq)]
pub struct RemoveInternalFacesOperator {
    /// Tolerance for geometric operations.
    pub tolerance: f64,
    /// Minimum face area threshold (faces below this are candidates for removal).
    pub min_face_area: f64,
    /// Whether to check for manifold connectivity before removal.
    pub check_manifold: bool,
    /// Whether to merge vertices after face removal.
    pub merge_vertices: bool,
    /// Whether to preserve faces that separate distinct material regions.
    pub preserve_material_boundaries: bool,
}

impl Default for RemoveInternalFacesOperator {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_ABS,
            min_face_area: TOLERANCE_LINEAR_ULTRA_STRICT,
            check_manifold: true,
            merge_vertices: true,
            preserve_material_boundaries: true,
        }
    }
}

impl RemoveInternalFacesOperator {
    /// Create a new RemoveInternalFacesOperator with specified tolerance.
    pub fn new(tolerance: f64) -> Self {
        Self {
            tolerance,
            ..Default::default()
        }
    }
}

/// Operator that performs comprehensive geometry healing.
///
/// This operator combines multiple repair operations into a single,
/// configurable healing pass. It can perform:
/// - Face orientation fixes
/// - SameParameter/SameRange repairs
/// - Wire closure verification
/// - Degenerate geometry removal
/// - Tolerance propagation
///
/// The sequence of operations is configurable, allowing customization
/// for different use cases (import cleanup, boolean post-processing, etc.).
///
/// Analogous to OCCT `ShapeFix_Shape` which orchestrates multiple
/// ShapeFix operations in a configurable sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct HealGeometryOperator {
    /// Tolerance for all geometric operations.
    pub tolerance: f64,
    /// Maximum number of healing passes.
    pub max_passes: usize,
    /// Whether to fix face orientation.
    pub fix_face_orientation: bool,
    /// Whether to fix SameParameter issues.
    pub fix_same_parameter: bool,
    /// Whether to fix SameRange issues.
    pub fix_same_range: bool,
    /// Whether to fix wire gaps.
    pub fix_wire_gaps: bool,
    /// Whether to remove degenerate faces.
    pub remove_degenerate_faces: bool,
    /// Whether to propagate tolerances.
    pub propagate_tolerances: bool,
    /// Whether to recompute face normals.
    pub recompute_normals: bool,
    /// Whether to fix UV bounds violations.
    pub fix_uv_bounds: bool,
    /// Whether to remove small edges.
    pub remove_small_edges: bool,
    /// Minimum edge length threshold for removal.
    pub min_edge_length: f64,
    /// Custom sequence of operations (if empty, uses default order).
    pub custom_sequence: Vec<HealGeometryStep>,
}

/// Step in the HealGeometry operator sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HealGeometryStep {
    /// Fix face orientation.
    FixFaceOrientation,
    /// Fix SameParameter issues.
    FixSameParameter,
    /// Fix SameRange issues.
    FixSameRange,
    /// Fix wire gaps.
    FixWireGaps,
    /// Remove degenerate faces.
    RemoveDegenerateFaces,
    /// Propagate tolerances.
    PropagateTolerances,
    /// Recompute face normals.
    RecomputeNormals,
    /// Fix UV bounds violations.
    FixUvBounds,
    /// Remove small edges.
    RemoveSmallEdges,
}

impl Default for HealGeometryOperator {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_ABS,
            max_passes: 3,
            fix_face_orientation: true,
            fix_same_parameter: true,
            fix_same_range: true,
            fix_wire_gaps: true,
            remove_degenerate_faces: true,
            propagate_tolerances: true,
            recompute_normals: true,
            fix_uv_bounds: true,
            remove_small_edges: false,
            min_edge_length: TOLERANCE_MESH_LEGACY,
            custom_sequence: Vec::new(),
        }
    }
}

impl HealGeometryOperator {
    /// Create a new HealGeometryOperator with specified tolerance.
    pub fn new(tolerance: f64) -> Self {
        Self {
            tolerance,
            ..Default::default()
        }
    }

    /// Create a minimal HealGeometryOperator for quick fixes.
    pub fn minimal(tolerance: f64) -> Self {
        Self {
            tolerance,
            max_passes: 1,
            fix_face_orientation: true,
            fix_same_parameter: true,
            fix_same_range: true,
            fix_wire_gaps: false,
            remove_degenerate_faces: false,
            propagate_tolerances: false,
            recompute_normals: true,
            fix_uv_bounds: false,
            remove_small_edges: false,
            min_edge_length: TOLERANCE_MESH_LEGACY,
            custom_sequence: Vec::new(),
        }
    }

    /// Create an aggressive HealGeometryOperator for thorough cleanup.
    pub fn aggressive(tolerance: f64) -> Self {
        Self {
            tolerance,
            max_passes: 5,
            fix_face_orientation: true,
            fix_same_parameter: true,
            fix_same_range: true,
            fix_wire_gaps: true,
            remove_degenerate_faces: true,
            propagate_tolerances: true,
            recompute_normals: true,
            fix_uv_bounds: true,
            remove_small_edges: true,
            min_edge_length: tolerance,
            custom_sequence: Vec::new(),
        }
    }

    /// Get the sequence of healing steps to execute.
    pub fn get_sequence(&self) -> Vec<HealGeometryStep> {
        if !self.custom_sequence.is_empty() {
            return self.custom_sequence.clone();
        }

        let mut steps = Vec::new();
        if self.recompute_normals {
            steps.push(HealGeometryStep::RecomputeNormals);
        }
        if self.fix_same_range {
            steps.push(HealGeometryStep::FixSameRange);
        }
        if self.fix_same_parameter {
            steps.push(HealGeometryStep::FixSameParameter);
        }
        if self.fix_face_orientation {
            steps.push(HealGeometryStep::FixFaceOrientation);
        }
        if self.fix_wire_gaps {
            steps.push(HealGeometryStep::FixWireGaps);
        }
        if self.fix_uv_bounds {
            steps.push(HealGeometryStep::FixUvBounds);
        }
        if self.remove_degenerate_faces {
            steps.push(HealGeometryStep::RemoveDegenerateFaces);
        }
        if self.remove_small_edges {
            steps.push(HealGeometryStep::RemoveSmallEdges);
        }
        if self.propagate_tolerances {
            steps.push(HealGeometryStep::PropagateTolerances);
        }
        steps
    }
}

/// Configuration parameters for individual healing operators.
#[derive(Debug, Clone)]
pub struct OperatorParams {
    /// Tolerance threshold for geometric operations.
    pub tolerance: f64,
    /// Area threshold for FixSmallAreaFaces.
    pub min_face_area: f64,
    /// Aspect ratio threshold for FixSliverFaces.
    pub max_sliver_aspect_ratio: f64,
    /// Whether to allow removal of internal faces.
    pub allow_internal_face_removal: bool,
    /// Parameters for SplitAngle operator.
    pub split_angle: SplitAngleOperator,
    /// Parameters for SplitContinuity operator.
    pub split_continuity: SplitContinuityOperator,
    /// Parameters for ConvertToBSpline operator.
    pub convert_to_bspline: ConvertToBSplineOperator,
    /// Parameters for SurfaceToBezier operator.
    pub surface_to_bezier: SurfaceToBezierOperator,
    /// Parameters for ScaleShape operator.
    pub scale_shape: ScaleShapeOperator,
    /// Parameters for DirectFaces operator.
    pub direct_faces: DirectFacesOperator,
    /// Parameters for SameParameter operator.
    pub same_parameter: SameParameterOperator,
    /// Parameters for RemoveInternalFaces operator.
    pub remove_internal_faces: RemoveInternalFacesOperator,
    /// Parameters for HealGeometry operator.
    pub heal_geometry: HealGeometryOperator,
}

impl Default for OperatorParams {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_ABS,
            min_face_area: TOLERANCE_LINEAR_ULTRA_STRICT,
            max_sliver_aspect_ratio: 100.0,
            allow_internal_face_removal: true,
            split_angle: SplitAngleOperator::default(),
            split_continuity: SplitContinuityOperator::default(),
            convert_to_bspline: ConvertToBSplineOperator::default(),
            surface_to_bezier: SurfaceToBezierOperator::default(),
            scale_shape: ScaleShapeOperator::default(),
            direct_faces: DirectFacesOperator::default(),
            same_parameter: SameParameterOperator::default(),
            remove_internal_faces: RemoveInternalFacesOperator::default(),
            heal_geometry: HealGeometryOperator::default(),
        }
    }
}

/// Report for one SameRange/SameParameter consistency pass.
#[derive(Debug, Clone, Default)]
pub struct ParametricConsistencyReport {
    pub same_range_fixed: usize,
    pub same_parameter_fixed: usize,
}

/// Report for a single healing operator execution.
#[derive(Debug, Clone)]
pub struct OperatorReport {
    /// The operator that was executed.
    pub operator: HealingOperator,
    /// Number of entities modified/removed.
    pub modifications: usize,
    /// Number of issues fixed by this operator.
    pub issues_fixed: usize,
    /// Whether the operator made any changes.
    pub changed: bool,
    /// Human-readable description of changes.
    pub description: String,
}

impl Default for OperatorReport {
    fn default() -> Self {
        Self {
            operator: HealingOperator::Repair,
            modifications: 0,
            issues_fixed: 0,
            changed: false,
            description: String::new(),
        }
    }
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Operator Result Aggregation, Rollback, and Progress Callbacks
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Aggregated results from a healing pipeline execution.
///
/// This struct collects results from multiple operator executions and
/// provides summary statistics and analysis capabilities.
#[derive(Debug, Clone, Default)]
pub struct OperatorResultAggregation {
    /// Individual operator results.
    pub results: Vec<OperatorResult>,
    /// Total number of operators executed (not skipped).
    pub total_executed: usize,
    /// Total number of operators skipped.
    pub total_skipped: usize,
    /// Total number of modifications across all operators.
    pub total_modifications: usize,
    /// Total number of issues fixed across all operators.
    pub total_issues_fixed: usize,
    /// Total execution time in seconds.
    pub total_elapsed_seconds: f64,
    /// Number of operators that made changes.
    pub operators_with_changes: usize,
    /// Number of operators that failed.
    pub operators_failed: usize,
    /// Whether rollback was triggered.
    pub rollback_triggered: bool,
    /// Reason for rollback (if triggered).
    pub rollback_reason: Option<String>,
}

impl OperatorResultAggregation {
    /// Create a new empty aggregation.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an operator result to the aggregation.
    pub fn add_result(&mut self, result: OperatorResult) {
        if result.skipped {
            self.total_skipped += 1;
        } else {
            self.total_executed += 1;
            self.total_modifications += result.modifications;
            self.total_issues_fixed += result.issues_fixed;
            self.total_elapsed_seconds += result.elapsed_seconds;
            if result.changed {
                self.operators_with_changes += 1;
            }
        }
        self.results.push(result);
    }

    /// Check if any operator made changes.
    pub fn has_changes(&self) -> bool {
        self.operators_with_changes > 0
    }

    /// Get success rate (executed operators that made changes).
    pub fn change_rate(&self) -> f64 {
        if self.total_executed == 0 {
            return 0.0;
        }
        self.operators_with_changes as f64 / self.total_executed as f64
    }

    /// Get the result for a specific operator index.
    pub fn get_result(&self, idx: usize) -> Option<&OperatorResult> {
        self.results.get(idx)
    }

    /// Find operators that made changes.
    pub fn operators_with_changes_iter(&self) -> impl Iterator<Item = &OperatorResult> {
        self.results.iter().filter(|r| r.changed)
    }

    /// Generate a summary string.
    pub fn summary(&self) -> String {
        if self.results.is_empty() {
            return "No operators executed".to_string();
        }

        let mut parts = Vec::new();
        parts.push(format!("{} executed", self.total_executed));
        if self.total_skipped > 0 {
            parts.push(format!("{} skipped", self.total_skipped));
        }
        parts.push(format!("{} modifications", self.total_modifications));
        parts.push(format!("{} issues fixed", self.total_issues_fixed));
        parts.push(format!("{:.3}s", self.total_elapsed_seconds));

        if self.rollback_triggered {
            parts.push("ROLLBACK".to_string());
        }

        parts.join(", ")
    }
}

/// Snapshot of BRep state for potential rollback.
///
/// This struct stores a clone of the BRep at a specific point in the
/// operator pipeline, allowing rollback to that state if needed.
#[derive(Debug, Clone)]
pub struct BRepSnapshot {
    /// The BRep state.
    pub brep: BRep,
    /// Operator index at which this snapshot was taken.
    pub operator_index: usize,
    /// Label for this snapshot.
    pub label: String,
    /// Timestamp when snapshot was created.
    pub timestamp_seconds: f64,
}

impl BRepSnapshot {
    /// Create a new snapshot.
    pub fn new(brep: &BRep, operator_index: usize, label: impl Into<String>, elapsed_seconds: f64) -> Self {
        Self {
            brep: brep.clone(),
            operator_index,
            label: label.into(),
            timestamp_seconds: elapsed_seconds,
        }
    }
}

/// Configuration for rollback behavior.
#[derive(Debug, Clone)]
pub struct RollbackConfig {
    /// Whether rollback is enabled.
    pub enabled: bool,
    /// Maximum number of issues that trigger rollback (0 = no auto-rollback).
    pub max_issues_threshold: usize,
    /// Whether to rollback on operator failure.
    pub rollback_on_failure: bool,
    /// Whether to rollback if issue count increases.
    pub rollback_on_regression: bool,
    /// Operator indices at which to create snapshots (for potential rollback).
    pub snapshot_indices: Vec<usize>,
    /// Whether to create snapshots before each operator.
    pub snapshot_before_each: bool,
    /// Maximum number of snapshots to retain.
    pub max_snapshots: usize,
}

impl Default for RollbackConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            max_issues_threshold: 0,
            rollback_on_failure: true,
            rollback_on_regression: true,
            snapshot_indices: Vec::new(),
            snapshot_before_each: false,
            max_snapshots: 10,
        }
    }
}

impl RollbackConfig {
    /// Create a rollback config that never rolls back.
    pub fn disabled() -> Self {
        Self {
            enabled: false,
            ..Default::default()
        }
    }

    /// Create a rollback config that snapshots at specific indices.
    pub fn with_snapshots(indices: Vec<usize>) -> Self {
        Self {
            snapshot_indices: indices,
            ..Default::default()
        }
    }
}

/// Progress callback for operator execution.
///
/// This trait allows external code to monitor the progress of a healing
/// pipeline execution and potentially cancel it.
pub trait ProgressCallback: Send + Sync {
    /// Called before an operator is executed.
    fn on_operator_start(&self, operator_index: usize, operator: &HealingOperator);

    /// Called after an operator completes.
    fn on_operator_complete(&self, operator_index: usize, result: &OperatorResult);

    /// Called when progress is made (0.0 to 1.0).
    fn on_progress(&self, progress: f64, message: &str);

    /// Called when an error occurs.
    fn on_error(&self, operator_index: usize, error: &str);

    /// Check if execution should be cancelled.
    fn is_cancelled(&self) -> bool;
}

/// A simple progress callback that tracks execution state.
#[derive(Debug, Default)]
pub struct SimpleProgressCallback {
    /// Current operator index.
    pub current_operator: usize,
    /// Total number of operators.
    pub total_operators: usize,
    /// Whether cancellation was requested.
    pub cancelled: std::sync::atomic::AtomicBool,
    /// Last progress message.
    pub last_message: String,
}

impl Clone for SimpleProgressCallback {
    fn clone(&self) -> Self {
        Self {
            current_operator: self.current_operator,
            total_operators: self.total_operators,
            cancelled: std::sync::atomic::AtomicBool::new(
                self.cancelled.load(std::sync::atomic::Ordering::SeqCst)
            ),
            last_message: self.last_message.clone(),
        }
    }
}

impl SimpleProgressCallback {
    /// Create a new simple progress callback.
    pub fn new(total_operators: usize) -> Self {
        Self {
            total_operators,
            ..Default::default()
        }
    }

    /// Request cancellation.
    pub fn cancel(&self) {
        self.cancelled.store(true, std::sync::atomic::Ordering::SeqCst);
    }

    /// Get progress as a fraction (0.0 to 1.0).
    pub fn progress(&self) -> f64 {
        if self.total_operators == 0 {
            return 1.0;
        }
        self.current_operator as f64 / self.total_operators as f64
    }
}

impl ProgressCallback for SimpleProgressCallback {
    fn on_operator_start(&self, operator_index: usize, _operator: &HealingOperator) {
        // Note: In a single-threaded context, we can't mutate, but this is for demonstration
        // In practice, this would use interior mutability (e.g., Mutex)
        let _ = operator_index;
    }

    fn on_operator_complete(&self, operator_index: usize, _result: &OperatorResult) {
        let _ = operator_index;
    }

    fn on_progress(&self, progress: f64, message: &str) {
        let _ = (progress, message);
    }

    fn on_error(&self, operator_index: usize, error: &str) {
        let _ = (operator_index, error);
    }

    fn is_cancelled(&self) -> bool {
        self.cancelled.load(std::sync::atomic::Ordering::SeqCst)
    }
}

/// Report from an operator pipeline execution with rollback support.
#[derive(Debug, Clone)]
pub struct PipelineExecutionReport {
    /// Aggregated results.
    pub aggregation: OperatorResultAggregation,
    /// Snapshots taken during execution.
    pub snapshots: Vec<BRepSnapshot>,
    /// Final BRep state.
    pub final_brep: BRep,
    /// Whether the pipeline completed successfully.
    pub completed: bool,
    /// Reason for failure (if not completed).
    pub failure_reason: Option<String>,
    /// Index to which rollback occurred (if any).
    pub rollback_index: Option<usize>,
}

impl PipelineExecutionReport {
    /// Check if the pipeline made any changes.
    pub fn has_changes(&self) -> bool {
        self.aggregation.has_changes()
    }

    /// Get a snapshot by index.
    pub fn get_snapshot(&self, operator_index: usize) -> Option<&BRepSnapshot> {
        self.snapshots.iter().find(|s| s.operator_index == operator_index)
    }

    /// Generate a summary.
    pub fn summary(&self) -> String {
        let status = if self.completed {
            "Completed"
        } else if let Some(ref reason) = self.failure_reason {
            reason
        } else {
            "Unknown status"
        };

        let rollback_info = if let Some(idx) = self.rollback_index {
            format!(" (rolled back to operator {})", idx)
        } else {
            String::new()
        };

        format!("{}: {}{}", status, self.aggregation.summary(), rollback_info)
    }
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Operator Chaining Improvements
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Condition for conditional operator execution.
#[derive(Debug, Clone, PartialEq)]
pub enum OperatorCondition {
    /// Always execute the operator.
    Always,
    /// Execute only if the shape has checker issues.
    OnlyIfIssues,
    /// Execute only if the shape is checker-clean.
    OnlyIfClean,
    /// Execute only if a specific issue type is present.
    OnlyIfIssueType(CheckIssuePredicate),
    /// Execute only if a previous operator made changes.
    OnlyIfPreviousChanged(usize),
    /// Execute only if a previous operator did NOT make changes.
    OnlyIfPreviousUnchanged(usize),
    /// Execute only if the number of issues exceeds a threshold.
    OnlyIfIssueCountAbove(usize),
    /// Execute only if the number of issues is below a threshold.
    OnlyIfIssueCountBelow(usize),
}

/// Predicate for checking specific issue types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckIssuePredicate {
    /// Any open wire issue.
    OpenWire,
    /// Any zero normal issue.
    ZeroNormal,
    /// Any degenerate face issue.
    DegenerateFace,
    /// Any non-manifold edge issue.
    NonManifoldEdge,
    /// Any self-intersection issue.
    SelfIntersection,
    /// Any geometric self-intersection.
    GeometricSelfIntersection,
}

impl OperatorCondition {
    /// Evaluate whether the condition is met.
    pub fn evaluate(&self, _brep: &BRep, report: &HealingReport, previous_results: &[OperatorResult]) -> bool {
        match self {
            OperatorCondition::Always => true,
            OperatorCondition::OnlyIfIssues => !report.final_result.is_valid(),
            OperatorCondition::OnlyIfClean => report.final_result.is_valid(),
            OperatorCondition::OnlyIfIssueType(pred) => {
                report.final_result.issues.iter().any(|issue| pred.matches(issue))
            }
            OperatorCondition::OnlyIfPreviousChanged(idx) => {
                previous_results.get(*idx).map(|r| r.changed).unwrap_or(false)
            }
            OperatorCondition::OnlyIfPreviousUnchanged(idx) => {
                previous_results.get(*idx).map(|r| !r.changed).unwrap_or(true)
            }
            OperatorCondition::OnlyIfIssueCountAbove(threshold) => {
                report.final_result.issues.len() > *threshold
            }
            OperatorCondition::OnlyIfIssueCountBelow(threshold) => {
                report.final_result.issues.len() < *threshold
            }
        }
    }
}

impl CheckIssuePredicate {
    fn matches(&self, issue: &CheckIssue) -> bool {
        match self {
            CheckIssuePredicate::OpenWire => matches!(issue, CheckIssue::OpenWire { .. }),
            CheckIssuePredicate::ZeroNormal => matches!(issue, CheckIssue::ZeroNormal { .. }),
            CheckIssuePredicate::DegenerateFace => matches!(issue, CheckIssue::DegenerateFace { .. }),
            CheckIssuePredicate::NonManifoldEdge => matches!(issue, CheckIssue::NonManifoldEdge { .. }),
            CheckIssuePredicate::SelfIntersection => matches!(issue, CheckIssue::SelfIntersectingWire { .. }),
            CheckIssuePredicate::GeometricSelfIntersection => matches!(issue, CheckIssue::GeometricSelfIntersection { .. }),
        }
    }
}

/// Result from executing a single operator in a chain.
#[derive(Debug, Clone)]
pub struct OperatorResult {
    /// The operator that was executed.
    pub operator: HealingOperator,
    /// Whether the operator made any changes.
    pub changed: bool,
    /// Number of entities modified.
    pub modifications: usize,
    /// Number of issues fixed.
    pub issues_fixed: usize,
    /// Description of changes.
    pub description: String,
    /// Execution time in seconds.
    pub elapsed_seconds: f64,
    /// Whether the operator was skipped due to a condition.
    pub skipped: bool,
    /// Reason for skipping (if skipped).
    pub skip_reason: Option<String>,
}

impl Default for OperatorResult {
    fn default() -> Self {
        Self {
            operator: HealingOperator::Repair,
            changed: false,
            modifications: 0,
            issues_fixed: 0,
            description: String::new(),
            elapsed_seconds: 0.0,
            skipped: false,
            skip_reason: None,
        }
    }
}

/// An operator with optional execution conditions and dependencies.
#[derive(Debug, Clone)]
pub struct HealingOperatorWithCondition {
    /// The operator to execute.
    pub operator: HealingOperator,
    /// Optional condition for execution.
    pub condition: Option<OperatorCondition>,
    /// Dependencies on other operators (indices in the chain).
    pub dependencies: Vec<usize>,
    /// Whether to skip this operator if dependencies failed.
    pub skip_on_dependency_failure: bool,
    /// Optional label for debugging/logging.
    pub label: Option<String>,
}

impl HealingOperatorWithCondition {
    /// Create a new operator that always executes.
    pub fn new(operator: HealingOperator) -> Self {
        Self {
            operator,
            condition: None,
            dependencies: Vec::new(),
            skip_on_dependency_failure: true,
            label: None,
        }
    }

    /// Create an operator with a condition.
    pub fn with_condition(operator: HealingOperator, condition: OperatorCondition) -> Self {
        Self {
            operator,
            condition: Some(condition),
            dependencies: Vec::new(),
            skip_on_dependency_failure: true,
            label: None,
        }
    }

    /// Add a dependency on another operator.
    pub fn depends_on(mut self, operator_idx: usize) -> Self {
        self.dependencies.push(operator_idx);
        self
    }

    /// Set the label for this operator.
    pub fn with_label(mut self, label: impl Into<String>) -> Self {
        self.label = Some(label.into());
        self
    }
}

impl From<HealingOperator> for HealingOperatorWithCondition {
    fn from(operator: HealingOperator) -> Self {
        Self::new(operator)
    }
}

/// Configuration for advanced operator chaining.
#[derive(Debug, Clone)]
pub struct OperatorChainConfig {
    /// Operators with conditions and dependencies.
    pub operators: Vec<HealingOperatorWithCondition>,
    /// Stop processing if the shape becomes checker-clean.
    pub stop_on_clean: bool,
    /// Maximum number of iterations (0 = run once).
    pub max_iterations: usize,
    /// Base tolerance for operations.
    pub base_tolerance: f64,
    /// Tolerance growth factor per iteration.
    pub tolerance_growth: f64,
    /// Maximum tolerance cap.
    pub tolerance_cap: f64,
    /// Healing options for internal passes.
    pub healing_options: HealingOptions,
    /// Operator parameters.
    pub operator_params: OperatorParams,
    /// Whether to collect detailed timing information.
    pub collect_timing: bool,
}

impl Default for OperatorChainConfig {
    fn default() -> Self {
        Self {
            operators: vec![
                HealingOperatorWithCondition::new(HealingOperator::ParametricConsistency),
                HealingOperatorWithCondition::new(HealingOperator::Repair),
                HealingOperatorWithCondition::new(HealingOperator::StopIfClean),
            ],
            stop_on_clean: true,
            max_iterations: 1,
            base_tolerance: TOLERANCE_ABS,
            tolerance_growth: 1.0,
            tolerance_cap: TOLERANCE_ADAPTIVE_MAX,
            healing_options: HealingOptions::default(),
            operator_params: OperatorParams::default(),
            collect_timing: true,
        }
    }
}

impl OperatorChainConfig {
    /// Create a preset for mesh preparation (split angles, convert to BSpline).
    pub fn mesh_prep_preset() -> Self {
        Self {
            operators: vec![
                HealingOperatorWithCondition::new(HealingOperator::SplitAngle(SplitAngleOperator {
                    max_angle: std::f64::consts::PI / 4.0, // 45 degrees
                    ..Default::default()
                })),
                HealingOperatorWithCondition::new(HealingOperator::ConvertToBSpline(ConvertToBSplineOperator {
                    convert_elementary: true,
                    convert_planes: false,
                    ..Default::default()
                })),
                HealingOperatorWithCondition::new(HealingOperator::ParametricConsistency),
                HealingOperatorWithCondition::new(HealingOperator::Repair),
            ],
            stop_on_clean: true,
            max_iterations: 2,
            base_tolerance: TOLERANCE_ABS,
            tolerance_growth: 1.5,
            tolerance_cap: TOLERANCE_ADAPTIVE_MAX,
            healing_options: HealingOptions::default(),
            operator_params: OperatorParams::default(),
            collect_timing: true,
        }
    }

    /// Create a preset for export preparation (convert to Bezier).
    pub fn export_prep_preset() -> Self {
        Self {
            operators: vec![
                HealingOperatorWithCondition::new(HealingOperator::ConvertToBSpline(ConvertToBSplineOperator::default())),
                HealingOperatorWithCondition::new(HealingOperator::SurfaceToBezier(SurfaceToBezierOperator::default())),
                HealingOperatorWithCondition::new(HealingOperator::ParametricConsistency),
                HealingOperatorWithCondition::new(HealingOperator::Repair),
            ],
            stop_on_clean: true,
            max_iterations: 1,
            base_tolerance: TOLERANCE_ABS,
            tolerance_growth: 1.0,
            tolerance_cap: TOLERANCE_ADAPTIVE_MAX,
            healing_options: HealingOptions::default(),
            operator_params: OperatorParams::default(),
            collect_timing: true,
        }
    }

    /// Create a preset for scaling operations.
    pub fn scale_preset(scale: f64) -> Self {
        Self {
            operators: vec![
                HealingOperatorWithCondition::new(HealingOperator::ScaleShape(ScaleShapeOperator::uniform(scale))),
                HealingOperatorWithCondition::new(HealingOperator::PropagateTolerances),
                HealingOperatorWithCondition::new(HealingOperator::Repair),
            ],
            stop_on_clean: true,
            max_iterations: 1,
            base_tolerance: TOLERANCE_ABS * scale,
            tolerance_growth: 1.0,
            tolerance_cap: TOLERANCE_ADAPTIVE_MAX * scale,
            healing_options: HealingOptions::default(),
            operator_params: OperatorParams::default(),
            collect_timing: true,
        }
    }
}

/// Extended report from running an advanced operator chain.
#[derive(Debug, Clone)]
pub struct OperatorChainReport {
    /// Results from each operator execution.
    pub operator_results: Vec<OperatorResult>,
    /// Initial check result.
    pub initial: CheckResult,
    /// Final check result.
    pub final_result: CheckResult,
    /// Initial issue stats.
    pub initial_stats: HealingIssueStats,
    /// Final issue stats.
    pub final_stats: HealingIssueStats,
    /// Total execution time in seconds.
    pub total_elapsed_seconds: f64,
    /// Number of operators executed (not skipped).
    pub operators_executed: usize,
    /// Number of operators skipped.
    pub operators_skipped: usize,
    /// Whether the shape is now clean.
    pub is_clean: bool,
    /// Summary description.
    pub summary: String,
}

/// Report for a single stage in the ShapeProcess pipeline.
#[derive(Debug, Clone)]
pub struct StageReport {
    /// The stage type.
    pub stage: HealingStage,
    /// Zero-based pass index (for multi-pass stages).
    pub pass_index: Option<usize>,
    /// Issue count before this stage.
    pub issue_count_before: usize,
    /// Issue count after this stage.
    pub issue_count_after: usize,
    /// Reports from individual operators executed in this stage.
    pub operator_reports: Vec<OperatorReport>,
    /// Wall-clock time for this stage (seconds).
    pub elapsed_seconds: f64,
}

impl StageReport {
    pub fn issues_fixed(&self) -> usize {
        self.issue_count_before.saturating_sub(self.issue_count_after)
    }

    pub fn is_improved(&self) -> bool {
        self.issue_count_after < self.issue_count_before
    }
}

/// Overall statistics for a ShapeProcess run.
#[derive(Debug, Clone, Default)]
pub struct ShapeProcessStats {
    /// Total number of operators executed.
    pub operators_executed: usize,
    /// Total number of modifications made.
    pub total_modifications: usize,
    /// Total number of issues fixed.
    pub total_issues_fixed: usize,
    /// Number of stages executed.
    pub stages_executed: usize,
    /// Total wall-clock time (seconds).
    pub total_elapsed_seconds: f64,
    /// Number of iterations (when max_iterations > 1).
    pub iterations: usize,
    /// Whether the process converged early (shape became clean).
    pub converged_early: bool,
    /// Final shape is checker-clean.
    pub is_clean: bool,
}

/// Complete report from a ShapeProcess pipeline run.
#[derive(Debug, Clone)]
pub struct ShapeProcessReport {
    /// Initial check result.
    pub initial: CheckResult,
    /// Final check result.
    pub final_result: CheckResult,
    /// Structured issue counts before processing.
    pub initial_stats: HealingIssueStats,
    /// Structured issue counts after processing.
    pub final_stats: HealingIssueStats,
    /// Per-stage reports.
    pub stages: Vec<StageReport>,
    /// Overall statistics.
    pub stats: ShapeProcessStats,
    /// Configuration used for this run.
    pub config_summary: String,
}

impl ShapeProcessReport {
    pub fn initial_issue_count(&self) -> usize {
        self.initial.issues.len()
    }

    pub fn final_issue_count(&self) -> usize {
        self.final_result.issues.len()
    }

    pub fn issues_fixed(&self) -> usize {
        self.initial_issue_count().saturating_sub(self.final_issue_count())
    }

    pub fn is_improved(&self) -> bool {
        self.final_issue_count() < self.initial_issue_count()
    }

    pub fn is_clean(&self) -> bool {
        self.final_result.is_valid()
    }

    pub fn summary(&self) -> String {
        if self.is_clean() {
            format!(
                "ShapeProcess: Clean result after {} operators, {} modifications in {:.3}s",
                self.stats.operators_executed,
                self.stats.total_modifications,
                self.stats.total_elapsed_seconds
            )
        } else {
            format!(
                "ShapeProcess: {} 鈫?{} issues ({} fixed) after {} operators in {:.3}s",
                self.initial_issue_count(),
                self.final_issue_count(),
                self.issues_fixed(),
                self.stats.operators_executed,
                self.stats.total_elapsed_seconds
            )
        }
    }
}

/// Configuration for the ShapeProcess pipeline.
///
/// This struct provides OCCT ShapeProcess-like configuration for running
/// a customizable sequence of healing operations on a BRep.
#[derive(Debug, Clone)]
pub struct ShapeProcessConfig {
    /// Sequence of healing operators to execute.
    pub operators: Vec<HealingOperator>,
    /// Stop processing if the shape becomes checker-clean.
    pub stop_on_clean: bool,
    /// Maximum number of iterations (0 = run once).
    pub max_iterations: usize,
    /// Tolerance growth factor per iteration.
    pub tolerance_growth: f64,
    /// Maximum tolerance cap.
    pub tolerance_cap: f64,
    /// Base tolerance for operations.
    pub base_tolerance: f64,
    /// Parameters for individual operators.
    pub operator_params: OperatorParams,
    /// Healing options for internal passes.
    pub healing_options: HealingOptions,
}

impl Default for ShapeProcessConfig {
    fn default() -> Self {
        Self {
            operators: vec![
                HealingOperator::ParametricConsistency,
                HealingOperator::Repair,
                HealingOperator::StopIfClean,
            ],
            stop_on_clean: true,
            max_iterations: 1,
            tolerance_growth: 1.0,
            tolerance_cap: TOLERANCE_ADAPTIVE_MAX,
            base_tolerance: TOLERANCE_ABS,
            operator_params: OperatorParams::default(),
            healing_options: HealingOptions::default(),
        }
    }
}
include!("e1.rs");
include!("e2.rs");
include!("heal_lib.rs");

