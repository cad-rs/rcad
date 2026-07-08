// OCCT-style port: many helpers and data paths are staged for parity; keep CI output clean.
#![allow(dead_code, reason = "staged OCCT-parity helpers and tables")]
#![allow(
    private_interfaces,
    reason = "internal fillet/tcol types re-exported at crate root"
)]
#![allow(
    unreachable_patterns,
    reason = "defensive matches over evolving geometry enums"
)]
#![allow(clippy::all)]

pub use brep_graph::{
    BRepGraphHistory, NamedGraph, NodeKind, TopoGraph, TopoGraphHistory, TopoGraphHistoryEvent,
    TopoGraphValidationIssue, TopoNode,
};
pub mod bnd_lib;
pub mod bnd_lib_2d;
pub mod boolean;
mod boolean_unit_octant;
pub mod bopds;
pub mod brep_algo;
pub mod brep_algo_api;
pub mod brep_bnd;
pub mod brep_check;
pub mod brep_check_parallel;
pub mod brep_graph;
pub mod brep_lib;
pub mod brep_repair;
pub mod brep_tools;
pub mod brep_top_adaptor;
pub mod bspline_approx_interp;
pub mod bspline_edit;
pub mod boptools;
pub mod bopalgo;
pub mod builder;
pub mod debug_trace;
pub mod pipeline_dump;
pub mod bvh;
pub mod classify;
pub mod defeature;
pub mod draft;
// pub mod ds_to_brep; // disabled during OCCT-alignment migration
pub mod features;
pub mod geom_convert;
pub mod geom_lib;
pub mod geom_populate;
pub mod gluer;
pub mod healing;
pub mod history;
pub mod hlr;
pub mod imprint;
pub mod non_manifold;
pub mod shape_algo;
pub mod shape_analysis;
pub mod shape_build;
pub mod shape_construct;
pub mod shape_custom;
pub mod shape_extend;
pub use brep_feat::{
    BRepFeatError, DraftFeatureParams, FeatureParams, FuseMode, GrooveParams, RibParams,
    apply_depouille, apply_draft_feature, make_drafted_prism, make_groove, make_linear_rib as make_linear_rib_feat,
    make_loft_feature, make_pipe_feature, make_prism_feature, make_revol_feature, make_rib,
    make_through_groove,
};
pub use defeature::{
    BlendFeature, ConicalFeature, CylindricalFeature, DefeaturingError, DefeaturingOptions,
    DefeaturingOptionsEnhanced, DefeaturingReport, DefeaturingReportEnhanced, FeatureGroup,
    HolePattern, HolePatternType, PocketFeature, SlotFeature, defeature_brep,
    defeature_brep_enhanced, detect_blend_features, detect_conical_features,
    detect_connected_feature_groups, detect_cylindrical_features, detect_hole_patterns,
    detect_pocket_features, detect_slot_features, identify_small_faces,
};
pub use features::{
    FeatureError, SplitShapeError, extrude_polygon_solid, make_cylindrical_hole, make_draft_prism,
    make_linear_rib, make_prism, make_revolution, make_revolution_rib, revolve_polygon_solid,
    split_face_by_wire,
};
pub use geom2d_api::{
    circle_through_three_points, circles_tangent_to_circle_and_line_through_point,
    circles_tangent_to_circle_and_two_lines, circles_tangent_to_circle_through_points,
    circles_tangent_to_line_through_points, circles_tangent_to_three_circles,
    circles_tangent_to_three_lines, circles_tangent_to_two_circles_and_line,
    circles_tangent_to_two_circles_through_point, circles_tangent_to_two_lines_through_point,
};
// pub mod adaptor3d; -- removed (triggers rustc ICE on 1.94.1)
pub mod approx_int;
pub mod adv_app2_var;
pub mod app_cont;
pub mod array;
pub mod blend;
mod bop_occt_union;
pub mod brep_adaptor;
pub mod brep_feat;
pub mod brep_int_curve_surface;
pub mod brep_mesh;
pub mod brep_offset;
pub mod brep_proj;
pub mod cells_builder;
pub mod chamfer;
pub mod elc_lib;
pub mod els_lib;
pub mod extrema;
pub mod fillet;
pub mod gcpnts;
pub mod gc_make;
pub mod geom2d_api;
pub mod int_ana;
pub mod inttools;
pub mod law;
pub mod lprop_cur_and_inf;
pub mod maker_volume;
pub mod math_utils;
pub mod medial_axis;
pub mod offset;
pub mod offset_prism;
pub mod pave_filler;
pub mod point_cloud;
pub mod projection;
pub mod section;
pub mod splitter;
pub mod sweep;
pub mod tcol_std;
pub mod thicken;
pub mod tkgeombase_algo;
pub mod tkgeombase_gtests;
pub mod tolerance;
pub use tolerance::TOLERANCE_MESH_LEGACY;
use crate::tolerance::*;

pub mod top_loc;
pub mod triangulate;

use serde::Serialize;

pub use approx_int::{
    ApproxOptions, ApproxResult, IntersectionApproximator, IntersectionSample,
    adjust_same_parameter, approximate_2d_curve, approximate_2d_curve_with_ctrl,
    approximate_intersection, approximate_polyline, compute_same_parameter,
    compute_same_parameter_bspline, sample_curve_segment, sample_intersection_points,
    sample_with_adaptive_density,
};
pub use brep_algo::{
    BRepAlgoError, OrientationIssue as BRepAlgoOrientationIssue, check_orientation,
    evaluate_edge_tangent, evaluate_face_normal, evaluate_vertex_normal, find_connected_components,
    fix_orientation, is_valid_brep, max_edge_length, max_face_area, min_face_area, total_edge_length,
    propagate_edge_tolerances, propagate_face_tolerances, reverse_face, total_surface_area,
    total_volume,
};
pub use brep_bnd::{
    BoundingBox, add_brep_to_bbox, add_edge_to_bbox, add_face_to_bbox, add_vertex_to_bbox,
    curve_bounds, curve_bounds_default, curve_bounds_with_range, surface_bounds,
    surface_bounds_with_domain,
};
pub use brep_lib::{
    BRepLibError, EdgeData, FaceData, FittedSurfaceType, FoundSurface, add_edge_with_curve,
    add_face_with_surface, compute_edge_bounds, compute_face_bounds, faces_share_surface,
    find_surface_through_edges, find_surface_through_points, make_edge_from_curve,
    make_face_from_surface, make_wire_from_edges, sort_faces_by_area, sort_faces_by_bounding_box,
    sort_faces_by_distance,
};
pub use brep_tools::{
    BRepToolsError, ShapeType, bounding_box, count_edges, count_faces, count_shells,
    count_vertices, count_wires, extract_shells, extract_solids, get_curve, get_edge_range,
    get_edge_tolerance, get_face_tolerance, get_inner_wires, get_outer_wire, get_pcurve,
    get_shape_type, get_surface, get_vertex_tolerance, is_closed, is_edge_degenerate,
    mirror_shape, n_ary_partition, read_brep_from_file, read_brep_from_string, rotate_shape,
    scale_shape, transform_shape, write_brep_to_file, write_brep_to_string,
};
pub use brep_top_adaptor::{
    EdgeAdaptor, EdgeExplorer, FaceAdaptor, FaceExplorer, OrientedEdge, ShapeIterator,
    VertexAdaptor, VertexExplorer, WireExplorer, edges_of_face, edges_of_vertex, face_count,
    faces_of_edge, faces_of_vertex, shell_count, vertices_of_edge, wire_count,
};
pub use bspline_edit::{
    move_bspline2_point, move_bspline2_tangent, move_bspline3_point, move_bspline3_tangent,
};
pub use bvh::{Aabb, Bvh, BvhStats};
pub use elc_lib::{
    bspline_derivative,
    // BSpline utilities
    bspline_point_at,
    circle_binormal_at,
    circle_derivative,
    circle_normal_at,
    circle_parameter,
    // Circle utilities
    circle_point_at,
    circle_tangent_at,
    ellipse_derivative,
    ellipse_parameter,
    // Ellipse utilities
    ellipse_point_at,
    hyperbola_derivative,
    // Hyperbola utilities
    hyperbola_point_at,
    line_closest_point,
    line_distance_to_point,
    line_parameter,
    // Line utilities
    line_point_at,
    parabola_derivative,
    // Parabola utilities
    parabola_point_at,
};
pub use els_lib::{
    bspline_surface_derivatives,
    bspline_surface_normal,
    // BSplineSurface utilities
    bspline_surface_point_at,
    cone_normal,
    cone_parameters,
    // Cone utilities
    cone_point_at,
    cylinder_normal,
    cylinder_parameters,
    // Cylinder utilities
    cylinder_point_at,
    cylinder_tangent_u,
    cylinder_tangent_v,
    plane_normal,
    plane_parameters,
    // Plane utilities
    plane_point_at,
    plane_tangent_u,
    plane_tangent_v,
    sphere_normal,
    sphere_parameters,
    // Sphere utilities
    sphere_point_at,
    sphere_tangent_u,
    sphere_tangent_v,
    torus_normal,
    torus_parameters,
    // Torus utilities
    torus_point_at,
    torus_tangent_u,
    torus_tangent_v,
};
pub use extrema::{
    closest_point_on_curve, closest_point_on_surface, distance_brep_brep, distance_curve_curve,
    distance_curve_surface, distance_point_curve, distance_point_point, distance_point_surface,
    distance_surface_surface, find_closest_points, find_furthest_points, find_supporting_edge,
    find_supporting_face,
};
pub use gcpnts::{
    adaptive_sample_curve, arc_length, point_at_arc_length, points_at_equal_arc_length,
    quasi_uniform, sample_surface_adaptive, sample_surface_grid, sample_surface_uniform,
    sample_u_isolines, sample_v_isolines, sampled_points_bounds, tangential_deflection,
    total_arc_length, uniform_abscissa, uniform_abscissa_points, uniform_deflection,
};
pub use geom_convert::{
    ConvertParams,
    approx_curve_to_bspline,
    approx_surface_to_bspline,
    bspline_surface_to_bezier,
    // BSpline operations
    bspline_to_bezier,
    circle_to_bspline,
    cone_to_bspline,
    curve_to_bspline,
    cylinder_to_bspline,
    ellipse_to_bspline,
    // Curve conversions
    line_to_bspline,
    // Surface conversions
    plane_to_bspline,
    sphere_to_bspline,
    surface_to_bspline,
    torus_to_bspline,
};
pub use geom_lib::{
    // Continuity checking
    check_curve_continuity,
    check_surface_continuity,
    // Normal estimation
    estimate_normal,
    estimate_normal_by_neighbors,
    // Closure checking
    is_curve_closed,
    is_surface_u_closed,
    is_surface_v_closed,
    // Degeneracy removal
    remove_degenerate_curve_sections,
    // Curve tools
    reverse_curve,
    // Surface tools
    reverse_surface_u,
    reverse_surface_v,
    transform_curve,
    transform_surface,
    trim_curve,
    trim_surface,
};
pub use geom2d_api::{
    Curve2dIntersection, curve2d_angle_at, curve2d_curvature_at, distance_between_curves2d,
    distance_point_to_curve2d, intersect_curves2d, points_to_bspline2d,
    points_to_bspline2d_interpolate, project_point_on_curve2d,
};
pub use top_loc::{
    Datum, Location, LocationManager, apply_location_to_shape, apply_location_to_shape_owned,
};

use rcad_kernel::BRep;
use rcad_kernel::topods;

// pub use adaptor3d::{Curve3dAdaptor, CurveOnSurfaceAdaptor, HSurfaceAdaptor, SurfaceAdaptor}; // removed
pub use array::{
    CircularPatternParams, LinearPatternParams, PatternError, circular_pattern, linear_pattern,
};
pub use blend::{
    BlendBoundary, BlendContinuity, BlendError, BlendMode, BlendParams, BlendQuality, BlendResult,
    RadiusLaw, SurfaceCurvePair, apply_blend_to_edge, blend_edge_to_face, blend_two_surfaces,
    blend_vertex, compute_blend_boundary_curves, compute_guide_curves, compute_pipe_blend,
    compute_rolling_ball_blend, compute_ruled_blend, compute_spine_curve,
};
pub use boolean::{
    BooleanAttemptDiagnostic, BooleanDiagnosticReport, BooleanFailureClass, FailureAnalyzer,
    FinalSuccessfulConfig, RecoveryStrategy, RetryPolicy, RetryPolicyBuilder,
};
pub use brep_algo_api::{
    // BRepAlgoAPI-style high-level boolean API
    BRepAlgoAPI_Common,
    BRepAlgoAPI_Cut,
    BRepAlgoAPI_Fuse,
    BRepAlgoAPI_Section,
    BRepHistory,
    BooleanApiOptions,
};
pub use brep_check::{
    CheckIssue,
    CheckResult,
    // Comprehensive check
    ComprehensiveCheckResult,
    EulerAnalysis,
    FaceSurfaceConsistencyDiagnosis,
    // Geometry validation (OCCT BRepCheck_Analyzer equivalent)
    GeometryValidationReport,
    OrientationIssue,
    OrientationReport,
    QualityMetricsConfig,
    // Quality metrics
    QualityMetricsReport,
    RicherValidityReport,
    SameParameterDiagnosis,
    SameRangeDiagnosis,
    ShellTopologyReport,
    SmallFeatureType,
    // Surface UV analysis (ShapeAnalysis_Surface equivalent)
    SurfaceAnalysisReport as SurfaceUvAnalysisReport,
    SuspectEdge,
    SuspectFaceSurfaceEdge,
    SuspectSameRangeEdge,
    // Tolerance checking
    ToleranceValidationReport,
    // Topology validation
    TopologyValidationReport,
    UvBoundsViolation,
    WireAnalysisReport,
    WireIssueReport,
    // Wire quality metrics (ShapeAnalysis_Wire enhancement)
    WireQualityMetrics,
    WireQualityReport,
    analyze_quality_metrics,
    analyze_shell_topology,
    analyze_surface_uv_consistency,
    analyze_wire_issues,
    analyze_wire_quality,
    brep_check_analyze,
    check_brep,
    check_comprehensive,
    check_curve_surface_consistency,
    check_edge_tolerance,
    check_orientation_consistency,
    check_tolerance_consistency,
    check_vertex_tolerance,
    diagnose_face_surface_consistency,
    diagnose_same_parameter,
    diagnose_same_range,
    euler_analysis,
    richer_validity_analysis,
    validate_nested_wires,
    validate_shell_orientation,
    validate_solid_closure,
    validate_wire_orientation,
};
pub use brep_check_parallel::{
    ParallelCheckIssue, ParallelCheckOptions, ParallelCheckResult, ParallelCheckStats,
    check_many_parallel, check_many_parallel_with_options, check_parallel,
    check_parallel_with_batch_size, check_parallel_with_options, check_parallel_with_stats,
};
pub use brep_int_curve_surface::{
    CurveBRepIntersection, CurveFaceIntersection, RayHit, intersect_curve_with_brep,
    intersect_curve_with_face, intersect_line_with_brep, intersect_line_with_face,
    is_point_inside_by_ray, ray_cast, shoot_ray,
};
pub use brep_mesh::{
    BRepMesh, Mesh, MeshParams, discretize_edge, discretize_edge_on_surface, mesh_aspect_ratio,
    mesh_brep as brep_mesh_brep, mesh_face, mesh_max_edge_length, mesh_min_angle, refine_mesh,
};
pub use brep_offset::{
    BRepOffsetOptions, EvolvedResult, MakeEvolved, MakeOffset, MakeOffsetShape, MakePipeShell,
    MakeThickSolid, OffsetMode, PipeShellResult, ThickSolidResult, WireOffsetResult, make_evolved,
    make_hollow_solid, make_pipe_shell, make_thick_solid, offset_shape_with_join,
    offset_shape_with_options, offset_wire,
};
pub use brep_proj::{BrepProjOptions, brep_proj_cylindrical};
pub use brep_repair::{
    AdaptiveToleranceConfig,
    AdaptiveToleranceMergeReport,
    // Comprehensive Tolerance Propagation (new)
    BooleanOpTypeForTolerance,
    ConflictResolutionPolicy,
    // Enhanced edge sewing and adaptive tolerance
    EdgeSewConfig,
    EdgeSewReport,
    EdgeValenceInfo,
    EnhancedEdgeSewReport,
    GapInfo,
    GapRepairFailureReason,
    // Enhanced Internal Face Detection and Removal
    InternalFaceDetectionConfig,
    InternalFaceDetectionReport,
    InternalFaceRemovalValidation,
    MakeConnectedReport,
    // MakeConnectedStrategy for configurable connectivity repair
    MakeConnectedStrategy,
    ManifoldRepairResult,
    NonManifoldEdgeInfo,
    PostBooleanRemovalConfig,
    PostBooleanRemovalReport,
    PostBooleanToleranceConfig,
    PostBooleanToleranceReport,
    PostSewToleranceConfig,
    PostSewToleranceReport,
    RepairReport,
    SeedDetectionConfig,
    SeedDetectionResult,
    // Seed detection for scoped make-connected
    SeedDetectionStrategy,
    ShellClosureResult,
    // Enhanced Shell Repair (ShapeFix_Shell extensions)
    ShellOrientationReport,
    ShellValidationReport,
    ToleranceAnalysisReport,
    ToleranceConsistencyReport,
    ToleranceFix,
    ToleranceFlowDirection,
    TolerancePropagationConfig,
    TolerancePropagationEngine,
    TolerancePropagationReport,
    ToleranceRule,
    ToleranceStats,
    ToleranceViolation,
    ToleranceViolationType,
    UnrepairedGap,
    UvBoundsRepairReport,
    // UV Gap Repair
    UvGapRepairConfig,
    UvGapRepairReport,
    VertexValenceInfo,
    WireGapRepairReport,
    analyze_tolerance_consistency,
    analyze_tolerances,
    apply_tolerance_fixes,
    detect_internal_faces,
    detect_internal_faces_with_config,
    detect_seeds_for_scoped_cleanup,
    fix_all_uv_gaps,
    fix_edge_pcurve_uv_bounds,
    fix_face_orientation,
    fix_same_parameter,
    fix_same_parameter_with_scan,
    fix_same_range_with_scan,
    fix_shell_orientation_advanced,
    fix_uv_bounds_violations,
    fix_uv_gaps,
    fix_wire_gaps,
    fix_wire_orientation,
    limit_tolerances,
    make_connected_baseline,
    make_connected_enhanced,
    make_connected_iterative,
    make_connected_iterative_scoped_with_growth_cap,
    make_connected_iterative_with_growth,
    make_connected_iterative_with_growth_cap,
    make_connected_scoped_auto,
    make_connected_with_strategy,
    merge_adjacent_faces_after_removal,
    merge_close_vertices,
    merge_close_vertices_topods,
    merge_vertices_adaptive,
    propagate_tolerances,
    propagate_tolerances_post_boolean,
    propagate_tolerances_post_boolean_op,
    propagate_tolerances_post_boolean_op_with_config,
    propagate_tolerances_post_sew,
    propagate_tolerances_post_sew_with_config,
    recompute_face_normals,
    remove_degenerate_faces,
    remove_internal_faces_post_boolean,
    remove_internal_faces_post_boolean_with_config,
    remove_small_edges,
    repair,
    repair_non_manifold_edges,
    repair_shell_closure,
    sew_close_edges,
    sew_edges_enhanced,
    validate_internal_face_removal,
    validate_shell_topology,
};
pub use builder::{
    BooleanError,
    BooleanOpType,
    // Glue path enhancement types
    GlueConfig,
    GlueFaceCache,
    GlueFacePair,
    apply_glue_optimization,
    compute_adaptive_glue_tolerance,
    detect_glue_faces,
};
pub use cells_builder::{CellExpr, CellsBuilder, CellsBuilderError};
pub use chamfer::{
    ChamferError, ChamferMode, ChamferParams, ChamferResult, ChamferWarning,
    compute_chamfer_curves, compute_chamfer_surface, make_chamfer_all_edges, make_chamfer_angle,
    make_chamfer_asymmetric, make_chamfer_edge, trim_adjacent_faces,
};
pub use fillet::{
    FilletContinuity, FilletError, FilletMode, FilletParams, FilletResult, VariableRadiusPoint,
    blend_adjacent_faces, compute_fillet_curves, compute_rollball_surface, make_fillet_all_edges,
    make_fillet_edge, make_fillet_edge_with_params, make_variable_fillet,
};
pub use gluer::{
    EdgeOrigin as GluerEdgeOrigin, FaceOrigin as GluerFaceOrigin, Gluer, GluerError, GluerHistory,
    GluerMode, GluerOptions, GluerResult, InterfaceInfo, VertexOrigin as GluerVertexOrigin,
    detect_interface, detect_interface_bvh, glue_at_interface, glue_shapes,
};
pub use healing::{
    BRepSnapshot,
    ComprehensiveDiagnosis,
    ComprehensiveHealingReport,
    // New ShapeProcess operators
    DirectFacesOperator,
    HealGeometryOperator,
    HealGeometryStep,
    HealingIssueStats,
    HealingMode,
    HealingOperator,
    HealingOptions,
    HealingReport,
    HealingStage,
    HealingStageReport,
    MakeConnectedPrepassMode,
    OperatorParams,
    OperatorReport,
    // Operator result aggregation and rollback
    OperatorResultAggregation,
    ParametricConsistencyReport,
    PipelineExecutionReport,
    // Progress callbacks
    ProgressCallback,
    RemoveInternalFacesOperator,
    RollbackConfig,
    SameParameterOperator,
    ShapeProcessConfig,
    ShapeProcessReport,
    ShapeProcessStats,
    SimpleProgressCallback,
    SolidFixReport,
    StageReport,
    WireFixReport,
    WireIssueLocation,
    analyze_and_heal,
    diagnose_all,
    // ShapeFix_Solid and ShapeFix_Wire equivalents
    fix_solid,
    fix_wire,
    heal,
    heal_comprehensive,
    run_healing_operator_chain,
    run_healing_pipeline_with_rollback,
    run_shape_process,
};
pub use history::{
    BooleanHistory, BooleanNamingPropagationReport, BooleanOperationType, ChainStatistics,
    DeletionReason, DeletionRecord, EdgeOrigin, EntityType, FaceOrigin, GenerationCause,
    GenerationRecord, HistoryChain, HistoryStatistics, HistoryTracker, InputSource,
    ModificationRecord, ModificationType, ShellOrigin, SolidOrigin, VertexOrigin,
};
pub use hlr::{
    AssemblyHlrResult, ComponentHlr, CurveHint, HlrCamera, HlrOptions, HlrResult, HlrSegment,
    SegmentType, SilhouetteCurve3, compute_hlr, compute_hlr_with_options,
    extract_silhouette_curves, hlr_assembly, hlr_to_svg,
};
pub use imprint::{
    Gap, GapOverlapReport, ImprintResult, Overlap, detect_gaps_overlaps, imprint_shape,
    min_distance,
};
pub use int_ana::{
    // Cylinder-Cylinder intersection (IntAna_IntCylCyl)
    CylCylResult,
    // Line-Surface intersections (IntAna_IntLinPln, IntAna_IntLinCyl, etc.)
    LinPlnIntersection,
    PlnConResult,
    PlnCylResult,
    // Plane-Surface intersections (IntAna_IntPlnPln, IntAna_IntPlnCyl, etc.)
    PlnPlnResult,
    PlnSphResult,
    intersect_cylinder_cylinder,
    intersect_line_cone,
    intersect_line_cylinder,
    intersect_line_plane,
    intersect_line_sphere,
    intersect_line_torus,
    intersect_plane_cone_intana,
    intersect_plane_cylinder_intana,
    intersect_plane_plane_intana,
    intersect_plane_sphere_intana,
};
pub use inttools::{
    ASPECT_RATIO_THRESHOLD,
    ASPECT_RATIO_VERY_HIGH,
    // Extreme geometry handling
    AspectRatioAdaptiveTolerance,
    DegenerateGeometryHandler,
    DegenerateType,
    ExtremeGeometryAnalysis,
    ExtremeGeometryAnalysisOptions,
    HighAspectRatioEdge,
    HighAspectRatioFace,
    NearDegenerateGeometry,
    NearTangentConfig,
    NearTangentHandler,
    NearTangentSeverity,
    SIZE_RATIO_THRESHOLD,
    SizeDifferenceAnalysis,
    SizeDifferenceHandler,
    SurfaceCurve,
    SurfaceIntersectionResult,
    SurfaceSurfaceIntersection,
    analyze_extreme_geometry,
    analyze_size_difference,
    detect_high_aspect_ratio_edges,
    detect_near_degenerate_geometry,
    detect_near_tangent_configurations,
    intersect_surfaces,
    intersect_surfaces_with_density,
    intersect_surfaces_with_density_tol,
    intersect_surfaces_with_tolerance,
};
pub use law::{
    BSplineLaw, CompositeLaw, ConstantLaw, InterpolateLaw, LawFunction, LinearLaw, SineLaw,
    SmoothStepLaw, sine_law, smooth_step_law,
};
pub use maker_volume::{
    MakerVolume, MakerVolumeError, MakerVolumeSelection, make_solid_from_cell_indices,
    make_solid_from_region, make_solid_from_region_with_history,
};
pub use math_utils::{
    bisection,
    determinant_3x3,
    // Eigenvalue/Matrix
    eigenvalues_2x2,
    eigenvalues_3x3,
    gaussian_quadrature,
    golden_section_max,
    // Optimization
    golden_section_min,
    inverse_3x3,
    // Multi-dimensional Newton
    newton_2d,
    newton_3d,
    // Root finding
    newton_raphson,
    secant,
    // Integration
    simpson_integrate,
    solve_cubic,
    // Polynomial solvers
    solve_linear,
    solve_quadratic,
    solve_quartic,
};
pub use medial_axis::{
    MedialAxis2d, MedialAxisOptions, MedialBranch2d, MedialEdge, MedialFace, MedialPoint2d,
    MedialSurface, MedialVertex, MidSurfaceResult, ThicknessMap, ThicknessSample, ThicknessStats,
    ThinRegion, VoronoiDiagram2d, VoronoiEdge2d, VoronoiVertex2d, WallThicknessResult,
    cluster_medial_vertices, compute_mat_2d, compute_medial_axis_2d, compute_medial_surface,
    compute_mid_surface, compute_thickness_map, compute_voronoi_2d, compute_wall_thickness,
    detect_thin_regions, find_max_inscribed_circle, generate_rib_paths, point_in_polygon_2d,
};
pub use non_manifold::{
    EdgeSplitReport, MakeManifoldOptions, MakeManifoldReport, MergeShellsOptions,
    MergeShellsResult, NonManifoldReport, NonManifoldTraversal, analyze_non_manifold,
    boundary_edges, is_manifold, make_manifold, make_manifold_with_options,
    merge_shells_at_interface, multi_face_edges, non_manifold_edges, non_manifold_vertices,
    orphan_edges, split_non_manifold_edges,
};
pub use offset::{
    JoinType, OffsetError, OffsetOptions, OffsetQuality, OffsetResult, VariableThickness,
    detect_self_intersection, hollow_solid, hollow_solid_with_options, offset_shape, offset_shell,
    offset_shell_with_options, offset_solid, offset_solid_with_options, offset_surface,
};
pub use point_cloud::{
    Dimensionality, FittedCylinder, FittedPlane, FittedPolygon, FittedSphere, OutlierPoint,
    PointCloud, PointCloudAnalysis, SamplingStrategy, analyze_point_cloud, compute_inertia,
    compute_pca, detect_outliers, estimate_dimensionality, estimate_normals,
    extract_points_from_brep_mesh, extract_points_from_brep_vertices, extract_points_from_mesh,
    fit_cylinder, fit_plane, fit_polygon, fit_sphere, remove_outliers,
    sample_points_from_brep_surfaces, simplify_point_cloud,
};
pub use projection::{
    PointBRepProjection, PointCurveProjection, PointSurfaceProjection, ProjectionDirection,
    ProjectionOptions, SilhouetteResult, compute_all_curve_surface_projections,
    compute_contour_edges, compute_silhouette_curves, directional_project_curve_on_surface,
    normal_project_curve_on_surface, project_curve_on_surface, project_point_on_brep,
    project_point_on_curve, project_point_on_curve_with_options, project_point_on_surface,
    project_point_on_surface_with_options, project_surface_on_surface, project_wire_on_face,
    project_wire_on_surface,
};
pub use section::{
    SectionCurve, brep_section, brep_triangle_soup, intersect_triangle_soups,
    intersect_triangle_soups_adaptive, intersect_triangle_soups_eps,
    intersect_triangle_soups_for_brep_tolerance, section, section_curves, section_polylines,
};
pub use shape_algo::{
    // Algorithm container
    AlgoContainer,
    // Geometry extraction structures
    BoxGeometry,
    ConeGeometry,
    CylinderGeometry,
    ShapeAlgorithm,
    SphereGeometry,
    TorusGeometry,
    // Geometry extraction functions
    get_box_geometry,
    get_cone_geometry,
    get_cylinder_geometry,
    get_sphere_geometry,
    get_torus_geometry,
    // Primitive detection
    is_box,
    is_cone,
    is_cylinder,
    is_sphere,
    is_torus,
};
pub use shape_analysis::{
    // Full BRep analysis
    BRepAnalysisReport,
    ContinuityLevel,
    // Curve analysis (ShapeAnalysis_Curve)
    CurveAnalysisReport,
    CurveSelfIntersection,
    // Face analysis (ShapeAnalysis_Face)
    FaceAnalysisReport,
    OverTrimmedRegion,
    ParamRangeIssue,
    SeamEdgeIssue,
    SingularPoint,
    SingularPointKind,
    // Surface analysis (ShapeAnalysis_Surface)
    SurfaceAnalysisReport as ShapeAnalysisSurfaceReport,
    // Enhanced ShapeAnalysis_Surface equivalent
    SurfaceBoundsAnalysis,
    SurfaceDeviation,
    SurfaceDeviationViolation,
    SurfaceWireIssue,
    SurfaceWireIssueKind,
    UnderTrimmedRegion,
    UvConsistencyReport as FaceUvConsistencyReport,
    UvFlipIssue,
    UvFlipType,
    UvInconsistency,
    UvInconsistencyKind,
    // Wire analysis (ShapeAnalysis_Wire)
    WireAnalysisReport as ShapeAnalysisWireReport,
    WireGap,
    WireSelfIntersection,
    analyze_brep,
    analyze_curve,
    analyze_face,
    analyze_surface,
    analyze_surface_bounds_for_face,
    analyze_wire,
    check_face_uv_consistency_by_idx,
    check_face_wires,
    check_uv_consistency,
    compute_surface_deviation,
    detect_surface_self_intersection,
};
pub use shape_build::{
    BRepBuilder, BuildError, BuildFace, BuildShell, BuildSolid, BuildVertex, BuildWire, Rebuild,
    validate_shell_closed, validate_solid_valid, validate_wire_closed,
};
pub use shape_construct::{
    // BSpline construction
    construct_bspline_curve,
    construct_bspline_surface,
    construct_circle_center_normal,
    construct_circle_from_3_points,
    construct_circle_wire,
    construct_cone_from_axis,
    construct_cylinder_from_axis,
    construct_ellipse_from_points,
    construct_face_from_boundary,
    // Curve construction
    construct_line,
    // Face construction
    construct_planar_face_from_wire,
    // Surface construction
    construct_plane_from_3_points,
    construct_plane_from_point_normal,
    // Wire construction
    construct_polygon_wire,
    construct_sphere_from_center_radius,
    construct_torus_from_center_radii,
};
pub use shape_custom::{
    BSplineSimplifyOptions, ConversionReport, GeometryRestrictions, SimplificationResult,
    convert_to_bspline, curve_degree, ensure_bspline_curve, ensure_bspline_surface,
    is_bspline_curve, is_bspline_surface, restrict_geometry, restrict_to_bspline,
    simplify_bspline_curve, simplify_bspline_surface, surface_degrees,
};
pub use shape_extend::{
    // ShapeExtend_CompositeSurface
    CompositeSurface,
    // ShapeExtend_BasicMsgRegistrator
    MessageRegistrator,
    MessageSeverity,
    ShapeContextMessage,
    // ShapeExtend_Explorer
    ShapeExplorer,
    ShapeMessage,
    // ShapeExtend_MsgRegistrator
    ShapeMessageRegistrator,
    // ShapeExtend_WireData
    WireData,
};
pub use sweep::{
    CornerMode, Law, PiecewiseLinearLaw, SweepError, SweepHistory, SweepMode, SweepOptions,
    handle_pipe_corners, linear_law_sweep, linear_sweep, linear_sweep_face, linear_sweep_wire,
    linear_sweep_with_history, linear_sweep_with_options, pipe_sweep, pipe_sweep_wire,
    pipe_sweep_with_history, pipe_sweep_with_options, pipe_with_rotation, rotational_sweep,
    rotational_sweep_face, rotational_sweep_wire, rotational_sweep_with_history,
    rotational_sweep_with_options, variable_section_sweep,
};
pub use thicken::{ThickeningResult, thicken_shell};
pub use triangulate::{
    AdaptiveSubdivider, BoundarySensitiveTessellator, FeatureEdge, IncrementalMesher, MeshDelta,
    MeshQualityMetrics, MeshSimplifier, SurfaceMesh, TessellationParams, compute_mesh_quality,
    mesh_brep, triangulate_surface,
};
include!("lib_inline/inc_1.rs");
include!("lib_inline/inc_2.rs");
include!("lib_inline/inc_3.rs");
include!("lib_inline/inc_4.rs");
include!("lib_inline/inc_5.rs");
include!("lib_inline/inc_6.rs");
include!("lib_inline/test_wrapper.rs");
