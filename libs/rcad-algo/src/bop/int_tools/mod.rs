pub mod classify_lin2d;
pub mod cone_cone;
pub mod context;
pub mod coplanar;
pub mod curve_range;
pub mod curve_surface;
pub mod cylinder_cone;
pub mod cylinder_cylinder;
pub mod cylinder_torus;
pub mod edge_edge;
pub mod edge_face;
pub mod ellipse_intersection;
pub mod extreme_geometry;
pub mod face_face;
pub mod fclass2d;
pub mod geom_abs_surface_type;
pub mod hyperbola_intersection;
pub mod imp_prm;
pub mod int_ana_quad_quad_geo;
pub mod int_patch_aline_to_wline;
pub mod int_patch_imp_imp_intersection;
pub mod int_patch_intersection;
pub mod int_patch_line;
pub mod int_patch_line_constructor;
pub mod int_patch_point;
pub mod int_patch_special_points;
pub mod int_patch_type;
pub mod int_patch_wline_tool;
pub mod int_surf_quadric;
pub mod intss;
pub mod marching;
pub mod parabola_intersection;
pub mod pcurve_derive;
pub mod plane_cone;
pub mod plane_cylinder;
pub mod plane_plane;
pub mod plane_sphere;
pub mod plane_torus;
pub mod sphere_cone;
pub mod sphere_cylinder;
pub mod sphere_torus;
pub mod torus_cone;
pub mod torus_torus;
pub mod vertex_ops;

pub use cylinder_torus::{
    CylinderTorusResult, intersect_cylinder_torus, intersect_cylinder_torus_with_tolerance,
};
pub use extreme_geometry::{
    ASPECT_RATIO_THRESHOLD, ASPECT_RATIO_VERY_HIGH, AspectRatioAdaptiveTolerance,
    DegenerateGeometryHandler, DegenerateType, ExtremeGeometryAnalysis,
    ExtremeGeometryAnalysisOptions, HighAspectRatioEdge, HighAspectRatioFace,
    NearDegenerateGeometry, NearTangentConfig, NearTangentHandler, NearTangentSeverity,
    SIZE_RATIO_THRESHOLD, SizeDifferenceAnalysis, SizeDifferenceHandler, analyze_extreme_geometry,
    analyze_size_difference, detect_high_aspect_ratio_edges, detect_near_degenerate_geometry,
    detect_near_tangent_configurations,
};
pub use intss::{
    SurfaceCurve, SurfaceIntersectionResult, SurfaceSurfaceIntersection,
    convert_polylines_to_bsplines, intersect_surfaces, intersect_surfaces_with_density,
    intersect_surfaces_with_density_tol, intersect_surfaces_with_tolerance, polyline_to_bspline,
};
pub use plane_torus::{
    PlaneTorusResult, intersect_plane_torus, intersect_plane_torus_with_tolerance,
};

/// Shared chord-error tolerance for adaptive refinement of analytic
/// intersection curves (`1e-6`). Equivalent to [`TOLERANCE_MESH_LEGACY`];
/// semantically distinct — chord tolerance in surface-surface intersection
/// curve refinement, not mesh merging.
pub const CHORD_TOLERANCE: f64 = crate::tolerance::TOLERANCE_MESH_LEGACY;

/// Maximum refinement depth for chord-error adaptive subdivision.
pub const CHORD_REFINE_DEPTH: usize = 2;

pub mod bean_face_intersector;
pub mod shrunk_range;
