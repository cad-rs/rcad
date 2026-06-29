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
pub mod extreme_geometry;
pub mod intss;
pub mod marching;
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
pub mod ellipse_intersection;
pub mod fclass2d;
pub mod hyperbola_intersection;
pub mod parabola_intersection;
pub mod vertex_ops;

pub use intss::{
    SurfaceCurve, SurfaceIntersectionResult, SurfaceSurfaceIntersection, convert_polylines_to_bsplines,
    intersect_surfaces, intersect_surfaces_with_density, intersect_surfaces_with_density_tol,
    intersect_surfaces_with_tolerance, polyline_to_bspline,
};
pub use extreme_geometry::{
    AspectRatioAdaptiveTolerance, DegenerateGeometryHandler, DegenerateType,
    HighAspectRatioEdge, HighAspectRatioFace, NearDegenerateGeometry,
    NearTangentConfig, NearTangentHandler, NearTangentSeverity,
    SizeDifferenceAnalysis, SizeDifferenceHandler,
    ExtremeGeometryAnalysis, ExtremeGeometryAnalysisOptions,
    analyze_extreme_geometry, analyze_size_difference,
    detect_high_aspect_ratio_edges, detect_near_degenerate_geometry,
    detect_near_tangent_configurations,
    ASPECT_RATIO_THRESHOLD, ASPECT_RATIO_VERY_HIGH, SIZE_RATIO_THRESHOLD,
};
pub use plane_torus::{PlaneTorusResult, intersect_plane_torus, intersect_plane_torus_with_tolerance};
pub use cylinder_torus::{CylinderTorusResult, intersect_cylinder_torus, intersect_cylinder_torus_with_tolerance};

pub mod shrunk_range;
