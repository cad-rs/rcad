// OCCT BRepExtrema — Minimum distance computation between shapes.
//
// OCCT ref: TKTopAlgo/BRepExtrema/BRepExtrema_DistShapeShape
//
// rcad: delegates to rcad-kernel::base::extrema.

pub mod dist_shape_shape;
pub use dist_shape_shape::{min_distance_edge_segments, min_distance_edge_vertex};
