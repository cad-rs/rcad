//! Triangulation — delegates to rcad-kernel face triangulation.
pub fn triangulate_polygon(_verts: &[glam::DVec2]) -> Vec<[u32; 3]> {
    Vec::new()
}
pub fn triangulate_polygon_with_holes(
    _outer: &[glam::DVec2], _holes: &[&[glam::DVec2]],
) -> Vec<[u32; 3]> {
    Vec::new()
}
