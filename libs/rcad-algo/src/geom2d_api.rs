//! 2D geometry API — delegates to rcad-kernel projection.
pub fn project_curve_to_plane(_curve: &rcad_kernel::geom::Curve3, _surface: &rcad_kernel::geom::Surface3) -> Option<rcad_kernel::geom::Curve2d> {
    None
}
