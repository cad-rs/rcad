//! Adaptor3d_HVertex (TKG3d) — a domain vertex of a TopolTool.
//!
//! 1:1 translation of OCCT `Adaptor3d_HVertex.hxx` (L24-52) + `.cxx`
//! (L24-62).

use glam::DVec2;
use rcad_kernel::math::el::elclib_line_parameter_2d;
use rcad_kernel::topods::Orientation;

use crate::geomalgo::geom2d_int::Curve2dAdaptor;

/// OCCT Adaptor3d_HVertex — a point vertex of the surface domain with an
/// orientation and a 2D parametric resolution.
#[derive(Debug, Clone, Copy)]
pub struct HVertex {
    my_pnt: DVec2,
    my_tol: f64,
    my_ori: Orientation,
}

impl HVertex {
    /// OCCT Adaptor3d_HVertex() (cxx L24-27) — empty.
    pub fn new() -> Self {
        HVertex {
            my_pnt: DVec2::ZERO,
            my_tol: 0.0,
            my_ori: Orientation::Forward,
        }
    }

    /// OCCT Adaptor3d_HVertex(P, Ori, Resolution) (cxx L29-36).
    pub fn new_with(p: DVec2, ori: Orientation, resolution: f64) -> Self {
        HVertex {
            my_pnt: p,
            my_tol: resolution,
            my_ori: ori,
        }
    }

    /// OCCT Value() (cxx L38-41).
    pub fn value(&self) -> DVec2 {
        self.my_pnt
    }

    /// OCCT Parameter(C) (cxx L43-46) — ElCLib::Parameter(C->Line(),
    /// myPnt).
    pub fn parameter(&self, c: &dyn Curve2dAdaptor) -> f64 {
        let l = c.line();
        elclib_line_parameter_2d(self.my_pnt, l.origin, l.direction)
    }

    /// OCCT Resolution(C) (cxx L48-51) — the stored parametric resolution.
    pub fn resolution(&self, _c: &dyn Curve2dAdaptor) -> f64 {
        self.my_tol
    }

    /// OCCT Orientation() (cxx L53-56).
    pub fn orientation(&self) -> Orientation {
        self.my_ori
    }

    /// OCCT IsSame(Other) (cxx L58-62) — distance <= Precision::Confusion().
    pub fn is_same(&self, other: &HVertex) -> bool {
        self.my_pnt.distance(other.value()) <= rcad_kernel::precision::CONFUSION
    }
}

impl Default for HVertex {
    fn default() -> Self {
        HVertex::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topalgo::adaptor2d::line2d::Line2dAdaptor;

    /// OCCT anchor: the parameter of a vertex on a restriction line is the
    /// projected arc-length parameter (ElCLib::Parameter).
    #[test]
    fn hvertex_parameter_anchor() {
        let l = Line2dAdaptor::new_pnt_dir(DVec2::new(0.0, 2.0), DVec2::X, -1.0, 5.0);
        let v = HVertex::new_with(DVec2::new(3.0, 2.0), Orientation::Forward, 1.0e-8);
        assert!((v.parameter(&l) - 3.0).abs() < 1e-12);
        assert!((v.resolution(&l) - 1.0e-8).abs() < 1e-18);
        assert!(v.is_same(&HVertex::new_with(
            DVec2::new(3.0, 2.0 + 1e-9),
            Orientation::Reversed,
            1e-8
        )));
        assert!(!v.is_same(&HVertex::new_with(
            DVec2::new(3.5, 2.0),
            Orientation::Forward,
            1e-8
        )));
    }
}
