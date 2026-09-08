//! GAP carrier: OCCT IntCurveSurface_HInter over a generic
//! (Adaptor3d_Curve, GeomAdaptor_Surface) pair as consumed by
//! GeomFill_LocationGuide / GeomFill_GuideTrihedronPlan.
//!
//! rcad hosts the IntCurveSurface engine itself
//! (`crate::geomalgo::int_curve_surface::HInter`), but driving it needs the
//! production `HCurveTool` / `HSurfaceTool` marker implementations for the
//! rcad `Curve3` / `Surface3` enums — infrastructure belonging to the
//! IntCurveSurface host-tool batch, not to the GeomFill sweep closure batch
//! (plan §9 E0 R4 staged).  Until those markers land, this carrier preserves
//! the OCCT failure path (panic at Perform) at the exact call sites.
//!
//! OCCT anchor: IntCurveSurface_HInter.hxx / IntCurveSurface_HInter.cxx
//! (Perform / NbPoints / Point) and
//! IntCurveSurface_IntersectionPoint.hxx (Pnt / W).

use glam::DVec3;

use rcad_kernel::geom::{Curve3, Surface3};

/// OCCT IntCurveSurface_IntersectionPoint — GAP carrier (only the Pnt()
/// and W() accessors consumed by the GeomFill sweep closure).
#[derive(Debug, Clone, Default)]
pub struct IntersectionPointGap {
    pnt: DVec3,
    w: f64,
}

impl IntersectionPointGap {
    /// OCCT IntCurveSurface_IntersectionPoint::Pnt().
    pub fn pnt(&self) -> DVec3 {
        self.pnt
    }

    /// OCCT IntCurveSurface_IntersectionPoint::W().
    pub fn w(&self) -> f64 {
        self.w
    }
}

/// OCCT IntCurveSurface_HInter — GAP carrier (see module doc).
#[derive(Debug, Clone, Default)]
pub struct IntCurveSurfaceHInter {
    _private: (),
}

impl IntCurveSurfaceHInter {
    /// OCCT IntCurveSurface_HInter::Perform(HCurve, Surface).
    pub fn perform(&mut self, _curve: &Curve3, _surface: &Surface3) {
        panic!(
            "GAP: IntCurveSurface_HInter host tools (HCurveTool/HSurfaceTool over Curve3/\
             Surface3) are not translated — see file header"
        )
    }

    /// OCCT IntCurveSurface_HInter::NbPoints().
    pub fn nb_points(&self) -> usize {
        panic!(
            "GAP: IntCurveSurface_HInter host tools (HCurveTool/HSurfaceTool over Curve3/\
             Surface3) are not translated — see file header"
        )
    }

    /// OCCT IntCurveSurface_HInter::Point(j).
    pub fn point(&self, _j: usize) -> IntersectionPointGap {
        panic!(
            "GAP: IntCurveSurface_HInter host tools (HCurveTool/HSurfaceTool over Curve3/\
             Surface3) are not translated — see file header"
        )
    }
}
