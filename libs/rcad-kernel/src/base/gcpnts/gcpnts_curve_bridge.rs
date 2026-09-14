//! Architecture glue (no direct OCCT counterpart): the adapter that rides an
//! `Extrema_CurveTool`-faced curve onto the GCPnts template interface.
//!
//! OCCT instantiates the GCPnts templates with `TheCurve = Adaptor3d_Curve`
//! — the same `Adaptor3d_Curve&` the `Extrema_CurveTool` statics forward to
//! (Extrema_CurveTool.hxx L38-142), so
//! `Extrema_CurveTool::DeflCurvIntervals` (cxx L83) constructs
//! `GCPnts_TangentialDeflection aPntGen(C, ...)` directly.  The rcad kernel
//! carries the two faces on separate traits: the Extrema statics on
//! `ExtremaCurveTool` (held as the `&dyn ExtPCurveTool` trait object the
//! Extrema algorithms store) and the GCPnts template surface on
//! [`GCPntsCurve`](super::gcpnts_curve::GCPntsCurve).  [`CurveToolAsGCPnts`]
//! restores the OCCT identity: it implements the GCPnts surface by routing
//! every query to its OCCT-named `Extrema_CurveTool` static.

use glam::DVec3;

use crate::base::extrema_ext_pc::ExtPCurveTool;
use crate::base::proj_lib::CurveType;
use crate::math::GeomAbsShape;

use super::gcpnts_curve::GCPntsCurve;

/// The `Adaptor3d_Curve` flavor of the GCPnts template interface, served over
/// an `Extrema_CurveTool` facade object (the `&dyn ExtPCurveTool` the Extrema
/// algorithms hold; consumed at `Extrema_CurveTool::DeflCurvIntervals` cxx
/// L83).
pub struct CurveToolAsGCPnts<'a> {
    /// The `Extrema_CurveTool` facade (Extrema_CurveTool.hxx L38-142).
    pub the_c: &'a dyn ExtPCurveTool,
}

impl GCPntsCurve for CurveToolAsGCPnts<'_> {
    /// OCCT Extrema_CurveTool::FirstParameter (hxx L43).
    fn first_parameter(&self) -> f64 {
        self.the_c.first_parameter()
    }

    /// OCCT Extrema_CurveTool::LastParameter (hxx L45).
    fn last_parameter(&self) -> f64 {
        self.the_c.last_parameter()
    }

    /// OCCT Extrema_CurveTool::D0 (hxx L84-87).
    fn d0(&self, u: f64) -> DVec3 {
        self.the_c.d0(u)
    }

    /// OCCT Extrema_CurveTool::D1 (hxx L89-92).
    fn d1(&self, u: f64) -> (DVec3, DVec3) {
        self.the_c.d1(u)
    }

    /// OCCT Extrema_CurveTool::D2 (hxx L94-101).
    fn d2(&self, u: f64) -> (DVec3, DVec3, DVec3) {
        self.the_c.d2(u)
    }

    /// OCCT Extrema_CurveTool::GetType (hxx L80).
    fn get_type(&self) -> CurveType {
        self.the_c.get_type()
    }

    /// OCCT Extrema_CurveTool::NbIntervals (hxx L51-54) at GeomAbs_CN.
    fn nb_intervals_cn(&self) -> usize {
        self.the_c.nb_intervals(GeomAbsShape::CN) as usize
    }

    /// OCCT Extrema_CurveTool::Intervals (hxx L59-64) at GeomAbs_CN.
    fn intervals_cn(&self) -> Vec<f64> {
        self.the_c.intervals(GeomAbsShape::CN)
    }

    /// OCCT Extrema_CurveTool::Circle (hxx L120) — `theC.Circle().Radius()`;
    /// raises Standard_NoSuchObject for a non-circle adaptor.
    fn circle_radius(&self) -> f64 {
        self.the_c.circle().radius
    }

    /// OCCT Extrema_CurveTool::Degree (hxx L128).
    fn curve_degree(&self) -> i32 {
        self.the_c.curve_degree()
    }

    /// OCCT Extrema_CurveTool::NbPoles (hxx L132).
    fn nb_poles(&self) -> i32 {
        self.the_c.nb_poles()
    }

    /// OCCT Extrema_CurveTool::IsRational (hxx L130).
    fn is_rational(&self) -> bool {
        self.the_c.is_rational()
    }

    /// OCCT Extrema_CurveTool::DN (hxx L113-116) at order 1.
    fn dn1(&self, u: f64) -> DVec3 {
        self.the_c.dn(u, 1)
    }
}
