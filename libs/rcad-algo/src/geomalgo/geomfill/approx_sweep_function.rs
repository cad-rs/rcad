//! OCCT Approx_SweepFunction (TKGeomBase/Approx) — 1:1 port of
//! Approx_SweepFunction.hxx + Approx_SweepFunction.cxx (whole file L24-84):
//! the abstract sweep-function interface consumed by AppBlend_AppSurf.
//!
//! Architecture mapping: the OCCT abstract class is a Rust trait; the
//! non-pure virtual bodies keep their OCCT defaults (Standard_NotImplemented
//! raises become panics carrying the OCCT message).

use glam::{DVec2, DVec3};

use rcad_kernel::math::GeomAbsShape;

/// OCCT Approx_SweepFunction.
pub trait ApproxSweepFunction {
    /// OCCT D0 — pure virtual: compute the section at Param into
    /// Poles/Poles2d/Weigths (First/Last are the current bounds).
    fn d0(
        &self,
        param: f64,
        first: f64,
        last: f64,
        poles: &mut [DVec3],
        poles2d: &mut [DVec2],
        weigths: &mut [f64],
    ) -> bool;

    /// OCCT D1 (Approx_SweepFunction.cxx L24-36) — default raises
    /// Standard_NotImplemented.
    #[allow(clippy::too_many_arguments)]
    fn d1(
        &self,
        _param: f64,
        _first: f64,
        _last: f64,
        _poles: &mut [DVec3],
        _dpoles: &mut [DVec3],
        _poles2d: &mut [DVec2],
        _dpoles2d: &mut [DVec2],
        _weigths: &mut [f64],
        _dweigths: &mut [f64],
    ) -> bool {
        panic!("Standard_NotImplemented: Approx_SweepFunction::D1");
    }

    /// OCCT D2 (L41-59) — default raises Standard_NotImplemented.
    #[allow(clippy::too_many_arguments)]
    fn d2(
        &self,
        _param: f64,
        _first: f64,
        _last: f64,
        _poles: &mut [DVec3],
        _dpoles: &mut [DVec3],
        _d2poles: &mut [DVec3],
        _poles2d: &mut [DVec2],
        _dpoles2d: &mut [DVec2],
        _d2poles2d: &mut [DVec2],
        _weigths: &mut [f64],
        _dweigths: &mut [f64],
        _d2weigths: &mut [f64],
    ) -> bool {
        panic!("Standard_NotImplemented: Approx_SweepFunction::D2");
    }

    /// OCCT Nb2dCurves — pure virtual.
    fn nb_2d_curves(&self) -> usize;

    /// OCCT SectionShape — pure virtual.
    fn section_shape(&self, nb_poles: &mut usize, nb_knots: &mut usize, degree: &mut usize);

    /// OCCT Knots — pure virtual.
    fn knots(&self, t_knots: &mut [f64]);

    /// OCCT Mults — pure virtual.
    fn mults(&self, t_mults: &mut [i32]);

    /// OCCT IsRational — pure virtual.
    fn is_rational(&self) -> bool;

    /// OCCT NbIntervals — pure virtual.
    fn nb_intervals(&self, s: GeomAbsShape) -> usize;

    /// OCCT Intervals — pure virtual.
    fn intervals(&self, t: &mut [f64], s: GeomAbsShape);

    /// OCCT SetInterval — pure virtual.
    fn set_interval(&mut self, first: f64, last: f64);

    /// OCCT Resolution (L63-68) — default raises Standard_NotImplemented.
    fn resolution(&self, _index: usize, _tol: f64, _tolu: &mut f64, _tolv: &mut f64) {
        panic!("Standard_NotImplemented: Approx_SweepFunction::Resolution");
    }

    /// OCCT GetTolerance — pure virtual.
    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, angle_tol: f64, tol3d: &mut [f64]);

    /// OCCT SetTolerance — pure virtual.
    fn set_tolerance(&mut self, tol3d: f64, tol2d: f64);

    /// OCCT BarycentreOfSurf (L70-73) — default raises
    /// Standard_NotImplemented.
    fn barycentre_of_surf(&self) -> DVec3 {
        panic!("Standard_NotImplemented: Approx_SweepFunction::BarycentreOfSurf");
    }

    /// OCCT MaximalSection (L75-78) — default raises
    /// Standard_NotImplemented.
    fn maximal_section(&self) -> f64 {
        panic!("Standard_NotImplemented: Approx_SweepFunction::MaximalSection()");
    }

    /// OCCT GetMinimalWeight (L80-84) — default raises
    /// Standard_NotImplemented.
    fn get_minimal_weight(&self, _weigths: &mut [f64]) {
        panic!("Standard_NotImplemented: Approx_SweepFunction::GetMinimalWeight");
    }
}
