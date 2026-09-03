//! OCCT GeomFill_SectionLaw (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_SectionLaw.hxx + GeomFill_SectionLaw.cxx (whole file): the base
//! trait for section laws along a swept path.

use glam::DVec3;

use rcad_kernel::geom::{BSplineSurface, Curve3};
use rcad_kernel::math::GeomAbsShape;

/// OCCT GeomFill_SectionLaw.
pub trait SectionLaw {
    /// OCCT D0 — pure virtual: the section at Param.
    fn d0(&self, param: f64, poles: &mut [DVec3], weights: &mut [f64]) -> bool;

    /// OCCT D1 — default throws Standard_NotImplemented.
    fn d1(
        &self,
        _param: f64,
        _poles: &mut [DVec3],
        _dpoles: &mut [DVec3],
        _weights: &mut [f64],
        _dweights: &mut [f64],
    ) -> bool {
        panic!("GeomFill_SectionLaw::D1");
    }

    /// OCCT D2 — default throws Standard_NotImplemented.
    #[allow(clippy::too_many_arguments)]
    fn d2(
        &self,
        _param: f64,
        _poles: &mut [DVec3],
        _dpoles: &mut [DVec3],
        _d2poles: &mut [DVec3],
        _weights: &mut [f64],
        _dweights: &mut [f64],
        _d2weights: &mut [f64],
    ) -> bool {
        panic!("GeomFill_SectionLaw::D2");
    }

    /// OCCT BSplineSurface() — default: null surface.
    fn bspline_surface(&self) -> Option<&BSplineSurface> {
        None
    }

    /// OCCT SectionShape — pure virtual.
    fn section_shape(&self, nb_poles: &mut usize, nb_knots: &mut usize, degree: &mut usize);

    /// OCCT Knots — pure virtual.
    fn knots(&self, t_knots: &mut [f64]);

    /// OCCT Mults — pure virtual.
    fn mults(&self, t_mults: &mut [i32]);

    /// OCCT IsRational — pure virtual.
    fn is_rational(&self) -> bool;

    /// OCCT IsUPeriodic — pure virtual.
    fn is_u_periodic(&self) -> bool;

    /// OCCT IsVPeriodic — pure virtual.
    fn is_v_periodic(&self) -> bool;

    /// OCCT NbIntervals — pure virtual.
    fn nb_intervals(&self, s: GeomAbsShape) -> usize;

    /// OCCT Intervals — pure virtual.
    fn intervals(&self, t: &mut [f64], s: GeomAbsShape);

    /// OCCT SetInterval — pure virtual.
    fn set_interval(&mut self, first: f64, last: f64);

    /// OCCT GetInterval — pure virtual.
    fn get_interval(&self, first: &mut f64, last: &mut f64);

    /// OCCT GetDomain — pure virtual.
    fn get_domain(&self, first: &mut f64, last: &mut f64);

    /// OCCT GetTolerance — pure virtual.
    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, angle_tol: f64, tol3d: &mut [f64]);

    /// OCCT SetTolerance — "Ne fait Rien".
    fn set_tolerance(&mut self, _tol3d: f64, _tol2d: f64) {}

    /// OCCT BarycentreOfSurf — default throws Standard_NotImplemented.
    fn barycentre_of_surf(&self) -> DVec3 {
        panic!("GeomFill_SectionLaw::BarycentreOfSurf");
    }

    /// OCCT MaximalSection — pure virtual.
    fn maximal_section(&self) -> f64;

    /// OCCT GetMinimalWeight — default throws Standard_NotImplemented.
    fn get_minimal_weight(&self, _weights: &mut [f64]) {
        panic!("GeomFill_SectionLaw::GetMinimalWeight");
    }

    /// OCCT IsConstant — default: Error = 0, false.
    fn is_constant(&self, error: &mut f64) -> bool {
        *error = 0.0;
        false
    }

    /// OCCT ConstantSection — default throws Standard_DomainError.
    fn constant_section(&self) -> Curve3 {
        panic!("Standard_DomainError: GeomFill_SectionLaw::ConstantSection");
    }

    /// OCCT IsConicalLaw — default: Error = 0, false.
    fn is_conical_law(&self, error: &mut f64) -> bool {
        *error = 0.0;
        false
    }

    /// OCCT CirclSection — default throws Standard_DomainError.
    fn circl_section(&self, _param: f64) -> Curve3 {
        panic!("Standard_DomainError: GeomFill_SectionLaw::CirclSection");
    }
}
