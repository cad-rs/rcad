//! OCCT GeomFill_LocationLaw (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_LocationLaw.hxx + GeomFill_LocationLaw.cxx (whole file L21-115):
//! the abstract base trait for location laws along a swept path.
//!
//! Architecture mapping: the OCCT abstract class is a Rust trait; the
//! non-pure virtual bodies keep their OCCT defaults (Standard_NotImplemented
//! raises become panics carrying the OCCT message).

use glam::{DVec2, DVec3};

use rcad_kernel::geom::Curve3;
use rcad_kernel::math::GeomAbsShape;

use super::gp_mat::GpMat;
use super::trihedron_law::PipeError;

/// OCCT GeomFill_LocationLaw.
pub trait LocationLaw {
    /// OCCT SetCurve — pure virtual.
    fn set_curve(&mut self, c: Curve3) -> bool;

    /// OCCT GetCurve — pure virtual.
    fn get_curve(&self) -> Option<Curve3>;

    /// OCCT SetTrsf — pure virtual.
    fn set_trsf(&mut self, transfo: GpMat);

    /// OCCT Copy — pure virtual.
    fn copy_law(&self) -> Box<dyn LocationLaw>;

    /// OCCT D0(Param, M, V) — pure virtual: transform an point P in MP+V.
    fn d0(&self, param: f64, m: &mut GpMat, v: &mut DVec3) -> bool;

    /// OCCT D0(Param, M, V, Pnts2d) — the 2d-array overload.
    fn d0_2d(&self, param: f64, m: &mut GpMat, v: &mut DVec3, pnts2d: &mut [DVec2]) -> bool;

    /// OCCT D1 (GeomFill_LocationLaw.cxx L21-31) — default raises
    /// Standard_NotImplemented.
    #[allow(clippy::too_many_arguments)]
    fn d1(
        &self,
        _param: f64,
        _m: &mut GpMat,
        _v: &mut DVec3,
        _dm: &mut GpMat,
        _dv: &mut DVec3,
        _pnts2d: &mut [DVec2],
        _vecs2d: &mut [DVec2],
    ) -> bool {
        panic!("Standard_NotImplemented: GeomFill_LocationLaw::D1");
    }

    /// OCCT D2 (GeomFill_LocationLaw.cxx L33-43) — default raises
    /// Standard_NotImplemented.
    #[allow(clippy::too_many_arguments)]
    fn d2(
        &self,
        _param: f64,
        _m: &mut GpMat,
        _v: &mut DVec3,
        _dm: &mut GpMat,
        _dv: &mut DVec3,
        _d2m: &mut GpMat,
        _d2v: &mut DVec3,
        _pnts2d: &mut [DVec2],
        _d1vecs2d: &mut [DVec2],
        _d2vecs2d: &mut [DVec2],
    ) -> bool {
        panic!("Standard_NotImplemented: GeomFill_LocationLaw::D2");
    }

    /// OCCT Nb2dCurves (GeomFill_LocationLaw.cxx L45-58).
    fn nb_2d_curves(&self) -> usize {
        let mut n = self.trace_number();
        if self.has_first_restriction() {
            n += 1;
        }
        if self.has_last_restriction() {
            n += 1;
        }
        n
    }

    /// OCCT HasFirstRestriction (L60-63) — default false.
    fn has_first_restriction(&self) -> bool {
        false
    }

    /// OCCT HasLastRestriction (L65-68) — default false.
    fn has_last_restriction(&self) -> bool {
        false
    }

    /// OCCT TraceNumber (L70-73) — default 0.
    fn trace_number(&self) -> usize {
        0
    }

    /// OCCT ErrorStatus (L77-80) — default GeomFill_PipeOk.
    fn error_status(&self) -> PipeError {
        PipeError::PipeOk
    }

    /// OCCT NbIntervals — pure virtual.
    fn nb_intervals(&self, s: GeomAbsShape) -> usize;

    /// OCCT Intervals — pure virtual.
    fn intervals(&self, t: &mut Vec<f64>, s: GeomAbsShape);

    /// OCCT SetInterval — pure virtual.
    fn set_interval(&mut self, first: f64, last: f64);

    /// OCCT GetInterval — pure virtual.
    fn get_interval(&self, first: &mut f64, last: &mut f64);

    /// OCCT GetDomain — pure virtual.
    fn get_domain(&self, first: &mut f64, last: &mut f64);

    /// OCCT Resolution (L82-90) — default raises Standard_NotImplemented.
    fn resolution(&self, _index: usize, _tol: f64, _tolu: &mut f64, _tolv: &mut f64) {
        panic!("Standard_NotImplemented: GeomFill_LocationLaw::Resolution");
    }

    /// OCCT SetTolerance (L92-95) — "Ne fait rien !!".
    fn set_tolerance(&mut self, _tol3d: f64, _tol2d: f64) {}

    /// OCCT GetMaximalNorm — pure virtual.
    fn get_maximal_norm(&self) -> f64;

    /// OCCT GetAverageLaw — pure virtual.
    fn get_average_law(&self, am: &mut GpMat, av: &mut DVec3);

    /// OCCT IsTranslation (L97-100) — default false.
    fn is_translation(&self, _error: &mut f64) -> bool {
        false
    }

    /// OCCT IsRotation (L102-105) — default false.
    fn is_rotation(&self, _error: &mut f64) -> bool {
        false
    }

    /// OCCT Rotation (L107-114) — default raises Standard_NotImplemented
    /// (the OCCT message keeps the SectionLaw label — a literal quirk).
    fn rotation(&self, _centre: &mut DVec3) {
        panic!("Standard_NotImplemented: GeomFill_SectionLaw::Rotation");
    }
}
