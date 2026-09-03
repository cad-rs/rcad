//! OCCT GeomFill_PipeError + GeomFill_TrihedronLaw (TKGeomAlgo/GeomFill) —
//! 1:1 port of GeomFill_PipeError.hxx L20-25, GeomFill_TrihedronLaw.hxx
//! L34-115 and GeomFill_TrihedronLaw.cxx (whole file).
//!
//! Architecture mapping: `Adaptor3d_Curve` -> rcad `Curve3` (enum dispatch
//! replaces virtual dispatch).  `myCurve->Trim(First, Last, 0)` ->
//! `Curve3::Trimmed(TrimmedCurve3)`, which keeps the base parametrization
//! exactly like OCCT `Trim(..., Adjust = False)`.  The OCCT inheritance is
//! expressed as a Rust trait: pure virtuals (= 0) are required methods,
//! the Standard_Transient base members are exposed through the embedded
//! [`TrihedronLawBase`].

use glam::DVec3;

use rcad_kernel::geom::{Curve3, CurveEval, TrimmedCurve3};
use rcad_kernel::math::GeomAbsShape;

/// OCCT GeomFill_PipeError (GeomFill_PipeError.hxx L20-25).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PipeError {
    PipeOk,
    PipeNotOk,
    PlaneNotIntersectGuide,
    ImpossibleContact,
}

/// The OCCT protected base-class members (GeomFill_TrihedronLaw.hxx
/// L112-114): `myCurve` and `myTrimmed`.  Concrete laws embed this struct.
#[derive(Debug, Clone, Default)]
pub struct TrihedronLawBase {
    pub(crate) my_curve: Option<Curve3>,
    pub(crate) my_trimmed: Option<Curve3>,
}

/// OCCT Adaptor3d_Curve::FirstParameter.
pub(crate) fn curve_first_parameter(c: &Curve3) -> f64 {
    match c {
        Curve3::Trimmed(tc) => tc.first,
        other => other.default_domain()[0],
    }
}

/// OCCT Adaptor3d_Curve::LastParameter.
pub(crate) fn curve_last_parameter(c: &Curve3) -> f64 {
    match c {
        Curve3::Trimmed(tc) => tc.last,
        other => other.default_domain()[1],
    }
}

/// OCCT GeomFill_TrihedronLaw.
pub trait TrihedronLaw {
    /// OCCT protected member `myCurve`.
    fn my_curve(&self) -> &Option<Curve3>;
    /// OCCT protected member `myTrimmed`.
    fn my_trimmed(&self) -> &Option<Curve3>;
    /// OCCT protected member `myCurve` (write access for the base class).
    fn set_my_curve(&mut self, c: Curve3);
    /// OCCT protected member `myTrimmed` (write access for the base class).
    fn set_my_trimmed(&mut self, c: Option<Curve3>);

    /// OCCT SetCurve (GeomFill_TrihedronLaw.cxx L29-35) — the default body;
    /// overriding implementations reuse
    /// [`trihedron_law_base_set_curve`].
    fn set_curve(&mut self, c: Curve3) -> bool {
        self.set_my_curve(c.clone());
        // myTrimmed = myCurve;
        self.set_my_trimmed(Some(c));
        true
    }

    /// OCCT Copy() — pure virtual.
    fn copy_law(&self) -> Box<dyn TrihedronLaw>;

    /// OCCT ErrorStatus() — "Returns PipeOk (default implementation)".
    fn error_status(&self) -> PipeError {
        PipeError::PipeOk
    }

    /// OCCT D0() — pure virtual: compute Trihedron on curve at Param.
    fn d0(
        &self,
        param: f64,
        tangent: &mut DVec3,
        normal: &mut DVec3,
        binormal: &mut DVec3,
    ) -> bool;

    /// OCCT D1() — default throws Standard_NotImplemented.
    fn d1(
        &self,
        _param: f64,
        _tangent: &mut DVec3,
        _dtangent: &mut DVec3,
        _normal: &mut DVec3,
        _dnormal: &mut DVec3,
        _binormal: &mut DVec3,
        _dbinormal: &mut DVec3,
    ) -> bool {
        // OCCT message text kept verbatim (the base D1 raises with the
        // D2 label — a literal in GeomFill_TrihedronLaw.cxx).
        panic!(" GeomFill_TrihedronLaw::D2");
    }

    /// OCCT D2() — default throws Standard_NotImplemented.
    fn d2(
        &self,
        _param: f64,
        _tangent: &mut DVec3,
        _dtangent: &mut DVec3,
        _d2tangent: &mut DVec3,
        _normal: &mut DVec3,
        _dnormal: &mut DVec3,
        _d2normal: &mut DVec3,
        _binormal: &mut DVec3,
        _dbinormal: &mut DVec3,
        _d2binormal: &mut DVec3,
    ) -> bool {
        panic!(" GeomFill_TrihedronLaw::D2");
    }

    /// OCCT NbIntervals() — pure virtual.
    fn nb_intervals(&self, s: GeomAbsShape) -> usize;

    /// OCCT Intervals() — pure virtual.
    fn intervals(&self, t: &mut Vec<f64>, s: GeomAbsShape);

    /// OCCT SetInterval (GeomFill_TrihedronLaw.cxx): myTrimmed =
    /// myCurve->Trim(First, Last, 0).
    fn set_interval(&mut self, first: f64, last: f64) {
        let curve = self
            .my_curve()
            .clone()
            .expect("GeomFill_TrihedronLaw::SetInterval with null curve");
        self.set_my_trimmed(Some(Curve3::Trimmed(TrimmedCurve3::new(curve, first, last))));
    }

    /// OCCT GetInterval (GeomFill_TrihedronLaw.cxx).
    fn get_interval(&self, first: &mut f64, last: &mut f64) {
        let trimmed = self.my_trimmed().as_ref().expect("null myTrimmed");
        *first = curve_first_parameter(trimmed);
        *last = curve_last_parameter(trimmed);
    }

    /// OCCT GetAverageLaw() — pure virtual.
    fn get_average_law(&self, atangent: &mut DVec3, anormal: &mut DVec3, abinormal: &mut DVec3);

    /// OCCT IsConstant() — "Return False by Default".
    fn is_constant(&self) -> bool {
        false
    }

    /// OCCT IsOnlyBy3dCurve() — "Return False by Default".
    fn is_only_by3d_curve(&self) -> bool {
        false
    }
}


/// OCCT GeomFill_TrihedronLaw::SetCurve base implementation — a Rust trait
/// default body cannot be invoked explicitly when a type overrides the
/// method (the fully-qualified call resolves to the override), so the base
/// logic lives in a free function.
pub fn trihedron_law_base_set_curve(law: &mut dyn TrihedronLaw, c: Curve3) -> bool {
    law.set_my_curve(c.clone());
    // myTrimmed = myCurve;
    law.set_my_trimmed(Some(c));
    true
}

/// OCCT GeomFill_TrihedronLaw::SetInterval base implementation — same
/// explicit-base-call pattern as [`trihedron_law_base_set_curve`].
pub fn trihedron_law_base_set_interval(law: &mut dyn TrihedronLaw, first: f64, last: f64) {
    let curve = law
        .my_curve()
        .clone()
        .expect("GeomFill_TrihedronLaw::SetInterval with null curve");
    law.set_my_trimmed(Some(Curve3::Trimmed(TrimmedCurve3::new(curve, first, last))));
}
