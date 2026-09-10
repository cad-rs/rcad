//! OCCT GeomFill_TrihedronWithGuide (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_TrihedronWithGuide.hxx + GeomFill_TrihedronWithGuide.cxx (whole
//! file L22-33).  Architecture mapping: the OCCT abstract base is expressed
//! as a Rust trait extending [`TrihedronLaw`]; the protected members
//! (myGuide / myTrimG / myCurPointOnGuide) live in the embedded
//! [`TrihedronWithGuideBase`].  `myCurPointOnGuide` is a D0/D1/D2 output
//! cache written through `&self`, hence the `Cell`.

use std::cell::Cell;

use glam::DVec3;

use rcad_kernel::geom::Curve3;

use super::trihedron_law::TrihedronLaw;

/// OCCT GeomFill_TrihedronWithGuide protected members
/// (GeomFill_TrihedronWithGuide.hxx L47-50).
#[derive(Debug, Clone)]
pub struct TrihedronWithGuideBase {
    /// OCCT handle(Adaptor3d_Curve) myGuide.
    pub(crate) my_guide: Option<Curve3>,
    /// OCCT handle(Adaptor3d_Curve) myTrimG.
    pub(crate) my_trim_g: Option<Curve3>,
    /// OCCT gp_Pnt myCurPointOnGuide (D0/D1/D2 output cache).
    pub(crate) my_cur_point_on_guide: Cell<DVec3>,
}

/// OCCT GeomFill_TrihedronWithGuide — the abstract base trait.
pub trait TrihedronWithGuide: TrihedronLaw {
    /// OCCT Guide() — the guide curve (pure virtual in OCCT).
    fn guide(&self) -> Option<Curve3>;

    /// OCCT Origine(Param1, Param2) — pure virtual.
    fn origine(&mut self, param1: f64, param2: f64);

    /// OCCT Copy() viewed through the concrete GeomFill_TrihedronWithGuide
    /// down-cast — the form consumed by GeomFill_LocationGuide::Copy
    /// (GeomFill_LocationGuide.cxx L520).
    fn copy_with_guide(&self) -> Box<dyn TrihedronWithGuide>;

    /// OCCT CurrentPointOnGuide (GeomFill_TrihedronWithGuide.cxx L27-32) —
    /// the current point on guide found by D0, D1 or D2.
    fn current_point_on_guide(&self) -> DVec3 {
        self.with_guide_base().my_cur_point_on_guide.get()
    }

    /// Access to the protected member block.
    fn with_guide_base(&self) -> &TrihedronWithGuideBase;
    /// Mutable access to the protected member block.
    fn with_guide_base_mut(&mut self) -> &mut TrihedronWithGuideBase;
}
