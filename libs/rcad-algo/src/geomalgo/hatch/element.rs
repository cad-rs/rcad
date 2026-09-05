//! OCCT Geom2dHatch_Element (TKGeomAlgo/Geom2dHatch).
//!
//! Geom2dHatch_Element.hxx L27-53 + .cxx L23-72 — one hatching boundary
//! element: a bounded 2D curve (OCCT Geom2dAdaptor_Curve, rcad [`Curve2d`])
//! and its orientation in the boundary.

use rcad_kernel::geom::Curve2d;
use rcad_kernel::topods::Orientation;

/// OCCT Geom2dHatch_Element — an element of the hatched face boundary.
#[derive(Debug, Clone)]
pub struct HatchElement {
    /// OCCT Geom2dAdaptor_Curve myCurve.
    my_curve: Curve2d,
    /// OCCT TopAbs_Orientation myOrientation.
    my_orientation: Orientation,
}

impl HatchElement {
    /// OCCT Geom2dHatch_Element() (cxx L23) — the default element.  OCCT
    /// leaves the curve adaptor with a null curve handle (using it is
    /// undefined); rcad initializes a neutral degenerate line so the value
    /// stays deterministic.
    pub fn empty() -> Self {
        HatchElement {
            my_curve: Curve2d::Line(rcad_kernel::geom::Line2d::new(glam::DVec2::ZERO, glam::DVec2::X)),
            my_orientation: Orientation::Forward,
        }
    }

    /// OCCT Geom2dHatch_Element(Curve, Orientation = TopAbs_FORWARD)
    /// (cxx L27-32) — creates an element.
    pub fn new(curve: Curve2d, orientation: Orientation) -> Self {
        HatchElement {
            my_curve: curve,
            my_orientation: orientation,
        }
    }

    /// OCCT Curve (cxx L39-42) — returns the curve associated to the element.
    pub fn curve(&self) -> &Curve2d {
        &self.my_curve
    }

    /// OCCT ChangeCurve (cxx L49-52).
    pub fn change_curve(&mut self) -> &mut Curve2d {
        &mut self.my_curve
    }

    /// OCCT Orientation(Orientation) setter (cxx L59-62).
    pub fn set_orientation(&mut self, orientation: Orientation) {
        self.my_orientation = orientation;
    }

    /// OCCT Orientation() getter (cxx L69-72).
    pub fn orientation(&self) -> Orientation {
        self.my_orientation
    }
}

impl Default for HatchElement {
    fn default() -> Self {
        Self::empty()
    }
}
