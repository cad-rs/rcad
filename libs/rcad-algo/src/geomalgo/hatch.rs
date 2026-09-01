//! OCCT Geom2dHatch package — the 2D hatching elements container and the
//! hatch intersector local geometry.
//!
//! 1:1 translations:
//!   - Geom2dHatch_Element.hxx/.cxx — an element (2D curve + orientation).
//!   - Geom2dHatch_Elements.hxx/.cxx — the element map with the wire/edge
//!     traversal (InitWires/InitEdges/CurrentEdge/MoreEdges/NextEdge).
//!   - Geom2dHatch_Intersector.cxx L74-101 — LocalGeometry (tangent, normal,
//!     curvature of a 2D curve at a parameter via GeomLProp_CLProps2d).

use glam::DVec2;
use rcad_kernel::base::geom_lprop::ClProps2d;
use rcad_kernel::core::precision::PCONFUSION;
use rcad_kernel::geom::{Curve2d, Curve2dEval};
use rcad_kernel::topods::Orientation;

/// OCCT Geom2dHatch_Element — a hatching element: the 2D curve and its
/// orientation.
#[derive(Debug, Clone)]
pub struct HatchElement {
    curve: Curve2d,
    orientation: Orientation,
}

impl HatchElement {
    /// OCCT Geom2dHatch_Element(Curve, Orientation = TopAbs_FORWARD)
    /// (Geom2dHatch_Element.cxx L29-34).
    pub fn new(curve: Curve2d, orientation: Orientation) -> Self {
        HatchElement {
            curve,
            orientation,
        }
    }

    /// OCCT Geom2dHatch_Element::Curve().
    pub fn curve(&self) -> &Curve2d {
        &self.curve
    }

    /// OCCT Geom2dHatch_Element::Orientation().
    pub fn orientation(&self) -> Orientation {
        self.orientation
    }

    /// OCCT Geom2dHatch_Element::Orientation(Orientation).
    pub fn set_orientation(&mut self, orientation: Orientation) {
        self.orientation = orientation;
    }
}

/// OCCT Geom2dHatch_Elements — a data map of hatching elements with the
/// wire/edge traversal state (Geom2dHatch_Elements.cxx whole).
#[derive(Debug, Clone)]
pub struct HatchElements {
    map: Vec<(usize, HatchElement)>,
    // Traversal state.
    num_wire: usize,
    num_edge: usize,
    cur_edge: usize,
    cur_edge_par: f64,
    iter_pos: usize,
}

/// OCCT Geom2dHatch_Elements.cxx L26-28 — the probing parameter constants.
const PROBING_START: f64 = 0.123;
const PROBING_END: f64 = 0.8;
const PROBING_STEP: f64 = 0.2111;

impl HatchElements {
    /// OCCT Geom2dHatch_Elements() (L41-47).
    pub fn new() -> Self {
        HatchElements {
            map: Vec::new(),
            num_wire: 0,
            num_edge: 0,
            cur_edge: 1,
            cur_edge_par: PROBING_START,
            iter_pos: 0,
        }
    }

    /// OCCT Geom2dHatch_Elements::Clear (L49-52).
    pub fn clear(&mut self) {
        self.map.clear();
    }

    /// OCCT Geom2dHatch_Elements::IsBound (L54-57).
    pub fn is_bound(&self, k: usize) -> bool {
        self.map.iter().any(|(key, _)| *key == k)
    }

    /// OCCT Geom2dHatch_Elements::UnBind (L59-62).
    pub fn unbind(&mut self, k: usize) -> bool {
        if let Some(pos) = self.map.iter().position(|(key, _)| *key == k) {
            self.map.remove(pos);
            true
        } else {
            false
        }
    }

    /// OCCT Geom2dHatch_Elements::Bind (L64-67) — returns false when the key
    /// is already bound.
    pub fn bind(&mut self, k: usize, element: HatchElement) -> bool {
        if self.is_bound(k) {
            return false;
        }
        self.map.push((k, element));
        true
    }

    /// OCCT Geom2dHatch_Elements::Find (L69-72) — panics when not bound.
    pub fn find(&self, k: usize) -> &HatchElement {
        self.map
            .iter()
            .find(|(key, _)| *key == k)
            .map(|(_, e)| e)
            .expect("Geom2dHatch_Elements::Find - not bound")
    }

    /// OCCT Geom2dHatch_Elements::InitWires (L201-204).
    pub fn init_wires(&mut self) {
        self.num_wire = 0;
    }

    /// OCCT Geom2dHatch_Elements::MoreWires (L239-242) — a single wire.
    pub fn more_wires(&self) -> bool {
        self.num_wire == 0
    }

    /// OCCT Geom2dHatch_Elements::NextWire (L246-249).
    pub fn next_wire(&mut self) {
        self.num_wire += 1;
    }

    /// OCCT Geom2dHatch_Elements::InitEdges (L215-219).
    pub fn init_edges(&mut self) {
        self.num_edge = 0;
        self.iter_pos = 0;
    }

    /// OCCT Geom2dHatch_Elements::MoreEdges (L253-256).
    pub fn more_edges(&self) -> bool {
        self.iter_pos < self.map.len()
    }

    /// OCCT Geom2dHatch_Elements::NextEdge (L260-262).
    pub fn next_edge(&mut self) {
        self.iter_pos += 1;
    }

    /// OCCT Geom2dHatch_Elements::CurrentEdge (L230-235) — the current edge
    /// curve and orientation.
    pub fn current_edge(&self) -> (Curve2d, Orientation) {
        let (_, element) = &self.map[self.iter_pos];
        (element.curve.clone(), element.orientation)
    }
}

impl Default for HatchElements {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT Geom2dHatch_Intersector — the hatch line/edge intersector.  The GTests
/// exercise only the local geometry computation.
#[derive(Debug, Clone)]
pub struct HatchIntersector {
    confusion_tolerance: f64,
    tangency_tolerance: f64,
}

impl HatchIntersector {
    /// OCCT Geom2dHatch_Intersector() (Geom2dHatch_Intersector.cxx L31-35).
    pub fn new() -> Self {
        HatchIntersector {
            confusion_tolerance: 0.0,
            tangency_tolerance: 0.0,
        }
    }

    /// OCCT Geom2dHatch_Intersector(Confusion, Tangency)
    /// (Geom2dHatch_Intersector.lxx L19-24).
    pub fn with_tolerances(confusion: f64, tangency: f64) -> Self {
        HatchIntersector {
            confusion_tolerance: confusion,
            tangency_tolerance: tangency,
        }
    }

    /// OCCT Geom2dHatch_Intersector::ConfusionTolerance (lxx L31-34).
    pub fn confusion_tolerance(&self) -> f64 {
        self.confusion_tolerance
    }

    /// OCCT Geom2dHatch_Intersector::TangencyTolerance (lxx L56-59).
    pub fn tangency_tolerance(&self) -> f64 {
        self.tangency_tolerance
    }

    /// OCCT Geom2dHatch_Intersector::LocalGeometry (Geom2dHatch_Intersector.cxx
    /// L74-101) — tangent, normal and curvature of the 2D curve at U via
    /// GeomLProp_CLProps2d.  Returns (Tang, Norm, C).
    pub fn local_geometry(&self, curve: &Curve2d, u: f64) -> (DVec2, DVec2, f64) {
        let mut prop = ClProps2d::with_param(curve, u, 2, PCONFUSION);

        let mut c = 0.0;
        let tang = if prop.is_tangent_defined() {
            c = prop.curvature();
            prop.tangent().unwrap_or(DVec2::new(1.0, 0.0))
        } else {
            DVec2::new(1.0, 0.0)
        };

        let norm = if c > PCONFUSION && c < f64::MAX {
            prop.normal().unwrap_or(DVec2::new(tang.y, -tang.x))
        } else {
            DVec2::new(tang.y, -tang.x)
        };

        (tang, norm, c)
    }
}

impl Default for HatchIntersector {
    fn default() -> Self {
        Self::new()
    }
}
