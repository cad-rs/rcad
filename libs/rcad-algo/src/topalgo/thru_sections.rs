//! OCCT BRepOffsetAPI_ThruSections (TKOffset) — loft between section wires.
//!
//! OCCT BRepOffsetAPI_ThruSections builds a shell/solid interpolating the
//! given section wires (BRepOffsetAPI_ThruSections.cxx, via BRepFill).
//! The rcad port keeps the OCCT method surface (new / add_wire / build /
//! is_done / shape); the actual loft (BRepFill) is a later work item — the
//! current implementation only records the wires and reports NotDone, so the
//! ported GTests stay #[ignore]d until the BRepFill port lands.

use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topo::topods::BRep;

/// OCCT BRepOffsetAPI_ThruSections — section-wire loft builder.
#[derive(Debug, Default)]
pub struct ThruSections {
    /// OCCT isSolid: whether to build a solid (vs an open shell).
    pub is_solid: bool,
    /// OCCT ruled: ruled (vs smooth) interpolation.
    pub ruled: bool,
    /// OCCT pres3d: 3D tolerance for the loft.
    pub pres3d: f64,
    wires: Vec<Shape>,
    done: bool,
    shape: Option<Shape>,
}

impl ThruSections {
    /// OCCT BRepOffsetAPI_ThruSections(isSolid, ruled, pres3d).
    pub fn new(is_solid: bool, ruled: bool, pres3d: f64) -> Self {
        ThruSections {
            is_solid,
            ruled,
            pres3d,
            wires: Vec::new(),
            done: false,
            shape: None,
        }
    }

    /// OCCT BRepOffsetAPI_ThruSections::AddWire(wire).
    pub fn add_wire(&mut self, wire: Shape) {
        self.wires.push(wire);
    }

    /// OCCT BRepOffsetAPI_ThruSections::Build().
    ///
    /// NOT YET IMPLEMENTED: the BRepFill loft (BRepOffsetAPI_ThruSections.cxx)
    /// is a later work item.  Until then the builder stays NotDone and
    /// `shape` is None; the ported GTests are #[ignore]d accordingly.
    pub fn build(&mut self, _brep: &mut BRep) {
        self.done = false;
        self.shape = None;
    }

    /// OCCT BRepOffsetAPI_ThruSections::IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT BRepOffsetAPI_ThruSections::Shape().
    pub fn shape(&self) -> Option<Shape> {
        self.shape.clone()
    }
}
