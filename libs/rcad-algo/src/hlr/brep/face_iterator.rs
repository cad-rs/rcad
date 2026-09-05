//! OCCT HLRBRep_FaceIterator (TKHLR HLRBRep package).
//!
//! 1:1 translation of `HLRBRep_FaceIterator.hxx` (L30-77) + `.cxx`
//! (L26-56) + `.lxx` (L22-96): an exploration iterator over the wires and
//! edges of a face.
//!
//! Documented deferrals:
//! - OCCT `InitEdge(HLRBRep_FaceData& fd)` reads `fd.Wires()`; the
//!   HLRBRep_FaceData structure lands with Stage 3f, so the rcad iterator
//!   binds the [`WiresBlock`] (what `fd.Wires()` returns) at construction.
//! - The cached `myEdges` handle is re-derived per access through
//!   `myWires->Wire(iWire)`: `WiresBlock::wire` needs `&mut`, and the
//!   accessed block is invariantly `Wire(iWire)` between `NextEdge` calls,
//!   so the lookups are equivalent.  The per-access accessors therefore
//!   take `&mut self` where OCCT has const methods.

use rcad_kernel::topods::Orientation;

use crate::hlr::algo::edges_block::EdgesBlock;
use crate::hlr::algo::wires_block::WiresBlock;

/// OCCT HLRBRep_FaceIterator.
pub struct FaceIterator<'a> {
    i_wire: i32,
    nb_wires: i32,
    i_edge: i32,
    nb_edges: i32,
    my_wires: &'a mut WiresBlock,
}

impl<'a> FaceIterator<'a> {
    /// OCCT HLRBRep_FaceIterator() — cxx L26 (`= default`; the ints and
    /// handles stay uninitialized/null).  The rcad iterator binds the
    /// WiresBlock reference at construction (see the module deferrals);
    /// the counters keep the neutral default 0.
    pub fn new(wires: &'a mut WiresBlock) -> Self {
        FaceIterator {
            i_wire: 0,
            nb_wires: 0,
            i_edge: 0,
            nb_edges: 0,
            my_wires: wires,
        }
    }

    /// OCCT InitEdge(HLRBRep_FaceData& fd) — cxx L30-39: begin an
    /// exploration of the edges of the face `fd` (the `fd.Wires()` read is
    /// the construction-time binding).
    pub fn init_edge(&mut self) {
        self.i_wire = 0;
        self.nb_wires = self.my_wires.nb_wires() as i32;

        self.i_edge = 0;
        self.nb_edges = 0;
        self.next_edge();
    }

    /// OCCT NextEdge — cxx L43-56.
    pub fn next_edge(&mut self) {
        self.i_edge += 1;
        if self.i_edge > self.nb_edges {
            self.i_wire += 1;
            if self.i_wire <= self.nb_wires {
                self.i_edge = 1;
                // myEdges = myWires->Wire(iWire) — fused into the per-access
                // lookups (see the module deferrals).
                self.nb_edges = self.my_wires.wire(self.i_wire as usize).nb_edges() as i32;
            }
        }
    }

    /// OCCT MoreEdge — lxx L22-25.
    pub fn more_edge(&self) -> bool {
        self.i_wire <= self.nb_wires
    }

    /// OCCT BeginningOfWire — lxx L29-32: true if the current edge is the
    /// first of a wire.
    pub fn beginning_of_wire(&self) -> bool {
        self.i_edge == 1
    }

    /// OCCT EndOfWire — lxx L36-39: true if the current edge is the last of
    /// a wire.
    pub fn end_of_wire(&self) -> bool {
        self.i_edge == self.nb_edges
    }

    /// OCCT SkipWire — lxx L43-47: skip the current wire in the
    /// exploration.
    pub fn skip_wire(&mut self) {
        self.i_edge = self.nb_edges;
        self.next_edge();
    }

    /// OCCT Wire — lxx L51-54: the edges of the current wire.
    pub fn wire(&mut self) -> &mut EdgesBlock {
        self.my_wires.wire(self.i_wire as usize)
    }

    /// OCCT Edge — lxx L58-61.
    pub fn edge(&mut self) -> i32 {
        self.my_wires.wire(self.i_wire as usize).edge(self.i_edge as usize)
    }

    /// OCCT Orientation — lxx L65-68.
    pub fn orientation(&mut self) -> Orientation {
        self.my_wires
            .wire(self.i_wire as usize)
            .orientation(self.i_edge as usize)
    }

    /// OCCT OutLine — lxx L72-75.
    pub fn out_line(&mut self) -> bool {
        self.my_wires.wire(self.i_wire as usize).out_line(self.i_edge as usize)
    }

    /// OCCT Internal — lxx L79-82.
    pub fn internal(&mut self) -> bool {
        self.my_wires
            .wire(self.i_wire as usize)
            .internal(self.i_edge as usize)
    }

    /// OCCT Double — lxx L86-89.
    pub fn double(&mut self) -> bool {
        self.my_wires.wire(self.i_wire as usize).double(self.i_edge as usize)
    }

    /// OCCT IsoLine — lxx L93-96.
    pub fn iso_line(&mut self) -> bool {
        self.my_wires
            .wire(self.i_wire as usize)
            .iso_line(self.i_edge as usize)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A face with two wires: wire 1 carries edges 1,2; wire 2 carries
    /// edge 3 (HLRAlgo_EdgesBlock 1-based indices).
    fn two_wire_block() -> WiresBlock {
        let mut wb = WiresBlock::new(2);
        let mut w1 = EdgesBlock::new(2);
        w1.set_edge(1, 10);
        w1.set_edge(2, 11);
        w1.set_orientation(1, Orientation::Forward);
        w1.set_orientation(2, Orientation::Reversed);
        w1.set_out_line(1, true);
        w1.set_internal(2, true);
        let mut w2 = EdgesBlock::new(1);
        w2.set_edge(1, 12);
        w2.set_orientation(1, Orientation::Internal);
        w2.set_double(1, true);
        w2.set_iso_line(1, true);
        wb.set(1, w1);
        wb.set(2, w2);
        wb
    }

    /// OCCT InitEdge / NextEdge / MoreEdge / BeginningOfWire / EndOfWire
    /// (cxx L30-56, lxx L22-39): the iteration walks wire 1 edges 10, 11
    /// then wire 2 edge 12.
    #[test]
    fn face_iterator_edge_walk_order() {
        let mut wb = two_wire_block();
        let mut it = FaceIterator::new(&mut wb);
        it.init_edge();
        // Edge 1 of wire 1.
        assert!(it.more_edge());
        assert!(it.beginning_of_wire());
        assert!(!it.end_of_wire());
        assert_eq!(it.edge(), 10);
        assert_eq!(it.orientation(), Orientation::Forward);
        assert!(it.out_line());
        assert!(!it.internal());
        it.next_edge();
        // Edge 2 of wire 1 (last of the wire).
        assert!(it.more_edge());
        assert!(!it.beginning_of_wire());
        assert!(it.end_of_wire());
        assert_eq!(it.edge(), 11);
        assert_eq!(it.orientation(), Orientation::Reversed);
        assert!(it.internal());
        assert!(!it.out_line());
        it.next_edge();
        // Edge 1 of wire 2.
        assert!(it.more_edge());
        assert!(it.beginning_of_wire());
        assert!(it.end_of_wire());
        assert_eq!(it.edge(), 12);
        assert_eq!(it.orientation(), Orientation::Internal);
        assert!(it.double());
        assert!(it.iso_line());
        it.next_edge();
        assert!(!it.more_edge());
    }

    /// OCCT SkipWire (lxx L43-47): from the first edge of wire 1 the
    /// exploration continues at the first edge of wire 2.
    #[test]
    fn face_iterator_skip_wire() {
        let mut wb = two_wire_block();
        let mut it = FaceIterator::new(&mut wb);
        it.init_edge();
        assert_eq!(it.edge(), 10);
        it.skip_wire();
        assert!(it.more_edge());
        assert!(it.beginning_of_wire());
        assert_eq!(it.edge(), 12);
        // Wire() returns the current wire block.
        assert_eq!(it.wire().nb_edges(), 1);
        it.next_edge();
        assert!(!it.more_edge());
    }
}
