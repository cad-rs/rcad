//! OCCT HLRBRep_FaceIterator (TKHLR HLRBRep package).
//!
//! 1:1 translation of `HLRBRep_FaceIterator.hxx` (L30-77) + `.cxx`
//! (L26-56) + `.lxx` (L22-96): an exploration iterator over the wires and
//! edges of a face.
//!
//! Architecture note (the HLRBRep_Data value-member form): OCCT
//! `myFaceItr1` / `myFaceItr2` are default-constructed members of
//! HLRBRep_Data — the `myWires` handle stays null until `InitEdge(fd)`
//! binds `fd.Wires()`.  The rcad iterator mirrors that with a raw
//! `*mut WiresBlock` (the `HLRBRep_Surface::myProj` /
//! `Intersector::mySurface` raw-handle precedent): the constructor keeps
//! the OCCT null-handle default and [`FaceIterator::init_edge`] performs
//! the `myWires = fd.Wires()` binding.  The per-access accessors keep the
//! OCCT const receivers (the raw handle needs no `&mut`).

use std::marker::PhantomData;
use std::sync::Arc;

use rcad_kernel::topods::Orientation;

use crate::hlr::algo::edges_block::EdgesBlock;
use crate::hlr::algo::wires_block::WiresBlock;
use crate::hlr::brep::face_data::FaceData;

/// OCCT HLRBRep_FaceIterator.
pub struct FaceIterator<'a> {
    i_wire: i32,
    nb_wires: i32,
    i_edge: i32,
    nb_edges: i32,
    /// OCCT `occ::handle<HLRAlgo_WiresBlock> myWires` — null until InitEdge
    /// (the default-constructed member handle).
    my_wires: *mut WiresBlock,
    /// The lifetime parameter of the OCCT value-member form (unused by the
    /// raw handle; kept so `FaceIterator<'a>` stays the Data member type).
    marker: PhantomData<&'a mut WiresBlock>,
}

impl<'a> FaceIterator<'a> {
    /// OCCT HLRBRep_FaceIterator() — cxx L26 (`= default`; the counters
    /// stay uninitialized, the handle is null).  `wires` is the OCCT null
    /// handle stand-in: pass `std::ptr::null_mut()` for the default state;
    /// [`FaceIterator::init_edge`] binds the real `fd.Wires()`.
    pub fn new(wires: *mut WiresBlock) -> Self {
        FaceIterator {
            i_wire: 0,
            nb_wires: 0,
            i_edge: 0,
            nb_edges: 0,
            my_wires: wires,
            marker: PhantomData,
        }
    }

    /// OCCT InitEdge(HLRBRep_FaceData& fd) — cxx L30-39: begin an
    /// exploration of the edges of the face `fd` (`myWires = fd.Wires()`).
    pub fn init_edge(&mut self, fd: *mut FaceData<'_>) {
        self.i_wire = 0;
        // myWires = fd.Wires();
        self.my_wires = unsafe { (*fd).wires() };
        self.nb_wires = unsafe { (*self.my_wires).nb_wires() } as i32;

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
                // myEdges = myWires->Wire(iWire).
                self.nb_edges =
                    unsafe { (*self.my_wires).wire(self.i_wire as usize) }.nb_edges() as i32;
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

    /// OCCT Wire — lxx L51-54: the edges of the current wire (the shared
    /// handle is handed out mutable, as the OCCT callers mutate the block
    /// through it).
    pub fn wire(&self) -> &mut EdgesBlock {
        unsafe { (*self.my_wires).wire(self.i_wire as usize) }
    }

    /// OCCT Edge — lxx L58-61.
    pub fn edge(&self) -> i32 {
        unsafe { (*self.my_wires).wire(self.i_wire as usize) }.edge(self.i_edge as usize)
    }

    /// OCCT Orientation — lxx L65-68.
    pub fn orientation(&self) -> Orientation {
        unsafe { (*self.my_wires).wire(self.i_wire as usize) }.orientation(self.i_edge as usize)
    }

    /// OCCT OutLine — lxx L72-75.
    pub fn out_line(&self) -> bool {
        unsafe { (*self.my_wires).wire(self.i_wire as usize) }.out_line(self.i_edge as usize)
    }

    /// OCCT Internal — lxx L79-82.
    pub fn internal(&self) -> bool {
        unsafe { (*self.my_wires).wire(self.i_wire as usize) }.internal(self.i_edge as usize)
    }

    /// OCCT Double — lxx L86-89.
    pub fn double(&self) -> bool {
        unsafe { (*self.my_wires).wire(self.i_wire as usize) }.double(self.i_edge as usize)
    }

    /// OCCT IsoLine — lxx L93-96.
    pub fn iso_line(&self) -> bool {
        unsafe { (*self.my_wires).wire(self.i_wire as usize) }.iso_line(self.i_edge as usize)
    }

    /// The test-only InitEdge form: binds the WiresBlock directly (the
    /// rcad FaceData-less construction of the iterator anchors).
    #[cfg(test)]
    fn init_edge_wires(&mut self, wires: *mut WiresBlock) {
        self.i_wire = 0;
        self.my_wires = wires;
        self.nb_wires = unsafe { (*self.my_wires).nb_wires() } as i32;

        self.i_edge = 0;
        self.nb_edges = 0;
        self.next_edge();
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
        let mut it = FaceIterator::new(std::ptr::null_mut());
        it.init_edge_wires(&mut wb);
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
        let mut it = FaceIterator::new(std::ptr::null_mut());
        it.init_edge_wires(&mut wb);
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

    /// The default-constructed iterator (cxx L26 `= default`): the handle
    /// is null and the counters sit at their neutral 0 (the OCCT members
    /// are uninitialized; the exploration is only defined after InitEdge —
    /// BeginningOfWire (iEdge == 1) is false).
    #[test]
    fn face_iterator_default_state() {
        let it = FaceIterator::<'_>::new(std::ptr::null_mut());
        assert!(!it.beginning_of_wire());
    }
}
