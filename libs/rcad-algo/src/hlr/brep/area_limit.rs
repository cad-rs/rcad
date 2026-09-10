//! OCCT HLRBRep_AreaLimit (TKHLR HLRBRep package).
//!
//! 1:1 translation of `HLRBRep_AreaLimit.hxx` (L30-87) + `.cxx` (L25-153).
//! The private nested class AreaLimit represents a vertex on the Edge with
//! the state on the left and the right.
//!
//! `occ::handle<HLRBRep_AreaLimit>` maps to [`SharedAreaLimit`]
//! (`Arc<AreaLimit>`); the OCCT mutations through a shared handle
//! (StateBefore/EdgeAfter setters from HLRBRep_EdgeBuilder, Previous/Next
//! linking, Clear) map to interior mutability — `Cell` for the Copy `State`
//! fields, `RefCell` for the two links.

use std::cell::{Cell, RefCell};
use std::sync::Arc;

use rcad_kernel::topods::State;

use crate::hlr::algo::intersection::Intersection;

/// OCCT `occ::handle<HLRBRep_AreaLimit>` — a shared, nullable handle.
pub type SharedAreaLimit = Arc<AreaLimit>;

/// OCCT HLRBRep_AreaLimit.
pub struct AreaLimit {
    my_vertex: Intersection,
    my_boundary: bool,
    my_interference: bool,
    my_state_before: Cell<State>,
    my_state_after: Cell<State>,
    my_edge_before: Cell<State>,
    my_edge_after: Cell<State>,
    my_previous: RefCell<Option<SharedAreaLimit>>,
    my_next: RefCell<Option<SharedAreaLimit>>,
}

impl AreaLimit {
    /// OCCT HLRBRep_AreaLimit(V, Boundary, Interference, StateBefore,
    /// StateAfter, EdgeBefore, EdgeAfter) — cxx L25-40 (the previous and
    /// next field are set to NULL).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        v: Intersection,
        boundary: bool,
        interference: bool,
        state_before: State,
        state_after: State,
        edge_before: State,
        edge_after: State,
    ) -> Self {
        AreaLimit {
            my_vertex: v,
            my_boundary: boundary,
            my_interference: interference,
            my_state_before: Cell::new(state_before),
            my_state_after: Cell::new(state_after),
            my_edge_before: Cell::new(edge_before),
            my_edge_after: Cell::new(edge_after),
            my_previous: RefCell::new(None),
            my_next: RefCell::new(None),
        }
    }

    /// OCCT StateBefore(St) — cxx L44-47.
    pub fn set_state_before(&self, st: State) {
        self.my_state_before.set(st);
    }

    /// OCCT StateAfter(St) — cxx L51-54.
    pub fn set_state_after(&self, st: State) {
        self.my_state_after.set(st);
    }

    /// OCCT EdgeBefore(St) — cxx L58-61.
    pub fn set_edge_before(&self, st: State) {
        self.my_edge_before.set(st);
    }

    /// OCCT EdgeAfter(St) — cxx L65-68.
    pub fn set_edge_after(&self, st: State) {
        self.my_edge_after.set(st);
    }

    /// OCCT Previous(P) — cxx L72-75.
    pub fn set_previous(&self, p: Option<SharedAreaLimit>) {
        *self.my_previous.borrow_mut() = p;
    }

    /// OCCT Next(N) — cxx L79-82.
    pub fn set_next(&self, n: Option<SharedAreaLimit>) {
        *self.my_next.borrow_mut() = n;
    }

    /// OCCT Vertex() — cxx L86-89.
    pub fn vertex(&self) -> &Intersection {
        &self.my_vertex
    }

    /// OCCT IsBoundary() — cxx L93-96.
    pub fn is_boundary(&self) -> bool {
        self.my_boundary
    }

    /// OCCT IsInterference() — cxx L100-103.
    pub fn is_interference(&self) -> bool {
        self.my_interference
    }

    /// OCCT StateBefore() — cxx L107-110.
    pub fn state_before(&self) -> State {
        self.my_state_before.get()
    }

    /// OCCT StateAfter() — cxx L114-117.
    pub fn state_after(&self) -> State {
        self.my_state_after.get()
    }

    /// OCCT EdgeBefore() — cxx L121-124.
    pub fn edge_before(&self) -> State {
        self.my_edge_before.get()
    }

    /// OCCT EdgeAfter() — cxx L128-131.
    pub fn edge_after(&self) -> State {
        self.my_edge_after.get()
    }

    /// OCCT Previous() — cxx L135-139.
    pub fn previous(&self) -> Option<SharedAreaLimit> {
        self.my_previous.borrow().clone()
    }

    /// OCCT Next() — cxx L142-145.
    pub fn next(&self) -> Option<SharedAreaLimit> {
        self.my_next.borrow().clone()
    }

    /// OCCT Clear() — cxx L149-153 (myPrevious.Nullify(); myNext.Nullify()).
    pub fn clear(&self) {
        *self.my_previous.borrow_mut() = None;
        *self.my_next.borrow_mut() = None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::topods::Orientation;

    /// OCCT ctor + the state-machine fields (hxx L35-41, cxx L25-40): the
    /// before/after states round-trip and the boundary/interference flags
    /// are stored verbatim.
    #[test]
    fn area_limit_state_machine_fields() {
        let v = Intersection::from_parts(
            Orientation::Forward,
            1,
            0,
            5,
            0.25,
            1.0e-4,
            State::On,
        );
        let lim = AreaLimit::new(
            v,
            true,
            false,
            State::Out,
            State::In,
            State::Unknown,
            State::Unknown,
        );
        assert_eq!(lim.vertex().parameter(), 0.25);
        assert_eq!(lim.vertex().index(), 5);
        assert!(lim.is_boundary());
        assert!(!lim.is_interference());
        assert_eq!(lim.state_before(), State::Out);
        assert_eq!(lim.state_after(), State::In);
        assert_eq!(lim.edge_before(), State::Unknown);
        assert_eq!(lim.edge_after(), State::Unknown);
        assert!(lim.previous().is_none());
        assert!(lim.next().is_none());
    }

    /// OCCT setters + Previous/Next linking + Clear (cxx L44-82, L149-153).
    #[test]
    fn area_limit_linking_and_clear() {
        let a = Arc::new(AreaLimit::new(
            Intersection::new(),
            true,
            false,
            State::Unknown,
            State::Unknown,
            State::Unknown,
            State::Unknown,
        ));
        let b = Arc::new(AreaLimit::new(
            Intersection::new(),
            false,
            true,
            State::In,
            State::In,
            State::In,
            State::In,
        ));
        a.set_next(Some(b.clone()));
        b.set_previous(Some(a.clone()));
        // Shared-handle identity: the returned handle is the same object.
        let n = a.next().expect("next");
        assert!(Arc::ptr_eq(&n, &b));
        let p = b.previous().expect("previous");
        assert!(Arc::ptr_eq(&p, &a));
        // State setters mutate through the shared handle.
        n.set_state_before(State::On);
        n.set_edge_after(State::Out);
        assert_eq!(b.state_before(), State::On);
        assert_eq!(b.edge_after(), State::Out);
        // Clear nullifies both links.
        b.clear();
        assert!(b.previous().is_none());
        assert!(b.next().is_none());
        // a's link to b is a strong handle; drop it to release.
        assert!(a.next().is_some());
    }
}
