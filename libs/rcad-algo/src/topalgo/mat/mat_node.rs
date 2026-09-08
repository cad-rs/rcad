//! OCCT MAT_Node — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT/
//!         MAT_Node.hxx L17-80, MAT_Node.cxx L25-151

use std::sync::{Arc, RwLock};

use rcad_kernel::core::precision::INFINITE_VALUE; // OCCT Precision::Infinite()

use super::mat_arc::HandleMatArc;
use super::mat_basic_elt::HandleMatBasicElt;
use super::MatSide;

/// OCCT `occ::handle<MAT_Node>`
pub type HandleMatNode = Arc<RwLock<MatNode>>;

/// Node of Graph.
///
/// OCCT MAT_Node.hxx L32-78. OCCT stores `aLinkedArc` as a raw `void*`;
/// here a strong handle (see mod.rs handle mapping).
pub struct MatNode {
    node_index: i32,
    geom_index: i32,
    a_linked_arc: Option<HandleMatArc>, // OCCT: void* aLinkedArc
    distance: f64,
}

impl MatNode {
    /// OCCT MAT_Node.cxx L25-33
    pub fn new(geom_index: i32, linked_arc: &HandleMatArc, distance: f64) -> Self {
        MatNode {
            node_index: 0,
            geom_index,
            a_linked_arc: Some(linked_arc.clone()),
            distance,
        }
    }

    /// Returns the index associated of the geometric
    /// representation of <me>.
    ///
    /// OCCT MAT_Node.cxx L37-40
    pub fn geom_index(&self) -> i32 {
        self.geom_index
    }

    /// Returns the index associated of the node
    ///
    /// OCCT MAT_Node.cxx L44-47
    pub fn index(&self) -> i32 {
        self.node_index
    }

    /// Returns in <S> the Arcs linked to <me>.
    ///
    /// OCCT MAT_Node.cxx L51-68. OCCT L54 creates `occ::handle<MAT_Node>
    /// Me = this`; Rust has no implicit handle to self, so the caller
    /// passes the outer node handle as <me> (architecture difference).
    pub fn linked_arcs(&self, me: &HandleMatNode, s: &mut Vec<HandleMatArc>) {
        s.clear();
        let la = self.a_linked_arc.clone().expect("MAT_Node::LinkedArcs");

        s.push(la.clone());

        if la.read().unwrap().has_neighbour(me, MatSide::Left) {
            let mut ca = la.read().unwrap().neighbour(me, MatSide::Left).expect("MAT_Node::LinkedArcs");
            while !Arc::ptr_eq(&ca, &la) {
                s.push(ca.clone());
                let next = ca.read().unwrap().neighbour(me, MatSide::Left).expect("MAT_Node::LinkedArcs");
                ca = next;
            }
        }
    }

    /// Returns in <S> the BasicElts equidistant
    /// to <me>.
    ///
    /// OCCT MAT_Node.cxx L72-107
    pub fn near_elts(&self, me: &HandleMatNode, s: &mut Vec<HandleMatBasicElt>) {
        s.clear();

        let la = self.a_linked_arc.clone().expect("MAT_Node::NearElts");

        s.push(la.read().unwrap().first_element().expect("MAT_Node::NearElts"));
        s.push(la.read().unwrap().second_element().expect("MAT_Node::NearElts"));

        if la.read().unwrap().has_neighbour(me, MatSide::Left) {
            let mut ca = la.read().unwrap().neighbour(me, MatSide::Left).expect("MAT_Node::NearElts");
            let mut pair = false;

            //---------------------------------------------------------
            // Recuperation des deux elements separes pour un arc sur
            // deux.
            //---------------------------------------------------------

            while !Arc::ptr_eq(&ca, &la) {
                if pair {
                    let ca_ref = ca.read().unwrap();
                    s.push(ca_ref.first_element().expect("MAT_Node::NearElts"));
                    s.push(ca_ref.second_element().expect("MAT_Node::NearElts"));
                } else {
                    pair = true;
                }
                let next = ca.read().unwrap().neighbour(me, MatSide::Left).expect("MAT_Node::NearElts");
                ca = next;
            }
        }
    }

    /// OCCT MAT_Node.cxx L111-114
    pub fn distance(&self) -> f64 {
        self.distance
    }

    /// Returns True if <me> is a pending Node.
    /// (ie : the number of Arc Linked = 1)
    ///
    /// OCCT MAT_Node.cxx L118-122. OCCT L120 creates `Me = this`; the
    /// caller passes the outer handle (see linked_arcs).
    pub fn pending_node(&self, me: &HandleMatNode) -> bool {
        let la = self.a_linked_arc.clone().expect("MAT_Node::PendingNode");
        let pending = !la.read().unwrap().has_neighbour(me, MatSide::Left);
        pending
    }

    /// Returns True if <me> belongs to the figure.
    ///
    /// OCCT MAT_Node.cxx L126-129
    pub fn on_basic_elt(&self) -> bool {
        self.distance() == 0.0
    }

    /// Returns True if the distance of <me> is Infinite
    ///
    /// OCCT MAT_Node.cxx L133-136
    pub fn infinite(&self) -> bool {
        self.distance() == INFINITE_VALUE
    }

    /// Set the index associated of the node
    ///
    /// OCCT MAT_Node.cxx L147-150
    pub fn set_index(&mut self, an_index: i32) {
        self.node_index = an_index;
    }

    /// OCCT MAT_Node.cxx L140-143
    pub fn set_linked_arc(&mut self, linked_arc: &HandleMatArc) {
        self.a_linked_arc = Some(linked_arc.clone());
    }
}
