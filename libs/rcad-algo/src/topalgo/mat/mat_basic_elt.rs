//! OCCT MAT_BasicElt — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT/
//!         MAT_BasicElt.hxx L17-67, MAT_BasicElt.cxx L27-89

use std::sync::{Arc, RwLock};

use super::mat_arc::HandleMatArc;

/// OCCT `occ::handle<MAT_BasicElt>`
pub type HandleMatBasicElt = Arc<RwLock<MatBasicElt>>;

/// A BasicELt is associated to each elementary
/// constituent of the figure.
///
/// OCCT MAT_BasicElt.hxx L29-65. OCCT stores `startLeftArc`/`endLeftArc`
/// as raw `void*`; here strong handles (see mod.rs handle mapping).
pub struct MatBasicElt {
    start_left_arc: Option<HandleMatArc>, // OCCT: void* startLeftArc
    end_left_arc: Option<HandleMatArc>,   // OCCT: void* endLeftArc
    index: i32,
    geom_index: i32,
}

impl MatBasicElt {
    /// Constructor, <anInteger> is the <index> of <me>.
    ///
    /// OCCT MAT_BasicElt.cxx L27-33
    pub fn new(an_integer: i32) -> Self {
        MatBasicElt {
            start_left_arc: None,
            end_left_arc: None,
            index: an_integer,
            geom_index: 0,
        }
    }

    /// Return <startArcLeft> or <startArcRight> corresponding
    /// to <aSide>.
    ///
    /// OCCT MAT_BasicElt.cxx L37-40
    pub fn start_arc(&self) -> Option<HandleMatArc> {
        self.start_left_arc.clone()
    }

    /// Return <endArcLeft> or <endArcRight> corresponding
    /// to <aSide>.
    ///
    /// OCCT MAT_BasicElt.cxx L44-47
    pub fn end_arc(&self) -> Option<HandleMatArc> {
        self.end_left_arc.clone()
    }

    /// Return the <index> of <me> in Graph.TheBasicElts.
    ///
    /// OCCT MAT_BasicElt.cxx L51-54
    pub fn index(&self) -> i32 {
        self.index
    }

    /// Return the <GeomIndex> of <me>.
    ///
    /// OCCT MAT_BasicElt.cxx L58-61
    pub fn geom_index(&self) -> i32 {
        self.geom_index
    }

    /// OCCT MAT_BasicElt.cxx L65-68.
    /// OCCT signature takes a handle which may be null (set from
    /// MAT_Graph::FusionOfBasicElts with Elt2->EndArc()); hence Option.
    pub fn set_start_arc(&mut self, an_arc: Option<&HandleMatArc>) {
        self.start_left_arc = an_arc.cloned();
    }

    /// OCCT MAT_BasicElt.cxx L72-75 (see set_start_arc for the Option param)
    pub fn set_end_arc(&mut self, an_arc: Option<&HandleMatArc>) {
        self.end_left_arc = an_arc.cloned();
    }

    /// OCCT MAT_BasicElt.cxx L79-82
    pub fn set_index(&mut self, an_integer: i32) {
        self.index = an_integer;
    }

    /// OCCT MAT_BasicElt.cxx L86-89
    pub fn set_geom_index(&mut self, an_integer: i32) {
        self.geom_index = an_integer;
    }
}
