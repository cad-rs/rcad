//! OCCT MAT_Arc — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT/
//!         MAT_Arc.hxx L17-110, MAT_Arc.cxx L26-271

use std::sync::{Arc, RwLock};

use super::mat_basic_elt::HandleMatBasicElt;
use super::mat_node::HandleMatNode;
use super::MatSide;

/// OCCT `occ::handle<MAT_Arc>`
pub type HandleMatArc = Arc<RwLock<MatArc>>;

/// An Arc is associated to each Bisecting of the mat.
///
/// OCCT MAT_Arc.hxx L30-108. OCCT stores the four neighbour arcs as raw
/// `void*` to avoid handle cycles; here strong handles (see mod.rs handle
/// mapping).
pub struct MatArc {
    arc_index: i32,
    geom_index: i32,
    first_element: Option<HandleMatBasicElt>,
    second_element: Option<HandleMatBasicElt>,
    first_node: Option<HandleMatNode>,
    second_node: Option<HandleMatNode>,
    first_arc_left: Option<HandleMatArc>,   // OCCT: void* firstArcLeft
    first_arc_right: Option<HandleMatArc>,  // OCCT: void* firstArcRight
    second_arc_right: Option<HandleMatArc>, // OCCT: void* secondArcRight
    second_arc_left: Option<HandleMatArc>,  // OCCT: void* secondArcLeft
}

impl MatArc {
    /// OCCT MAT_Arc.cxx L26-39
    pub fn new(
        arc_index: i32,
        geom_index: i32,
        first_element: &HandleMatBasicElt,
        second_element: &HandleMatBasicElt,
    ) -> Self {
        MatArc {
            arc_index,
            geom_index,
            first_arc_left: None,
            first_arc_right: None,
            second_arc_right: None,
            second_arc_left: None,
            first_element: Some(first_element.clone()),
            second_element: Some(second_element.clone()),
            first_node: None,
            second_node: None,
        }
    }

    /// Returns the index of <me> in Graph.theArcs.
    ///
    /// OCCT MAT_Arc.cxx L43-46
    pub fn index(&self) -> i32 {
        self.arc_index
    }

    /// Returns the index associated of the geometric
    /// representation of <me>.
    ///
    /// OCCT MAT_Arc.cxx L50-53
    pub fn geom_index(&self) -> i32 {
        self.geom_index
    }

    /// Returns one of the BasicElt equidistant from <me>.
    ///
    /// OCCT MAT_Arc.cxx L57-60
    pub fn first_element(&self) -> Option<HandleMatBasicElt> {
        self.first_element.clone()
    }

    /// Returns the other BasicElt equidistant from <me>.
    ///
    /// OCCT MAT_Arc.cxx L64-67
    pub fn second_element(&self) -> Option<HandleMatBasicElt> {
        self.second_element.clone()
    }

    /// Returns one Node extremity of <me>.
    ///
    /// OCCT MAT_Arc.cxx L71-74
    pub fn first_node(&self) -> Option<HandleMatNode> {
        self.first_node.clone()
    }

    /// Returns the other Node extremity of <me>.
    ///
    /// OCCT MAT_Arc.cxx L78-81
    pub fn second_node(&self) -> Option<HandleMatNode> {
        self.second_node.clone()
    }

    /// An Arc has two Node, if <aNode> equals one
    /// Returns the other.
    ///
    /// if <aNode> is not oh <me>
    ///
    /// OCCT MAT_Arc.cxx L85-99
    pub fn the_other_node(&self, a_node: &HandleMatNode) -> Option<HandleMatNode> {
        if self
            .first_node
            .as_ref()
            .map(|n| Arc::ptr_eq(n, a_node))
            .unwrap_or(false)
        {
            return self.second_node.clone();
        } else if self
            .second_node
            .as_ref()
            .map(|n| Arc::ptr_eq(n, a_node))
            .unwrap_or(false)
        {
            return self.first_node.clone();
        } else {
            // throws Standard_DomainError("MAT_Arc::TheOtherNode")
            panic!("MAT_Arc::TheOtherNode");
        }
    }

    /// Returns True if there is an arc linked to
    /// the Node <aNode> located on the side <aSide> of <me>;
    /// if <aNode> is not on <me>
    ///
    /// OCCT MAT_Arc.cxx L103-132
    pub fn has_neighbour(&self, a_node: &HandleMatNode, a_side: MatSide) -> bool {
        if a_side == MatSide::Left {
            //    if (aNode == FirstNode())  return (!firstArcLeft  == NULL);
            if self
                .first_node
                .as_ref()
                .map(|n| Arc::ptr_eq(n, a_node))
                .unwrap_or(false)
            {
                return self.first_arc_left.is_some();
            }
            //    if (aNode == SecondNode()) return (!secondArcLeft == NULL);
            if self
                .second_node
                .as_ref()
                .map(|n| Arc::ptr_eq(n, a_node))
                .unwrap_or(false)
            {
                return self.second_arc_left.is_some();
            }
        } else {
            //    if (aNode == FirstNode())  return (!firstArcRight  == NULL);
            if self
                .first_node
                .as_ref()
                .map(|n| Arc::ptr_eq(n, a_node))
                .unwrap_or(false)
            {
                return self.first_arc_right.is_some();
            }
            //    if (aNode == SecondNode()) return (!secondArcRight == NULL);
            if self
                .second_node
                .as_ref()
                .map(|n| Arc::ptr_eq(n, a_node))
                .unwrap_or(false)
            {
                return self.second_arc_right.is_some();
            }
        }
        // throws Standard_DomainError("MAT_Arc::HasNeighbour")
        panic!("MAT_Arc::HasNeighbour");
    }

    /// Returns the first arc linked to the Node <aNode>
    /// located on the side <aSide> of <me>;
    /// if HasNeighbour() returns FALSE.
    ///
    /// OCCT MAT_Arc.cxx L136-163
    pub fn neighbour(&self, a_node: &HandleMatNode, a_side: MatSide) -> Option<HandleMatArc> {
        if a_side == MatSide::Left {
            if self
                .first_node
                .as_ref()
                .map(|n| Arc::ptr_eq(n, a_node))
                .unwrap_or(false)
            {
                return self.first_arc_left.clone();
            }
            if self
                .second_node
                .as_ref()
                .map(|n| Arc::ptr_eq(n, a_node))
                .unwrap_or(false)
            {
                return self.second_arc_left.clone();
            }
        } else {
            if self
                .first_node
                .as_ref()
                .map(|n| Arc::ptr_eq(n, a_node))
                .unwrap_or(false)
            {
                return self.first_arc_right.clone();
            }
            if self
                .second_node
                .as_ref()
                .map(|n| Arc::ptr_eq(n, a_node))
                .unwrap_or(false)
            {
                return self.second_arc_right.clone();
            }
        }
        // throws Standard_DomainError("MAT_Arc::Neighbour")
        panic!("MAT_Arc::Neighbour");
    }

    /// OCCT MAT_Arc.cxx L167-170
    pub fn set_index(&mut self, an_integer: i32) {
        self.arc_index = an_integer;
    }

    /// OCCT MAT_Arc.cxx L174-177
    pub fn set_geom_index(&mut self, an_integer: i32) {
        self.geom_index = an_integer;
    }

    /// OCCT MAT_Arc.cxx L181-184
    pub fn set_first_element(&mut self, a_basic_elt: &HandleMatBasicElt) {
        self.first_element = Some(a_basic_elt.clone());
    }

    /// OCCT MAT_Arc.cxx L188-191
    pub fn set_second_element(&mut self, a_basic_elt: &HandleMatBasicElt) {
        self.second_element = Some(a_basic_elt.clone());
    }

    /// OCCT MAT_Arc.cxx L195-198
    pub fn set_first_node(&mut self, a_node: &HandleMatNode) {
        self.first_node = Some(a_node.clone());
    }

    /// OCCT MAT_Arc.cxx L202-205
    pub fn set_second_node(&mut self, a_node: &HandleMatNode) {
        self.second_node = Some(a_node.clone());
    }

    /// OCCT MAT_Arc.cxx L209-219.
    /// OCCT signature takes a handle which may be null (set from
    /// MAT_Graph::FusionOfArcs with EmptyArc); hence Option.
    pub fn set_first_arc(&mut self, a_side: MatSide, an_arc: Option<&HandleMatArc>) {
        if a_side == MatSide::Left {
            self.first_arc_left = an_arc.cloned();
        } else {
            self.first_arc_right = an_arc.cloned();
        }
    }

    /// OCCT MAT_Arc.cxx L223-233 (see set_first_arc for the Option param)
    pub fn set_second_arc(&mut self, a_side: MatSide, an_arc: Option<&HandleMatArc>) {
        if a_side == MatSide::Left {
            self.second_arc_left = an_arc.cloned();
        } else {
            self.second_arc_right = an_arc.cloned();
        }
    }

    /// OCCT MAT_Arc.cxx L237-271
    pub fn set_neighbour(&mut self, a_side: MatSide, a_node: &HandleMatNode, an_arc: &HandleMatArc) {
        if a_side == MatSide::Left {
            if self
                .first_node
                .as_ref()
                .map(|n| Arc::ptr_eq(n, a_node))
                .unwrap_or(false)
            {
                self.first_arc_left = Some(an_arc.clone());
            } else if self
                .second_node
                .as_ref()
                .map(|n| Arc::ptr_eq(n, a_node))
                .unwrap_or(false)
            {
                self.second_arc_left = Some(an_arc.clone());
            } else {
                // throws Standard_DomainError("MAT_Arc::SetNeighbour")
                panic!("MAT_Arc::SetNeighbour");
            }
        } else {
            if self
                .first_node
                .as_ref()
                .map(|n| Arc::ptr_eq(n, a_node))
                .unwrap_or(false)
            {
                self.first_arc_right = Some(an_arc.clone());
            } else if self
                .second_node
                .as_ref()
                .map(|n| Arc::ptr_eq(n, a_node))
                .unwrap_or(false)
            {
                self.second_arc_right = Some(an_arc.clone());
            } else {
                // throws Standard_DomainError("MAT_Arc::SetNeighbour")
                panic!("MAT_Arc::SetNeighbour");
            }
        }
    }
}
