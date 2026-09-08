//! OCCT MAT_Zone — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT/
//!         MAT_Zone.hxx L17-73, MAT_Zone.cxx L27-179

use std::sync::{Arc, RwLock};

use super::mat_arc::HandleMatArc;
use super::mat_basic_elt::HandleMatBasicElt;
use super::mat_node::HandleMatNode;
use super::MatSide;

/// OCCT `occ::handle<MAT_Zone>`
pub type HandleMatZone = Arc<RwLock<MatZone>>;

/// Definition of Zone of Proximity of a BasicElt :
/// ----------------------------------------------
/// A Zone of proximity is the set of the points which are
/// more near from the BasicElt than any other.
///
/// OCCT MAT_Zone.hxx L36-70
pub struct MatZone {
    frontier: Vec<HandleMatArc>,
    limited: bool,
}

impl MatZone {
    /// OCCT MAT_Zone.cxx L27-30 (default constructor)
    pub fn new() -> Self {
        MatZone {
            frontier: Vec::new(),
            limited: true,
        }
    }

    /// Compute the frontier of the Zone of proximity.
    ///
    /// OCCT MAT_Zone.cxx L34-37 (constructor overload)
    pub fn new_from_basic_elt(a_basic_elt: &HandleMatBasicElt) -> Self {
        let mut zone = MatZone::new();
        zone.perform(a_basic_elt);
        zone
    }

    /// Compute the frontier of the Zone of proximity.
    ///
    /// OCCT MAT_Zone.cxx L41-117
    pub fn perform(&mut self, a_basic_elt: &HandleMatBasicElt) {
        let mut next_node: Option<HandleMatNode>;
        let start_node: Option<HandleMatNode>;
        let mut current_arc: Option<HandleMatArc>;

        self.limited = true;
        self.frontier.clear();
        // ------------------------------------------------------------------------
        // Si le premier arc correspondant a la zone est Null => Sequence vide.
        // ------------------------------------------------------------------------
        if a_basic_elt.read().unwrap().end_arc().is_none() {
            return;
        }

        // ----------------------------
        // Angle rentrant => Zone Vide.
        // ----------------------------
        // if(aBasicElt->EndArc() == aBasicElt->StartArc()) return;

        // --------------------------------
        // Initialisation de la frontier.
        // --------------------------------
        current_arc = a_basic_elt.read().unwrap().end_arc();
        self.frontier.push(current_arc.clone().unwrap());

        // --------------------------------------------------------------------------
        // Determination du premier noeud qui permet de construire la zone en tournant
        // surla gauche.
        // --------------------------------------------------------------------------
        next_node = self.node_for_turn(current_arc.as_ref().unwrap(), a_basic_elt, MatSide::Left);
        start_node = current_arc
            .as_ref()
            .unwrap()
            .read()
            .unwrap()
            .the_other_node(next_node.as_ref().expect("MAT_Zone::Perform"));

        // -------------------------------------------------------------------------
        // Exploration du Graph toujours sur les arcs voisins a gauche jusqu'a
        // - retour sur la Figure .
        // - l acces a un noeud infini .
        // (Ces deux  cas correspondent a des noeuds pendants.)
        // - retour sur l arc de depart si le basicElt est ferme.
        // -------------------------------------------------------------------------

        while !next_node
            .as_ref()
            .expect("MAT_Zone::Perform")
            .read()
            .unwrap()
            .pending_node(next_node.as_ref().expect("MAT_Zone::Perform"))
            && !super::same_handle(&next_node, &start_node)
        {
            let next = current_arc
                .as_ref()
                .unwrap()
                .read()
                .unwrap()
                .neighbour(next_node.as_ref().expect("MAT_Zone::Perform"), MatSide::Left)
                .expect("MAT_Zone::Perform");
            current_arc = Some(next);
            self.frontier.push(current_arc.clone().unwrap());
            next_node = current_arc
                .as_ref()
                .unwrap()
                .read()
                .unwrap()
                .the_other_node(next_node.as_ref().expect("MAT_Zone::Perform"));
        }

        // -----------------------------------------------------------------------
        // Si NextNode est a l infini : exploration du graph a partir du StartArc
        //   sur <aBasicElt>.
        //   exploration sur les arcs voisins a droite.
        // Sinon => Fin.
        // -----------------------------------------------------------------------

        if next_node
            .as_ref()
            .expect("MAT_Zone::Perform")
            .read()
            .unwrap()
            .infinite()
        {
            self.limited = false;
            current_arc = a_basic_elt.read().unwrap().start_arc();
            self.frontier.push(current_arc.clone().unwrap());
            // --------------------------------------------------------------------------
            // Determination du premier noeud qui permet de construire la zone en
            // tournan surla droite.
            // --------------------------------------------------------------------------
            next_node = self.node_for_turn(current_arc.as_ref().unwrap(), a_basic_elt, MatSide::Right);

            // -----------------------------------------------------
            // Cette branche est aussi terminee par un noeud infini.
            // -----------------------------------------------------
            while !next_node
                .as_ref()
                .expect("MAT_Zone::Perform")
                .read()
                .unwrap()
                .infinite()
            {
                let next = current_arc
                    .as_ref()
                    .unwrap()
                    .read()
                    .unwrap()
                    .neighbour(next_node.as_ref().expect("MAT_Zone::Perform"), MatSide::Right)
                    .expect("MAT_Zone::Perform");
                current_arc = Some(next);
                self.frontier.push(current_arc.clone().unwrap());
                next_node = current_arc
                    .as_ref()
                    .unwrap()
                    .read()
                    .unwrap()
                    .the_other_node(next_node.as_ref().expect("MAT_Zone::Perform"));
            }
        }
    }

    /// Return the number Of Arcs On the frontier of <me>.
    ///
    /// OCCT MAT_Zone.cxx L121-124
    pub fn number_of_arcs(&self) -> i32 {
        self.frontier.len() as i32
    }

    /// Return the Arc number <Index> on the frontier.
    /// of <me>.
    ///
    /// OCCT MAT_Zone.cxx L128-131
    pub fn arc_on_frontier(&self, index: i32) -> HandleMatArc {
        self.frontier[(index - 1) as usize].clone()
    }

    /// Return TRUE if <me> is not empty .
    ///
    /// OCCT MAT_Zone.cxx L135-138
    pub fn no_empty_zone(&self) -> bool {
        !self.frontier.is_empty()
    }

    /// Return TRUE if <me> is Limited.
    ///
    /// OCCT MAT_Zone.cxx L142-145
    pub fn limited(&self) -> bool {
        self.limited
    }

    /// OCCT MAT_Zone.cxx L149-179 (private in OCCT)
    fn node_for_turn(
        &self,
        an_arc: &HandleMatArc,
        a_be: &HandleMatBasicElt,
        a_side: MatSide,
    ) -> Option<HandleMatNode> {
        let mut neighbour_arc: Option<HandleMatArc>;
        let mut node_sol: Option<HandleMatNode>;

        {
            let arc = an_arc.read().unwrap();
            node_sol = arc.first_node();
            neighbour_arc = match node_sol.as_ref() {
                Some(n) => arc.neighbour(n, a_side),
                None => None,
            };
            if neighbour_arc.is_none() {
                node_sol = arc.second_node();
                neighbour_arc = match node_sol.as_ref() {
                    Some(n) => arc.neighbour(n, a_side),
                    None => None,
                };
            }
        }
        if neighbour_arc.is_none() {
            return node_sol;
        }
        if super::same_handle(
            &neighbour_arc.as_ref().unwrap().read().unwrap().first_element(),
            &Some(a_be.clone()),
        ) {
            return node_sol;
        } else if super::same_handle(
            &neighbour_arc.as_ref().unwrap().read().unwrap().second_element(),
            &Some(a_be.clone()),
        ) {
            return node_sol;
        } else {
            return an_arc.read().unwrap().the_other_node(node_sol.as_ref().expect("MAT_Zone::NodeForTurn"));
        }
    }
}

impl Default for MatZone {
    fn default() -> Self {
        Self::new()
    }
}
