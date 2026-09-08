//! OCCT MAT_Graph — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT/
//!         MAT_Graph.hxx L17-131, MAT_Graph.cxx L35-581
//!         (including the file-static MakeArc, translated as a module
//!         private free function).

use std::collections::HashMap;
use std::sync::{Arc, RwLock};

use rcad_kernel::core::precision::INFINITE_VALUE; // OCCT Precision::Infinite()

use super::mat_arc::{HandleMatArc, MatArc};
use super::mat_basic_elt::{HandleMatBasicElt, MatBasicElt};
use super::mat_bisector::HandleMatBisector;
use super::mat_list_of_bisector_0::HandleMatListOfBisector;
use super::mat_node::{HandleMatNode, MatNode};
use super::mat_zone::MatZone;
use super::MatSide;

/// OCCT `occ::handle<MAT_Graph>`
pub type HandleMatGraph = Arc<RwLock<MatGraph>>;

/// The Class Graph permits the exploration of the
/// Bisector Locus.
///
/// OCCT MAT_Graph.hxx L36-129. OCCT NCollection_DataMap<int, handle> maps
/// to HashMap<i32, Handle> (architecture layer).
pub struct MatGraph {
    the_arcs: HashMap<i32, HandleMatArc>,
    the_basic_elts: HashMap<i32, HandleMatBasicElt>,
    the_nodes: HashMap<i32, HandleMatNode>,
    number_of_arcs: i32,
    number_of_nodes: i32,
    number_of_basic_elts: i32,
    number_of_infinite_nodes: i32,
}

impl MatGraph {
    /// Empty constructor.
    ///
    /// OCCT MAT_Graph.cxx L44-50
    pub fn new() -> Self {
        MatGraph {
            the_arcs: HashMap::new(),
            the_basic_elts: HashMap::new(),
            the_nodes: HashMap::new(),
            number_of_arcs: 0,
            number_of_nodes: 0,
            number_of_basic_elts: 0,
            number_of_infinite_nodes: 0,
        }
    }

    /// function : Perform
    /// purpose  : Creation du graphe contenant le resultat.
    ///
    /// OCCT MAT_Graph.cxx L56-170
    pub fn perform(
        &mut self,
        semi_infinite: bool,
        the_roots: &HandleMatListOfBisector,
        nb_basic_elts: i32,
        nb_arcs: i32,
    ) {
        let mut nb_roots: i32;
        let mut first_arc: Option<HandleMatArc> = None;
        let mut current_arc: Option<HandleMatArc> = None;
        let mut extremite: Option<HandleMatNode> = None;
        let mut ind_tab_arcs: i32 = 1;
        let mut ind_tab_nodes: i32;
        let dist_ext: f64;
        let ind_ext: i32;
        let mut previous_arc: Option<HandleMatArc> = current_arc.clone();

        //------------------------
        // Construction du graphe.
        //------------------------

        if semi_infinite {
            nb_roots = the_roots.read().unwrap().number();
            self.number_of_infinite_nodes = nb_roots;
        } else {
            nb_roots = 1;
            self.number_of_infinite_nodes = 0;
        }

        self.number_of_arcs = nb_arcs;
        self.number_of_basic_elts = nb_basic_elts;
        self.number_of_nodes = nb_roots + nb_arcs;
        ind_tab_nodes = self.number_of_nodes;

        //---------------------------
        //... Creation des BasicElts.
        //---------------------------
        for i in 1..=nb_basic_elts {
            self.the_basic_elts
                .insert(i, Arc::new(RwLock::new(MatBasicElt::new(i))));
            self.the_basic_elts
                .get(&i)
                .expect("MAT_Graph::Perform")
                .write()
                .unwrap()
                .set_geom_index(i);
        }

        //--------------------------------------------------------------------
        // ... Creation des ARCS et des NODES.
        //     Construction des arbres d arcs a partir des <Bisector> racines.
        //--------------------------------------------------------------------

        if semi_infinite {
            // Plusieurs points d entree a l infini.
            //--------------------------------------
            the_roots.write().unwrap().first();

            while the_roots.read().unwrap().more() {
                current_arc = Some(make_arc(
                    the_roots
                        .read()
                        .unwrap()
                        .current()
                        .as_ref()
                        .expect("MAT_Graph::Perform"),
                    &self.the_basic_elts,
                    &mut self.the_arcs,
                    &mut ind_tab_arcs,
                ));
                extremite = Some(Arc::new(RwLock::new(MatNode::new(
                    0,
                    current_arc.as_ref().expect("MAT_Graph::Perform"),
                    INFINITE_VALUE,
                ))));
                extremite
                    .as_ref()
                    .expect("MAT_Graph::Perform")
                    .write()
                    .unwrap()
                    .set_index(ind_tab_nodes);
                current_arc
                    .as_ref()
                    .expect("MAT_Graph::Perform")
                    .write()
                    .unwrap()
                    .set_second_node(extremite.as_ref().expect("MAT_Graph::Perform"));
                self.the_nodes
                    .insert(ind_tab_nodes, extremite.clone().expect("MAT_Graph::Perform"));
                the_roots.write().unwrap().next();
                ind_tab_nodes -= 1;
            }
        } else {
            // -----------------------------------------------
            // Un seul point d entree .
            // Creation d un premier ARC et du NODE racine.
            // -----------------------------------------------
            nb_roots = 1;
            the_roots.write().unwrap().first();
            current_arc = Some(make_arc(
                the_roots
                    .read()
                    .unwrap()
                    .current()
                    .as_ref()
                    .expect("MAT_Graph::Perform"),
                &self.the_basic_elts,
                &mut self.the_arcs,
                &mut ind_tab_arcs,
            ));
            dist_ext = the_roots
                .read()
                .unwrap()
                .current()
                .as_ref()
                .expect("MAT_Graph::Perform")
                .read()
                .unwrap()
                .first_edge()
                .expect("MAT_Graph::Perform")
                .read()
                .unwrap()
                .distance();
            ind_ext = the_roots
                .read()
                .unwrap()
                .current()
                .as_ref()
                .expect("MAT_Graph::Perform")
                .read()
                .unwrap()
                .end_point();

            extremite = Some(Arc::new(RwLock::new(MatNode::new(
                ind_ext,
                current_arc.as_ref().expect("MAT_Graph::Perform"),
                dist_ext,
            ))));
            extremite
                .as_ref()
                .expect("MAT_Graph::Perform")
                .write()
                .unwrap()
                .set_index(ind_tab_nodes);
            current_arc
                .as_ref()
                .expect("MAT_Graph::Perform")
                .write()
                .unwrap()
                .set_second_node(extremite.as_ref().expect("MAT_Graph::Perform"));
            self.the_nodes
                .insert(ind_tab_nodes, extremite.clone().expect("MAT_Graph::Perform"));
            ind_tab_nodes -= 1;

            // -----------------------------------------------------------
            // ...Creation des ARCs issues de la racine.
            //    Codage des voisinages sur ces arcs et mise a jour de la
            //    sequence des arcs issue du Node racine.
            // -----------------------------------------------------------
            first_arc = current_arc.clone();
            previous_arc = first_arc.clone();
            the_roots.write().unwrap().next();

            while the_roots.read().unwrap().more() {
                current_arc = Some(make_arc(
                    the_roots
                        .read()
                        .unwrap()
                        .current()
                        .as_ref()
                        .expect("MAT_Graph::Perform"),
                    &self.the_basic_elts,
                    &mut self.the_arcs,
                    &mut ind_tab_arcs,
                ));
                current_arc
                    .as_ref()
                    .expect("MAT_Graph::Perform")
                    .write()
                    .unwrap()
                    .set_second_node(extremite.as_ref().expect("MAT_Graph::Perform"));
                current_arc
                    .as_ref()
                    .expect("MAT_Graph::Perform")
                    .write()
                    .unwrap()
                    .set_neighbour(
                        MatSide::Left,
                        extremite.as_ref().expect("MAT_Graph::Perform"),
                        previous_arc.as_ref().expect("MAT_Graph::Perform"),
                    );
                previous_arc
                    .as_ref()
                    .expect("MAT_Graph::Perform")
                    .write()
                    .unwrap()
                    .set_neighbour(
                        MatSide::Right,
                        extremite.as_ref().expect("MAT_Graph::Perform"),
                        current_arc.as_ref().expect("MAT_Graph::Perform"),
                    );

                previous_arc = current_arc.clone();
                the_roots.write().unwrap().next();
            }
            first_arc
                .as_ref()
                .expect("MAT_Graph::Perform")
                .write()
                .unwrap()
                .set_neighbour(
                    MatSide::Left,
                    extremite.as_ref().expect("MAT_Graph::Perform"),
                    current_arc.as_ref().expect("MAT_Graph::Perform"),
                );
            current_arc
                .as_ref()
                .expect("MAT_Graph::Perform")
                .write()
                .unwrap()
                .set_neighbour(
                    MatSide::Right,
                    extremite.as_ref().expect("MAT_Graph::Perform"),
                    first_arc.as_ref().expect("MAT_Graph::Perform"),
                );
        }

        // ----------------------------------------------------
        // Les sequence des Arcs des Nodes racines sont a jour.
        // Mise a jour des sequences des autres Nodes.
        // ----------------------------------------------------
        self.update_nodes(&mut ind_tab_nodes);
    }

    /// Return the Arc of index <Index> in <theArcs>.
    ///
    /// OCCT MAT_Graph.cxx L174-177
    pub fn arc(&self, index: i32) -> HandleMatArc {
        self.the_arcs
            .get(&index)
            .expect("MAT_Graph::Arc")
            .clone()
    }

    /// Return the BasicElt of index <Index> in <theBasicElts>.
    ///
    /// OCCT MAT_Graph.cxx L181-184
    pub fn basic_elt(&self, index: i32) -> HandleMatBasicElt {
        self.the_basic_elts
            .get(&index)
            .expect("MAT_Graph::BasicElt")
            .clone()
    }

    /// Return the Node of index <Index> in <theNodes>.
    ///
    /// OCCT MAT_Graph.cxx L188-191
    pub fn node(&self, index: i32) -> HandleMatNode {
        self.the_nodes
            .get(&index)
            .expect("MAT_Graph::Node")
            .clone()
    }

    /// Return the number of arcs of <me>.
    ///
    /// OCCT MAT_Graph.cxx L195-198
    pub fn number_of_arcs(&self) -> i32 {
        self.number_of_arcs
    }

    /// Return the number of nodes of <me>.
    ///
    /// OCCT MAT_Graph.cxx L202-205
    pub fn number_of_nodes(&self) -> i32 {
        self.number_of_nodes
    }

    /// Return the number of basic elements of <me>.
    ///
    /// OCCT MAT_Graph.cxx L216-219
    pub fn number_of_basic_elts(&self) -> i32 {
        self.number_of_basic_elts
    }

    /// Return the number of infinites nodes of <me>.
    ///
    /// OCCT MAT_Graph.cxx L209-212
    pub fn number_of_infinite_nodes(&self) -> i32 {
        self.number_of_infinite_nodes
    }

    /// Merge two BasicElts. The End of the BasicElt Elt1
    /// of IndexElt1 becomes The End of the BasicElt Elt2
    /// of IndexElt2. Elt2 is replaced in the arcs by
    /// Elt1, Elt2 is eliminated.
    ///
    /// <MergeArc1> is True if the fusion of the BasicElts =>
    /// a fusion of two Arcs which separated the same elements.
    /// In this case <GeomIndexArc1> and <GeomIndexArc2> are the
    /// Geometric Index of this arcs.
    ///
    /// If the BasicElt corresponds to a close line,
    /// the StartArc and the EndArc of Elt1 can separate the same
    /// elements.
    /// In this case there is a fusion of this arcs, <MergeArc2>
    /// is true and <GeomIndexArc3> and <GeomIndexArc4> are the
    /// Geometric Index of this arcs.
    ///
    /// OCCT MAT_Graph.cxx L223-321
    pub fn fusion_of_basic_elts(
        &mut self,
        index_elt1: i32,
        index_elt2: i32,
        merge_arc1: &mut bool,
        i_geom_arc1: &mut i32,
        i_geom_arc2: &mut i32,
        merge_arc2: &mut bool,
        i_geom_arc3: &mut i32,
        i_geom_arc4: &mut i32,
    ) {
        let elt1 = self
            .the_basic_elts
            .get(&index_elt1)
            .expect("MAT_Graph::FusionOfBasicElts")
            .clone();
        let elt2 = self
            .the_basic_elts
            .get(&index_elt2)
            .expect("MAT_Graph::FusionOfBasicElts")
            .clone();

        if super::same_handle(&Some(elt1.clone()), &Some(elt2.clone())) {
            return;
        }

        let zone2 = Arc::new(RwLock::new(MatZone::new_from_basic_elt(&elt2)));

        //--------------------------------------------------------------------
        // Les arcs de la zone de Elt2 ne separent plus Elt2 et qq chose mais
        // Elt1 et qq chose.
        //--------------------------------------------------------------------
        for i in 1..=zone2.read().unwrap().number_of_arcs() {
            let arc_on_frontier = zone2.read().unwrap().arc_on_frontier(i);
            if super::same_handle(
                &arc_on_frontier.read().unwrap().first_element(),
                &Some(elt2.clone()),
            ) {
                let idx = arc_on_frontier.read().unwrap().index();
                self.the_arcs
                    .get(&idx)
                    .expect("MAT_Graph::FusionOfBasicElts")
                    .write()
                    .unwrap()
                    .set_first_element(&elt1);
            } else {
                let idx = arc_on_frontier.read().unwrap().index();
                self.the_arcs
                    .get(&idx)
                    .expect("MAT_Graph::FusionOfBasicElts")
                    .write()
                    .unwrap()
                    .set_second_element(&elt1);
            }
        }

        //-------------------------------------------------------------------
        // le EndArc de Elt1 et le StartArc de Elt2 peuvent separes les memes
        // elements de base => Fusion des deux arcs et mise a jour des noeuds.
        //-------------------------------------------------------------------
        let mut ea1 = elt1.read().unwrap().end_arc();
        let sa2 = elt2.read().unwrap().start_arc();

        let mut e1 = ea1
            .as_ref()
            .expect("MAT_Graph::FusionOfBasicElts")
            .read()
            .unwrap()
            .first_element();
        let mut e2 = ea1
            .as_ref()
            .expect("MAT_Graph::FusionOfBasicElts")
            .read()
            .unwrap()
            .second_element();
        let e3 = sa2
            .as_ref()
            .expect("MAT_Graph::FusionOfBasicElts")
            .read()
            .unwrap()
            .first_element();
        let e4 = sa2
            .as_ref()
            .expect("MAT_Graph::FusionOfBasicElts")
            .read()
            .unwrap()
            .second_element();
        *merge_arc1 = false;

        if (super::same_handle(&e1, &e3) || super::same_handle(&e1, &e4))
            && (super::same_handle(&e2, &e3) || super::same_handle(&e2, &e4))
        {
            let ea1_index = ea1
                .as_ref()
                .expect("MAT_Graph::FusionOfBasicElts")
                .read()
                .unwrap()
                .index();
            let sa2_index = sa2
                .as_ref()
                .expect("MAT_Graph::FusionOfBasicElts")
                .read()
                .unwrap()
                .index();
            let arc1 = self
                .the_arcs
                .get(&ea1_index)
                .expect("MAT_Graph::FusionOfBasicElts")
                .clone();
            let arc2 = self
                .the_arcs
                .get(&sa2_index)
                .expect("MAT_Graph::FusionOfBasicElts")
                .clone();
            self.fusion_of_arcs(&arc1, &arc2);
            *merge_arc1 = true;
            *i_geom_arc1 = ea1
                .as_ref()
                .expect("MAT_Graph::FusionOfBasicElts")
                .read()
                .unwrap()
                .geom_index();
            *i_geom_arc2 = sa2
                .as_ref()
                .expect("MAT_Graph::FusionOfBasicElts")
                .read()
                .unwrap()
                .geom_index();
        }

        //-------------------------------------------------
        // La fin de Elt1 devient la fin de Elt2.
        //-------------------------------------------------
        elt1.write().unwrap().set_end_arc(elt2.read().unwrap().end_arc().as_ref());

        //-------------------------------------------------------------------
        // le EndArc de Elt1 et le StartArc de Elt1 peuvent separer les memes
        // elements de base.
        // si les noeuds des arcs ne sont pas sur le contour
        //    => fusion des arcs.(contour ferme compose d un seul BasicElt)
        // sinon rien            (contour ferme compose de deux BasicElts)
        //-------------------------------------------------------------------
        let sa1 = elt1.read().unwrap().start_arc();
        ea1 = elt1.read().unwrap().end_arc();

        if !super::same_handle(&ea1, &sa1) {
            e1 = ea1
                .as_ref()
                .expect("MAT_Graph::FusionOfBasicElts")
                .read()
                .unwrap()
                .first_element();
            e2 = ea1
                .as_ref()
                .expect("MAT_Graph::FusionOfBasicElts")
                .read()
                .unwrap()
                .second_element();
            let e3 = sa1
                .as_ref()
                .expect("MAT_Graph::FusionOfBasicElts")
                .read()
                .unwrap()
                .first_element();
            let e4 = sa1
                .as_ref()
                .expect("MAT_Graph::FusionOfBasicElts")
                .read()
                .unwrap()
                .second_element();

            let on_fig = {
                let ea1_ref = ea1.as_ref().expect("MAT_Graph::FusionOfBasicElts");
                let ea1_node1 = ea1_ref.read().unwrap().first_node();
                let ea1_node2 = ea1_ref.read().unwrap().second_node();
                let sa1_ref = sa1.as_ref().expect("MAT_Graph::FusionOfBasicElts");
                let sa1_node1 = sa1_ref.read().unwrap().first_node();
                let sa1_node2 = sa1_ref.read().unwrap().second_node();
                ea1_node1
                    .expect("MAT_Graph::FusionOfBasicElts")
                    .read()
                    .unwrap()
                    .on_basic_elt()
                    || ea1_node2
                        .expect("MAT_Graph::FusionOfBasicElts")
                        .read()
                        .unwrap()
                        .on_basic_elt()
                    || sa1_node1
                        .expect("MAT_Graph::FusionOfBasicElts")
                        .read()
                        .unwrap()
                        .on_basic_elt()
                    || sa1_node2
                        .expect("MAT_Graph::FusionOfBasicElts")
                        .read()
                        .unwrap()
                        .on_basic_elt()
            };

            *merge_arc2 = false;

            if (super::same_handle(&e1, &e3) || super::same_handle(&e1, &e4))
                && (super::same_handle(&e2, &e3) || super::same_handle(&e2, &e4))
                && !on_fig
            {
                let ea1_index = ea1
                    .as_ref()
                    .expect("MAT_Graph::FusionOfBasicElts")
                    .read()
                    .unwrap()
                    .index();
                let sa1_index = sa1
                    .as_ref()
                    .expect("MAT_Graph::FusionOfBasicElts")
                    .read()
                    .unwrap()
                    .index();
                let arc1 = self
                    .the_arcs
                    .get(&ea1_index)
                    .expect("MAT_Graph::FusionOfBasicElts")
                    .clone();
                let arc2 = self
                    .the_arcs
                    .get(&sa1_index)
                    .expect("MAT_Graph::FusionOfBasicElts")
                    .clone();
                self.fusion_of_arcs(&arc1, &arc2);
                *merge_arc2 = true;
                *i_geom_arc3 = ea1
                    .as_ref()
                    .expect("MAT_Graph::FusionOfBasicElts")
                    .read()
                    .unwrap()
                    .geom_index();
                *i_geom_arc4 = sa1
                    .as_ref()
                    .expect("MAT_Graph::FusionOfBasicElts")
                    .read()
                    .unwrap()
                    .geom_index();
            }
        }

        //----------------------------------------------------
        // un element de base a ete elimine.
        //----------------------------------------------------
        let elt2_index = elt2.read().unwrap().index();
        self.the_basic_elts.remove(&elt2_index);
        self.number_of_basic_elts -= 1;
    }

    /// function : FusionOfArcs
    /// purpose  : Fusion de deux arcs separant les memes elements.
    ///            l <Arc1> ira du Second noeud de <Arc2> au second Noeud de <Arc1>.
    ///
    /// OCCT MAT_Graph.cxx L328-403 (private in OCCT)
    fn fusion_of_arcs(&mut self, arc1: &HandleMatArc, arc2: &HandleMatArc) {
        let old_node1 = arc1.read().unwrap().first_node();
        let old_node2 = arc2.read().unwrap().first_node();

        let arc2_second_node = arc2.read().unwrap().second_node();
        arc1.write().unwrap().set_first_node(
            arc2_second_node
                .as_ref()
                .expect("MAT_Graph::FusionOfArcs"),
        );

        //--------------------------------------------------------------------
        // Mise a jour des voisinages autour du nouveau premier noeud de Arc1.
        //--------------------------------------------------------------------
        if !arc2_second_node
            .as_ref()
            .expect("MAT_Graph::FusionOfArcs")
            .read()
            .unwrap()
            .infinite()
        {
            let l_neighbour = arc2
                .read()
                .unwrap()
                .neighbour(
                    arc2_second_node.as_ref().expect("MAT_Graph::FusionOfArcs"),
                    MatSide::Left,
                )
                .expect("MAT_Graph::FusionOfArcs");
            let r_neighbour = arc2
                .read()
                .unwrap()
                .neighbour(
                    arc2_second_node.as_ref().expect("MAT_Graph::FusionOfArcs"),
                    MatSide::Right,
                )
                .expect("MAT_Graph::FusionOfArcs");

            arc1.write().unwrap().set_first_arc(MatSide::Left, Some(&l_neighbour));
            arc1.write()
                .unwrap()
                .set_first_arc(MatSide::Right, Some(&r_neighbour));
            let l_neighbour_index = l_neighbour.read().unwrap().index();
            self.the_arcs
                .get(&l_neighbour_index)
                .expect("MAT_Graph::FusionOfArcs")
                .write()
                .unwrap()
                .set_neighbour(
                    MatSide::Right,
                    arc2_second_node.as_ref().expect("MAT_Graph::FusionOfArcs"),
                    arc1,
                );
            let r_neighbour_index = r_neighbour.read().unwrap().index();
            self.the_arcs
                .get(&r_neighbour_index)
                .expect("MAT_Graph::FusionOfArcs")
                .write()
                .unwrap()
                .set_neighbour(
                    MatSide::Left,
                    arc2_second_node.as_ref().expect("MAT_Graph::FusionOfArcs"),
                    arc1,
                );
        } else {
            let empty_arc: Option<HandleMatArc> = None;
            arc1.write().unwrap().set_first_arc(MatSide::Left, empty_arc.as_ref());
            arc1.write()
                .unwrap()
                .set_first_arc(MatSide::Right, empty_arc.as_ref());
        }

        //-------------------------------------------------------------------
        // Mise a jour du premier noeud Arc1.
        //-----------------------------------------------------------------
        let arc1_first_node = arc1
            .read()
            .unwrap()
            .first_node()
            .expect("MAT_Graph::FusionOfArcs");
        arc1_first_node.write().unwrap().set_linked_arc(arc1);

        //------------------------------------
        // Elimination de Arc2 et des OldNode
        //------------------------------------
        let old_node1_index = old_node1
            .as_ref()
            .expect("MAT_Graph::FusionOfArcs")
            .read()
            .unwrap()
            .index();
        if self.the_nodes.contains_key(&old_node1_index) {
            self.the_nodes.remove(&old_node1_index);
            self.number_of_nodes -= 1;
        }
        let old_node2_index = old_node2
            .as_ref()
            .expect("MAT_Graph::FusionOfArcs")
            .read()
            .unwrap()
            .index();
        if self.the_nodes.contains_key(&old_node2_index) {
            self.the_nodes.remove(&old_node2_index);
            self.number_of_nodes -= 1;
        }

        // Note: the Arc2 is actually a reference to a handle contained in theArcs map;
        // it is necessary to create copy of that handle and use only it to access
        // that object, since the handle contained in the map is destroyed by UnBind()
        let an_arc2 = arc2.clone();
        let arc2_index = arc2.read().unwrap().index();
        self.the_arcs.remove(&arc2_index);
        self.number_of_arcs -= 1;

        for i in 1..=2 {
            let be: HandleMatBasicElt;
            if i == 1 {
                let fe = an_arc2
                    .read()
                    .unwrap()
                    .first_element()
                    .expect("MAT_Graph::FusionOfArcs");
                let fe_index = fe.read().unwrap().index();
                be = self
                    .the_basic_elts
                    .get(&fe_index)
                    .expect("MAT_Graph::FusionOfArcs")
                    .clone();
            } else {
                let se = an_arc2
                    .read()
                    .unwrap()
                    .second_element()
                    .expect("MAT_Graph::FusionOfArcs");
                let se_index = se.read().unwrap().index();
                be = self
                    .the_basic_elts
                    .get(&se_index)
                    .expect("MAT_Graph::FusionOfArcs")
                    .clone();
            }

            if super::same_handle(&be.read().unwrap().start_arc(), &Some(an_arc2.clone())) {
                be.write().unwrap().set_start_arc(Some(arc1));
            }
            if super::same_handle(&be.read().unwrap().end_arc(), &Some(an_arc2.clone())) {
                be.write().unwrap().set_end_arc(Some(arc1));
            }
        }
    }

    /// function : CompactArcs
    /// purpose  : Decalage des Arcs pour boucher les trous.
    ///
    /// OCCT MAT_Graph.cxx L409-433
    pub fn compact_arcs(&mut self) {
        let mut i_find: i32 = 0;
        let mut i: i32 = 1;
        let mut ya_decalage = false;

        while i_find < self.number_of_arcs {
            if !self.the_arcs.contains_key(&i) {
                ya_decalage = true;
            } else {
                i_find += 1;
                if ya_decalage {
                    let an_arc = self
                        .the_arcs
                        .get(&i)
                        .expect("MAT_Graph::CompactArcs")
                        .clone();
                    an_arc.write().unwrap().set_index(i_find);
                    self.the_arcs.insert(i_find, an_arc);
                    self.the_arcs.remove(&i);
                }
            }
            i += 1;
        }
    }

    /// function : CompactNodes
    /// purpose  : Decalage des Nodes pour boucher les trous.
    ///
    /// OCCT MAT_Graph.cxx L439-463
    pub fn compact_nodes(&mut self) {
        let mut i_find: i32 = 0;
        let mut i: i32 = 1;
        let mut ya_decalage = false;

        while i_find < self.number_of_nodes {
            if !self.the_nodes.contains_key(&i) {
                ya_decalage = true;
            } else {
                i_find += 1;
                if ya_decalage {
                    let a_node = self
                        .the_nodes
                        .get(&i)
                        .expect("MAT_Graph::CompactNodes")
                        .clone();
                    a_node.write().unwrap().set_index(i_find);
                    self.the_nodes.insert(i_find, a_node);
                    self.the_nodes.remove(&i);
                }
            }
            i += 1;
        }
    }

    /// OCCT MAT_Graph.cxx L467-475
    pub fn change_basic_elts(&mut self, new_map: &HashMap<i32, HandleMatBasicElt>) {
        self.the_basic_elts = new_map.clone();
        for (key, value) in self.the_basic_elts.iter() {
            value.write().unwrap().set_index(*key);
        }
    }

    /// OCCT MAT_Graph.cxx L479-482
    pub fn change_basic_elt(&self, index: i32) -> HandleMatBasicElt {
        self.the_basic_elts
            .get(&index)
            .expect("MAT_Graph::ChangeBasicElt")
            .clone()
    }

    /// function : UpDateNodes
    /// purpose  : Mise a jour des sequence d'ARC de chaque FirstNode de chaque arc.
    ///            et stockage de chaque noeud dans la table des noeuds.
    ///            Les noeuds racines sont traites dans PERFORM.
    ///
    /// OCCT MAT_Graph.cxx L490-504 (private in OCCT)
    fn update_nodes(&mut self, ind_tab_nodes: &mut i32) {
        for i in 1..=self.number_of_arcs {
            let bout = self
                .the_arcs
                .get(&i)
                .expect("MAT_Graph::UpDateNodes")
                .read()
                .unwrap()
                .first_node()
                .expect("MAT_Graph::UpDateNodes");
            self.the_nodes.insert(*ind_tab_nodes, bout.clone());
            bout.write().unwrap().set_index(*ind_tab_nodes);
            *ind_tab_nodes -= 1;
            bout.write().unwrap().set_linked_arc(
                self.the_arcs
                    .get(&i)
                    .expect("MAT_Graph::UpDateNodes"),
            );
        }
    }
}

/// function : MakeArc
/// purpose  : Creation des <ARCS> en parcourant l'arbre issue de <aBisector>.
///
/// OCCT MAT_Graph.cxx L35-39 (forward declaration) + L510-581 (file-static
/// function, translated as a module private free function).
fn make_arc(
    a_bisector: &HandleMatBisector,
    the_basic_elts: &HashMap<i32, HandleMatBasicElt>,
    the_arcs: &mut HashMap<i32, HandleMatArc>,
    ind_tab_arcs: &mut i32,
) -> HandleMatArc {
    let current_arc: HandleMatArc;
    let mut prev_arc: Option<HandleMatArc>;
    let mut next_arc: Option<HandleMatArc> = None;
    let extremite: HandleMatNode;
    let bisector_list: HandleMatListOfBisector;
    let mut dist_ext: f64;

    // The OCCT_DEBUG_Graph std::cout blocks (L523-526, L536-539, L572-575)
    // are conditional debug output, not compiled in release OCCT.

    // OCCT L528-531 evaluates aBisector->FirstEdge()->EdgeNumber() and
    // aBisector->SecondEdge()->EdgeNumber() inline in the MAT_Arc
    // constructor call; the two edge numbers are bound here for the
    // borrow checker (identical evaluation order).
    let first_edge_number = a_bisector
        .read()
        .unwrap()
        .first_edge()
        .expect("MAT_Graph::MakeArc")
        .read()
        .unwrap()
        .edge_number();
    let second_edge_number = a_bisector
        .read()
        .unwrap()
        .second_edge()
        .expect("MAT_Graph::MakeArc")
        .read()
        .unwrap()
        .edge_number();

    current_arc = Arc::new(RwLock::new(MatArc::new(
        *ind_tab_arcs,
        a_bisector.read().unwrap().bisector_number(),
        the_basic_elts
            .get(&first_edge_number)
            .expect("MAT_Graph::MakeArc"),
        the_basic_elts
            .get(&second_edge_number)
            .expect("MAT_Graph::MakeArc"),
    )));
    dist_ext = a_bisector.read().unwrap().dist_issue_point();
    if dist_ext == INFINITE_VALUE {
        dist_ext = 1.0;
    }

    extremite = Arc::new(RwLock::new(MatNode::new(
        a_bisector.read().unwrap().issue_point(),
        &current_arc,
        dist_ext,
    )));

    current_arc
        .write()
        .unwrap()
        .set_first_node(&extremite);
    bisector_list = a_bisector.read().unwrap().list();
    bisector_list.write().unwrap().first();

    if !bisector_list.read().unwrap().more() {
        // -------------------
        // Arc sur le contour.
        // -------------------
        the_basic_elts
            .get(&second_edge_number)
            .expect("MAT_Graph::MakeArc")
            .write()
            .unwrap()
            .set_start_arc(Some(&current_arc));
        the_basic_elts
            .get(&first_edge_number)
            .expect("MAT_Graph::MakeArc")
            .write()
            .unwrap()
            .set_end_arc(Some(&current_arc));
    } else {
        prev_arc = Some(current_arc.clone());

        while bisector_list.read().unwrap().more() {
            let item = bisector_list.read().unwrap().current();
            next_arc = Some(make_arc(
                item.as_ref().expect("MAT_Graph::MakeArc"),
                the_basic_elts,
                the_arcs,
                ind_tab_arcs,
            ));
            next_arc
                .as_ref()
                .expect("MAT_Graph::MakeArc")
                .write()
                .unwrap()
                .set_second_node(&extremite);
            next_arc
                .as_ref()
                .expect("MAT_Graph::MakeArc")
                .write()
                .unwrap()
                .set_neighbour(
                    MatSide::Left,
                    &extremite,
                    prev_arc.as_ref().expect("MAT_Graph::MakeArc"),
                );
            prev_arc
                .as_ref()
                .expect("MAT_Graph::MakeArc")
                .write()
                .unwrap()
                .set_neighbour(
                    MatSide::Right,
                    &extremite,
                    next_arc.as_ref().expect("MAT_Graph::MakeArc"),
                );
            prev_arc = next_arc.clone();
            bisector_list.write().unwrap().next();
        }
        current_arc.write().unwrap().set_neighbour(
            MatSide::Left,
            &extremite,
            next_arc.as_ref().expect("MAT_Graph::MakeArc"),
        );
        next_arc
            .as_ref()
            .expect("MAT_Graph::MakeArc")
            .write()
            .unwrap()
            .set_neighbour(
                MatSide::Right,
                &extremite,
                &current_arc,
            );
    }

    current_arc.write().unwrap().set_index(*ind_tab_arcs);
    the_arcs.insert(*ind_tab_arcs, current_arc.clone());
    *ind_tab_arcs = *ind_tab_arcs + 1;

    current_arc
}

impl Default for MatGraph {
    fn default() -> Self {
        Self::new()
    }
}
