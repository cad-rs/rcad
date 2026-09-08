//! OCCT MAT2d_Mat2d — the generic algorithm of computation of the
//! bisecting locus.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT2d/
//!         MAT2d_Mat2d.hxx L17-105, MAT2d_Mat2d.cxx L30-1910
//!
//! File split (AGENTS.md <2000-line rule): CreateMat lives here
//! (MAT2d_Mat2d.cxx L854-1627); CreateMatOpen lives in mat2d_mat2d_b.rs
//! (L121-852).

use std::collections::HashMap;
use std::sync::{Arc, RwLock};

use crate::topalgo::mat::{
    HandleMatBisector, HandleMatEdge, HandleMatListOfBisector, HandleMatListOfEdge, MatBisector,
    MatEdge, MatListOfBisector, MatListOfEdge,
};
use crate::topalgo::mat2d::mat2d_tool2d::Mat2dTool2d;

/// OCCT Precision::Infinite().
pub(crate) const INFINITE: f64 = rcad_kernel::precision::INFINITE_VALUE;

/// OCCT MAT2d_Mat2d (MAT2d_Mat2d.hxx L36-102).
pub struct Mat2dMat2d {
    /// OCCT bool myIsOpenResult.
    my_is_open_result: bool,
    /// OCCT int thenumberofbisectors.
    thenumberofbisectors: i32,
    /// OCCT int thenumberofedges.
    thenumberofedges: i32,
    /// OCCT bool semiInfinite.
    semi_infinite: bool,
    /// OCCT occ::handle<MAT_ListOfEdge> theedgelist.
    theedgelist: Option<HandleMatListOfEdge>,
    /// OCCT occ::handle<MAT_ListOfEdge> RemovedEdgesList.
    removed_edges_list: Option<HandleMatListOfEdge>,
    /// OCCT NCollection_DataMap<int, int> typeofbisectortoremove.
    typeofbisectortoremove: HashMap<i32, i32>,
    /// OCCT NCollection_DataMap<int, handle<MAT_Bisector>> bisectoronetoremove.
    bisectoronetoremove: HashMap<i32, HandleMatBisector>,
    /// OCCT NCollection_DataMap<int, handle<MAT_Bisector>> bisectortwotoremove.
    bisectortwotoremove: HashMap<i32, HandleMatBisector>,
    /// OCCT NCollection_DataMap<int, handle<MAT_Bisector>> bisectormap.
    bisectormap: HashMap<i32, HandleMatBisector>,
    /// OCCT occ::handle<MAT_ListOfBisector> roots.
    roots: Option<HandleMatListOfBisector>,
    /// OCCT bool isDone.
    is_done: bool,
}

impl Mat2dMat2d {
    /// OCCT MAT2d_Mat2d::MAT2d_Mat2d(IsOpenResult) (cxx L30-37).
    pub fn new(is_open_result: bool) -> Self {
        Mat2dMat2d {
            my_is_open_result: is_open_result,
            thenumberofbisectors: 0,
            thenumberofedges: 0,
            semi_infinite: false,
            theedgelist: None,
            removed_edges_list: None,
            typeofbisectortoremove: HashMap::new(),
            bisectoronetoremove: HashMap::new(),
            bisectortwotoremove: HashMap::new(),
            bisectormap: HashMap::new(),
            roots: None,
            is_done: false,
        }
    }

    /// OCCT MAT2d_Mat2d::CreateMat(MAT2d_Tool2d& atool) (cxx L854-1627) —
    /// algorithm of computation of the bisecting locus.
    pub fn create_mat(&mut self, atool: &mut Mat2dTool2d) {
        let mut interrupt = false;

        let mut edgetoremove: Option<HandleMatEdge>;
        let mut previousedge: Option<HandleMatEdge>;
        let mut currentedge: Option<HandleMatEdge>;

        let mut noofbisectorstoremove: i32;
        let mut firstbisector: Option<HandleMatBisector>;
        let mut secondbisector: Option<HandleMatBisector>;
        let mut edge: Option<HandleMatEdge>;
        let mut intersectionpoint: i32 = 0;
        let mut beginbisector: i32;
        let mut noofbisectors: i32;

        let mut nb_iter_bis: i32 = 0;
        let even_nb_iter_bis: i32 = 10;
        // OCCT NCollection_Array1<int> EdgeNumbers(1, EvenNbIterBis + 1),
        // Init(-1) — slot 0 unused.
        let mut edge_numbers = [-1i32; 12];
        let mut to_nullify_noofbisectorstoremove = false;

        let mut currentbisectorlist: Option<HandleMatListOfBisector>;

        let mut bisectortoremove: Option<HandleMatBisector>;
        let mut currentbisector: Option<HandleMatBisector>;
        // OCCT locals lastbisector / previousbisector are declared but never
        // used on this path (cxx L881-882).
        let _lastbisector: Option<HandleMatBisector> = None;
        let _previousbisector: Option<HandleMatBisector> = None;

        let noofedges;
        let number_max_of_ite;
        let toleranceofconfusion;

        noofedges = atool.number_of_items();
        toleranceofconfusion = atool.tolerance_of_confusion();
        number_max_of_ite = noofedges * noofedges;

        // OCCT NCollection_Array1<int> firstarea/lastarea/noofarea(0,
        // noofedges).
        let mut firstarea = vec![0i32; (noofedges + 1) as usize];
        let mut lastarea = vec![0i32; (noofedges + 1) as usize];
        let mut noofarea = vec![0i32; (noofedges + 1) as usize];

        let mut parama = [0i32; 2];
        let mut paramb = [0i32; 2];
        //
        // Patch to prevent infinite loop because of bad geometry (cxx
        // L899-902).
        let mut a_nb_of_narea1: i32 = 0;
        let mut a_pref_narea: i32 = 0;
        let a_nb_max_narea1: i32 = 10;
        let mut a_nb_elts = [0i32; 2];
        let mut a_count_elts = [0i32; 2];
        let mut is_break = false;

        // -----------------------------------------
        // Initialisation et remise a zero des maps.
        // -----------------------------------------
        self.bisectoronetoremove.clear();
        self.bisectortwotoremove.clear();
        self.typeofbisectortoremove.clear();
        self.bisectormap.clear();

        self.is_done = true;
        noofbisectors = noofedges;
        beginbisector = 0;

        // --------------------------------------------------------------------
        // Construction de <theedgelist> un edge correspond a un element simple
        // du contour.  (cxx L916-931)
        // --------------------------------------------------------------------
        self.theedgelist = Some(Arc::new(RwLock::new(MatListOfEdge::new())));
        self.removed_edges_list = Some(Arc::new(RwLock::new(MatListOfEdge::new())));

        {
            let theedgelist = self.theedgelist.as_ref().unwrap();
            let mut list = theedgelist.write().unwrap();
            for i in 0..noofedges {
                let mut new_edge = MatEdge::new();
                new_edge.set_edge_number(i + 1);
                new_edge.set_distance(-1.0);
                list.back_add(&Arc::new(RwLock::new(new_edge)));
            }

            // OCCT theedgelist->Loop().
            list.r#loop();
        }

        //---------------------------------------------------
        // Initialisation des bissectrices issues du contour.  (cxx L933-952)
        //---------------------------------------------------
        {
            let theedgelist = self.theedgelist.as_ref().unwrap();
            let mut list = theedgelist.write().unwrap();
            let mut dist: f64 = 0.0;
            list.first();

            for i in 0..list.number() {
                let mut new_bisector = MatBisector::new();
                new_bisector.set_index_number(i);
                let current = list.current().expect("MAT_ListOfEdge::Current");
                new_bisector.set_first_edge(&current);
                let edge_number = current.read().unwrap().edge_number();
                new_bisector
                    .set_first_vector(atool.tangent_before(edge_number, self.my_is_open_result));
                list.next();
                let current = list.current().expect("MAT_ListOfEdge::Current");
                new_bisector.set_second_edge(&current);
                let edge_number = current.read().unwrap().edge_number();
                new_bisector.set_issue_point(atool.first_point(edge_number, &mut dist));
                new_bisector.set_dist_issue_point(dist);
                new_bisector.set_second_vector(atool.tangent_after(edge_number, self.my_is_open_result));
                self.bisectormap.insert(i, Arc::new(RwLock::new(new_bisector)));
            }
        }

        //----------------------------------------------------
        // Affectation a chaque edge de ses deux bissectrices.  (cxx L954-964)
        //----------------------------------------------------
        {
            let theedgelist = self.theedgelist.as_ref().unwrap();
            let mut list = theedgelist.write().unwrap();
            list.first();

            for i in 0..list.number() {
                let current = list.current().expect("MAT_ListOfEdge::Current");
                let first = self.bisectormap[&((i - 1 + noofbisectors) % noofbisectors)].clone();
                let second = self.bisectormap[&i].clone();
                let mut cm = current.write().unwrap();
                cm.set_first_bisector(&first);
                cm.set_second_bisector(&second);
                list.next();
            }
        }

        //===========================================================================
        //                         Boucle Principale   (etape 2)
        //===========================================================================
        let mut number_of_ite: i32 = 0;

        loop {
            // ------------------------------------------------------------------
            //  Creation des geometries des bissectrices via le tool. (etape 2.1)
            // -------------------------------------------------------------------
            let a_nb_bis = noofbisectors - beginbisector;
            for i in beginbisector..noofbisectors {
                let bis = self.bisectormap[&i].clone();
                atool.create_bisector(&bis);
                self.thenumberofbisectors += 1;
            }

            // Patch to prevent infinit loop because of bad geometry
            // (cxx L992-1032).
            if a_nb_bis == 1 {
                if a_pref_narea == 1 {
                    a_nb_of_narea1 += 1;
                    let begin_map = self.bisectormap[&beginbisector].clone();
                    let begin = begin_map.read().unwrap();
                    let edge1number = begin.first_edge().unwrap().read().unwrap().edge_number();
                    let edge2number = begin.second_edge().unwrap().read().unwrap().edge_number();
                    if a_nb_elts[0] == edge1number {
                        a_count_elts[0] += 1;
                    } else {
                        a_count_elts[0] = 0;
                        a_nb_elts[0] = edge1number;
                    }
                    if a_nb_elts[1] == edge2number {
                        a_count_elts[1] += 1;
                    } else {
                        a_count_elts[1] = 0;
                        a_nb_elts[1] = edge2number;
                    }
                    if a_nb_of_narea1 >= a_nb_max_narea1
                        && (a_count_elts[0] >= a_nb_max_narea1
                            || a_count_elts[1] >= a_nb_max_narea1)
                    {
                        is_break = true;
                    }
                } else {
                    a_nb_of_narea1 = 0;
                    a_count_elts[0] = 0;
                    a_count_elts[1] = 0;
                }
            }
            a_pref_narea = a_nb_bis;

            // ---------------------------------------------
            //  Condition de sortie de la boucle principale.
            // ---------------------------------------------

            //  Modified by Sergey KHROMOV - Fri Nov 17 10:28:28 2000 Begin
            {
                let theedgelist = self.theedgelist.as_ref().unwrap();
                if theedgelist.read().unwrap().number() < 3 {
                    break;
                }
            }
            //  Modified by Sergey KHROMOV - Fri Nov 17 10:28:37 2000 End

            //---------------------------------------------------
            // loop 2: While there are bisectors to remove.
            //---------------------------------------------------
            loop {
                nb_iter_bis += 1;

                noofbisectorstoremove = 0;
                {
                    let theedgelist = self.theedgelist.as_ref().unwrap();
                    theedgelist.write().unwrap().first();
                }

                //--------------------------------------------------------------
                // Calcul des intersections des bisectrices voisines.(etape 2.2)
                //--------------------------------------------------------------

                if nb_iter_bis <= even_nb_iter_bis + 1 {
                    let n = {
                        let theedgelist = self.theedgelist.as_ref().unwrap();
                        theedgelist.read().unwrap().number()
                    };
                    edge_numbers[nb_iter_bis as usize] = n;
                } else {
                    for k in 1..=even_nb_iter_bis {
                        edge_numbers[k as usize] = edge_numbers[(k + 1) as usize];
                    }
                    let n = {
                        let theedgelist = self.theedgelist.as_ref().unwrap();
                        theedgelist.read().unwrap().number()
                    };
                    edge_numbers[(even_nb_iter_bis + 1) as usize] = n;
                }
                if edge_numbers[(even_nb_iter_bis + 1) as usize] == edge_numbers[1] {
                    to_nullify_noofbisectorstoremove = true;
                }

                {
                    let theedgelist = self.theedgelist.as_ref().unwrap().clone();
                    let mut list = theedgelist.write().unwrap();
                    let number = list.number();
                    for _i in 0..number {
                        edge = list.current();
                        let edge_h = edge.clone().expect("MAT_ListOfEdge::Current");
                        let (fb, sb, distance) = {
                            let em = edge_h.read().unwrap();
                            (em.first_bisector(), em.second_bisector(), em.distance())
                        };
                        firstbisector = fb;
                        secondbisector = sb;
                        if distance == -1.0 {
                            let fb = firstbisector.clone().expect("MAT_Bisector");
                            let sb = secondbisector.clone().expect("MAT_Bisector");
                            let dist =
                                atool.intersect_bisector(&fb, &sb, &mut intersectionpoint);
                            edge_h.write().unwrap().set_distance(dist);
                            edge_h
                                .write()
                                .unwrap()
                                .set_intersection_point(intersectionpoint);

                            if dist == INFINITE {
                                let idx_ok = {
                                    let f = fb.read().unwrap().index_number();
                                    let s = sb.read().unwrap().index_number();
                                    f >= beginbisector || s >= beginbisector
                                };
                                if idx_ok {
                                    self.intersect(atool, 0, &mut noofbisectorstoremove, &fb, &sb);
                                }
                            } else {
                                let idx_f = fb.read().unwrap().index_number();
                                if idx_f >= beginbisector {
                                    self.intersect(atool, 1, &mut noofbisectorstoremove, &fb, &sb);
                                }
                                let idx_s = sb.read().unwrap().index_number();
                                if idx_s >= beginbisector {
                                    self.intersect(atool, 2, &mut noofbisectorstoremove, &fb, &sb);
                                }
                            }
                        }
                        list.next();
                    }
                }

                //-------------------------------
                // Test de sortie de la boucle 2.
                //-------------------------------

                if to_nullify_noofbisectorstoremove {
                    noofbisectorstoremove = 0;
                }
                if noofbisectorstoremove == 0 {
                    break;
                }

                //---------------------------------------------------
                // Annulation des bissectrices a effacer. (etape 2.4)
                //---------------------------------------------------

                for i in 0..noofbisectorstoremove {
                    bisectortoremove = Some(self.bisectoronetoremove[&i].clone());

                    //---------------------------------------------------------------
                    // Destruction des bisectrices descendantes de <bisectortoremove>
                    // On descend dans l arbre jusqu a ce qu on atteigne
                    // <bisectortwotoremove(i)>.
                    //---------------------------------------------------------------

                    loop {
                        // ----------------------------------
                        // Annulation de <bisectortoremove>.
                        // ----------------------------------
                        self.thenumberofbisectors -= 1;
                        {
                            let btr = bisectortoremove.clone().unwrap();
                            currentbisectorlist = Some(btr.read().unwrap().list());
                        }
                        {
                            let cbl = currentbisectorlist.clone().unwrap();
                            let mut cl = cbl.write().unwrap();
                            cl.first();
                        }
                        currentbisector = {
                            let cbl = currentbisectorlist.clone().unwrap();
                            let cl = cbl.read().unwrap();
                            cl.first_item()
                        };
                        previousedge = {
                            let cb = currentbisector.clone().unwrap();
                            cb.read().unwrap().first_edge()
                        };
                        {
                            let theedgelist = self.theedgelist.as_ref().unwrap();
                            let mut list = theedgelist.write().unwrap();
                            list.init(previousedge.as_ref().expect("MAT_Edge"));
                        }
                        {
                            let pe = previousedge.clone().unwrap();
                            let mut pem = pe.write().unwrap();
                            pem.set_distance(-1.0);
                            let fb = pem.first_bisector().expect("MAT_Bisector");
                            fb.write().unwrap().set_second_parameter(INFINITE);
                            let sb = pem.second_bisector().expect("MAT_Bisector");
                            sb.write().unwrap().set_first_parameter(INFINITE);
                        }

                        //------------------------------------------
                        // Annulation des fils de <currentbisector>.
                        //------------------------------------------

                        loop {
                            let more = {
                                let cbl = currentbisectorlist.clone().unwrap();
                                let cl = cbl.read().unwrap();
                                cl.more()
                            };
                            if !more {
                                break;
                            }
                            currentbisector = {
                                let cbl = currentbisectorlist.clone().unwrap();
                                let cl = cbl.read().unwrap();
                                cl.current()
                            };
                            currentedge = {
                                let cb = currentbisector.clone().unwrap();
                                cb.read().unwrap().second_edge()
                            };

                            //---------------------------------------
                            // Reinsertion de l edge dans le contour.
                            //---------------------------------------
                            {
                                let theedgelist = self.theedgelist.as_ref().unwrap();
                                let mut list = theedgelist.write().unwrap();
                                list.link_after(currentedge.as_ref().expect("MAT_Edge"));
                                list.next();
                            }

                            {
                                let cb = currentbisector.clone().unwrap();
                                let ce = currentedge.clone().unwrap();
                                ce.write().unwrap().set_first_bisector(&cb);
                            }
                            {
                                let pe = previousedge.clone().unwrap();
                                let cb = currentbisector.clone().unwrap();
                                pe.write().unwrap().set_second_bisector(&cb);
                            }

                            //------------------------------------------------------
                            // Annulation de l intersection ie les fils qui
                            // ont generes l intersection sont prolonges a l infini.
                            //------------------------------------------------------

                            {
                                let cb = currentbisector.clone().unwrap();
                                let mut cbm = cb.write().unwrap();
                                cbm.set_first_parameter(INFINITE);
                                cbm.set_second_parameter(INFINITE);
                            }

                            atool.trim_bisector(currentbisector.as_ref().unwrap());

                            {
                                let ce = currentedge.clone().unwrap();
                                let mut cem = ce.write().unwrap();
                                cem.set_distance(-1.0);
                                let fb = cem.first_bisector().expect("MAT_Bisector");
                                fb.write().unwrap().set_second_parameter(INFINITE);
                                let sb = cem.second_bisector().expect("MAT_Bisector");
                                sb.write().unwrap().set_first_parameter(INFINITE);
                            }

                            previousedge = currentedge.clone();
                            {
                                let cbl = currentbisectorlist.clone().unwrap();
                                let mut cl = cbl.write().unwrap();
                                cl.next();
                            }
                        }

                        {
                            let theedgelist = self.theedgelist.as_ref().unwrap();
                            let removed = self.removed_edges_list.as_ref().unwrap();
                            let current = {
                                let list = theedgelist.read().unwrap();
                                list.current().expect("MAT_ListOfEdge::Current")
                            };
                            removed.write().unwrap().back_add(&current);
                            let mut list = theedgelist.write().unwrap();
                            list.unlink();
                        }

                        //-----------------------------------------------------------
                        // Test de sortie de la boucle d annulation des bissectrices.
                        //-----------------------------------------------------------

                        {
                            let btr = bisectortoremove.clone().unwrap();
                            let two = self.bisectortwotoremove[&i].clone();
                            let exit = btr.read().unwrap().bisector_number()
                                == two.read().unwrap().bisector_number();
                            if exit {
                                break;
                            }
                        }

                        //-----------------------
                        // Descente dans l arbre.
                        //-----------------------

                        if self.typeofbisectortoremove[&i] == 1 {
                            let next = {
                                let btr = bisectortoremove.clone().unwrap();
                                btr.read().unwrap().first_bisector()
                            };
                            bisectortoremove = next;
                        } else {
                            let next = {
                                let btr = bisectortoremove.clone().unwrap();
                                btr.read().unwrap().last_bisector()
                            };
                            bisectortoremove = next;
                        }
                    }
                }

                // Fin Boucle 2 — the OCCT `for(;;)` re-runs the
                // intersection pass after the annulation phase.
            }

            // ----------------------------------------------------------------------
            // Analyse des parametres des intersections sur les bisectrices de chaque
            // edge et determination des portions de contour a supprimees. (etape 2.5)
            // ----------------------------------------------------------------------

            {
                let theedgelist = self.theedgelist.as_ref().unwrap();
                theedgelist.write().unwrap().first();
            }

            currentbisector = {
                let theedgelist = self.theedgelist.as_ref().unwrap();
                let list = theedgelist.read().unwrap();
                let current = list.current().expect("MAT_ListOfEdge::Current");
                current.read().unwrap().first_bisector()
            };
            // OCCT L1246-1284: the parameter-analysis chain for parama[0] /
            // paramb[0].
            {
                let cb = currentbisector.clone().unwrap();
                let cbm = cb.read().unwrap();
                if cbm.first_parameter() == INFINITE && cbm.second_parameter() == INFINITE {
                    parama[0] = -1;
                    paramb[0] = -1;
                } else if cbm.first_parameter() == INFINITE {
                    parama[0] = -1;
                    paramb[0] = 1;
                } else if cbm.second_parameter() == INFINITE {
                    paramb[0] = -1;
                    parama[0] = 1;
                } else if atool.distance(&cb, cbm.first_parameter(), cbm.second_parameter())
                    > toleranceofconfusion
                {
                    if (cbm.first_parameter() - cbm.second_parameter()) * cbm.sense() > 0.0 {
                        parama[0] = -1;
                        paramb[0] = 1;
                    } else {
                        paramb[0] = -1;
                        parama[0] = 1;
                    }
                } else {
                    parama[0] = 1;
                    paramb[0] = 1;
                }
            }

            let mut narea: i32 = -1;

            {
                let theedgelist = self.theedgelist.as_ref().unwrap();
                let mut list = theedgelist.write().unwrap();
                let number = list.number();
                for _i in 0..number {
                    currentbisector = {
                        let current = list.current().expect("MAT_ListOfEdge::Current");
                        current.read().unwrap().second_bisector()
                    };
                    // OCCT L1291-1329: the chain for parama[1] / paramb[1].
                    {
                        let cb = currentbisector.clone().unwrap();
                        let cbm = cb.read().unwrap();
                        if cbm.first_parameter() == INFINITE
                            && cbm.second_parameter() == INFINITE
                        {
                            parama[1] = -1;
                            paramb[1] = -1;
                        } else if cbm.first_parameter() == INFINITE {
                            parama[1] = -1;
                            paramb[1] = 1;
                        } else if cbm.second_parameter() == INFINITE {
                            paramb[1] = -1;
                            parama[1] = 1;
                        } else if atool.distance(
                            &cb,
                            cbm.first_parameter(),
                            cbm.second_parameter(),
                        ) > toleranceofconfusion
                        {
                            if (cbm.first_parameter() - cbm.second_parameter()) * cbm.sense()
                                > 0.0
                            {
                                parama[1] = -1;
                                paramb[1] = 1;
                            } else {
                                paramb[1] = -1;
                                parama[1] = 1;
                            }
                        } else {
                            parama[1] = 1;
                            paramb[1] = 1;
                        }
                    }

                    //-----------------------------------------------------------------
                    // Test si l edge est a enlever du contour
                    // Construction des portions de contour a eliminer.
                    //-----------------------------------------------------------------

                    if paramb[0] > 0 && parama[1] > 0 {
                        let index = list.index();
                        if narea < 0 {
                            narea += 1;
                            firstarea[narea as usize] = index;
                            lastarea[narea as usize] = firstarea[narea as usize];
                            noofarea[narea as usize] = 1;
                        } else if index == lastarea[narea as usize] + 1 {
                            lastarea[narea as usize] += 1;
                            noofarea[narea as usize] += 1;
                        } else {
                            narea += 1;
                            firstarea[narea as usize] = index;
                            lastarea[narea as usize] = firstarea[narea as usize];
                            noofarea[narea as usize] = 1;
                        }
                    }
                    parama[0] = parama[1];
                    paramb[0] = paramb[1];
                    list.next();
                }
            }

            let mut compact: i32 = 0;
            let number = {
                let theedgelist = self.theedgelist.as_ref().unwrap();
                theedgelist.read().unwrap().number()
            };
            if narea > 0 && lastarea[narea as usize] == number && firstarea[0] == 1 {
                firstarea[0] = firstarea[narea as usize];
                noofarea[0] += noofarea[narea as usize];
                compact = noofarea[narea as usize];
                narea -= 1;
            }

            narea += 1;

            //------------------------------------------------------------------
            // Sortie de la boucle principale si il n y a pas d edge a eliminer.
            // (etape 2.6)
            //------------------------------------------------------------------
            //
            // Patch to break infinite loop.
            if narea == 1 && is_break {
                narea = 0;
            }
            //
            if narea == 0 {
                interrupt = true;
                break;
            }

            //----------------------------------------------------------------
            // Elimination des edges a enlever du contour
            // => Mise a jour du nouveau contour.
            // => Creation des bissectrices entre les nouvelles edges voisines.
            //----------------------------------------------------------------

            beginbisector = noofbisectors;
            let mut shift: i32 = 0;
            let mut all: i32 = 0;
            if narea == 1 && noofarea[0] == number {
                all = 1;
            }

            for i in 0..narea {
                if i == 1 {
                    shift -= compact;
                }
                {
                    let theedgelist = self.theedgelist.as_ref().unwrap();
                    theedgelist.write().unwrap().first();
                }
                edgetoremove = {
                    let theedgelist = self.theedgelist.as_ref().unwrap();
                    let mut list = theedgelist.write().unwrap();
                    list.brackets(firstarea[i as usize] - shift)
                };

                {
                    let etm = edgetoremove.clone().unwrap();
                    let em = etm.read().unwrap();
                    let fb = em.first_bisector().expect("MAT_Bisector");
                    drop(em);
                    let ipoint = {
                        let em = etm.read().unwrap();
                        em.intersection_point()
                    };
                    fb.write().unwrap().set_end_point(ipoint);
                    let second = fb.read().unwrap().second_parameter();
                    fb.write().unwrap().set_first_parameter(second);
                    atool.trim_bisector(&fb);
                }

                {
                    let mut new_bisector = MatBisector::new();
                    new_bisector.set_index_number(noofbisectors);
                    let etm = edgetoremove.clone().unwrap();
                    let em = etm.read().unwrap();
                    new_bisector.set_dist_issue_point(em.distance());
                    new_bisector.set_issue_point(em.intersection_point());
                    let fb = em.first_bisector().expect("MAT_Bisector");
                    drop(em);
                    let theedgelist = self.theedgelist.as_ref().unwrap();
                    let previous = {
                        let list = theedgelist.read().unwrap();
                        list.previous_item()
                    };
                    new_bisector.set_first_edge(&previous.expect("MAT_ListOfEdge::PreviousItem"));
                    new_bisector.add_bisector(&fb);
                    self.bisectormap
                        .insert(noofbisectors, Arc::new(RwLock::new(new_bisector)));
                }

                for j in 0..noofarea[i as usize] {
                    {
                        let theedgelist = self.theedgelist.as_ref().unwrap();
                        let removed = self.removed_edges_list.as_ref().unwrap();
                        let current = {
                            let list = theedgelist.read().unwrap();
                            list.current().expect("MAT_ListOfEdge::Current")
                        };
                        removed.write().unwrap().back_add(&current);
                        let mut list = theedgelist.write().unwrap();
                        list.unlink();
                        list.next();
                    }
                    shift += 1;

                    let etm = edgetoremove.clone().unwrap();
                    let sb = {
                        let em = etm.read().unwrap();
                        em.second_bisector().expect("MAT_Bisector")
                    };
                    if all == 0 || j + 1 != noofarea[i as usize] {
                        self.bisectormap[&noofbisectors]
                            .read()
                            .unwrap()
                            .add_bisector(&sb);
                    }
                    {
                        let em = etm.read().unwrap();
                        let ipoint = em.intersection_point();
                        sb.write().unwrap().set_end_point(ipoint);
                        let first = sb.read().unwrap().first_parameter();
                        sb.write().unwrap().set_second_parameter(first);
                    }
                    atool.trim_bisector(&sb);
                    edgetoremove = {
                        let theedgelist = self.theedgelist.as_ref().unwrap();
                        let list = theedgelist.read().unwrap();
                        list.current()
                    };
                }
                {
                    let current = edgetoremove.clone().expect("MAT_ListOfEdge::Current");
                    let new_bis = self.bisectormap[&noofbisectors].clone();
                    new_bis.write().unwrap().set_second_edge(&current);
                }

                {
                    let theedgelist = self.theedgelist.as_ref().unwrap();
                    let list = theedgelist.read().unwrap();
                    let previous = list.previous_item().expect("MAT_ListOfEdge::PreviousItem");
                    let new_bis = self.bisectormap[&noofbisectors].clone();
                    previous.write().unwrap().set_second_bisector(&new_bis);
                    let current = list.current().expect("MAT_ListOfEdge::Current");
                    current.write().unwrap().set_first_bisector(&new_bis);
                }

                {
                    let new_bis = self.bisectormap[&noofbisectors].clone();
                    let bnum = new_bis.read().unwrap();
                    let fb = bnum.first_bisector().expect("MAT_Bisector");
                    let fb_number = fb.read().unwrap().bisector_number();
                    let v = atool.tangent(fb_number);
                    drop(bnum);
                    new_bis.write().unwrap().set_first_vector(v);
                }
                {
                    let new_bis = self.bisectormap[&noofbisectors].clone();
                    let bnum = new_bis.read().unwrap();
                    let lb = bnum.last_bisector().expect("MAT_Bisector");
                    let lb_number = lb.read().unwrap().bisector_number();
                    let v = atool.tangent(lb_number);
                    drop(bnum);
                    new_bis.write().unwrap().set_second_vector(v);
                }

                noofbisectors += 1;

                {
                    let theedgelist = self.theedgelist.as_ref().unwrap();
                    let list = theedgelist.read().unwrap();
                    let previous = list.previous_item().expect("MAT_ListOfEdge::PreviousItem");
                    let mut pm = previous.write().unwrap();
                    pm.set_distance(-1.0);
                    let pfb = pm.first_bisector().expect("MAT_Bisector");
                    pfb.write().unwrap().set_second_parameter(INFINITE);
                }
                {
                    let theedgelist = self.theedgelist.as_ref().unwrap();
                    let list = theedgelist.read().unwrap();
                    let current = list.current().expect("MAT_ListOfEdge::Current");
                    let mut cm = current.write().unwrap();
                    cm.set_distance(-1.0);
                    let csb = cm.second_bisector().expect("MAT_Bisector");
                    csb.write().unwrap().set_first_parameter(INFINITE);
                }
            }

            //-----------------------------------------------------------------------
            // Test sur le nombre d iterations (cxx L1505-1517).
            //-----------------------------------------------------------------------
            if number_of_ite > number_max_of_ite {
                self.is_done = false; // Echec calcul de la carte.
                break;
            }
            number_of_ite += 1;
        } //===============================================
         //            Fin Boucle Principale.
         //===============================================

        //----------
        // etape 3.  (cxx L1523-1627)
        //----------

        //----------------------------------------------
        // interupt = True => bissectrices semi_infinies.
        //----------------------------------------------

        if interrupt {
            self.semi_infinite = true;
        } else {
            self.semi_infinite = false;

            //------------------------------------------------------------------
            // Si le nombre d edge > 1 => le nombre d edge = 2
            //              (cf test sortie boucle principale)
            // Les deux dernieres bisectrices separent les memes edges .
            // Soit elles sont confondues si calcul a l interieur, soit elles
            // sont semi-Infinies (exemple : contour compose seulement de deux
            // arcs de cercles).
            //------------------------------------------------------------------

            {
                let theedgelist = self.theedgelist.as_ref().unwrap();
                let number = theedgelist.read().unwrap().number();
                if number > 1 {
                    // Now this branch is never reachable
                    // because the case edgenumber = 2 is processed in the main
                    // loop
                    let mut list = theedgelist.write().unwrap();
                    list.first();
                    edge = list.current();
                    let edge_h = edge.clone().expect("MAT_ListOfEdge::Current");
                    let (fb, sb) = {
                        let em = edge_h.read().unwrap();
                        (
                            em.first_bisector().expect("MAT_Bisector"),
                            em.second_bisector().expect("MAT_Bisector"),
                        )
                    };
                    if fb.read().unwrap().index_number() == noofbisectors - 1 {
                        //  Modified by skv - Tue Sep 13 12:13:28 2005 IDEM Begin
                        let issue = fb.read().unwrap().issue_point();
                        if atool.trim_bisector_point(&sb, issue) {
                            if sb.read().unwrap().end_point() == 0 {
                                sb.write().unwrap().set_end_point(issue);
                            }
                            self.bisectormap[&(noofbisectors - 1)]
                                .read()
                                .unwrap()
                                .add_bisector(&sb);
                        } else {
                            self.semi_infinite = true;
                        }
                        //  Modified by skv - Tue Sep 13 12:13:28 2005 IDEM End
                    } else {
                        //  Modified by skv - Tue Sep 13 12:13:28 2005 IDEM Begin
                        let issue = sb.read().unwrap().issue_point();
                        if atool.trim_bisector_point(&fb, issue) {
                            if fb.read().unwrap().end_point() == 0 {
                                fb.write().unwrap().set_end_point(issue);
                            }
                            self.bisectormap[&(noofbisectors - 1)]
                                .read()
                                .unwrap()
                                .add_bisector(&fb);
                        } else {
                            self.semi_infinite = true;
                        }
                        //  Modified by skv - Tue Sep 13 12:13:28 2005 IDEM End
                    }
                    if !self.semi_infinite {
                        self.thenumberofbisectors -= 1;
                        let bm = self.bisectormap[&(noofbisectors - 1)].clone();
                        let mut bmw = bm.write().unwrap();
                        bmw.set_second_edge(&edge_h);
                        bmw.set_bisector_number(-1);
                    }
                }
            }
        }
        if self.semi_infinite {
            beginbisector = noofbisectors;
            let theedgelist = self.theedgelist.as_ref().unwrap();
            let mut list = theedgelist.write().unwrap();
            list.first();
            for _i in 0..list.number() {
                edge = list.current();
                let edge_h = edge.clone().expect("MAT_ListOfEdge::Current");
                let sb = edge_h
                    .read()
                    .unwrap()
                    .second_bisector()
                    .expect("MAT_Bisector");
                self.bisectormap.insert(noofbisectors, sb);
                noofbisectors += 1;
                list.next();
            }
        }

        //---------------------------
        // Recuperations des racines.  (cxx L1608-1627)
        //---------------------------

        self.roots = Some(Arc::new(RwLock::new(MatListOfBisector::new())));

        {
            let last = self.bisectormap[&(noofbisectors - 1)].clone();
            if last.read().unwrap().bisector_number() == -1 {
                self.roots = Some(last.read().unwrap().list());
                let roots = self.roots.as_ref().unwrap();
                let mut rl = roots.write().unwrap();
                rl.first();
                let current = rl.current().expect("MAT_ListOfBisector::Current");
                let fe = current.read().unwrap().first_edge().expect("MAT_Edge");
                let dist_issue = last.read().unwrap().dist_issue_point();
                fe.write().unwrap().set_distance(dist_issue);
            } else {
                let roots = self.roots.as_ref().unwrap();
                let mut rl = roots.write().unwrap();
                for i in beginbisector..noofbisectors {
                    rl.back_add(&self.bisectormap[&i].clone());
                }
            }
        }
    }

    //========================================================================
    //  function : LoadBisectorsToRemove (cxx L1633-1699)
    //  purpose  : Chargement des bisectrices a effacer.
    //========================================================================
    fn load_bisectors_to_remove(
        &mut self,
        noofbisectorstoremove: &mut i32,
        distance1: f64,
        distance2: f64,
        firstbisectortoremove1: &HandleMatBisector,
        firstbisectortoremove2: &HandleMatBisector,
        lastbisectortoremove1: &HandleMatBisector,
        lastbisectortoremove2: &HandleMatBisector,
    ) {
        let mut firstbisectortoremove = [None, None];
        let mut lastbisectortoremove = [None, None];

        firstbisectortoremove[0] = Some(firstbisectortoremove1.clone());
        firstbisectortoremove[1] = Some(firstbisectortoremove2.clone());
        lastbisectortoremove[0] = Some(lastbisectortoremove1.clone());
        lastbisectortoremove[1] = Some(lastbisectortoremove2.clone());

        let index: i32 = if distance1 < INFINITE && distance2 == INFINITE {
            0
        } else if distance2 < INFINITE && distance1 == INFINITE {
            1
        } else {
            -1
        };

        if index != -1 {
            let index = index as usize;
            let mut found = *noofbisectorstoremove;
            for j in 0..*noofbisectorstoremove {
                let one = self.bisectoronetoremove[&j].clone();
                let first = firstbisectortoremove[index].clone().unwrap();
                if one.read().unwrap().bisector_number()
                    == first.read().unwrap().bisector_number()
                {
                    found = j;
                    let two = self.bisectortwotoremove[&j].clone();
                    let last = lastbisectortoremove[index].clone().unwrap();
                    if two.read().unwrap().bisector_number()
                        < last.read().unwrap().bisector_number()
                    {
                        found = -1;
                    }
                    break;
                }
            }

            if found != -1 {
                let found = found as usize;
                self.bisectoronetoremove
                    .insert(found as i32, firstbisectortoremove[index].clone().unwrap());
                self.bisectortwotoremove
                    .insert(found as i32, lastbisectortoremove[index].clone().unwrap());
                self.typeofbisectortoremove
                    .insert(found as i32, index as i32 + 1);

                if found as i32 == *noofbisectorstoremove {
                    *noofbisectorstoremove += 1;
                }
            }
        }
    }

    //========================================================================
    //  function : Intersect (cxx L1701-1825)
    //  purpose  : Si <aside=0> Intersection de <firstbisector> avec les
    //                 descendants de <secondbisector> les plus a gauche
    //                (ie secondbisector->FirstBisector()->FirstBisector...)
    //                          Intersection de <secondbisector> avec les
    //                 descendants de <firstbisector> les plus a droite
    //                (ie firstbisector->LastBisector()->LastBisector...)
    //             Si <aside=1> Intersection de <firstbisector> avec ses
    //                descendants les plus a gauche et les plus a droite.
    //             Si <aside=2> Intersection de <secondbisector> avec ses
    //                descendants les plus a gauche et les plus a droite.
    //========================================================================
    pub(crate) fn intersect(
        &mut self,
        atool: &mut Mat2dTool2d,
        aside: i32,
        noofbisectortoremove: &mut i32,
        firstbisector: &HandleMatBisector,
        secondbisector: &HandleMatBisector,
    ) {
        let mut distance = [INFINITE; 2];
        let mut intersectionpoint: i32 = 0;
        let mut firstbisectortoremove: [Option<HandleMatBisector>; 2] = [None, None];
        let mut lastbisectortoremove: [Option<HandleMatBisector>; 2] = [None, None];

        for bisectornumber in 0..2i32 {
            if aside == 0 {
                if bisectornumber == 0 {
                    firstbisectortoremove[bisectornumber as usize] = Some(secondbisector.clone());
                } else {
                    firstbisectortoremove[bisectornumber as usize] = Some(firstbisector.clone());
                }
            } else if aside == 1 {
                firstbisectortoremove[bisectornumber as usize] = Some(firstbisector.clone());
            } else {
                firstbisectortoremove[bisectornumber as usize] = Some(secondbisector.clone());
            }

            let mut lastbisector = firstbisectortoremove[bisectornumber as usize]
                .clone()
                .unwrap();

            let previousbisector: HandleMatBisector;
            if aside == 0 {
                previousbisector = firstbisectortoremove[bisectornumber as usize]
                    .clone()
                    .unwrap();
            } else {
                let is_empty = lastbisector
                    .read()
                    .unwrap()
                    .list()
                    .read()
                    .unwrap()
                    .is_empty();
                if is_empty {
                    continue;
                }

                if bisectornumber == 0 {
                    previousbisector = lastbisector
                        .read()
                        .unwrap()
                        .first_bisector()
                        .expect("MAT_Bisector::FirstBisector");
                } else {
                    previousbisector = lastbisector
                        .read()
                        .unwrap()
                        .last_bisector()
                        .expect("MAT_Bisector::LastBisector");
                }
            }

            let mut previousbisector = previousbisector;
            let mut distant = distance[bisectornumber as usize];
            loop {
                let is_empty = previousbisector
                    .read()
                    .unwrap()
                    .list()
                    .read()
                    .unwrap()
                    .is_empty();
                if is_empty {
                    break;
                }

                if bisectornumber == 0 {
                    let next = previousbisector
                        .read()
                        .unwrap()
                        .first_bisector()
                        .expect("MAT_Bisector::FirstBisector");
                    previousbisector = next;
                } else {
                    let next = previousbisector
                        .read()
                        .unwrap()
                        .last_bisector()
                        .expect("MAT_Bisector::LastBisector");
                    previousbisector = next;
                }

                if aside == 1 || (aside == 0 && bisectornumber == 0) {
                    let saveparameter = previousbisector.read().unwrap().first_parameter();
                    distant = atool
                        .intersect_bisector(firstbisector, &previousbisector, &mut intersectionpoint);
                    previousbisector
                        .write()
                        .unwrap()
                        .set_first_parameter(saveparameter);
                } else {
                    let saveparameter = previousbisector.read().unwrap().second_parameter();
                    distant = atool
                        .intersect_bisector(&previousbisector, secondbisector, &mut intersectionpoint);
                    previousbisector
                        .write()
                        .unwrap()
                        .set_second_parameter(saveparameter);
                }

                if distant < INFINITE {
                    distance[bisectornumber as usize] = distant;
                    lastbisectortoremove[bisectornumber as usize] = Some(lastbisector.clone());
                }

                lastbisector = previousbisector.clone();
            }
        }

        //---------------------------------------
        // Chargement des bissectrices a effacer.
        //---------------------------------------

        // OCCT passes the two C arrays; slots left null on this path are
        // never dereferenced (the corresponding distance stayed infinite).
        self.load_bisectors_to_remove(
            noofbisectortoremove,
            distance[0],
            distance[1],
            firstbisectortoremove[0].as_ref().expect("MAT_Bisector"),
            firstbisectortoremove[1].as_ref().expect("MAT_Bisector"),
            lastbisectortoremove[0].as_ref().expect("MAT_Bisector"),
            lastbisectortoremove[1].as_ref().expect("MAT_Bisector"),
        );
    }

    /// OCCT MAT2d_Mat2d::Init() (cxx L1829-1832) — initialize an iterator on
    /// the set of the roots of the trees of bisectors.
    pub fn init(&mut self) {
        let roots = self.roots.as_ref().unwrap();
        roots.write().unwrap().first();
    }

    /// OCCT MAT2d_Mat2d::More() const (cxx L1836-1839) — false if there is
    /// no more root.
    pub fn more(&self) -> bool {
        let roots = self.roots.as_ref().unwrap();
        roots.read().unwrap().more()
    }

    /// OCCT MAT2d_Mat2d::Next() (cxx L1843-1846).
    pub fn next(&mut self) {
        let roots = self.roots.as_ref().unwrap();
        roots.write().unwrap().next();
    }

    /// OCCT MAT2d_Mat2d::Bisector() const (cxx L1850-1853) — the current
    /// root.
    pub fn bisector(&self) -> HandleMatBisector {
        let roots = self.roots.as_ref().unwrap();
        roots
            .read()
            .unwrap()
            .current()
            .expect("MAT_ListOfBisector::Current")
    }

    /// OCCT MAT2d_Mat2d::NumberOfBisectors() const (cxx L1857-1860).
    pub fn number_of_bisectors(&self) -> i32 {
        self.thenumberofbisectors
    }

    /// OCCT MAT2d_Mat2d::SemiInfinite() const (cxx L1864-1867).
    pub fn semi_infinite(&self) -> bool {
        self.semi_infinite
    }

    /// OCCT MAT2d_Mat2d::IsDone() const (cxx L1871-1874).
    pub fn is_done(&self) -> bool {
        self.is_done
    }
}

/// OCCT ~MAT2d_Mat2d (cxx L1878-1910) — the handle cycles are broken
/// explicitly.  Pre-written surface (MAT package): MatBisector::
/// clear_first_edge()/clear_second_edge(), MatEdge::clear_first_bisector()/
/// clear_second_bisector() — reconciled at joint compile.
impl Drop for Mat2dMat2d {
    fn drop(&mut self) {
        for a_bisector in self.bisectormap.values() {
            let mut b = a_bisector.write().unwrap();
            b.clear_first_edge();
            b.clear_second_edge();
        }

        if let Some(theedgelist) = &self.theedgelist {
            let mut list = theedgelist.write().unwrap();
            list.first();
            for _ in 0..list.number() {
                let an_edge = list.current().expect("MAT_ListOfEdge::Current");
                let mut em = an_edge.write().unwrap();
                em.clear_first_bisector();
                em.clear_second_bisector();
                list.next();
            }
        }
        if let Some(removed_edges_list) = &self.removed_edges_list {
            let mut list = removed_edges_list.write().unwrap();
            list.first();
            for _ in 0..list.number() {
                let an_edge = list.current().expect("MAT_ListOfEdge::Current");
                let mut em = an_edge.write().unwrap();
                em.clear_first_bisector();
                em.clear_second_bisector();
                list.next();
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Field accessors shared with mat2d_mat2d_b.rs (CreateMatOpen).  The struct
// fields are private to this module; these helpers expose exactly the member
// operations OCCT performs on them.
// ---------------------------------------------------------------------------
impl Mat2dMat2d {
    pub(crate) fn theedgelist_handle(&self) -> &HandleMatListOfEdge {
        self.theedgelist.as_ref().unwrap()
    }

    pub(crate) fn removed_edges_list_handle(&self) -> &HandleMatListOfEdge {
        self.removed_edges_list.as_ref().unwrap()
    }

    pub(crate) fn roots_handle(&self) -> &HandleMatListOfBisector {
        self.roots.as_ref().unwrap()
    }

    pub(crate) fn my_is_open_result(&self) -> bool {
        self.my_is_open_result
    }

    pub(crate) fn bisectormap_get(&self, i: i32) -> HandleMatBisector {
        self.bisectormap[&i].clone()
    }

    pub(crate) fn bisectormap_insert(&mut self, i: i32, b: HandleMatBisector) {
        self.bisectormap.insert(i, b);
    }

    pub(crate) fn bisectormap_clear(&mut self) {
        self.bisectormap.clear();
    }

    pub(crate) fn bisectoronetoremove_get(&self, i: i32) -> HandleMatBisector {
        self.bisectoronetoremove[&i].clone()
    }

    pub(crate) fn bisectoronetoremove_clear(&mut self) {
        self.bisectoronetoremove.clear();
    }

    pub(crate) fn bisectortwotoremove_get(&self, i: i32) -> HandleMatBisector {
        self.bisectortwotoremove[&i].clone()
    }

    pub(crate) fn bisectortwotoremove_clear(&mut self) {
        self.bisectortwotoremove.clear();
    }

    pub(crate) fn typeofbisectortoremove_get(&self, i: i32) -> i32 {
        self.typeofbisectortoremove[&i]
    }

    pub(crate) fn typeofbisectortoremove_clear(&mut self) {
        self.typeofbisectortoremove.clear();
    }

    pub(crate) fn theedgelist_set(&mut self, l: Option<HandleMatListOfEdge>) {
        self.theedgelist = l;
    }

    pub(crate) fn removed_edges_list_set(&mut self, l: Option<HandleMatListOfEdge>) {
        self.removed_edges_list = l;
    }

    pub(crate) fn roots_set(&mut self, l: Option<HandleMatListOfBisector>) {
        self.roots = l;
    }

    pub(crate) fn is_done_set(&mut self, v: bool) {
        self.is_done = v;
    }

    pub(crate) fn semi_infinite_flag(&self) -> bool {
        self.semi_infinite
    }

    pub(crate) fn semi_infinite_set(&mut self, v: bool) {
        self.semi_infinite = v;
    }

    pub(crate) fn thenumberofbisectors_inc(&mut self) {
        self.thenumberofbisectors += 1;
    }

    pub(crate) fn thenumberofbisectors_dec(&mut self) {
        self.thenumberofbisectors -= 1;
    }
}
