//! OCCT MAT2d_Mat2d::CreateMatOpen — the bisecting locus computation for
//! open wires (companion of mat2d_mat2d.rs, see its module docs for the
//! file split rationale).
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT2d/
//!         MAT2d_Mat2d.cxx L121-852

use std::sync::{Arc, RwLock};

use crate::topalgo::mat::{MatBisector, MatEdge};
use crate::topalgo::mat2d::mat2d_mat2d::{Mat2dMat2d, INFINITE};
use crate::topalgo::mat2d::mat2d_tool2d::Mat2dTool2d;

impl Mat2dMat2d {
    /// OCCT MAT2d_Mat2d::CreateMatOpen(MAT2d_Tool2d& atool)
    /// (MAT2d_Mat2d.cxx L121-852) — algorithm of computation of the
    /// bisecting locus for open wire.
    pub fn create_mat_open(&mut self, atool: &mut Mat2dTool2d) {
        let mut interrupt = false;

        let mut edgetoremove: Option<Arc<RwLock<MatEdge>>>;
        let mut previousedge: Option<Arc<RwLock<MatEdge>>>;
        let mut currentedge: Option<Arc<RwLock<MatEdge>>>;

        let mut noofbisectorstoremove: i32;
        let mut firstbisector: Option<Arc<RwLock<MatBisector>>>;
        let mut secondbisector: Option<Arc<RwLock<MatBisector>>>;
        let mut edge: Option<Arc<RwLock<MatEdge>>>;
        let mut intersectionpoint: i32 = 0;
        let mut beginbisector: i32;
        let mut noofbisectors: i32;

        let mut nb_iter_bis: i32 = 0;
        let even_nb_iter_bis: i32 = 10;
        // OCCT NCollection_Array1<int> EdgeNumbers(1, EvenNbIterBis + 1),
        // Init(-1) — slot 0 unused.
        let mut edge_numbers = [-1i32; 12];
        let mut to_nullify_noofbisectorstoremove = false;

        let mut currentbisectorlist: Option<Arc<RwLock<crate::topalgo::mat::MatListOfBisector>>>;

        let mut bisectortoremove: Option<Arc<RwLock<MatBisector>>>;
        let mut currentbisector: Option<Arc<RwLock<MatBisector>>>;
        // OCCT locals lastbisector / previousbisector are declared but never
        // used on this path (cxx L148-149).
        let _lastbisector: Option<Arc<RwLock<MatBisector>>> = None;
        let _previousbisector: Option<Arc<RwLock<MatBisector>>> = None;

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

        // -----------------------------------------
        // Initialisation et remise a zero des maps.
        // -----------------------------------------
        self.bisectoronetoremove_clear();
        self.bisectortwotoremove_clear();
        self.typeofbisectortoremove_clear();
        self.bisectormap_clear();

        self.is_done_set(true);
        noofbisectors = noofedges - 1;
        beginbisector = 0;

        // --------------------------------------------------------------------
        // Construction de <theedgelist> un edge correspond a un element simple
        // du contour.  (cxx L179-194)
        // --------------------------------------------------------------------
        self.theedgelist_set(Some(Arc::new(RwLock::new(
            crate::topalgo::mat::MatListOfEdge::new(),
        ))));
        self.removed_edges_list_set(Some(Arc::new(RwLock::new(
            crate::topalgo::mat::MatListOfEdge::new(),
        ))));

        {
            let theedgelist = self.theedgelist_handle().clone();
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
        // Initialisation des bissectrices issues du contour.  (cxx L196-215)
        //---------------------------------------------------
        {
            let theedgelist = self.theedgelist_handle().clone();
            let mut list = theedgelist.write().unwrap();
            let mut dist: f64 = 0.0;
            list.first();

            for i in 0..list.number() - 1 {
                let mut new_bisector = MatBisector::new();
                new_bisector.set_index_number(i);
                let current = list.current().expect("MAT_ListOfEdge::Current");
                new_bisector.set_first_edge(&current);
                let edge_number = current.read().unwrap().edge_number();
                new_bisector
                    .set_first_vector(atool.tangent_before(edge_number, self.my_is_open_result()));
                list.next();
                let current = list.current().expect("MAT_ListOfEdge::Current");
                new_bisector.set_second_edge(&current);
                let edge_number = current.read().unwrap().edge_number();
                new_bisector.set_issue_point(atool.first_point(edge_number, &mut dist));
                new_bisector.set_dist_issue_point(dist);
                new_bisector
                    .set_second_vector(atool.tangent_after(edge_number, self.my_is_open_result()));
                self.bisectormap_insert(i, Arc::new(RwLock::new(new_bisector)));
            }
        }

        //----------------------------------------------------
        // Affectation a chaque edge de ses deux bissectrices.  (cxx L217-233)
        //----------------------------------------------------
        {
            let theedgelist = self.theedgelist_handle().clone();
            let mut list = theedgelist.write().unwrap();
            list.first();
            {
                let current = list.current().expect("MAT_ListOfEdge::Current");
                let first = self.bisectormap_get(0);
                let mut cm = current.write().unwrap();
                cm.set_first_bisector(&first);
                cm.set_second_bisector(&first);
            }
            list.next();

            for i in 1..list.number() - 1 {
                let current = list.current().expect("MAT_ListOfEdge::Current");
                let first = self.bisectormap_get(i - 1);
                let second = self.bisectormap_get(i);
                let mut cm = current.write().unwrap();
                cm.set_first_bisector(&first);
                cm.set_second_bisector(&second);
                list.next();
            }

            {
                let current = list.current().expect("MAT_ListOfEdge::Current");
                let number = list.number();
                let first = self.bisectormap_get(number - 2);
                let mut cm = current.write().unwrap();
                cm.set_first_bisector(&first);
                cm.set_second_bisector(&first);
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
            for i in beginbisector..noofbisectors {
                let bis = self.bisectormap_get(i);
                atool.create_bisector(&bis);
                self.thenumberofbisectors_inc();
            }

            // ---------------------------------------------
            //  Condition de sortie de la boucle principale.
            // ---------------------------------------------

            //  Modified by Sergey KHROMOV - Fri Nov 17 10:28:28 2000 Begin
            {
                let theedgelist = self.theedgelist_handle().clone();
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
                    let theedgelist = self.theedgelist_handle().clone();
                    let mut list = theedgelist.write().unwrap();
                    list.first();
                    list.next();
                }

                //--------------------------------------------------------------
                // Calcul des intersections des bisectrices voisines.(etape 2.2)
                //--------------------------------------------------------------

                if nb_iter_bis <= even_nb_iter_bis + 1 {
                    let n = {
                        let theedgelist = self.theedgelist_handle().clone();
                        theedgelist.read().unwrap().number()
                    };
                    edge_numbers[nb_iter_bis as usize] = n;
                } else {
                    for k in 1..=even_nb_iter_bis {
                        edge_numbers[k as usize] = edge_numbers[(k + 1) as usize];
                    }
                    let n = {
                        let theedgelist = self.theedgelist_handle().clone();
                        theedgelist.read().unwrap().number()
                    };
                    edge_numbers[(even_nb_iter_bis + 1) as usize] = n;
                }
                if edge_numbers[(even_nb_iter_bis + 1) as usize] == edge_numbers[1] {
                    to_nullify_noofbisectorstoremove = true;
                }

                {
                    let theedgelist = self.theedgelist_handle().clone();
                    let mut list = theedgelist.write().unwrap();
                    let number = list.number();
                    for _i in 1..number - 1 {
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
                    bisectortoremove = Some(self.bisectoronetoremove_get(i));

                    //---------------------------------------------------------------
                    // Destruction des bisectrices descendantes de <bisectortoremove>
                    // On descend dans l arbre jusqu a ce qu on atteigne
                    // <bisectortwotoremove(i)>.
                    //---------------------------------------------------------------

                    loop {
                        // ----------------------------------
                        // Annulation de <bisectortoremove>.
                        // ----------------------------------
                        self.thenumberofbisectors_dec();
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
                            let theedgelist = self.theedgelist_handle().clone();
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
                                let theedgelist = self.theedgelist_handle().clone();
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
                            let theedgelist = self.theedgelist_handle().clone();
                            let removed = self.removed_edges_list_handle().clone();
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
                            let two = self.bisectortwotoremove_get(i);
                            let exit = btr.read().unwrap().bisector_number()
                                == two.read().unwrap().bisector_number();
                            if exit {
                                break;
                            }
                        }

                        //-----------------------
                        // Descente dans l arbre.
                        //-----------------------

                        if self.typeofbisectortoremove_get(i) == 1 {
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
                let theedgelist = self.theedgelist_handle().clone();
                let mut list = theedgelist.write().unwrap();
                list.first();
                list.next();
            }

            currentbisector = {
                let theedgelist = self.theedgelist_handle().clone();
                let list = theedgelist.read().unwrap();
                let current = list.current().expect("MAT_ListOfEdge::Current");
                current.read().unwrap().first_bisector()
            };
            // OCCT L476-514: the parameter-analysis chain for parama[0] /
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
                let theedgelist = self.theedgelist_handle().clone();
                let mut list = theedgelist.write().unwrap();
                let number = list.number();
                for _i in 1..number - 1 {
                    currentbisector = {
                        let current = list.current().expect("MAT_ListOfEdge::Current");
                        current.read().unwrap().second_bisector()
                    };
                    // OCCT L520-559: the chain for parama[1] / paramb[1].
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
                            // OCCT firstarea(++narea) = Index.
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
                let theedgelist = self.theedgelist_handle().clone();
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
                    let theedgelist = self.theedgelist_handle().clone();
                    let mut list = theedgelist.write().unwrap();
                    list.first();
                    list.next();
                }
                edgetoremove = {
                    let theedgelist = self.theedgelist_handle().clone();
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
                    let theedgelist = self.theedgelist_handle().clone();
                    let previous = {
                        let list = theedgelist.read().unwrap();
                        list.previous_item()
                    };
                    new_bisector.set_first_edge(&previous.expect("MAT_ListOfEdge::PreviousItem"));
                    new_bisector.add_bisector(&fb);
                    self.bisectormap_insert(noofbisectors, Arc::new(RwLock::new(new_bisector)));
                }

                for j in 0..noofarea[i as usize] {
                    {
                        let theedgelist = self.theedgelist_handle().clone();
                        let removed = self.removed_edges_list_handle().clone();
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
                        self.bisectormap_get(noofbisectors)
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
                        let theedgelist = self.theedgelist_handle().clone();
                        let list = theedgelist.read().unwrap();
                        list.current()
                    };
                }
                {
                    let current = edgetoremove.clone().expect("MAT_ListOfEdge::Current");
                    let new_bis = self.bisectormap_get(noofbisectors);
                    new_bis.write().unwrap().set_second_edge(&current);
                }

                {
                    let theedgelist = self.theedgelist_handle().clone();
                    let list = theedgelist.read().unwrap();
                    let previous = list.previous_item().expect("MAT_ListOfEdge::PreviousItem");
                    let new_bis = self.bisectormap_get(noofbisectors);
                    previous.write().unwrap().set_second_bisector(&new_bis);
                    let current = list.current().expect("MAT_ListOfEdge::Current");
                    current.write().unwrap().set_first_bisector(&new_bis);
                }

                {
                    let new_bis = self.bisectormap_get(noofbisectors);
                    let bnum = new_bis.read().unwrap();
                    let fb = bnum.first_bisector().expect("MAT_Bisector");
                    let fb_number = fb.read().unwrap().bisector_number();
                    let v = atool.tangent(fb_number);
                    drop(bnum);
                    new_bis.write().unwrap().set_first_vector(v);
                }
                {
                    let new_bis = self.bisectormap_get(noofbisectors);
                    let bnum = new_bis.read().unwrap();
                    let lb = bnum.last_bisector().expect("MAT_Bisector");
                    let lb_number = lb.read().unwrap().bisector_number();
                    let v = atool.tangent(lb_number);
                    drop(bnum);
                    new_bis.write().unwrap().set_second_vector(v);
                }

                noofbisectors += 1;

                {
                    let theedgelist = self.theedgelist_handle().clone();
                    let list = theedgelist.read().unwrap();
                    let previous = list.previous_item().expect("MAT_ListOfEdge::PreviousItem");
                    let mut pm = previous.write().unwrap();
                    pm.set_distance(-1.0);
                    let pfb = pm.first_bisector().expect("MAT_Bisector");
                    pfb.write().unwrap().set_second_parameter(INFINITE);
                }
                {
                    let theedgelist = self.theedgelist_handle().clone();
                    let list = theedgelist.read().unwrap();
                    let current = list.current().expect("MAT_ListOfEdge::Current");
                    let mut cm = current.write().unwrap();
                    cm.set_distance(-1.0);
                    let csb = cm.second_bisector().expect("MAT_Bisector");
                    csb.write().unwrap().set_first_parameter(INFINITE);
                }
            }

            //-----------------------------------------------------------------------
            // Test sur le nombre d iterations (cxx L729-741).
            //-----------------------------------------------------------------------
            if number_of_ite > number_max_of_ite {
                self.is_done_set(false); // Echec calcul de la carte.
                break;
            }
            number_of_ite += 1;
        } //===============================================
         //            Fin Boucle Principale.
         //===============================================

        //----------
        // etape 3.  (cxx L747-852)
        //----------

        //----------------------------------------------
        // interupt = True => bissectrices semi_infinies.
        //----------------------------------------------

        if interrupt {
            self.semi_infinite_set(true);
        } else {
            self.semi_infinite_set(false);

            //------------------------------------------------------------------
            // Si le nombre d edge > 1 => le nombre d edge = 2
            //              (cf test sortie boucle principale)
            // Les deux dernieres bisectrices separent les memes edges .
            // Soit elles sont confondues si calcul a l interieur, soit elles
            // sont semi-Infinies (exemple : contour compose seulement de deux
            // arcs de cercles).
            //------------------------------------------------------------------

            {
                let theedgelist = self.theedgelist_handle().clone();
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
                            self.bisectormap_get(noofbisectors - 1)
                                .read()
                                .unwrap()
                                .add_bisector(&sb);
                        } else {
                            self.semi_infinite_set(true);
                        }
                        //  Modified by skv - Tue Sep 13 12:13:28 2005 IDEM End
                    } else {
                        //  Modified by skv - Tue Sep 13 12:13:28 2005 IDEM Begin
                        let issue = sb.read().unwrap().issue_point();
                        if atool.trim_bisector_point(&fb, issue) {
                            if fb.read().unwrap().end_point() == 0 {
                                fb.write().unwrap().set_end_point(issue);
                            }
                            self.bisectormap_get(noofbisectors - 1)
                                .read()
                                .unwrap()
                                .add_bisector(&fb);
                        } else {
                            self.semi_infinite_set(true);
                        }
                        //  Modified by skv - Tue Sep 13 12:13:28 2005 IDEM End
                    }
                    if !self.semi_infinite_flag() {
                        self.thenumberofbisectors_dec();
                        let bm = self.bisectormap_get(noofbisectors - 1);
                        let mut bmw = bm.write().unwrap();
                        bmw.set_second_edge(&edge_h);
                        bmw.set_bisector_number(-1);
                    }
                }
            }
        }
        if self.semi_infinite_flag() {
            beginbisector = noofbisectors;
            let theedgelist = self.theedgelist_handle().clone();
            let mut list = theedgelist.write().unwrap();
            list.first();
            for _i in 1..list.number() {
                edge = list.current();
                let edge_h = edge.clone().expect("MAT_ListOfEdge::Current");
                let sb = edge_h
                    .read()
                    .unwrap()
                    .second_bisector()
                    .expect("MAT_Bisector");
                self.bisectormap_insert(noofbisectors, sb);
                noofbisectors += 1;
                list.next();
            }
        }

        //---------------------------
        // Recuperations des racines.  (cxx L833-852)
        //---------------------------

        self.roots_set(Some(Arc::new(RwLock::new(
            crate::topalgo::mat::MatListOfBisector::new(),
        ))));

        {
            let last = self.bisectormap_get(noofbisectors - 1);
            if last.read().unwrap().bisector_number() == -1 {
                let list = last.read().unwrap().list();
                self.roots_set(Some(list.clone()));
                let roots = self.roots_handle().clone();
                let mut rl = roots.write().unwrap();
                rl.first();
                let current = rl.current().expect("MAT_ListOfBisector::Current");
                let fe = current.read().unwrap().first_edge().expect("MAT_Edge");
                let dist_issue = last.read().unwrap().dist_issue_point();
                fe.write().unwrap().set_distance(dist_issue);
            } else {
                let roots = self.roots_handle().clone();
                let mut rl = roots.write().unwrap();
                for i in beginbisector..noofbisectors {
                    rl.back_add(&self.bisectormap_get(i));
                }
            }
        }
    }
}
