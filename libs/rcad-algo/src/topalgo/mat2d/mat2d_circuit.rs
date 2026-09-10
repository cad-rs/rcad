//! OCCT MAT2d_Circuit — constructs a circuit on a set of lines; EquiCircuit
//! gives a circuit passing by all the lines in a set and all the connexions
//! of the minipath associated.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT2d/
//!         MAT2d_Circuit.hxx L17-117, MAT2d_Circuit.cxx L34-939

use std::collections::HashMap;


use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Curve2dEval};

use crate::geomalgo::geom2d_int::{AdaptorOffsetCurve, Curve2dAdaptor, GInter};
use crate::topalgo::bisector::bisector::GeomAbsJoinType;
use crate::topalgo::mat2d::mat2d_connexion::HandleMat2dConnexion;
use crate::topalgo::mat2d::mat2d_mini_path::Mat2dMiniPath;
use crate::topalgo::mat2d::{Geom2dGeometry, Mat2dBiInt};

/// OCCT static double CrossProd(Geom1, Geom2, DotProd) (cxx L52-54,
/// body L860-872) — the vectorial and scalar products between the directions
/// of the tangents at the end of Geom1 and at the beginning of Geom2.
fn cross_prod(geom1: &Geom2dGeometry, geom2: &Geom2dGeometry, dot_prod: &mut f64) -> f64 {
    let curve1 = geom1.curve();
    let dir1 = curve1.dn(curve1.last_parameter(), 1).normalize_or_zero();
    let curve2 = geom2.curve();
    let dir2 = curve2.dn(curve2.first_parameter(), 1).normalize_or_zero();
    *dot_prod = dir1.dot(dir2);
    // OCCT Dir1 ^ Dir2 — the 2d cross product.
    dir1.x * dir2.y - dir1.y * dir2.x
}

/// OCCT static void SubSequence(S1, IF, IL, S2) (cxx L346-356).
fn sub_sequence(s1: &[Geom2dGeometry], i_f: i32, i_l: i32, s2: &mut Vec<Geom2dGeometry>) {
    s2.clear();
    for i in i_f..=i_l {
        s2.push(s1[(i - 1) as usize].clone());
    }
}

/// OCCT Geom2d_TrimmedCurve::ReversedParameter(U) (Geom2d_TrimmedCurve.cxx
/// L91-94: `return basisCurve->ReversedParameter(U)`) — the parameter of
/// the same point on the reversed curve.
fn trimmed_reversed_parameter(curve: &Curve2d, u: f64) -> f64 {
    curve.inner().reversed_parameter(u)
}

/// OCCT Geom2d_TrimmedCurve::Reverse() (Geom2d_TrimmedCurve.cxx L81-87):
/// U1 = basis.ReversedParameter(Last), U2 = basis.ReversedParameter(First),
/// basis.Reverse(), SetTrim(U1, U2).
fn trimmed_reverse(curve: &mut Curve2d) {
    if let Curve2d::Trimmed(tc) = curve {
        let u1 = tc.curve.reversed_parameter(tc.t_max);
        let u2 = tc.curve.reversed_parameter(tc.t_min);
        reverse_basis(&mut tc.curve);
        tc.t_min = u1;
        tc.t_max = u2;
    }
}

/// OCCT Geom2d_<Type>::Reverse() on the basis curve.  GAP note (dependency,
/// plan §0.6): the kernel Geom2d port currently carries the line and circle
/// reversals (Geom2d_Line.cxx / Geom2d_Circle.cxx Reverse); other basis
/// types keep the loud failure path until their kernel reversal lands.
fn reverse_basis(basis: &mut Curve2d) {
    match basis {
        // OCCT Geom2d_Line::Reverse: the direction is reversed.
        Curve2d::Line(l) => l.direction = -l.direction,
        // OCCT Geom2d_Circle::Reverse: gp_Circ2d::Reverse reverses the
        // Y direction of the local frame.
        Curve2d::Circle(c) => c.y_dir = -c.y_dir,
        _ => unimplemented!("GAP: kernel Curve2d::Reverse for this basis type"),
    }
}

/// OCCT MAT2d_Circuit (MAT2d_Circuit.hxx L39-114).
pub struct Mat2dCircuit {
    /// OCCT double direction.
    direction: f64,
    /// OCCT NCollection_Sequence<handle<Geom2d_Geometry>> geomElements.
    geom_elements: Vec<Geom2dGeometry>,
    /// OCCT NCollection_DataMap<int, handle<MAT2d_Connexion>> connexionMap.
    connexion_map: HashMap<i32, HandleMat2dConnexion>,
    /// OCCT NCollection_DataMap<MAT2d_BiInt, NCollection_Sequence<int>>
    /// linkRefEqui.
    link_ref_equi: HashMap<Mat2dBiInt, Vec<i32>>,
    /// OCCT NCollection_Sequence<int> linesLength.
    lines_length: Vec<i32>,
    /// OCCT GeomAbs_JoinType myJoinType.
    my_join_type: GeomAbsJoinType,
    /// OCCT bool myIsOpenResult.
    my_is_open_result: bool,
}

impl Mat2dCircuit {
    /// OCCT MAT2d_Circuit::MAT2d_Circuit(aJoinType, IsOpenResult)
    /// (cxx L58-63).
    pub fn new(a_join_type: GeomAbsJoinType, is_open_result: bool) -> Self {
        Mat2dCircuit {
            direction: 0.0,
            geom_elements: Vec::new(),
            connexion_map: HashMap::new(),
            link_ref_equi: HashMap::new(),
            lines_length: Vec::new(),
            my_join_type: a_join_type,
            my_is_open_result: is_open_result,
        }
    }

    /// OCCT MAT2d_Circuit::Perform(FigItem, IsClosed, IndRefLine, Trigo)
    /// (cxx L67-219).
    pub fn perform(
        &mut self,
        fig_item: &mut Vec<Vec<Geom2dGeometry>>,
        is_closed: &[bool],
        ind_ref_line: i32,
        trigo: bool,
    ) {
        let nb_lines = fig_item.len() as i32;
        // OCCT NCollection_Array1<bool> Open(1, NbLines) — 0-based here.
        let mut open = vec![false; nb_lines as usize];
        // OCCT NCollection_Sequence<handle> SVide — the empty connexion set.
        let s_vide: Vec<HandleMat2dConnexion> = Vec::new();
        // OCCT handle<MAT2d_Connexion> ConnexionNul — a null handle.
        let connexion_nul: Option<HandleMat2dConnexion> = None;

        if trigo {
            self.direction = 1.0;
        } else {
            self.direction = -1.0;
        }

        //---------------------
        // Reinitialisation SD.
        //---------------------
        self.geom_elements.clear();
        self.connexion_map.clear();
        self.link_ref_equi.clear();
        self.lines_length.clear();

        //----------------------------
        // Detection Lignes ouvertes.  (OCCT L99-122)
        //----------------------------
        for i in 1..=nb_lines {
            let line = &fig_item[(i - 1) as usize];
            // OCCT: Curve = down_cast<Trimmed>(FigItem.Value(i).First());
            // P1 = Curve->StartPoint(); ... P2 = Curve->EndPoint().
            let first = line.first().expect("MAT2d_Circuit::Perform");
            let last = line.last().expect("MAT2d_Circuit::Perform");
            let p1 = first.curve().value(first.curve().first_parameter());
            let p2 = last.curve().value(last.curve().last_parameter());
            // Modified by Sergey KHROMOV - Wed Mar 6 16:59:01 2002:
            if is_closed[(i - 1) as usize] {
                open[(i - 1) as usize] = false;
            } else if p1.abs_diff_eq(p2, rcad_kernel::precision::CONFUSION) {
                open[(i - 1) as usize] = false;
            } else {
                open[(i - 1) as usize] = true;
            }
        }

        //---------------------------------------------------------------
        // Insertion des cassures saillantes ou
        // ajout des extremites de chaque courbe si la ligne est ouverte.
        //---------------------------------------------------------------
        for i in 1..=nb_lines {
            if open[(i - 1) as usize] {
                self.init_open(&mut fig_item[(i - 1) as usize]);
                self.lines_length
                    .push(fig_item[(i - 1) as usize].len() as i32);
            } else {
                self.insert_corner(&mut fig_item[(i - 1) as usize]);
                self.lines_length
                    .push(fig_item[(i - 1) as usize].len() as i32);
            }
        }

        //---------------------------------
        // Une seule ligne => Rien a faire.  (OCCT L145-156)
        //---------------------------------
        if nb_lines == 1 {
            if open[0] {
                self.double_line(&mut fig_item[0], &s_vide, connexion_nul, self.direction);
                self.lines_length[0] = fig_item[0].len() as i32;
            }
            self.geom_elements = fig_item[0].clone();
            let len = self.geom_elements.len() as i32;
            self.update_link(1, 1, 1, len);
            self.lines_length.push(fig_item[0].len() as i32);
            return;
        }

        //------------------
        // Plusieurs lignes.  (OCCT L162-166)
        //------------------

        //---------------------------------------------------------
        // Calcul de l ensemble des connexions realisant le chemin.
        //---------------------------------------------------------
        let mut road = Mat2dMiniPath::new();
        road.perform(fig_item, ind_ref_line, trigo);

        //------------------------
        // Fermeture ligne ouverte.  (OCCT L171-194)
        //-------------------------
        for i in 1..=nb_lines {
            if open[(i - 1) as usize] {
                let cf: Option<HandleMat2dConnexion>;
                if road.is_root(i) {
                    cf = connexion_nul.clone();
                } else {
                    cf = Some(road.father(i));
                }
                if road.is_connexions_from(i) {
                    let from = road.connexions_from(i).clone();
                    self.double_line(&mut fig_item[(i - 1) as usize], &from, cf, self.direction);
                } else {
                    self.double_line(&mut fig_item[(i - 1) as usize], &s_vide, cf, self.direction);
                }
                self.lines_length[(i - 1) as usize] = fig_item[(i - 1) as usize].len() as i32;
            }
        }

        //------------------------
        // Construction du chemin.  (OCCT L199)
        //------------------------
        road.run_on_connexions();

        //-------------------------
        // Construction du Circuit.  (OCCT L218)
        //-------------------------
        self.construct_circuit(fig_item, ind_ref_line, &road);
    }

    //=======================================================================
    // function : IsSharpCorner (cxx L227-342)
    // purpose  : Return True si le point commun entre Geom1 et Geom2 est
    //            une cassure saillante par rapport a Direction.
    //=======================================================================
    fn is_sharp_corner(
        &self,
        geom1: &Geom2dGeometry,
        geom2: &Geom2dGeometry,
        direction: f64,
    ) -> bool {
        let mut dot_prod: f64 = 0.0;
        let mut pro_vec = cross_prod(geom1, geom2, &mut dot_prod);
        let mut nb_test: i32 = 1;
        let du: f64 = rcad_kernel::precision::CONFUSION;
        // OCCT: C1 = down_cast<Trimmed>(Geom1); C2 = down_cast<Trimmed>(Geom2).
        let c1 = geom1.curve();
        let c2 = geom2.curve();
        // Modified by Sergey KHROMOV - Thu Oct 24 19:02:46 2002:
        let mut tol_ang: f64 = 1.0e-8;

        if self.my_join_type == GeomAbsJoinType::Arc {
            while nb_test <= 10 {
                if pro_vec * direction < -tol_ang {
                    return true; // Saillant.
                }
                if pro_vec * direction > tol_ang {
                    return false; // Rentrant.
                } else {
                    if dot_prod > 0.0 {
                        return false; // Plat.
                    }
                    tol_ang = 1.0e-8;
                    let u1 = c1.last_parameter() - nb_test as f64 * du;
                    let u2 = c2.first_parameter() + nb_test as f64 * du;
                    // OCCT gp_Dir2d — normalized directions.
                    let dir1 = c1.dn(u1, 1).normalize_or_zero();
                    let dir2 = c2.dn(u2, 1).normalize_or_zero();
                    dot_prod = dir1.dot(dir2);
                    pro_vec = dir1.x * dir2.y - dir1.y * dir2.x;
                    nb_test += 1;
                }
            }

            // Rebroussement.  (OCCT L274-309)
            // on calcule des paralleles aux deux courbes du cote du domaine
            // de calcul.  Si pas d'intersection => saillant, sinon rentrant.
            let tol: f64 = rcad_kernel::precision::CONFUSION;
            let mil_c1 = (c1.last_parameter() + c1.first_parameter()) * 0.5;
            let mil_c2 = (c2.last_parameter() + c2.first_parameter()) * 0.5;
            let p = c1.value(c1.last_parameter());
            let p1 = c1.value(mil_c1);
            let p2 = c2.value(mil_c2);

            let mut d = p1.distance(p).min(p2.distance(p));
            d /= 10.0;

            if direction < 0.0 {
                d = -d;
            }

            // OCCT L295-300: HC1 = new Geom2dAdaptor_Curve(C1);
            // OC1 = Adaptor2d_OffsetCurve(HC1, D, MilC1, C1->LastParameter()).
            let oc1 = AdaptorOffsetCurve::new_bounded(c1.clone(), d, mil_c1, c1.last_parameter());
            let oc2 = AdaptorOffsetCurve::new_bounded(c2.clone(), d, c2.first_parameter(), mil_c2);
            let mut intersect = GInter::new();
            intersect.perform_cc(&oc1, &oc2, tol, tol);

            // OCCT: !Intersect.IsDone() || Intersect.IsEmpty() — the rcad
            // GInter delegates IsEmpty to the intersection base (no points
            // and no segments).
            let is_empty = intersect.nb_points() == 0 && intersect.nb_segments() == 0;
            return !intersect.is_done() || is_empty;
        } else if self.my_join_type == GeomAbsJoinType::Intersection {
            if pro_vec.abs() <= tol_ang && dot_prod < 0.0 {
                while nb_test <= 10 {
                    let u1 = c1.last_parameter() - nb_test as f64 * du;
                    let u2 = c2.first_parameter() + nb_test as f64 * du;
                    let dir1 = c1.dn(u1, 1).normalize_or_zero();
                    let dir2 = c2.dn(u2, 1).normalize_or_zero();
                    dot_prod = dir1.dot(dir2);
                    pro_vec = dir1.x * dir2.y - dir1.y * dir2.x;
                    if pro_vec * direction < -tol_ang {
                        return true; // Saillant.
                    }
                    if pro_vec * direction > tol_ang {
                        return false; // Rentrant.
                    }

                    nb_test += 1;
                }
                return false;
            } else {
                return false;
            }
        }
        false
    }

    /// OCCT MAT2d_Circuit::ConstructCircuit(FigItem, IndRefLine, Road)
    /// (cxx L360-480).
    fn construct_circuit(
        &mut self,
        fig_item: &Vec<Vec<Geom2dGeometry>>,
        ind_ref_line: i32,
        road: &Mat2dMiniPath,
    ) {
        let mut set_of_item: Vec<Geom2dGeometry> = Vec::new();
        let nb_connexions = road.path().len() as i32;

        //-----------------------------------------------------
        // Depart du premier element de la ligne de reference.  (OCCT L376-380)
        //-----------------------------------------------------
        let mut prev_c = road.path()[0].clone();
        sub_sequence(
            &fig_item[(ind_ref_line - 1) as usize],
            1,
            prev_c.read().unwrap().index_item_on_first(),
            &mut self.geom_elements,
        );
        self.update_link(
            1,
            ind_ref_line,
            1,
            prev_c.read().unwrap().index_item_on_first(),
        );
        self.connexion_map
            .insert(self.geom_elements.len() as i32 + 1, prev_c.clone());
        let mut i_last_item = self.geom_elements.len() as i32;

        //-----------------------------------------------------------------------
        // Ajout des portion de lignes delimites par deux connexions successives.
        //-----------------------------------------------------------------------
        for i in 2..=nb_connexions {
            let cur_c = road.path()[(i - 1) as usize].clone();
            if self.pass_by_last(&prev_c, &cur_c) {
                //------------------------------------------------------
                // La portion passe par le dernier element de la ligne.
                //------------------------------------------------------
                let cur = cur_c.read().unwrap();
                let ind_last = fig_item[(cur.index_first_line() - 1) as usize].len() as i32;
                sub_sequence(
                    &fig_item[(cur.index_first_line() - 1) as usize],
                    prev_c.read().unwrap().index_item_on_second(),
                    ind_last,
                    &mut set_of_item,
                );
                self.update_link(
                    i_last_item + 1,
                    cur.index_first_line(),
                    prev_c.read().unwrap().index_item_on_second(),
                    ind_last,
                );
                self.geom_elements.append(&mut set_of_item);
                i_last_item = self.geom_elements.len() as i32;

                if fig_item[(cur.index_first_line() - 1) as usize].len() > 1 {
                    sub_sequence(
                        &fig_item[(cur.index_first_line() - 1) as usize],
                        1,
                        cur.index_item_on_first(),
                        &mut set_of_item,
                    );
                    self.update_link(
                        i_last_item + 1,
                        cur.index_first_line(),
                        1,
                        cur.index_item_on_first(),
                    );
                    self.geom_elements.append(&mut set_of_item);
                    i_last_item = self.geom_elements.len() as i32;
                }
                self.connexion_map.insert(i_last_item + 1, cur_c.clone());
            } else {
                //------------------------------------------------------
                // La portion ne passe par le dernier element de la ligne.
                //------------------------------------------------------
                let cur = cur_c.read().unwrap();
                sub_sequence(
                    &fig_item[(cur.index_first_line() - 1) as usize],
                    prev_c.read().unwrap().index_item_on_second(),
                    cur.index_item_on_first(),
                    &mut set_of_item,
                );
                self.update_link(
                    i_last_item + 1,
                    cur.index_first_line(),
                    prev_c.read().unwrap().index_item_on_second(),
                    cur.index_item_on_first(),
                );
                self.geom_elements.append(&mut set_of_item);
                i_last_item = self.geom_elements.len() as i32;
                self.connexion_map.insert(i_last_item + 1, cur_c.clone());
            }
            // OCCT L433: PrevC = CurC (handle assignment).
            prev_c = cur_c.clone();
        }

        //-------------------------------------------------------------
        // Fermeture : de la derniere connexion au dernier element de la
        //             ligne de reference.  (OCCT L440-452)
        //-------------------------------------------------------------
        let ind_last = fig_item[(ind_ref_line - 1) as usize].len() as i32;
        if ind_last == 1 {
            self.connexion_map.insert(1, road.path()[(nb_connexions - 1) as usize].clone());
            self.connexion_map.remove(&(i_last_item + 1));
        } else {
            sub_sequence(
                &fig_item[(ind_ref_line - 1) as usize],
                prev_c.read().unwrap().index_item_on_second(),
                ind_last,
                &mut set_of_item,
            );
            self.update_link(
                i_last_item + 1,
                ind_ref_line,
                prev_c.read().unwrap().index_item_on_second(),
                ind_last,
            );
            self.geom_elements.append(&mut set_of_item);
        }

        //--------------------------------------
        // Tri des RefToEqui pour chaque element.  (OCCT L457-465)
        //--------------------------------------
        let keys: Vec<Mat2dBiInt> = self.link_ref_equi.keys().copied().collect();
        for key in keys {
            if self.link_ref_equi[&key].len() > 1 {
                self.sort_ref_to_equi(key);
            }
        }
    }

    /// OCCT MAT2d_Circuit::InitOpen(Line) const (cxx L484-503).
    fn init_open(&self, line: &mut Vec<Geom2dGeometry>) {
        // OCCT L489-492: insert the start point of the first curve and the
        // end point of the last curve.
        let curve = line.first().unwrap().curve();
        let start = curve.value(curve.first_parameter());
        line.insert(0, Geom2dGeometry::Point(start));
        let curve = line.last().unwrap().curve();
        let end = curve.value(curve.last_parameter());
        line.push(Geom2dGeometry::Point(end));

        // OCCT L494-502.
        let mut i: usize = 2;
        while i + 1 <= line.len() - 1 {
            let mut dot_prod: f64 = 0.0;
            let cross = cross_prod(&line[i - 1], &line[i], &mut dot_prod);
            if cross.abs() > 1.0e-8 || dot_prod < 0.0 {
                let curve = line[i - 1].curve();
                let end = curve.value(curve.last_parameter());
                line.insert(i, Geom2dGeometry::Point(end));
                i += 1;
            }
            i += 1;
        }
    }

    /// OCCT MAT2d_Circuit::DoubleLine(Line, ConnexionFrom, ConnexionFather,
    /// SideRef) const (cxx L507-670).
    fn double_line(
        &self,
        line: &mut Vec<Geom2dGeometry>,
        connexion_from: &Vec<HandleMat2dConnexion>,
        connexion_father: Option<HandleMat2dConnexion>,
        side_ref: f64,
    ) {
        let nb_items = line.len() as i32;

        //--------------------------
        // Completion de la ligne.  (OCCT L522-538)
        //--------------------------
        if !self.my_is_open_result {
            let mut i = nb_items - 1;
            while i > 1 {
                let value = line[(i - 1) as usize].clone();
                if value.is_point() {
                    line.push(value);
                } else {
                    // OCCT L533-535: Curve = down_cast<Trimmed>(Value(i)->Copy());
                    // Curve->Reverse(); Line.Append(Curve).
                    let mut curve = value.curve().clone();
                    trimmed_reverse(&mut curve);
                    line.push(Geom2dGeometry::Curve(curve));
                }
                i -= 1;
            }
        }

        //------------------------------------------
        // Repartition des connexions sur la ligne  (OCCT L543-578)
        //------------------------------------------
        // OCCT mutates the shared connexion handles inside ConnexionFrom; a
        // local mutable copy of the handle sequence models the OCCT
        // InsertAfter/Remove traffic on the sequence itself.
        let mut connexion_from = connexion_from.clone();
        let mut i_after = connexion_from.len() as i32;
        let nb_connexions = i_after;

        let mut i: i32 = 1;
        while i <= i_after {
            let cc = connexion_from[(i - 1) as usize].clone();
            let mut ind_cof = cc.read().unwrap().index_item_on_first();
            let type_is_point = line[(ind_cof - 1) as usize].is_point();

            if type_is_point {
                if ind_cof != nb_items && ind_cof != 1 {
                    let mut dot_prod: f64 = 0.0;
                    let pro_vec = cross_prod(
                        &line[(ind_cof - 2) as usize],
                        &line[ind_cof as usize],
                        &mut dot_prod,
                    );
                    if pro_vec * side_ref > 0.0 {
                        cc.write().unwrap().set_index_item_on_first(2 * nb_items - ind_cof);
                        connexion_from.insert(i_after as usize, cc.clone());
                        connexion_from.remove((i - 1) as usize);
                        i_after -= 1;
                        i -= 1;
                    }
                }
            } else if self.side(&cc, line) != side_ref {
                ind_cof = cc.read().unwrap().index_item_on_first();
                // OCCT L570-572: the reversed copy of the item.
                let curve = line[(ind_cof - 1) as usize].curve();
                let param = cc.read().unwrap().parameter_on_first();
                let reversed = trimmed_reversed_parameter(curve, param);
                let mut ccm = cc.write().unwrap();
                ccm.set_index_item_on_first(2 * nb_items - ind_cof);
                ccm.set_parameter_on_first(reversed);
                drop(ccm);
                connexion_from.insert(i_after as usize, cc.clone());
                connexion_from.remove((i - 1) as usize);
                i_after -= 1;
                i -= 1;
            }
            i += 1;
        }

        //---------------------------
        // Mise a jour connexion pere.  (OCCT L583-607)
        //---------------------------
        if let Some(father) = &connexion_father {
            let cc = father.read().unwrap().reverse();
            let mut ind_cof = cc.read().unwrap().index_item_on_first();
            let type_is_point = line[(ind_cof - 1) as usize].is_point();

            if type_is_point {
                if ind_cof != nb_items && ind_cof != 1 {
                    let mut dot_prod: f64 = 0.0;
                    let pro_vec = cross_prod(
                        &line[(ind_cof - 2) as usize],
                        &line[ind_cof as usize],
                        &mut dot_prod,
                    );
                    if pro_vec * side_ref > 0.0 {
                        father
                            .write()
                            .unwrap()
                            .set_index_item_on_second(2 * nb_items - ind_cof);
                    }
                }
            } else if self.side(&cc, line) != side_ref {
                ind_cof = cc.read().unwrap().index_item_on_first();
                let curve = line[(ind_cof - 1) as usize].curve();
                let param = father.read().unwrap().parameter_on_second();
                let reversed = trimmed_reversed_parameter(curve, param);
                let mut fm = father.write().unwrap();
                fm.set_index_item_on_second(2 * nb_items - ind_cof);
                fm.set_parameter_on_second(reversed);
            }
        }

        //-------------------------------------
        // Suppression des cassures rentrantes.  (OCCT L612-668)
        //-------------------------------------
        let mut ind_line: i32 = 1;
        let mut i_corres: i32 = 1;
        // OCCT NCollection_Array1<int> Corres(1, Line.Length()).
        let mut corres = vec![0i32; line.len()];

        // OCCT L616: while (Line.Value(IndLine) != Line.Last()) — the walk
        // stops on the last element; handles carry no identity here, so the
        // positional form is used.
        while (ind_line as usize) < line.len() {
            corres[(i_corres - 1) as usize] = ind_line;
            let type_is_point = line[(ind_line - 1) as usize].is_point();

            if type_is_point && i_corres != 1 && i_corres != nb_items {
                if !self.is_sharp_corner(
                    &line[(ind_line - 2) as usize],
                    &line[ind_line as usize],
                    side_ref,
                ) {
                    line.remove((ind_line - 1) as usize);
                    ind_line -= 1;
                    corres[(i_corres - 1) as usize] = 0;
                }
            }
            ind_line += 1;
            i_corres += 1;
        }
        corres[(i_corres - 1) as usize] = ind_line;

        if !self.my_is_open_result {
            for i in 1..2 * nb_items - 2 {
                if corres[(i - 1) as usize] == 0 {
                    corres[(i - 1) as usize] = corres[(2 * nb_items - i - 1) as usize];
                }
            }

            //----------------------------
            // Mise a jour des Connexions.  (OCCT L659-668)
            //----------------------------
            for i in 1..=nb_connexions {
                let cc = connexion_from[(i - 1) as usize].clone();
                let item = cc.read().unwrap().index_item_on_first();
                let mapped = corres[(item - 1) as usize];
                cc.write().unwrap().set_index_item_on_first(mapped);
            }

            if let Some(father) = &connexion_father {
                let item = father.read().unwrap().index_item_on_second();
                let mapped = corres[(item - 1) as usize];
                father.write().unwrap().set_index_item_on_second(mapped);
            }
        }
    }

    /// OCCT MAT2d_Circuit::InsertCorner(Line) const (cxx L674-702).
    fn insert_corner(&self, line: &mut Vec<Geom2dGeometry>) {
        let mut i: i32 = 1;
        while i <= line.len() as i32 {
            let isuiv = if i == line.len() as i32 { 1 } else { i + 1 };
            let insert = self.is_sharp_corner(
                &line[(i - 1) as usize],
                &line[(isuiv - 1) as usize],
                self.direction,
            );

            if insert {
                let curve = line[(isuiv - 1) as usize].curve();
                let start = curve.value(curve.first_parameter());
                line.insert(i as usize, Geom2dGeometry::Point(start));
                // OCCT L699: i++ (in addition to the loop increment).
                i += 1;
            }
            i += 1;
        }
    }

    /// OCCT MAT2d_Circuit::NumberOfItems() const (cxx L706-709).
    pub fn number_of_items(&self) -> i32 {
        self.geom_elements.len() as i32
    }

    /// OCCT MAT2d_Circuit::LineLength(I) const (cxx L713-716).
    pub fn line_length(&self, i: i32) -> i32 {
        self.lines_length[(i - 1) as usize]
    }

    /// OCCT MAT2d_Circuit::Value(Index) const (cxx L720-723) — the item at
    /// position Index.
    pub fn value(&self, index: i32) -> &Geom2dGeometry {
        &self.geom_elements[(index - 1) as usize]
    }

    /// OCCT MAT2d_Circuit::RefToEqui(IndLine, IndCurve) const (cxx L727-732)
    /// — the set of index of the items corresponding to the curve IndCurve
    /// on the line IndLine from the initial figure.
    pub fn ref_to_equi(&self, ind_line: i32, ind_curve: i32) -> &Vec<i32> {
        let key = Mat2dBiInt::new(ind_line, ind_curve);
        self.link_ref_equi
            .get(&key)
            .expect("MAT2d_Circuit::RefToEqui")
    }

    /// OCCT MAT2d_Circuit::SortRefToEqui(BiRef) (cxx L736-755).
    fn sort_ref_to_equi(&mut self, bi_ref: Mat2dBiInt) {
        let s = self.link_ref_equi[&bi_ref].clone();

        let mut i: usize = 1;
        while i <= s.len() {
            if !self.connexion_on(s[i - 1]) {
                break;
            }
            i += 1;
        }
        // OCCT L749-754: if (i > 1 && i <= S.Length()) { SFin = S;
        // SFin.Split(i, S); S.Append(SFin); }  — the first i-1 entries are
        // moved to the back (left rotation by i-1).
        if i > 1 && i <= s.len() {
            let s_ref = self.link_ref_equi.get_mut(&bi_ref).unwrap();
            let head: Vec<i32> = s_ref.drain(..i - 1).collect();
            s_ref.extend(head);
        }
    }

    /// OCCT MAT2d_Circuit::Connexion(I) const (cxx L759-762).
    pub fn connexion(&self, i: i32) -> HandleMat2dConnexion {
        self.connexion_map[&i].clone()
    }

    /// OCCT MAT2d_Circuit::ConnexionOn(I) const (cxx L766-769).
    pub fn connexion_on(&self, i: i32) -> bool {
        self.connexion_map.contains_key(&i)
    }

    /// OCCT MAT2d_Circuit::Side(C1, Line) const (cxx L773-790).
    fn side(&self, c1: &HandleMat2dConnexion, line: &Vec<Geom2dGeometry>) -> f64 {
        let c = c1.read().unwrap();
        // OCCT L778-779: gp_Vec2d Vect1(PointOnSecond() - PointOnFirst()).
        let vect1 = DVec2::new(
            c.point_on_second().x - c.point_on_first().x,
            c.point_on_second().y - c.point_on_first().y,
        );
        let curve = line[(c.index_item_on_first() - 1) as usize].curve();
        let vect2 = curve.dn(c.parameter_on_first(), 1);
        let cross = vect1.x * vect2.y - vect1.y * vect2.x;
        if cross > 0.0 {
            -1.0
        } else {
            1.0
        }
    }

    /// OCCT MAT2d_Circuit::PassByLast(C1, C2) const (cxx L794-825).
    fn pass_by_last(&self, c1: &HandleMat2dConnexion, c2: &HandleMat2dConnexion) -> bool {
        let c1 = c1.read().unwrap();
        let c2 = c2.read().unwrap();

        if c2.index_first_line() == c1.index_second_line() {
            if c2.index_item_on_first() < c1.index_item_on_second() {
                return true;
            } else if c2.index_item_on_first() == c1.index_item_on_second() {
                if c1.index_first_line() == c2.index_second_line() {
                    return true;
                }
                if c2.parameter_on_first() == c1.parameter_on_second() {
                    // OCCT L811-816.
                    let vect1 = DVec2::new(
                        c1.point_on_second().x - c1.point_on_first().x,
                        c1.point_on_second().y - c1.point_on_first().y,
                    );
                    let vect2 = DVec2::new(
                        c2.point_on_first().x - c2.point_on_second().x,
                        c2.point_on_first().y - c2.point_on_second().y,
                    );
                    let cross = vect1.x * vect2.y - vect1.y * vect2.x;
                    if cross * self.direction > 0.0 {
                        return true;
                    }
                } else if c2.parameter_on_first() < c1.parameter_on_second() {
                    return true;
                }
            }
        }
        false
    }

    /// OCCT MAT2d_Circuit::UpDateLink(IFirst, ILine, ICurveFirst, ICurveLast)
    /// (cxx L829-852).
    fn update_link(&mut self, i_first: i32, i_line: i32, i_curve_first: i32, i_curve_last: i32) {
        let mut i_equi = i_first;

        for i in i_curve_first..=i_curve_last {
            let key = Mat2dBiInt::new(i_line, i);
            let entry = self.link_ref_equi.entry(key).or_default();
            entry.push(i_equi);
            i_equi += 1;
        }
    }
}
