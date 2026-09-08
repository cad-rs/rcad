//! OCCT MAT2d_MiniPath — computes a path to link all the lines in a set of
//! lines; the path is described as a set of connexions.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT2d/
//!         MAT2d_MiniPath.hxx L17-115, MAT2d_MiniPath.cxx L33-438

use std::collections::HashMap;
use std::sync::{Arc, RwLock};

use glam::DVec2;
use rcad_kernel::base::extrema::{ExtPC2d, POnCurve2d};
use rcad_kernel::precision::CONFUSION;

use crate::geomalgo::geom2d_int::Curve2dAdaptor;
use crate::topalgo::mat2d::mat2d_connexion::{HandleMat2dConnexion, Mat2dConnexion};
use crate::topalgo::mat2d::Geom2dGeometry;

/// OCCT RealLast().
const REAL_LAST: f64 = f64::MAX;

/// GAP carrier (dependency of another package).
///
/// OCCT dependency: `Extrema_ExtCC2d` (TKGeomBase/Extrema,
/// Extrema_ExtCC2d.hxx/.cxx) — the extrema of distance between two 2d
/// curves.  The 2d curve/curve extrema are not translated yet (the kernel
/// only carries the 3d `Extrema_ExtCC`), so this carrier reproduces the
/// public interface consumed by MAT2d_MiniPath::MinimumL1L2 (cxx L406-423:
/// `Extrema_ExtCC2d Extremas(C1, C2)`, `IsParallel()`, `IsDone()`,
/// `NbExt()`, `SquareDistance(i)`, `Points(i, P1, P2)`) with the OCCT
/// *failure path*: `IsDone()` returns false, so the guard
/// `if (!Extremas.IsParallel() && Extremas.IsDone())` (cxx L407) skips the
/// block exactly as OCCT does when the computation fails.  When the
/// Extrema 2d port lands, its implementation replaces the body of
/// [`ExtCC2d::new`].
#[derive(Debug)]
pub struct ExtCC2d {
    done: bool,
    parallel: bool,
    points: Vec<(POnCurve2d, POnCurve2d)>,
    sq_distances: Vec<f64>,
}

impl ExtCC2d {
    /// OCCT Extrema_ExtCC2d(C1, C2) — GAP: see type docs.
    pub fn new(_c1: &rcad_kernel::geom::Curve2d, _c2: &rcad_kernel::geom::Curve2d) -> Self {
        ExtCC2d {
            done: false,
            parallel: false,
            points: Vec::new(),
            sq_distances: Vec::new(),
        }
    }

    /// OCCT Extrema_ExtCC2d::IsParallel().
    pub fn is_parallel(&self) -> bool {
        self.parallel
    }

    /// OCCT Extrema_ExtCC2d::IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Extrema_ExtCC2d::NbExt().
    pub fn nb_ext(&self) -> usize {
        self.sq_distances.len()
    }

    /// OCCT Extrema_ExtCC2d::SquareDistance(N) (1-based).
    pub fn square_distance(&self, n: usize) -> f64 {
        self.sq_distances[n - 1]
    }

    /// OCCT Extrema_ExtCC2d::Points(N, P1, P2) (1-based output form).
    pub fn points(&self, n: usize) -> (POnCurve2d, POnCurve2d) {
        self.points[n - 1].clone()
    }
}

/// OCCT MAT2d_MiniPath (MAT2d_MiniPath.hxx L44-114).
///
/// The set of connexions can be seen as an arbitrary Tree; the nodes of the
/// tree are the lines, the arcs are the connexions.  The children of a line
/// are ordered by the relation IsAfter defined on the connexions.
#[derive(Debug)]
pub struct Mat2dMiniPath {
    /// OCCT NCollection_DataMap<int, NCollection_Sequence<handle>>
    /// theConnexions.
    the_connexions: HashMap<i32, Vec<HandleMat2dConnexion>>,
    /// OCCT NCollection_DataMap<int, handle> theFather.
    the_father: HashMap<i32, HandleMat2dConnexion>,
    /// OCCT NCollection_Sequence<handle> thePath.
    the_path: Vec<HandleMat2dConnexion>,
    /// OCCT double theDirection.
    the_direction: f64,
    /// OCCT int indStart.
    ind_start: i32,
}

impl Default for Mat2dMiniPath {
    fn default() -> Self {
        Self::new()
    }
}

impl Mat2dMiniPath {
    /// OCCT MAT2d_MiniPath::MAT2d_MiniPath() (cxx L33-37).
    pub fn new() -> Self {
        Mat2dMiniPath {
            the_connexions: HashMap::new(),
            the_father: HashMap::new(),
            the_path: Vec::new(),
            the_direction: 1.0,
            ind_start: 0,
        }
    }

    /// OCCT MAT2d_MiniPath::Perform(Figure, IndStart, Sense) (cxx L45-130).
    ///
    /// Computes the path to link the lines in `figure`; the path starts on
    /// the line of index `ind_start`.  `sense` = true if the Circuit turns
    /// in the trigonometric sense.
    pub fn perform(&mut self, figure: &[Vec<Geom2dGeometry>], ind_start: i32, sense: bool) {
        let nb_lines = figure.len() as i32;
        // OCCT NCollection_Array2<handle> Connexion(1, NbLines, 1, NbLines)
        // — stored 0-based here, OCCT (i, j) -> [i-1][j-1].
        let mut connexion: Vec<Vec<Option<HandleMat2dConnexion>>> =
            vec![vec![None; nb_lines as usize]; nb_lines as usize];

        self.ind_start = ind_start;
        self.the_direction = 1.0;
        if sense {
            self.the_direction = -1.0;
        }

        //----------------------------------------------------------------------
        // Calcul des connexions qui realisent le minimum de distance entre les
        // differents elements de la figure.  (cxx L66-73)
        //----------------------------------------------------------------------
        for i in 1..nb_lines {
            for j in (i + 1)..=nb_lines {
                connexion[(i - 1) as usize][(j - 1) as usize] =
                    Some(self.minimum_l1_l2(figure, i, j));
                let rev = connexion[(i - 1) as usize][(j - 1) as usize]
                    .as_ref()
                    .unwrap()
                    .read()
                    .unwrap()
                    .reverse();
                connexion[(j - 1) as usize][(i - 1) as usize] = Some(rev);
            }
        }

        // OCCT L75-79: NCollection_Sequence<int> Set1, Set2; ...
        let mut set1: Vec<i32> = Vec::new();
        let mut set2: Vec<i32> = Vec::new();
        let mut dist_s1_s2: f64;
        let mut indice_line1: i32;
        let mut indice_line2: i32;
        let mut i_suiv: usize = 0;
        let mut min_on_set1: i32 = 0;
        let mut min_on_set2: i32 = 0;

        //---------------------------------------------------------------------------
        // - 0 Set1 est initialise avec la ligne de depart.
        //     Set2 contient toutes les autres.  (cxx L85-93)
        //---------------------------------------------------------------------------
        set1.push(ind_start);

        for i in 1..=nb_lines {
            if i != ind_start {
                set2.push(i);
            }
        }

        //---------------------------------------------------------------------------
        // - 1 Recherche de la connexion C la plus courte entre Set1 et Set2.
        // - 2 La ligne de Set2 realisant le minimum de distance est inseree dans
        //     Set1 et supprime dans Set2.
        // - 3 Insertion de la connexion dans l ensemble des connexions.
        // - 4 Si Set2 est non vide retour en 1.  (cxx L103-124)
        //---------------------------------------------------------------------------
        while !set2.is_empty() {
            dist_s1_s2 = REAL_LAST;
            for i in 1..=set1.len() {
                indice_line1 = set1[i - 1];
                for j in 1..=set2.len() {
                    indice_line2 = set2[j - 1];
                    let cij = connexion[(indice_line1 - 1) as usize][(indice_line2 - 1) as usize]
                        .as_ref()
                        .unwrap();
                    if cij.read().unwrap().distance() < dist_s1_s2 {
                        i_suiv = j;
                        dist_s1_s2 = cij.read().unwrap().distance();
                        min_on_set1 = indice_line1;
                        min_on_set2 = indice_line2;
                    }
                }
            }
            set1.push(set2[i_suiv - 1]);
            set2.remove(i_suiv - 1);
            let c = connexion[(min_on_set1 - 1) as usize][(min_on_set2 - 1) as usize]
                .clone()
                .unwrap();
            self.append(&c);
        }

        //----------------------------------------------------------------
        // Construction du chemin en parcourant l ensemble des connexions.
        // (cxx L129)
        //----------------------------------------------------------------
        self.run_on_connexions();
    }

    /// OCCT MAT2d_MiniPath::Append(C) (cxx L141-181) — insertion d une
    /// nouvelle connexion dans le chemin.
    ///
    /// Les connexions et les lignes constituent un arbre dont
    /// - les noeuds sont les lignes,
    /// - les connexions sont les branches.
    fn append(&mut self, c: &HandleMat2dConnexion) {
        let (index_first_line, index_second_line) = {
            let cc = c.read().unwrap();
            (cc.index_first_line(), cc.index_second_line())
        };

        if !self.the_connexions.contains_key(&index_first_line) {
            self.the_connexions
                .entry(index_first_line)
                .or_default()
                .push(c.clone());
            self.the_father.insert(index_second_line, c.clone());
            return;
        }

        let seq = self.the_connexions.get_mut(&index_first_line).unwrap();
        let mut index_after: usize = 0;
        let nb_connexions = seq.len();

        for i in 1..=nb_connexions {
            let cc = seq[i - 1].clone();
            if cc.read().unwrap().is_after(c, self.the_direction) {
                index_after = i;
                break;
            }
        }
        //----------------------------------------------------------------------
        // Insertion de <C> avant <IAfter>.
        // Si <IAfter> = 0 => Pas de connexions apres <C> => <C> est la
        // derniere.
        //----------------------------------------------------------------------
        if index_after == 0 {
            seq.push(c.clone());
        } else {
            seq.insert(index_after - 1, c.clone());
        }
        self.the_father.insert(index_second_line, c.clone());
    }

    /// OCCT MAT2d_MiniPath::Path() const (cxx L187-190) — the sequence of
    /// connexions corresponding to the path.
    pub fn path(&self) -> &Vec<HandleMat2dConnexion> {
        &self.the_path
    }

    /// OCCT MAT2d_MiniPath::IsConnexionsFrom(i) const (cxx L194-197) — true
    /// if there is one connexion which starts on line `i`.
    pub fn is_connexions_from(&self, i: i32) -> bool {
        self.the_connexions.contains_key(&i)
    }

    /// OCCT MAT2d_MiniPath::ConnexionsFrom(i) (cxx L203-206) — the
    /// connexions which start on line `i`.
    pub fn connexions_from(&mut self, i: i32) -> &mut Vec<HandleMat2dConnexion> {
        self.the_connexions.get_mut(&i).unwrap()
    }

    /// OCCT MAT2d_MiniPath::IsRoot(ILine) const (cxx L210-213) — true if the
    /// line `i_line` is the root.
    pub fn is_root(&self, i_line: i32) -> bool {
        i_line == self.ind_start
    }

    /// OCCT MAT2d_MiniPath::Father(ILine) (cxx L219-222) — the first
    /// connexion which arrives on line `i_line`.
    pub fn father(&mut self, i_line: i32) -> HandleMat2dConnexion {
        self.the_father.get(&i_line).unwrap().clone()
    }

    /// OCCT MAT2d_MiniPath::RunOnConnexions() (cxx L228-243) — construction
    /// de <thePath> en parcourant <theConnexions>.
    pub fn run_on_connexions(&mut self) {
        // Architecture note: OCCT iterates a const reference to
        // theConnexions(indStart) while appending to thePath (a distinct
        // member); the Vec clone is observationally identical (the tree is
        // not modified during the walk).
        let sc = self.the_connexions[&self.ind_start].clone();

        self.the_path.clear();

        for i in 1..=sc.len() {
            let c = sc[i - 1].clone();
            self.the_path.push(c.clone());
            Self::explo_sons(&self.the_connexions, self.the_direction, &mut self.the_path, &c);
            let rev = c.read().unwrap().reverse();
            self.the_path.push(rev);
        }
    }

    /// OCCT MAT2d_MiniPath::ExploSons(CResult, CRef) (cxx L247-287).
    fn explo_sons(
        connexions: &HashMap<i32, Vec<HandleMat2dConnexion>>,
        the_direction: f64,
        c_result: &mut Vec<HandleMat2dConnexion>,
        c_ref: &HandleMat2dConnexion,
    ) {
        let index = c_ref.read().unwrap().index_second_line();

        if !connexions.contains_key(&index) {
            return;
        }

        // Architecture note: Vec clone as in run_on_connexions.
        let sc = connexions[&index].clone();
        let crr = c_ref.read().unwrap().reverse();

        for i in 1..=sc.len() {
            let c = sc[i - 1].clone();
            if c.read().unwrap().is_after(&crr, the_direction) {
                c_result.push(c.clone());
                Self::explo_sons(connexions, the_direction, c_result, &c);
                let rev = c.read().unwrap().reverse();
                c_result.push(rev);
            }
        }

        for i in 1..=sc.len() {
            let c = sc[i - 1].clone();
            if !c.read().unwrap().is_after(&crr, the_direction) {
                c_result.push(c.clone());
                Self::explo_sons(connexions, the_direction, c_result, &c);
                let rev = c.read().unwrap().reverse();
                c_result.push(rev);
            } else {
                break;
            }
        }
    }

    /// OCCT MAT2d_MiniPath::MinimumL1L2(Figure, IL1, IL2) const
    /// (cxx L294-438) — the connexion realising the minimum of distance
    /// between the lines of index IL1 and IL2 in `figure`, oriented from
    /// IL1 to IL2.
    fn minimum_l1_l2(
        &self,
        figure: &[Vec<Geom2dGeometry>],
        il1: i32,
        il2: i32,
    ) -> HandleMat2dConnexion {
        let mut point_on_curv1 = POnCurve2d {
            param: 0.0,
            point: DVec2::ZERO,
        };
        let mut point_on_curv2 = POnCurve2d {
            param: 0.0,
            point: DVec2::ZERO,
        };
        let mut i_min_c1: i32 = 0;
        let mut i_min_c2: i32 = 0;
        let mut parameter_on_c1: f64 = 0.0;
        let mut parameter_on_c2: f64 = 0.0;
        let mut point1 = DVec2::ZERO;
        let mut point2 = DVec2::ZERO;
        let mut p1 = DVec2::ZERO;
        let mut p2 = DVec2::ZERO;
        // OCCT occ::handle<Geom2d_Curve> Item1, Item2 — the down-cast
        // curves; None models a null handle (the item is a point).
        let mut item1: Option<&rcad_kernel::geom::Curve2d> = None;
        let mut item2: Option<&rcad_kernel::geom::Curve2d> = None;

        // OCCT L1 = Figure.Value(IL1); L2 = Figure.Value(IL2).
        let l1 = &figure[(il1 - 1) as usize];
        let l2 = &figure[(il2 - 1) as usize];

        let mut dist_l1_l2_2 = REAL_LAST;

        //---------------------------------------------------------------------------
        // Calcul des extremas de distances entre les composants de L1 et de L2.
        //---------------------------------------------------------------------------

        for ic1 in 1..=l1.len() {
            // OCCT L320-328.
            let g1 = &l1[ic1 - 1];
            if !g1.is_point() {
                item1 = Some(g1.curve());
            } else {
                p1 = g1.pnt2d();
            }

            for ic2 in 1..=l2.len() {
                // OCCT L333-341.
                let g2 = &l2[ic2 - 1];
                if !g2.is_point() {
                    item2 = Some(g2.curve());
                } else {
                    p2 = g2.pnt2d();
                }

                // OCCT L343-357: both items are points.
                if g1.is_point() && g2.is_point() {
                    let dist_p1_p2_2 = p1.distance_squared(p2);
                    if dist_p1_p2_2 <= dist_l1_l2_2 {
                        dist_l1_l2_2 = dist_p1_p2_2;
                        i_min_c1 = ic1 as i32;
                        i_min_c2 = ic2 as i32;
                        point1 = p1;
                        point2 = p2;
                        parameter_on_c1 = 0.0;
                        parameter_on_c2 = 0.0;
                    }
                }
                // OCCT L358-379: Type1 is a point.
                else if g1.is_point() {
                    let c2 = item2.unwrap();
                    let mut extremas = ExtPC2d::new(
                        p1,
                        c2,
                        CONFUSION,
                        c2.first_parameter(),
                        c2.last_parameter(),
                    );
                    if extremas.is_done() {
                        for i in 1..=extremas.nb_ext() {
                            if extremas.square_distance(i) < dist_l1_l2_2 {
                                dist_l1_l2_2 = extremas.square_distance(i);
                                i_min_c1 = ic1 as i32;
                                i_min_c2 = ic2 as i32;
                                let pon = extremas.point(i);
                                point_on_curv2 = POnCurve2d {
                                    param: pon.param,
                                    point: pon.point,
                                };
                                parameter_on_c1 = 0.0;
                                parameter_on_c2 = point_on_curv2.param;
                                point1 = p1;
                                point2 = point_on_curv2.point;
                            }
                        }
                    }
                }
                // OCCT L380-401: Type2 is a point.
                else if g2.is_point() {
                    let c1 = item1.unwrap();
                    let mut extremas = ExtPC2d::new(
                        p2,
                        c1,
                        CONFUSION,
                        c1.first_parameter(),
                        c1.last_parameter(),
                    );
                    if extremas.is_done() {
                        for i in 1..=extremas.nb_ext() {
                            if extremas.square_distance(i) < dist_l1_l2_2 {
                                dist_l1_l2_2 = extremas.square_distance(i);
                                i_min_c1 = ic1 as i32;
                                i_min_c2 = ic2 as i32;
                                let pon = extremas.point(i);
                                point_on_curv1 = POnCurve2d {
                                    param: pon.param,
                                    point: pon.point,
                                };
                                parameter_on_c2 = 0.0;
                                parameter_on_c1 = point_on_curv1.param;
                                point1 = point_on_curv1.point;
                                point2 = p2;
                            }
                        }
                    }
                }
                // OCCT L402-424: two curves.
                else {
                    let c1 = item1.unwrap();
                    let c2 = item2.unwrap();
                    let extremas = ExtCC2d::new(c1, c2);
                    if !extremas.is_parallel() && extremas.is_done() {
                        for i in 1..=extremas.nb_ext() {
                            if extremas.square_distance(i) < dist_l1_l2_2 {
                                dist_l1_l2_2 = extremas.square_distance(i);
                                i_min_c1 = ic1 as i32;
                                i_min_c2 = ic2 as i32;
                                let (poc1, poc2) = extremas.points(i);
                                point_on_curv1 = poc1;
                                point_on_curv2 = poc2;
                                parameter_on_c1 = point_on_curv1.param;
                                parameter_on_c2 = point_on_curv2.param;
                                point1 = point_on_curv1.point;
                                point2 = point_on_curv2.point;
                            }
                        }
                    }
                }
            }
        }

        // OCCT L427-437.
        let connexion_l1_l2 = Mat2dConnexion::new_full(
            il1,
            il2,
            i_min_c1,
            i_min_c2,
            dist_l1_l2_2.sqrt(),
            parameter_on_c1,
            parameter_on_c2,
            point1,
            point2,
        );
        Arc::new(RwLock::new(connexion_l1_l2))
    }
}
