//! OCCT MAT2d_Tool2d — set of the methods useful for the MAT's computation;
//! Tool2d contains the geometry of the bisecting locus.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT2d/
//!         MAT2d_Tool2d.hxx L17-177, MAT2d_Tool2d.cxx L17-1481

use std::collections::HashMap;
use std::sync::{Arc, RwLock};

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Line2d, TrimmedCurve2};

use crate::fillet::chfi3d_builder_0::elclib_adjust_periodic;
use crate::geomalgo::geom2d_int::{curve2d_type_of, Curve2dAdaptor, Curve2dType, GInter};
use crate::geomalgo::int_res2d::Domain;
use crate::topalgo::bisector::bisector::GeomAbsJoinType;
use crate::topalgo::bisector::bisector_bisec::{BisectorBisec, HandleGeom2dTrimmedCurve};
use crate::topalgo::bisector::bisector_bisec_ana;
use crate::topalgo::bisector::bisector_bisec_cc;
use crate::topalgo::bisector::bisector_curve::{
    BisectorCurve, CurveKind, Geom2dCurveHandle, Geom2dCurveAdaptor,
};
use crate::topalgo::bisector::bisector_inter::BisectorInter;
use crate::topalgo::mat::HandleMatBisector;
use crate::topalgo::mat2d::{Geom2dGeometry, Mat2dCircuit, MatSide};

// ---------------------------------------------------------------------------
// Pre-written cross-package surface (Bisector package, signatures per .hxx,
// plan §0.6 — reconciled at joint compile time):
//   - BisectorBisec (Bisector_Bisec.hxx L17-94): value struct with the four
//     Perform overloads; Geom2d_Point arguments map to DVec2; Value() and
//     ChangeValue() return HandleGeom2dTrimmedCurve = Arc<RwLock<TrimmedCurve>>.
//   - BisectorInter (Bisector_Inter.hxx L17-47): new() + Perform(C1, D1, C2,
//     D2, TolConf, Tol, ComunElement) + the IntRes2d_Intersection surface.
//   - bisector_bisec_ana::geom2d_curve_of(&Arc<dyn BisectorCurve>) ->
//     Option<Curve2d>  — down_cast<Bisector_BisecAna>()->Geom2dCurve()
//     (MAT2d_Tool2d.cxx L1342).
//   - bisector_bisec_ana::set_trim_of(&Arc<dyn BisectorCurve>, uf, ul) —
//     down_cast<Bisector_BisecAna>()->SetTrim(uf, ul) (cxx L1245).
//   - bisector_bisec_cc::cast_from(Arc<dyn BisectorCurve>) ->
//     Option<Arc<BisectorBisecCC>> — down_cast<Bisector_BisecCC>
//     (cxx L1213); Bisector_BisecCC::Curve(i)/FirstParameter()/Parameter(P).
// ---------------------------------------------------------------------------

/// OCCT static double MAT2d_TOLCONF = 1.e-7 (cxx L79).
const MAT2D_TOLCONF: f64 = 1.0e-7;

/// OCCT Precision::Infinite().
const INFINITE: f64 = rcad_kernel::precision::INFINITE_VALUE;

/// OCCT static bool AreNeighbours(IEdge1, IEdge2, NbEdge) (cxx L71,
/// body L1271-1285) — TRUE if IEdge1 and IEdge2 are consecutive elements on
/// a closed contour of NbEdge elements.
fn are_neighbours(i_edge1: i32, i_edge2: i32, nb_edge: i32) -> bool {
    if (i_edge1 - i_edge2).abs() == 1 {
        true
    } else {
        (i_edge1 - i_edge2).abs() == nb_edge - 1
    }
}

/// OCCT static IntRes2d_Domain Domain(Bisector1, Tolerance) (cxx L66-67,
/// body L1331-1381).
fn domain(bisector1: &HandleGeom2dTrimmedCurve, tolerance: f64) -> Domain {
    let bis = bisector1.read().unwrap();
    let param1 = bis.first_parameter();
    let mut param2 = bis.last_parameter();
    if param2 > 10000.0 {
        param2 = 10000.0;
        // OCCT L1338-1344: Type1 = Type(Bisector1->BasisCurve())
        // (the trimmed unwrap of Type()); if Bisector_BisecAna:
        // BasisCurve = ->Geom2dCurve(); Type1 = its type.
        let mut basis_curve: Option<Curve2d> = None;
        if bis.basis().kind() == CurveKind::BisecAna {
            basis_curve = bisector_bisec_ana::geom2d_curve_of(bis.basis());
        }
        let limit: f64 = 50000.0;
        if let Some(bc) = &basis_curve {
            match bc {
                // OCCT L1349-1356 (Geom2d_Parabola).
                Curve2d::Parabola(parabola) => {
                    let focus = parabola.focal_param;
                    let val1 = (limit * focus).sqrt();
                    let val2 = (limit * limit).sqrt();
                    param2 = if val1 <= val2 { val1 } else { val2 };
                }
                // OCCT L1357-1367 (Geom2d_Hyperbola).
                Curve2d::Hyperbola(hyperbola) => {
                    let majr = hyperbola.semi_major;
                    let minr = hyperbola.semi_minor;
                    let valu1 = limit / majr;
                    let valu2 = limit / minr;
                    let val1 = (valu1 + (valu1 * valu1 - 1.0).sqrt()).ln();
                    let val2 = (valu2 + (valu2 * valu2 + 1.0).sqrt()).ln();
                    param2 = if val1 <= val2 { val1 } else { val2 };
                }
                _ => {}
            }
        }
    }

    // OCCT L1370-1375.
    let mut domain1 = Domain::bounded(
        bis.value(param1),
        param1,
        tolerance,
        bis.value(param2),
        param2,
        tolerance,
    );
    // OCCT L1376-1379.
    if bis.basis().is_periodic() {
        domain1.set_equivalent_parameters(0.0, 2.0 * std::f64::consts::PI);
    }
    domain1
}

/// OCCT static void SetTrim(Bisector_Bisec& Bis, handle<Geom2d_Curve> Line1)
/// (cxx L73, body L1289-1327) — restriction of the bisector by the
/// intersection point of smallest parameter.
fn set_trim(bis: &mut BisectorBisec, line1: &Curve2d) {
    let mut intersect = GInter::new();
    let tolerance = MAT2D_TOLCONF;
    // OCCT L1294: handle<Geom2d_TrimmedCurve> Bisector = Bis.ChangeValue().
    let bisector = bis.change_value();

    // OCCT L1296: IntRes2d_Domain Domain1 = Domain(Bisector, Tolerance).
    let domain1 = domain(&bisector, tolerance);

    let (ub1, ub2, first_point_bisector, bis_arc) = {
        let trimmed = bisector.read().unwrap();
        (
            trimmed.first_parameter(),
            trimmed.last_parameter(),
            trimmed.value(trimmed.first_parameter()),
            Arc::new(trimmed.clone()) as Arc<dyn BisectorCurve>,
        )
    };
    let mut u_trim = INFINITE;

    // OCCT L1303-1305: Geom2dAdaptor_Curve AdapBisector(Bisector);
    // AdapLine1(Line1); Intersect.Perform(AdapBisector, Domain1, AdapLine1,
    // Tolerance, Tolerance).  Architecture note: the trimmed bisector is
    // cloned into the adaptor (read-only use).
    let adap_bisector = Geom2dCurveAdaptor::new(bis_arc);
    intersect.perform_cd_c(&adap_bisector, &domain1, line1, tolerance, tolerance);

    // OCCT L1307-1318.
    if intersect.is_done() && !ginter_is_empty(&intersect) {
        for i in 1..=intersect.nb_points() {
            let p_int: DVec2 = intersect.point(i).value();
            let distance = first_point_bisector.distance(p_int);
            if distance > 10.0 * tolerance && intersect.point(i).param_on_first() < u_trim {
                u_trim = intersect.point(i).param_on_first();
            }
        }
    }
    // OCCT L1319-1326: restriction of the bisector by the intersection point
    // of smallest parameter.
    if u_trim < ub2 && u_trim > ub1 {
        bisector.write().unwrap().set_trim(ub1, u_trim);
    }
}

/// OCCT IntRes2d_Intersection::IsEmpty() on a GInter result — no points and
/// no segments (the rcad GInter type does not carry the base accessor).
fn ginter_is_empty(intersect: &GInter) -> bool {
    intersect.nb_points() == 0 && intersect.nb_segments() == 0
}

/// OCCT static bool CheckEnds(Elt, PCom, Distance, Tol) (cxx L74-77,
/// body L1385-1415) — true if one of the ends of the element is at
/// `distance` from PCom (within Tol).
fn check_ends(elt: &Geom2dGeometry, p_com: DVec2, distance: f64, tol: f64) -> bool {
    match elt {
        // OCCT L1394-1397.
        Geom2dGeometry::Point(_) => false,
        Geom2dGeometry::Curve(curve) => {
            // OCCT L1400-1404: Curve = down_cast<Trimmed>(Elt);
            // aPf = StartPoint(); aPl = EndPoint().
            let a_pf = curve.value(curve.first_parameter());
            let a_pl = curve.value(curve.last_parameter());
            let df = p_com.distance(a_pf);
            let dl = p_com.distance(a_pl);
            if (df - distance).abs() <= tol {
                return true;
            }
            if (dl - distance).abs() <= tol {
                return true;
            }
            false
        }
    }
}

/// OCCT MAT2d_Tool2d (MAT2d_Tool2d.hxx L41-176).
pub struct Mat2dTool2d {
    /// OCCT double theDirection.
    the_direction: f64,
    /// OCCT GeomAbs_JoinType theJoinType.
    the_join_type: GeomAbsJoinType,
    /// OCCT int theNumberOfBisectors.
    the_number_of_bisectors: i32,
    /// OCCT int theNumberOfPnts.
    the_number_of_pnts: i32,
    /// OCCT int theNumberOfVecs.
    the_number_of_vecs: i32,
    /// OCCT occ::handle<MAT2d_Circuit> theCircuit.
    the_circuit: Option<Arc<RwLock<Mat2dCircuit>>>,
    /// OCCT NCollection_DataMap<int, Bisector_Bisec> theGeomBisectors.
    the_geom_bisectors: HashMap<i32, BisectorBisec>,
    /// OCCT NCollection_DataMap<int, gp_Pnt2d> theGeomPnts.
    the_geom_pnts: HashMap<i32, DVec2>,
    /// OCCT NCollection_DataMap<int, gp_Vec2d> theGeomVecs.
    the_geom_vecs: HashMap<i32, DVec2>,
    /// OCCT NCollection_Sequence<int> theLinesLength.
    the_lines_length: Vec<i32>,
}

impl Mat2dTool2d {
    /// OCCT MAT2d_Tool2d::MAT2d_Tool2d() (cxx L83-90).
    pub fn new() -> Self {
        Mat2dTool2d {
            the_direction: 1.0,
            the_join_type: GeomAbsJoinType::Arc, // default
            the_number_of_bisectors: 0,
            the_number_of_vecs: 0,
            the_number_of_pnts: 0,
            the_circuit: None,
            the_geom_bisectors: HashMap::new(),
            the_geom_pnts: HashMap::new(),
            the_geom_vecs: HashMap::new(),
            the_lines_length: Vec::new(),
        }
    }

    /// OCCT MAT2d_Tool2d::InitItems(EquiCircuit) (cxx L94-105) — InitItems
    /// cuts the line in Items; the items are the geometric representations
    /// of the BasicElts from MAT.
    pub fn init_items(&mut self, equi_circuit: Arc<RwLock<Mat2dCircuit>>) {
        self.the_geom_bisectors.clear();
        self.the_geom_pnts.clear();
        self.the_geom_vecs.clear();
        self.the_lines_length.clear();
        self.the_number_of_bisectors = 0;
        self.the_number_of_vecs = 0;
        self.the_number_of_pnts = 0;

        self.the_circuit = Some(equi_circuit);
    }

    /// OCCT MAT2d_Tool2d::Sense(aside) (cxx L109-119) — aSide defines the
    /// side of the computation of the map.
    pub fn sense(&mut self, aside: MatSide) {
        if aside == MatSide::Left {
            self.the_direction = 1.0;
        } else {
            self.the_direction = -1.0;
        }
    }

    /// OCCT MAT2d_Tool2d::SetJoinType(aJoinType) (cxx L123-126).
    pub fn set_join_type(&mut self, a_join_type: GeomAbsJoinType) {
        self.the_join_type = a_join_type;
    }

    /// OCCT MAT2d_Tool2d::NumberOfItems() const (cxx L130-133).
    pub fn number_of_items(&self) -> i32 {
        self.the_circuit
            .as_ref()
            .unwrap()
            .read()
            .unwrap()
            .number_of_items()
    }

    /// OCCT MAT2d_Tool2d::ToleranceOfConfusion() const (cxx L137-140) —
    /// tolerance to test the confusion of two points.
    pub fn tolerance_of_confusion(&self) -> f64 {
        2.0 * MAT2D_TOLCONF
    }

    /// OCCT MAT2d_Tool2d::FirstPoint(anitem, dist) (cxx L144-174) — creates
    /// the point at the origin of the bisector between anitem and the
    /// previous item; dist is the distance from the FirstPoint to anitem.
    /// Returns the index of this point in theGeomPnts.
    pub fn first_point(&mut self, anitem: i32, dist: &mut f64) -> i32 {
        self.the_number_of_pnts += 1;

        let circuit = self.the_circuit.as_ref().unwrap().clone();
        let circuit = circuit.read().unwrap();

        if circuit.connexion_on(anitem) {
            let connexion = circuit.connexion(anitem);
            let connexion = connexion.read().unwrap();
            let p1 = connexion.point_on_first();
            let p2 = connexion.point_on_second();
            self.the_geom_pnts.insert(
                self.the_number_of_pnts,
                DVec2::new((p1.x + p2.x) * 0.5, (p1.y + p2.y) * 0.5),
            );
            *dist = p1.distance(p2) * 0.5;
            return self.the_number_of_pnts;
        }

        let value = circuit.value(anitem);
        *dist = 0.0;

        // OCCT L163-172.
        if !value.is_point() {
            let curve = value.curve();
            self.the_geom_pnts
                .insert(self.the_number_of_pnts, curve.value(curve.first_parameter()));
        } else {
            let point = value.pnt2d();
            self.the_geom_pnts.insert(self.the_number_of_pnts, point);
        }
        self.the_number_of_pnts
    }

    /// OCCT MAT2d_Tool2d::TangentBefore(anitem, IsOpenResult) (cxx L178-218)
    /// — creates the tangent at the end of the item anitem; returns the
    /// index of this vector in theGeomVecs.
    pub fn tangent_before(&mut self, anitem: i32, is_open_result: bool) -> i32 {
        self.the_number_of_vecs += 1;

        let circuit = self.the_circuit.as_ref().unwrap().clone();
        let circuit = circuit.read().unwrap();

        let item: i32;
        if !is_open_result {
            item = if anitem == circuit.number_of_items() {
                1
            } else {
                anitem + 1
            };
        } else {
            item = if anitem == circuit.number_of_items() {
                anitem - 1
            } else {
                anitem + 1
            };
        }
        if circuit.connexion_on(item) {
            let connexion = circuit.connexion(item);
            let connexion = connexion.read().unwrap();
            let p1 = connexion.point_on_first();
            let p2 = connexion.point_on_second();
            // OCCT L197: gp_Vec2d((x2 - x1), (y2 - y1)).
            self.the_geom_vecs
                .insert(self.the_number_of_vecs, DVec2::new(p2.x - p1.x, p2.y - p1.y));
            return self.the_number_of_vecs;
        }

        // OCCT L201-215.
        let value = circuit.value(anitem);
        if !value.is_point() {
            let curve = value.curve();
            self.the_geom_vecs
                .insert(self.the_number_of_vecs, curve.dn(curve.last_parameter(), 1));
        } else {
            let curve = circuit.value(item).curve();
            let param = if is_open_result && anitem == circuit.number_of_items() {
                curve.last_parameter()
            } else {
                curve.first_parameter()
            };
            self.the_geom_vecs
                .insert(self.the_number_of_vecs, curve.dn(param, 1));
        }

        self.the_number_of_vecs
    }

    /// OCCT MAT2d_Tool2d::TangentAfter(anitem, IsOpenResult) (cxx L222-262)
    /// — creates the reversed tangent at the origin of the item anitem.
    pub fn tangent_after(&mut self, anitem: i32, is_open_result: bool) -> i32 {
        self.the_number_of_vecs += 1;

        let circuit = self.the_circuit.as_ref().unwrap().clone();
        let circuit = circuit.read().unwrap();

        let thevector: DVec2;

        if circuit.connexion_on(anitem) {
            let connexion = circuit.connexion(anitem);
            let connexion = connexion.read().unwrap();
            let p1 = connexion.point_on_first();
            let p2 = connexion.point_on_second();
            // OCCT L234: gp_Vec2d((x1 - x2), (y1 - y2)).
            self.the_geom_vecs
                .insert(self.the_number_of_vecs, DVec2::new(p1.x - p2.x, p1.y - p2.y));
            return self.the_number_of_vecs;
        }

        // OCCT L238-259.
        let value = circuit.value(anitem);
        if !value.is_point() {
            let curve = value.curve();
            thevector = curve.dn(curve.first_parameter(), 1);
        } else {
            let item: i32;
            if !is_open_result {
                item = if anitem == 1 {
                    circuit.number_of_items()
                } else {
                    anitem - 1
                };
            } else {
                item = if anitem == 1 { 2 } else { anitem - 1 };
            }

            let curve = circuit.value(item).curve();
            let param = if is_open_result && anitem == 1 {
                curve.first_parameter()
            } else {
                curve.last_parameter()
            };
            thevector = curve.dn(param, 1);
        }
        // OCCT L260: thevector.Reversed().
        self.the_geom_vecs
            .insert(self.the_number_of_vecs, -thevector);
        self.the_number_of_vecs
    }

    /// OCCT MAT2d_Tool2d::Tangent(bisector) (cxx L266-272) — creates the
    /// tangent at the end of the bisector.
    pub fn tangent(&mut self, bisector: i32) -> i32 {
        self.the_number_of_vecs += 1;
        let trimmed = self.geom_bis(bisector).value();
        let trimmed = trimmed.read().unwrap();
        self.the_geom_vecs.insert(
            self.the_number_of_vecs,
            trimmed.dn(trimmed.last_parameter(), 1),
        );
        self.the_number_of_vecs
    }

    /// OCCT MAT2d_Tool2d::CreateBisector(abisector) (cxx L276-410) —
    /// creates the geometric bisector defined by abisector.
    pub fn create_bisector(&mut self, abisector: &HandleMatBisector) {
        let tolerance = MAT2D_TOLCONF;

        let (edge1number, edge2number) = {
            let abis = abisector.read().unwrap();
            let edge1 = abis.first_edge().expect("MAT_Bisector::FirstEdge");
            let edge2 = abis.second_edge().expect("MAT_Bisector::SecondEdge");
            (
                edge1.read().unwrap().edge_number(),
                edge2.read().unwrap().edge_number(),
            )
        };
        let mut ontheline = are_neighbours(edge1number, edge2number, self.number_of_items());
        let initial_neighbour = ontheline;

        let circuit = self.the_circuit.as_ref().unwrap().clone();
        let circuit = circuit.read().unwrap();

        if circuit.connexion_on(edge2number) {
            ontheline = false;
        }

        let elt1 = circuit.value(edge1number);
        let elt2 = circuit.value(edge2number);

        let type1_is_point = elt1.is_point();
        let type2_is_point = elt2.is_point();

        // OCCT L300-311: item1/item2 = down_cast<Geom2d_Curve>(elt) — a
        // kernel Geom2d curve as a Geom2d/Bisector curve handle.
        let item1: Option<Arc<dyn BisectorCurve>> = if !type1_is_point {
            Some(Arc::new(Geom2dCurveHandle::new(elt1.curve().clone())))
        } else {
            None
        };
        let item2: Option<Arc<dyn BisectorCurve>> = if !type2_is_point {
            Some(Arc::new(Geom2dCurveHandle::new(elt2.curve().clone())))
        } else {
            None
        };
        let point1 = if type1_is_point { Some(elt1.pnt2d()) } else { None };
        let point2 = if type2_is_point { Some(elt2.pnt2d()) } else { None };

        let mut bisector = BisectorBisec::new();

        // OCCT L332-382 — the four Perform overloads.
        if !type1_is_point && !type2_is_point {
            let issue = self.geom_pnt(abisector.read().unwrap().issue_point());
            let v1 = self.geom_vec(abisector.read().unwrap().first_vector());
            let v2 = self.geom_vec(abisector.read().unwrap().second_vector());
            bisector.perform_curve_curve(
                item1.as_ref().unwrap(),
                item2.as_ref().unwrap(),
                issue,
                v1,
                v2,
                self.the_direction,
                self.the_join_type,
                tolerance,
                ontheline,
            );
        } else if type1_is_point && type2_is_point {
            let issue = self.geom_pnt(abisector.read().unwrap().issue_point());
            let v1 = self.geom_vec(abisector.read().unwrap().first_vector());
            let v2 = self.geom_vec(abisector.read().unwrap().second_vector());
            bisector.perform_point_point(
                point1.unwrap(),
                point2.unwrap(),
                issue,
                v1,
                v2,
                self.the_direction,
                tolerance,
                ontheline,
            );
        } else if type1_is_point {
            let issue = self.geom_pnt(abisector.read().unwrap().issue_point());
            let v1 = self.geom_vec(abisector.read().unwrap().first_vector());
            let v2 = self.geom_vec(abisector.read().unwrap().second_vector());
            bisector.perform_point_curve(
                point1.unwrap(),
                item2.as_ref().unwrap(),
                issue,
                v1,
                v2,
                self.the_direction,
                tolerance,
                ontheline,
            );
        } else {
            let issue = self.geom_pnt(abisector.read().unwrap().issue_point());
            let v1 = self.geom_vec(abisector.read().unwrap().first_vector());
            let v2 = self.geom_vec(abisector.read().unwrap().second_vector());
            bisector.perform_curve_point(
                item1.as_ref().unwrap(),
                point2.unwrap(),
                issue,
                v1,
                v2,
                self.the_direction,
                tolerance,
                ontheline,
            );
        }

        //------------------------------
        // Restriction de la bisectrice.  (OCCT L384-388)
        //-----------------------------
        drop(circuit);
        self.trim_bisec(&mut bisector, edge1number, initial_neighbour, 1);
        self.trim_bisec(&mut bisector, edge2number, initial_neighbour, 2);

        self.the_number_of_bisectors += 1;
        self.the_geom_bisectors
            .insert(self.the_number_of_bisectors, bisector);

        let mut abis = abisector.write().unwrap();
        abis.set_bisector_number(self.the_number_of_bisectors);
        abis.set_sense(1.0);
    }

    /// OCCT MAT2d_Tool2d::TrimBisec(B1, IndexEdge, InitialNeighbour,
    /// StartOrEnd) const (cxx L419-486) — restriction de la bisectrice
    /// separant deux elements lies par une connexion ou l un au moins des
    /// elements est un cercle.
    fn trim_bisec(
        &self,
        b1: &mut BisectorBisec,
        index_edge: i32,
        initial_neighbour: bool,
        start_or_end: i32,
    ) {
        let circuit = self.the_circuit.as_ref().unwrap().clone();
        let circuit = circuit.read().unwrap();

        let i_next = if index_edge == circuit.number_of_items() {
            1
        } else {
            index_edge + 1
        };

        let edge_value = circuit.value(index_edge);

        if !edge_value.is_point() {
            if !initial_neighbour {
                // OCCT L439-440: the trimmed item and its basis curve.
                let curve = edge_value.curve().inner();

                let mut line1: Option<Curve2d> = None;
                let mut line2: Option<Curve2d> = None;

                //-------------------------------------------------------------------
                // si l edge est liee a sa voisine precedente par une connexion.
                //-------------------------------------------------------------------
                if circuit.connexion_on(index_edge) && start_or_end == 1 {
                    if let Curve2d::Circle(circle) = curve {
                        let ori = circle.center;
                        let p2 =
                            circuit.connexion(index_edge).read().unwrap().point_on_first();
                        line1 = Some(Curve2d::Line(Line2d::new(
                            ori,
                            DVec2::new(p2.x - ori.x, p2.y - ori.y),
                        )));
                    }
                }
                //-----------------------------------------------------------------------
                // Si l edge est liee a sa voisine suivante par une connexion.
                //-----------------------------------------------------------------------
                if circuit.connexion_on(i_next) && start_or_end == 2 {
                    if let Curve2d::Circle(circle) = curve {
                        let ori = circle.center;
                        let p2 = circuit.connexion(i_next).read().unwrap().point_on_second();
                        line2 = Some(Curve2d::Line(Line2d::new(
                            ori,
                            DVec2::new(p2.x - ori.x, p2.y - ori.y),
                        )));
                    }
                }
                if line1.is_none() && line2.is_none() {
                    return;
                }

                //-----------------------------------------------------------------------
                // Restriction de la bisectrice par les demi-droites liees aux
                // connexions si elles existent.
                //-----------------------------------------------------------------------
                if let Some(l1) = line1 {
                    // OCCT L476: new Geom2d_TrimmedCurve(Line1, 0., Infinite()).
                    let line = Curve2d::Trimmed(TrimmedCurve2 {
                        curve: Box::new(l1),
                        t_min: 0.0,
                        t_max: INFINITE,
                    });
                    set_trim(b1, &line);
                }
                if let Some(l2) = line2 {
                    let line = Curve2d::Trimmed(TrimmedCurve2 {
                        curve: Box::new(l2),
                        t_min: 0.0,
                        t_max: INFINITE,
                    });
                    set_trim(b1, &line);
                }
            }
        }
    }

    /// OCCT MAT2d_Tool2d::TrimBisector(abisector) (cxx L490-518) — trims the
    /// geometric bisector by the firstparameter of abisector; returns false
    /// if the parameter is out of the bisector, else true.
    pub fn trim_bisector(&mut self, abisector: &HandleMatBisector) -> bool {
        let mut param = abisector.read().unwrap().first_parameter();

        let bisector = self
            .change_geom_bis(abisector.read().unwrap().bisector_number())
            .value();

        {
            let bis = bisector.read().unwrap();
            if bis.basis().is_periodic() && param == INFINITE {
                param = bis.first_parameter() + 2.0 * std::f64::consts::PI;
            }
            if param > bis.basis().last_parameter() {
                param = bis.basis().last_parameter();
            }
            if bis.first_parameter() == param {
                return false;
            }
        }

        let first = bisector.read().unwrap().first_parameter();
        bisector.write().unwrap().set_trim(first, param);
        true
    }

    /// OCCT MAT2d_Tool2d::TrimBisector(abisector, apoint) (cxx L522-557) —
    /// trims the geometric bisector by the point of index apoint in
    /// theGeomPnts; returns false if the point is out of the bisector, else
    /// true.
    pub fn trim_bisector_point(&mut self, abisector: &HandleMatBisector, apoint: i32) -> bool {
        let bisector = self
            .change_geom_bis(abisector.read().unwrap().bisector_number())
            .value();

        // OCCT L528-531: Bis = down_cast<Bisector_Curve>(BasisCurve());
        // Param = Bis->Parameter(GeomPnt(apoint)) — the merged trait call is
        // the direct mapping of the down-cast (architecture note).
        let mut param = {
            let bis = bisector.read().unwrap();
            bis.basis().parameter(self.geom_pnt(apoint))
        };

        {
            let bis = bisector.read().unwrap();
            if bis.basis().is_periodic() && bis.first_parameter() > param {
                param += 2.0 * std::f64::consts::PI;
            }
            if bis.first_parameter() >= param {
                return false;
            }
            if bis.last_parameter() < param {
                return false;
            }
        }
        let first = bisector.read().unwrap().first_parameter();
        bisector.write().unwrap().set_trim(first, param);

        true
    }

    /// OCCT MAT2d_Tool2d::Projection(IEdge, PCom, Distance) const
    /// (cxx L561-652) — true if the point can be projected on the element
    /// IndexElt; Distance is the minimum of the distances to the projections.
    fn projection(&self, i_edge: i32, p_com: DVec2, distance: &mut f64) -> bool {
        let circuit = self.the_circuit.as_ref().unwrap().clone();
        let circuit = circuit.read().unwrap();

        let elt = circuit.value(i_edge);
        let eps = MAT2D_TOLCONF;

        if elt.is_point() {
            let p_edge = elt.pnt2d();
            *distance = p_com.distance(p_edge);
        } else {
            *distance = INFINITE;
            let curve = elt.curve();
            //-----------------------------------------------------------------------
            // Calcul des parametres MinMax sur l edge si celui ci est lie a ses
            // voisins par des connexions la courbe de calcul est limitee par
            // celles_ci.
            //-----------------------------------------------------------------------
            let mut param_min = curve.first_parameter();
            let mut param_max = curve.last_parameter();
            if circuit.connexion_on(i_edge) {
                param_min = circuit
                    .connexion(i_edge)
                    .read()
                    .unwrap()
                    .parameter_on_second();
            }
            let i_next = if i_edge == circuit.number_of_items() {
                1
            } else {
                i_edge + 1
            };
            if circuit.connexion_on(i_next) {
                param_max = circuit
                    .connexion(i_next)
                    .read()
                    .unwrap()
                    .parameter_on_first();
                // OCCT L594: Curve->BasisCurve()->IsPeriodic().
                if curve.inner().is_periodic() {
                    // OCCT L596: ElCLib::AdjustPeriodic(0., 2*M_PI, Eps,
                    // ParamMin, ParamMax).
                    let (a1, a2) = elclib_adjust_periodic(
                        0.0,
                        2.0 * std::f64::consts::PI,
                        eps,
                        param_min,
                        param_max,
                    );
                    param_min = a1;
                    param_max = a2;
                }
            }
            //---------------------------------------------------------------------
            // Constuction de la courbe pour les extremas et ajustement des bornes.
            //---------------------------------------------------------------------
            let type_c1 = curve2d_type_of(curve.inner());
            if type_c1 == Curve2dType::Circle {
                let r = match curve.inner() {
                    Curve2d::Circle(c) => c.radius,
                    _ => unreachable!(),
                };
                let mut eps_circ = 100.0 * eps;
                if r < 1.0 {
                    eps_circ = eps / r;
                }
                if (param_max - param_min + 2.0 * eps_circ) < 2.0 * std::f64::consts::PI {
                    param_max += eps_circ;
                    param_min -= eps_circ;
                }
            } else {
                param_max += eps;
                param_min -= eps;
            }
            //-----------------------------------------------------
            // Calcul des extremas et stockage minimum de distance.
            //-----------------------------------------------------
            let mut extremas = rcad_kernel::base::extrema::ExtPC2d::new(
                p_com,
                curve,
                rcad_kernel::precision::CONFUSION,
                param_min,
                param_max,
            );
            if extremas.is_done() {
                *distance = INFINITE;
                if extremas.nb_ext() < 1 {
                    return false;
                }
                for i in 1..=extremas.nb_ext() {
                    if extremas.square_distance(i) < *distance {
                        *distance = extremas.square_distance(i);
                    }
                }
                *distance = distance.sqrt();
            } else if type_c1 == Curve2dType::Circle {
                let r = match curve.inner() {
                    Curve2d::Circle(c) => c.radius,
                    _ => unreachable!(),
                };
                *distance = r;
            }
        }
        true
    }

    /// OCCT MAT2d_Tool2d::IsSameDistance(BisectorOne, BisectorTwo, PCom,
    /// Distance) const (cxx L656-784) — true if the point is equidistant to
    /// the elements separated by the two bisectors.
    fn is_same_distance(
        &self,
        bisector_one: &HandleMatBisector,
        bisector_two: &HandleMatBisector,
        p_com: DVec2,
        distance: &mut f64,
    ) -> bool {
        // OCCT L661: NCollection_Array1<double> Dist(1, 4) — 0-based here.
        let mut dist = [0.0f64; 4];
        let eps = 1.0e-7;

        let (i_edge1, i_edge2) = {
            let one = bisector_one.read().unwrap();
            (
                one.first_edge().unwrap().read().unwrap().edge_number(),
                one.second_edge().unwrap().read().unwrap().edge_number(),
            )
        };
        let (i_edge3, i_edge4) = {
            let two = bisector_two.read().unwrap();
            (
                two.first_edge().unwrap().read().unwrap().edge_number(),
                two.second_edge().unwrap().read().unwrap().edge_number(),
            )
        };

        let is_done1 = self.projection(i_edge1, p_com, &mut dist[0]);
        let is_done2 = self.projection(i_edge2, p_com, &mut dist[1]);

        // OCCT L673-696.
        if is_done1 {
            if !is_done2 {
                let circuit = self.the_circuit.as_ref().unwrap().clone();
                let circuit = circuit.read().unwrap();
                let elt = circuit.value(i_edge2);
                let tol = rcad_kernel::precision::CONFUSION.max(eps * dist[0]);
                if check_ends(elt, p_com, dist[0], tol) {
                    dist[1] = dist[0];
                }
            }
        } else if is_done2 {
            let circuit = self.the_circuit.as_ref().unwrap().clone();
            let circuit = circuit.read().unwrap();
            let elt = circuit.value(i_edge1);
            let tol = rcad_kernel::precision::CONFUSION.max(eps * dist[1]);
            if check_ends(elt, p_com, dist[1], tol) {
                dist[0] = dist[1];
            }
        }

        // OCCT L698-748.
        let mut is_done3 = true;
        let mut is_done4 = true;
        if i_edge3 == i_edge1 {
            dist[2] = dist[0];
        } else if i_edge3 == i_edge2 {
            dist[2] = dist[1];
        } else {
            is_done3 = self.projection(i_edge3, p_com, &mut dist[2]);
        }

        if i_edge4 == i_edge1 {
            dist[3] = dist[0];
        } else if i_edge4 == i_edge2 {
            dist[3] = dist[1];
        } else {
            is_done4 = self.projection(i_edge4, p_com, &mut dist[3]);
        }

        if is_done3 {
            if !is_done4 {
                let circuit = self.the_circuit.as_ref().unwrap().clone();
                let circuit = circuit.read().unwrap();
                let elt = circuit.value(i_edge4);
                let tol = rcad_kernel::precision::CONFUSION.max(eps * dist[2]);
                if check_ends(elt, p_com, dist[2], tol) {
                    dist[3] = dist[2];
                }
            }
        } else if is_done4 {
            let circuit = self.the_circuit.as_ref().unwrap().clone();
            let circuit = circuit.read().unwrap();
            let elt = circuit.value(i_edge3);
            let tol = rcad_kernel::precision::CONFUSION.max(eps * dist[3]);
            if check_ends(elt, p_com, dist[3], tol) {
                dist[2] = dist[3];
            }
        }

        // OCCT L758-783.
        let eps_dist = MAT2D_TOLCONF * 300.0;
        *distance = dist[0];
        if self.the_join_type == GeomAbsJoinType::Intersection
            && rcad_kernel::precision::is_infinite_value(*distance)
        {
            for item in dist.iter().skip(1) {
                if !rcad_kernel::precision::is_infinite_value(*item) {
                    *distance = *item;
                    break;
                }
            }
        }
        for item in dist.iter() {
            if self.the_join_type == GeomAbsJoinType::Intersection
                && rcad_kernel::precision::is_infinite_value(*item)
            {
                continue;
            }
            if (item - *distance).abs() > eps_dist {
                *distance = INFINITE;
                return false;
            }
        }
        true
    }

    /// OCCT MAT2d_Tool2d::IntersectBisector(BisectorOne, BisectorTwo, IntPnt)
    /// (cxx L788-1102) — computes the point of intersection between the two
    /// bisectors; if it exists IntPnt is its index in theGeomPnts and the
    /// distance of the point from the bisectors is returned, else RealLast.
    pub fn intersect_bisector(
        &mut self,
        bisector_one: &HandleMatBisector,
        bisector_two: &HandleMatBisector,
        int_pnt: &mut i32,
    ) -> f64 {
        let tolerance = MAT2D_TOLCONF;
        let mut param1: f64;
        let mut param2: f64;
        let mut parama: f64;
        let mut paramb: f64;
        let mut distance: f64 = 0.0;
        let mut point_solution = DVec2::ZERO;

        let bisector1 = self
            .change_geom_bis(bisector_one.read().unwrap().bisector_number())
            .value();
        let bisector2 = self
            .change_geom_bis(bisector_two.read().unwrap().bisector_number())
            .value();

        // (OCCT L805-808: null handles cannot occur once bound.)

        //-------------------------------------------------------------------------
        // Si les deux bissectrices separent des elements consecutifs et qu elles
        // sont issues des connexions C1 et C2.
        // Si C1 est la reverse de C2, alors les deux bissectrices sont issues
        // du meme point. Dans ce cas l intersection n est pas validee.
        //-------------------------------------------------------------------------
        let (is1, is2, if1, if2) = {
            let one = bisector_one.read().unwrap();
            let two = bisector_two.read().unwrap();
            (
                one.second_edge().unwrap().read().unwrap().edge_number(),
                two.second_edge().unwrap().read().unwrap().edge_number(),
                one.first_edge().unwrap().read().unwrap().edge_number(),
                two.first_edge().unwrap().read().unwrap().edge_number(),
            )
        };

        {
            let circuit = self.the_circuit.as_ref().unwrap().clone();
            let circuit = circuit.read().unwrap();
            if are_neighbours(if1, is1, self.number_of_items())
                && are_neighbours(if2, is2, self.number_of_items())
                && circuit.connexion_on(is2)
                && circuit.connexion_on(is1)
            {
                let c1 = circuit.connexion(is1);
                let c2 = circuit.connexion(is2);
                let (c1r, c2r) = (c1.read().unwrap(), c2.read().unwrap());
                if c2r.index_first_line() == c1r.index_second_line()
                    && c1r.index_first_line() == c2r.index_second_line()
                {
                    return INFINITE;
                }
            }
        }

        // -----------------------------------------
        // Construction des domaines d intersection.  (OCCT L834-847)
        // -----------------------------------------
        let domain1 = domain(&bisector1, tolerance);
        let domain2 = domain(&bisector2, tolerance);

        if domain1.last_parameter() - domain1.first_parameter() < tolerance {
            return INFINITE;
        }
        if domain2.last_parameter() - domain2.first_parameter() < tolerance {
            return INFINITE;
        }

        // -------------------------
        // Calcul de l intersection.  (OCCT L871-882)
        // -------------------------
        let mut intersect = BisectorInter::new();
        intersect.perform(
            self.geom_bis(bisector_one.read().unwrap().bisector_number()),
            &domain1,
            self.geom_bis(bisector_two.read().unwrap().bisector_number()),
            &domain2,
            tolerance,
            tolerance,
            true,
        );

        // -------------------------------------------------------------------------
        // Exploitation du resultat de l intersection et selection du point solution
        // equidistant des deux edges et le plus proche en parametre de l origine
        // des bissectrices.  (OCCT L887-906)
        // -------------------------------------------------------------------------
        if !intersect.is_done() {
            return INFINITE;
        }
        if intersect.is_empty() {
            return INFINITE;
        }

        let mut distance_mini = INFINITE;
        param1 = INFINITE;
        param2 = INFINITE;
        let mut solution_valide = false;

        // OCCT L908-969: segments.
        if intersect.nb_segments() >= 1 {
            let max_segment_length = 10.0 * tolerance;
            for i in 1..=intersect.nb_segments() {
                let segment = intersect.segment(i).clone();
                let mut point_retenu = false;
                let mut point_on_segment = DVec2::ZERO;
                //----------------------------------------------------------------
                // Si les segments sont petits, recherche des points sur le segment
                // equidistants des edges.
                //----------------------------------------------------------------
                if segment.has_first_point() && segment.has_last_point() {
                    let p1 = segment.first_point().value();
                    let p2 = segment.last_point().value();
                    let segment_length = p1.distance(p2);
                    if segment_length <= tolerance {
                        point_on_segment = p1;
                        if self.is_same_distance(
                            bisector_one,
                            bisector_two,
                            point_on_segment,
                            &mut distance,
                        ) {
                            point_retenu = true;
                        }
                    } else if segment_length <= max_segment_length {
                        // OCCT L937: gp_Dir2d Dir(P2.X()-P1.X(), P2.Y()-P1.Y()).
                        let dir = DVec2::new(p2.x - p1.x, p2.y - p1.y).normalize_or_zero();
                        let mut dist = 0.0f64;
                        while dist <= segment_length + tolerance {
                            // OCCT L941: P1.Translated(Dist * Dir).
                            point_on_segment = p1 + dist * dir;
                            if self.is_same_distance(
                                bisector_one,
                                bisector_two,
                                point_on_segment,
                                &mut distance,
                            ) {
                                point_retenu = true;
                                break;
                            }
                            dist += tolerance;
                        }
                    }
                }

                //----------------------------------------------------------------
                // Sauvegarde du point equidistant des edges de plus petit
                // parametre sur les bissectrices.
                //----------------------------------------------------------------
                if point_retenu {
                    // OCCT L958-959: down_cast<Bisector_Curve>(BasisCurve())
                    // ->Parameter(...) — merged trait call.
                    parama = bisector1.read().unwrap().basis().parameter(point_on_segment);
                    paramb = bisector2.read().unwrap().basis().parameter(point_on_segment);
                    if parama < param1 && paramb < param2 {
                        param1 = parama;
                        param2 = paramb;
                        distance_mini = distance;
                        point_solution = point_on_segment;
                        solution_valide = true;
                    }
                }
            }
        }

        // OCCT L972-998: points.
        if intersect.nb_points() != 1 {
            for i in 1..=intersect.nb_points() {
                let p_i = intersect.point(i).value();
                if self.is_same_distance(bisector_one, bisector_two, p_i, &mut distance)
                    && distance > tolerance
                {
                    parama = intersect.point(i).param_on_first();
                    paramb = intersect.point(i).param_on_second();
                    if parama < param1 && paramb < param2 {
                        param1 = parama;
                        param2 = paramb;
                        distance_mini = distance;
                        point_solution = p_i;
                        solution_valide = true;
                    }
                }
            }
        } else {
            point_solution = intersect.point(1).value();
            param1 = intersect.point(1).param_on_first();
            param2 = intersect.point(1).param_on_second();
            solution_valide = self.is_same_distance(
                bisector_one,
                bisector_two,
                point_solution,
                &mut distance_mini,
            );
        }

        if !solution_valide {
            return INFINITE;
        }
        self.the_number_of_pnts += 1;
        self.the_geom_pnts
            .insert(self.the_number_of_pnts, point_solution);
        *int_pnt = self.the_number_of_pnts;

        //-----------------------------------------------------------------------
        // Si le point d intersection est quasi confondue avec une des extremites
        // de l une ou l autre des bisectrices, l intersection n est pas validee.
        //
        // SAUF si une des bisectrices est issue d une connexion et que les
        // edges separes par les bissectrices sont des voisines sur le contour
        // initiales.  (OCCT L1008-1055)
        //-----------------------------------------------------------------------
        let (index_edge1, index_edge2, index_edge3, index_edge4) = {
            let one = bisector_one.read().unwrap();
            let two = bisector_two.read().unwrap();
            (
                one.first_edge().unwrap().read().unwrap().edge_number(),
                one.second_edge().unwrap().read().unwrap().edge_number(),
                two.first_edge().unwrap().read().unwrap().edge_number(),
                two.second_edge().unwrap().read().unwrap().edge_number(),
            )
        };
        let mut extremite_controle = true;

        {
            let circuit = self.the_circuit.as_ref().unwrap().clone();
            let circuit = circuit.read().unwrap();
            if circuit.connexion_on(index_edge2) {
                // --------------------------------------
                // BisectorOne est issue d une connexion.
                // --------------------------------------
                if are_neighbours(index_edge1, index_edge2, self.number_of_items())
                    && are_neighbours(index_edge3, index_edge4, self.number_of_items())
                    && index_edge2 == index_edge3
                {
                    extremite_controle = false;
                    param1 += tolerance;
                }
            }

            if circuit.connexion_on(index_edge4) {
                //--------------------------------------
                // BisectorTwo est issue d une connexion.
                //--------------------------------------
                if are_neighbours(index_edge1, index_edge2, self.number_of_items())
                    && are_neighbours(index_edge3, index_edge4, self.number_of_items())
                    && index_edge2 == index_edge3
                {
                    extremite_controle = false;
                    param2 += tolerance;
                }
            }
        }

        // OCCT L1063-1074 (StartPoint() = Value(FirstParameter())).
        if extremite_controle {
            {
                let b1 = bisector1.read().unwrap();
                let start = b1.value(b1.first_parameter());
                if start.distance(point_solution) < tolerance {
                    return INFINITE;
                }
            }
            {
                let b2 = bisector2.read().unwrap();
                let start = b2.value(b2.first_parameter());
                if start.distance(point_solution) < tolerance {
                    return INFINITE;
                }
            }
        }

        // OCCT L1076-1086.
        {
            let one = bisector_one.read().unwrap();
            if one.second_parameter() < INFINITE
                && one.second_parameter() < param1 * (1.0 - tolerance)
            {
                return INFINITE;
            }
        }
        {
            let two = bisector_two.read().unwrap();
            if two.first_parameter() < INFINITE && two.first_parameter() < param2 * (1.0 - tolerance)
            {
                return INFINITE;
            }
        }

        bisector_one.write().unwrap().set_second_parameter(param1);
        bisector_two.write().unwrap().set_first_parameter(param2);

        distance_mini
    }

    /// OCCT MAT2d_Tool2d::Distance(Bis, Param1, Param2) const (cxx
    /// L1106-1119) — the distance between the two points designed by their
    /// parameters on abisector.
    pub fn distance(&self, bis: &HandleMatBisector, param1: f64, param2: f64) -> f64 {
        let mut dist = INFINITE;

        if param1 != INFINITE && param2 != INFINITE {
            let trimmed = self.geom_bis(bis.read().unwrap().bisector_number()).value();
            let trimmed = trimmed.read().unwrap();
            let p1 = trimmed.value(param1);
            let p2 = trimmed.value(param2);
            dist = p1.distance(p2);
        }
        dist
    }

    /// OCCT MAT2d_Tool2d::Dump(bisector, erease) const (cxx L1123-1140) —
    /// the release build raises Standard_NotImplemented.
    pub fn dump(&self, _bisector: i32, _erease: i32) {
        panic!("Standard_NotImplemented: MAT2d_Tool2d::Dump");
    }

    /// OCCT MAT2d_Tool2d::GeomBis(Index) const (cxx L1144-1147) — the Bisec
    /// of index Index in theGeomBisectors.
    pub fn geom_bis(&self, index: i32) -> &BisectorBisec {
        self.the_geom_bisectors
            .get(&index)
            .expect("MAT2d_Tool2d::GeomBis")
    }

    /// OCCT MAT2d_Tool2d::ChangeGeomBis(Index) (cxx L1151-1154).
    pub fn change_geom_bis(&mut self, index: i32) -> &mut BisectorBisec {
        self.the_geom_bisectors
            .get_mut(&index)
            .expect("MAT2d_Tool2d::ChangeGeomBis")
    }

    /// OCCT MAT2d_Tool2d::GeomElt(Index) const (cxx L1158-1161) — the
    /// geometry of index Index in theGeomElts.
    pub fn geom_elt(&self, index: i32) -> Geom2dGeometry {
        self.the_circuit
            .as_ref()
            .unwrap()
            .read()
            .unwrap()
            .value(index)
            .clone()
    }

    /// OCCT MAT2d_Tool2d::GeomPnt(Index) const (cxx L1165-1168) — the point
    /// of index Index in theGeomPnts.
    pub fn geom_pnt(&self, index: i32) -> DVec2 {
        self.the_geom_pnts[&index]
    }

    /// OCCT MAT2d_Tool2d::GeomVec(Index) const (cxx L1172-1175) — the
    /// vector of index Index in theGeomVecs.
    pub fn geom_vec(&self, index: i32) -> DVec2 {
        self.the_geom_vecs[&index]
    }

    /// OCCT MAT2d_Tool2d::Circuit() const (cxx L1179-1182).
    pub fn circuit(&self) -> Arc<RwLock<Mat2dCircuit>> {
        self.the_circuit.as_ref().unwrap().clone()
    }

    /// OCCT MAT2d_Tool2d::BisecFusion(I1, I2) (cxx L1186-1249) — fusion of
    /// the two bisectors.
    pub fn bisec_fusion(&mut self, i1: i32, i2: i32) {
        let mut bisector1 = self.geom_bis(i1).value();
        let bisector2 = self.geom_bis(i2).value();

        let mut uf1;
        let mut ul1;

        {
            let b1 = bisector1.read().unwrap();
            uf1 = b1.first_parameter();
            ul1 = b1.last_parameter();
        }

        let type1 = bisector1.read().unwrap().basis().kind();

        if type1 == CurveKind::BisecCC {
            //------------------------------------------------------------------------------------
            // les bissectrice courbe/courbe sont construites avec un point de depart
            // elles ne peuvent pas etre trimes par un point se trouvant de l autre cote du
            // point de depart.
            // pour faire la fusion des deux bissectrices on reconstruit la bissectrice entre
            // les deux courbes avec comme point de depart le dernier point de la Bisector2.
            // on trime ensuite la courbe par le dernier point de Bisector1.
            //------------------------------------------------------------------------------------
            let tolerance = MAT2D_TOLCONF;
            let mut bis = BisectorBisec::new();
            let v_bid = DVec2::new(1.0, 0.0);
            let p2 = bisector2
                .read()
                .unwrap()
                .value(bisector2.read().unwrap().last_parameter());
            let p1 = bisector1
                .read()
                .unwrap()
                .value(bisector1.read().unwrap().last_parameter());
            // OCCT L1213: BCC1 = down_cast<Bisector_BisecCC>(BasisCurve()).
            let bcc1 = bisector_bisec_cc::cast_from(bisector1.read().unwrap().basis().clone())
                .expect("down_cast<Bisector_BisecCC>");

            bis.perform_curve_curve(
                &bcc1.curve(2),
                &bcc1.curve(1),
                p2,
                v_bid,
                v_bid,
                self.the_direction,
                self.the_join_type,
                tolerance,
                false,
            );

            bisector1 = bis.value();
            let bcc1 = bisector_bisec_cc::cast_from(bisector1.read().unwrap().basis().clone())
                .expect("down_cast<Bisector_BisecCC>");
            uf1 = bcc1.first_parameter();
            ul1 = bcc1.parameter(p1);
            bisector1.write().unwrap().set_trim(uf1, ul1);
            self.the_geom_bisectors.insert(i1, bis);
        } else {
            // OCCT L1234-1235.
            let du = bisector2.read().unwrap().last_parameter()
                - bisector2.read().unwrap().first_parameter();
            uf1 -= du;

            // OCCT L1237-1245: BAna = down_cast<Bisector_BisecAna>(...);
            // BAna->SetTrim(UF1, UL1).
            let basis = bisector1.read().unwrap().basis().clone();
            bisector_bisec_ana::set_trim_of(&basis, uf1, ul1);

            bisector1.write().unwrap().set_trim(uf1, ul1);
        }
    }
}

impl Default for Mat2dTool2d {
    fn default() -> Self {
        Self::new()
    }
}
