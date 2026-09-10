//! OCCT Geom2dHatch_Elements (TKGeomAlgo/Geom2dHatch).
//!
//! Geom2dHatch_Elements.hxx L33-93 + .cxx L26-263 — the data map of hatching
//! elements that plays the "face explorer" role of the TopClass_FaceClassifier
//! for the 2D hatching: it provides the CheckPoint/Reject/Segment probing of
//! the classifier and the wire/edge iteration (a single wire whose edges are
//! the bound elements).

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Line2d};
use rcad_kernel::precision::{is_negative_infinite_value, is_positive_infinite_value};
use rcad_kernel::topods::Orientation;

use crate::geomalgo::geom2d_int::Curve2dAdaptor;
use crate::geomalgo::hatch::element::HatchElement;

// OCCT Geom2dHatch_Elements.cxx L26-28 — the probing parameter constants.
const PROBING_START: f64 = 0.123;
const PROBING_END: f64 = 0.8;
const PROBING_STEP: f64 = 0.2111;

/// OCCT Geom2dHatch_Elements — the element map with the traversal state.
///
/// Architecture note: OCCT NCollection_DataMap<int, Geom2dHatch_Element>
/// becomes a `Vec<(usize, HatchElement)>` so that the DataMap iterator
/// (`Iter`) maps to a plain position index with deterministic order.
#[derive(Debug)]
pub struct HatchElements {
    /// OCCT NCollection_DataMap<int, Geom2dHatch_Element> myMap.
    my_map: Vec<(usize, HatchElement)>,
    /// OCCT NCollection_DataMap<int, Geom2dHatch_Element>::Iterator Iter —
    /// the current edge position (`my_map.len()` = past the end, i.e. the
    /// uninitialized/default iterator with More() == false).
    iter: usize,
    /// OCCT int NumWire.
    num_wire: i32,
    /// OCCT int NumEdge.
    num_edge: i32,
    /// OCCT int myCurEdge.
    my_cur_edge: usize,
    /// OCCT double myCurEdgePar.
    my_cur_edge_par: f64,
}

impl HatchElements {
    /// OCCT Geom2dHatch_Elements() (cxx L41-47).
    pub fn new() -> Self {
        HatchElements {
            my_map: Vec::new(),
            iter: usize::MAX, // default-constructed iterator (More() == false)
            num_wire: 0,
            num_edge: 0,
            my_cur_edge: 1,
            my_cur_edge_par: PROBING_START,
        }
    }

    /// OCCT Clear (cxx L49-52).
    pub fn clear(&mut self) {
        self.my_map.clear();
    }

    /// OCCT IsBound (cxx L54-57).
    pub fn is_bound(&self, k: usize) -> bool {
        self.my_map.iter().any(|(key, _)| *key == k)
    }

    /// OCCT UnBind (cxx L59-62).
    pub fn un_bind(&mut self, k: usize) -> bool {
        if let Some(pos) = self.my_map.iter().position(|(key, _)| *key == k) {
            self.my_map.remove(pos);
            true
        } else {
            false
        }
    }

    /// OCCT Bind (cxx L64-67) — NCollection_DataMap::Bind returns true when
    /// the key was not bound already; an already bound key gets its item
    /// overridden and returns false.
    pub fn bind(&mut self, k: usize, element: HatchElement) -> bool {
        if let Some(pos) = self.my_map.iter().position(|(key, _)| *key == k) {
            self.my_map[pos].1 = element;
            return false;
        }
        self.my_map.push((k, element));
        true
    }

    /// OCCT Find (cxx L69-72) — raises NoSuchObject when the key is unbound.
    pub fn find(&self, k: usize) -> &HatchElement {
        self.my_map
            .iter()
            .find(|(key, _)| *key == k)
            .map(|(_, e)| e)
            .expect("Geom2dHatch_Elements::Find - Standard_NoSuchObject")
    }

    /// OCCT ChangeFind (cxx L74-77) — raises NoSuchObject when unbound.
    pub fn change_find(&mut self, k: usize) -> &mut HatchElement {
        self.my_map
            .iter_mut()
            .find(|(key, _)| *key == k)
            .map(|(_, e)| e)
            .expect("Geom2dHatch_Elements::ChangeFind - Standard_NoSuchObject")
    }

    /// OCCT CheckPoint (cxx L81-84) — no-op check, always accepts the point.
    pub fn check_point(&mut self, p: &mut DVec2) -> bool {
        let _ = p;
        true
    }

    /// OCCT Reject (cxx L88-91) — no rejection.
    pub fn reject(&self, p: DVec2) -> bool {
        let _ = p;
        false
    }

    /// OCCT Segment (cxx L95-100) — starts the probing from the first edge.
    pub fn segment(&mut self, p: DVec2, l: &mut Line2d, par: &mut f64) -> bool {
        self.my_cur_edge = 1;
        self.my_cur_edge_par = PROBING_START;
        self.other_segment(p, l, par)
    }

    /// OCCT OtherSegment (cxx L104-197) — probes the edges from myCurEdge on,
    /// building a line from P to a probed point on an edge that is neither
    /// tangent to the edge nor passing through its endpoints.
    pub fn other_segment(&mut self, p: DVec2, l: &mut Line2d, par: &mut f64) -> bool {
        let mut entry = 0usize; // NCollection_DataMap iterator position
        let mut i = 1usize; // 1-based ordinal of the current entry
        while entry < self.my_map.len() {
            // for (Itertemp.Initialize(myMap), i = 1; Itertemp.More();
            //      Itertemp.Next(), i++)
            if i < self.my_cur_edge {
                // continue -> Itertemp.Next(), i++
                entry += 1;
                i += 1;
                continue;
            }

            // OCCT uses myMap.ChangeFind(Itertemp.Key()) / Item.ChangeCurve()
            // (mutable access) although the curve is only read below (D1 /
            // Value) — a C++ const-ness workaround.
            let item = &self.my_map[entry].1;
            let e = item.curve();
            let or = item.orientation();
            if or == Orientation::Forward || or == Orientation::Reversed {
                let mut a_f_par = e.first_parameter();
                let mut a_l_par = e.last_parameter();
                if is_negative_infinite_value(a_f_par) {
                    if is_positive_infinite_value(a_l_par) {
                        a_f_par = -1.;
                        a_l_par = 1.;
                    } else {
                        a_f_par = a_l_par - 1.;
                    }
                } else if is_positive_infinite_value(a_l_par) {
                    a_l_par = a_f_par + 1.;
                }

                // for (; myCurEdgePar < Probing_End; myCurEdgePar += Probing_Step)
                while self.my_cur_edge_par < PROBING_END {
                    let a_param =
                        self.my_cur_edge_par * a_f_par + (1. - self.my_cur_edge_par) * a_l_par;
                    // OCCT E.D1(aParam, aPOnC, aTanVec).
                    let (a_p_on_c, a_tan_vec) = e.d1(a_param);
                    // OCCT gp_Vec2d aLinVec(P, aPOnC).
                    let a_lin_vec = a_p_on_c - p;
                    *par = a_lin_vec.length_squared();
                    if *par > rcad_kernel::precision::square_p_confusion() {
                        // OCCT gp_Dir2d aLinDir(aLinVec) — normalizes.
                        let a_lin_dir = a_lin_vec.normalize();
                        let a_tan_mod = a_tan_vec.length_squared();
                        if a_tan_mod < rcad_kernel::precision::square_p_confusion() {
                            // continue -> myCurEdgePar += Probing_Step
                            self.my_cur_edge_par += PROBING_STEP;
                            continue;
                        }

                        let a_tan_vec = a_tan_vec / a_tan_mod.sqrt();
                        // OCCT aTanVec.Crossed(aLinDir).
                        let a_sin_a = a_tan_vec.x * a_lin_dir.y - a_tan_vec.y * a_lin_dir.x;
                        if a_sin_a.abs() < 0.001 {
                            // too small angle - line and edge may be considered
                            // as tangent which is bad for classifier
                            if self.my_cur_edge_par + PROBING_STEP < PROBING_END {
                                // continue -> myCurEdgePar += Probing_Step
                                self.my_cur_edge_par += PROBING_STEP;
                                continue;
                            }
                        }

                        *l = Line2d::new(p, a_lin_dir);

                        let mut a_p_on_c = e.value(a_f_par);
                        if line_square_distance(l, a_p_on_c)
                            > rcad_kernel::precision::square_p_confusion()
                        {
                            a_p_on_c = e.value(a_l_par);
                            if line_square_distance(l, a_p_on_c)
                                > rcad_kernel::precision::square_p_confusion()
                            {
                                self.my_cur_edge_par += PROBING_STEP;
                                if self.my_cur_edge_par >= PROBING_END {
                                    self.my_cur_edge += 1;
                                    self.my_cur_edge_par = PROBING_START;
                                }
                                *par = par.sqrt();
                                return true;
                            }
                        }
                    }
                    // loop increment: myCurEdgePar += Probing_Step
                    self.my_cur_edge_par += PROBING_STEP;
                }
            }
            self.my_cur_edge += 1;
            self.my_cur_edge_par = PROBING_START;
            entry += 1;
            i += 1;
        }

        *par = f64::MAX; // RealLast()
        *l = Line2d::new(p, DVec2::X); // gp_Dir2d(gp_Dir2d::D::X)
        false
    }

    /// OCCT InitWires (cxx L201-204) — a single wire.
    pub fn init_wires(&mut self) {
        self.num_wire = 0;
    }

    /// OCCT MoreWires (cxx L239-242).
    pub fn more_wires(&self) -> bool {
        self.num_wire == 0
    }

    /// OCCT NextWire (cxx L246-249).
    pub fn next_wire(&mut self) {
        self.num_wire += 1;
    }

    /// OCCT RejectWire (cxx L208-211) — no rejection.
    pub fn reject_wire(&self, l: &Line2d, par: f64) -> bool {
        let _ = (l, par);
        false
    }

    /// OCCT InitEdges (cxx L215-219) — resets the edge iterator over the map.
    pub fn init_edges(&mut self) {
        self.num_edge = 0;
        self.iter = 0; // Iter.Initialize(myMap)
    }

    /// OCCT MoreEdges (cxx L253-256).
    pub fn more_edges(&self) -> bool {
        self.iter < self.my_map.len()
    }

    /// OCCT NextEdge (cxx L260-263).
    pub fn next_edge(&mut self) {
        self.iter += 1;
    }

    /// OCCT RejectEdge (cxx L223-226) — no rejection.
    pub fn reject_edge(&self, l: &Line2d, par: f64) -> bool {
        let _ = (l, par);
        false
    }

    /// OCCT CurrentEdge (cxx L230-235) — the current edge curve and
    /// orientation (OCCT writes the out-parameters E and Or; rcad returns
    /// them as a pair).
    pub fn current_edge(&self) -> (Curve2d, Orientation) {
        let an_item = &self.my_map[self.iter].1;
        (an_item.curve().clone(), an_item.orientation())
    }
}

impl Clone for HatchElements {
    /// OCCT Geom2dHatch_Elements(const Geom2dHatch_Elements&) (cxx L30-39) —
    /// the "magic" copy constructor: it does NOT copy the source map; it
    /// creates an empty map with the traversal state zeroed (NumWire =
    /// NumEdge = myCurEdge = 0, myCurEdgePar = 0.0) and leaves the iterator
    /// uninitialized (OCC12627).
    fn clone(&self) -> Self {
        HatchElements {
            my_map: Vec::new(),
            iter: usize::MAX, // the DataMap iterator stays uninitialized
            num_wire: 0,
            num_edge: 0,
            my_cur_edge: 0,
            my_cur_edge_par: 0.0,
        }
    }
}

impl Default for HatchElements {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT gp_Lin2d::SquareDistance(gp_Pnt2d) — squared distance of a point to
/// the (unit-direction) line.
fn line_square_distance(l: &Line2d, p: DVec2) -> f64 {
    let d = p - l.origin;
    let t = d.dot(l.direction);
    let v = d - l.direction * t;
    v.length_squared()
}
