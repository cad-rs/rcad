//! OCCT ShapeAnalysis package class (TKShHealing): `ShapeAnalysis_BoxBndTree`
//! (`ShapeAnalysis_BoxBndTree.hxx` L23-110 + `.cxx` L1-195).
//!
//! The selector used to connect wires through the bounding-box tree: given
//! the first/last vertices (or points) of a growing wire, `Accept` decides
//! whether the candidate wire continues it (and at which end), driving the
//! ShapeExtend statuses DONE1..DONE4 / FAIL2.
//!
//! Architecture bridges (the W1-5 free_bounds.rs carrier conventions):
//! 1. `BRep` pool argument — `BRep_Tool::Pnt` reads the TShape graph through
//!    `rcad_kernel::BRep` (the edge.rs bridge #1).
//! 2. The OCCT selector stores the `NCollection_HArray1<TopoDS_Shape>` wire
//!    array handle; the rcad value model receives the array slice on
//!    `Accept` (bridge #4 of free_bounds.rs — the aliasing must stay
//!    visible at the call site).
//! 3. `Bnd_Box` -> `rcad_kernel::math::bnd::BndBox`.
//! 4. `NCollection_Map<int> myList` -> `HashSet<i32>`;
//!    `NCollection_Array1<int> myArrIndices(1, 2)` -> `[i32; 2]` indexed by
//!    the Accept-local `First = 1, Last = 2` enum.
//! 5. `ShapeAnalysis::FindBounds` -> the W2 1:1 static
//!    [`super::analysis::find_bounds`].
//!
//! GAP carriers (untranslated other-package dependencies):
//! - `NCollection_UBTree<int, Bnd_Box>` +
//!   `NCollection_UBTreeFiller<int, Bnd_Box>` (TKFoundation, untranslated):
//!   the local [`NCollectionUBTree`] / [`NCollectionUBTreeFiller`] re-hosts
//!   keep the OCCT Add/Fill/Select call form; `Select` walks the filled
//!   entries in insertion order invoking Reject/Accept — the documented
//!   reduction of the OCCT tree traversal (the traversal order is the
//!   observable difference; the selector semantics are unchanged).  GAP:
//!   closes with the TKFoundation NCollection_UBTree batch.

use std::collections::HashSet;

use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, TShape};

use super::analysis::find_bounds;
use crate::shhealing::shape_extend::status::{decode_status, encode_status, ShapeExtendStatus};

// OCCT Standard_Real.hxx RealSmall(): the smallest positive representable
// real (DBL_MIN).
const REAL_SMALL: f64 = f64::MIN_POSITIVE;

/// OCCT BRep_Tool::Pnt(V).
fn brep_tool_pnt(v: &Shape) -> glam::DVec3 {
    match v.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => glam::DVec3::ZERO,
    }
}

/// OCCT TopoDS_Shape::IsSame(S): same TShape and same Location.
fn shape_is_same(a: &Shape, b: &Shape) -> bool {
    a.is_same(b)
}

/// OCCT ShapeAnalysis_BoxBndTreeSelector (hxx L33-110 + cxx L29-195).
pub struct ShapeAnalysisBoxBndTreeSelector {
    /// OCCT `myFBox`.
    pub(crate) my_fbox: BndBox,
    /// OCCT `myLBox`.
    pub(crate) my_lbox: BndBox,
    /// OCCT `myShared`.
    pub(crate) my_shared: bool,
    /// OCCT `myNb`.
    pub(crate) my_nb: i32,
    /// OCCT `myFVertex`.
    pub(crate) my_fvertex: Shape,
    /// OCCT `myLVertex`.
    pub(crate) my_lvertex: Shape,
    /// OCCT `myFPnt`.
    pub(crate) my_f_pnt: glam::DVec3,
    /// OCCT `myLPnt`.
    pub(crate) my_l_pnt: glam::DVec3,
    /// OCCT `NCollection_Map<int> myList`.
    pub(crate) my_list: HashSet<i32>,
    /// OCCT `myTol`.
    pub(crate) my_tol: f64,
    /// OCCT `myMin3d`.
    pub(crate) my_min3d: f64,
    /// OCCT `NCollection_Array1<int> myArrIndices(1, 2)`.
    pub(crate) my_arr_indices: [i32; 2],
    /// OCCT `myStatus`.
    pub(crate) my_status: i32,
    /// OCCT `myStop` (the NCollection_UBTree::Selector traversal stop flag).
    pub(crate) my_stop: bool,
}

/// OCCT Accept's local enum { First = 1, Last = 2 } (cxx L51-55).
const SELECTOR_FIRST: usize = 1;
const SELECTOR_LAST: usize = 2;

impl ShapeAnalysisBoxBndTreeSelector {
    /// OCCT ShapeAnalysis_BoxBndTreeSelector(theSeq, theShared)
    /// (hxx L36-47).  Bridge #2: theSeq is passed on Accept instead of
    /// being stored.
    pub fn new(the_shared: bool) -> Self {
        ShapeAnalysisBoxBndTreeSelector {
            my_fbox: BndBox::new(),
            my_lbox: BndBox::new(),
            my_shared: the_shared,
            my_nb: 0,
            my_fvertex: Shape::null(),
            my_lvertex: Shape::null(),
            my_f_pnt: glam::DVec3::ZERO,
            my_l_pnt: glam::DVec3::ZERO,
            my_list: HashSet::new(),
            my_tol: 1e-7,
            my_min3d: 1e-7,
            my_arr_indices: [0; 2],
            my_status: encode_status(ShapeExtendStatus::Ok),
            my_stop: false,
        }
    }

    /// OCCT DefineBoxes (hxx L49-54).
    pub fn define_boxes(&mut self, the_fbox: &BndBox, the_lbox: &BndBox) {
        self.my_fbox = the_fbox.clone();
        self.my_lbox = the_lbox.clone();
        self.my_arr_indices = [0; 2]; // myArrIndices.Init(0)
    }

    /// OCCT DefineVertexes (hxx L56-61).
    pub fn define_vertexes(&mut self, the_vf: &Shape, the_vl: &Shape) {
        self.my_fvertex = the_vf.clone();
        self.my_lvertex = the_vl.clone();
        self.my_status = encode_status(ShapeExtendStatus::Ok);
    }

    /// OCCT DefinePnt (hxx L63-68).
    pub fn define_pnt(&mut self, the_f_pnt: glam::DVec3, the_l_pnt: glam::DVec3) {
        self.my_f_pnt = the_f_pnt;
        self.my_l_pnt = the_l_pnt;
        self.my_status = encode_status(ShapeExtendStatus::Ok);
    }

    /// OCCT GetNb (hxx L70).
    pub fn get_nb(&self) -> i32 {
        self.my_nb
    }

    /// OCCT SetNb (hxx L72).
    pub fn set_nb(&mut self, the_nb: i32) {
        self.my_nb = the_nb;
    }

    /// OCCT LoadList (hxx L74).
    pub fn load_list(&mut self, elem: i32) {
        self.my_list.insert(elem);
    }

    /// OCCT SetStop (hxx L76).
    pub fn set_stop(&mut self) {
        self.my_stop = false;
    }

    /// OCCT SetTolerance (hxx L78-83).
    pub fn set_tolerance(&mut self, the_tol: f64) {
        self.my_tol = the_tol;
        self.my_min3d = the_tol;
        self.my_status = encode_status(ShapeExtendStatus::Ok);
    }

    /// OCCT ContWire (hxx L85).
    pub fn cont_wire(&self, nb_wire: i32) -> bool {
        self.my_list.contains(&nb_wire)
    }

    /// OCCT LastCheckStatus (hxx L87-90).
    pub fn last_check_status(&self, the_status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status, the_status)
    }

    /// OCCT Reject (cxx L29-34).
    pub fn reject(&self, the_bnd: &BndBox) -> bool {
        let fch = self.my_fbox.is_out_box(the_bnd);
        let lch = self.my_lbox.is_out_box(the_bnd);
        fch && lch
    }

    /// OCCT Accept (cxx L38-195).  Bridge #1/#2: the pool and the aliased
    /// arrwires array are passed (the OCCT selector reads them through
    /// stored handles).  The unused_assignments allowance keeps the OCCT
    /// L168-169 `dm1 = dm2` statement, which is never re-read in OCCT
    /// either.
    #[allow(unused_assignments)]
    pub fn accept(&mut self, brep: &mut BRep, arrwires: &[Shape], the_obj: i32) -> bool {
        // OCCT L40-44: the range check (Standard_NoSuchObject; rcad panics).
        if the_obj < 1 || the_obj > arrwires.len() as i32 {
            panic!("ShapeAnalysis_BoxBndTreeSelector::Accept : no such object for current index");
        }
        let mut is_accept = false;
        // OCCT L46-49: myList.Contains(theObj).
        if self.my_list.contains(&the_obj) {
            return false;
        }

        // OCCT L57-59: W = TopoDS::Wire(mySeq->Value(theObj));
        // ShapeAnalysis::FindBounds(W, V1, V2) — the W2 static outputs the
        // null shape for the absent bound (the rcad Option None -> null).
        let w = &arrwires[(the_obj - 1) as usize];
        let mut v1_opt: Option<Shape> = None;
        let mut v2_opt: Option<Shape> = None;
        find_bounds(brep, w, &mut v1_opt, &mut v2_opt);
        let v1 = v1_opt.unwrap_or_else(Shape::null);
        let v2 = v2_opt.unwrap_or_else(Shape::null);

        if self.my_shared {
            // OCCT L60-98: the IsSame chain.
            if shape_is_same(&self.my_lvertex, &v1) {
                self.my_status = encode_status(ShapeExtendStatus::Done1);
                is_accept = true;
                self.my_arr_indices[SELECTOR_LAST - 1] = the_obj;
            } else if shape_is_same(&self.my_lvertex, &v2) {
                self.my_status = encode_status(ShapeExtendStatus::Done2);
                is_accept = true;
                self.my_arr_indices[SELECTOR_LAST - 1] = the_obj;
            } else if shape_is_same(&self.my_fvertex, &v2) {
                self.my_status = encode_status(ShapeExtendStatus::Done3);
                is_accept = true;
                self.my_arr_indices[SELECTOR_FIRST - 1] = the_obj;
            } else if shape_is_same(&self.my_fvertex, &v1) {
                self.my_status = encode_status(ShapeExtendStatus::Done4);
                is_accept = true;
                self.my_arr_indices[SELECTOR_FIRST - 1] = the_obj;
            } else {
                self.my_status = encode_status(ShapeExtendStatus::Fail2);
            }

            // OCCT L100-112.
            if is_accept {
                self.set_nb(the_obj);
                if self.my_arr_indices[SELECTOR_LAST - 1] != 0 {
                    self.my_stop = true;
                }
                return true;
            } else {
                self.my_stop = false;
            }
        } else {
            // OCCT L117-118: p1 = BRep_Tool::Pnt(V1); p2 = BRep_Tool::Pnt(V2).
            let p1 = brep_tool_pnt(&v1);
            let p2 = brep_tool_pnt(&v2);

            // OCCT L120-136: the four distances and the res1/res2 picks.
            let tailhead = p1.distance(self.my_l_pnt);
            let tailtail = p2.distance(self.my_l_pnt);
            let headhead = p1.distance(self.my_f_pnt);
            let headtail = p2.distance(self.my_f_pnt);
            let mut dm1 = tailhead;
            let mut dm2 = headtail;
            let mut res1 = 0;
            let mut res2 = 0;
            if tailhead > tailtail {
                res1 = 1;
                dm1 = tailtail;
            }
            if headtail > headhead {
                res2 = 1;
                dm2 = headhead;
            }
            let mut result = res1;
            let min3d = dm1.min(dm2);
            // OCCT L140-143.
            if min3d > self.my_min3d {
                return false;
            }

            // OCCT L145-151.
            let min_ind = if dm1 > dm2 { SELECTOR_FIRST } else { SELECTOR_LAST };
            let max_ind = if dm1 > dm2 { SELECTOR_LAST } else { SELECTOR_FIRST };
            self.my_arr_indices[min_ind - 1] = the_obj;
            if (min3d - self.my_min3d) > REAL_SMALL {
                self.my_arr_indices[max_ind - 1] = 0;
            }

            // OCCT L153-158.
            self.my_min3d = min3d;
            if min3d > self.my_tol {
                self.my_status = encode_status(ShapeExtendStatus::Fail2);
                return false;
            }

            // OCCT L160-161.
            let an_obj = if self.my_arr_indices[SELECTOR_LAST - 1] != 0 {
                self.my_arr_indices[SELECTOR_LAST - 1]
            } else {
                self.my_arr_indices[SELECTOR_FIRST - 1]
            };
            self.set_nb(an_obj);

            // OCCT L163-166.
            if min3d == 0.0 && min_ind == SELECTOR_LAST {
                self.my_stop = true;
            }

            // OCCT L168-190: the status/result switch.
            if dm1 > dm2 {
                dm1 = dm2;
                result = res2 + 2;
            }
            if an_obj == the_obj {
                match result {
                    0 => self.my_status = encode_status(ShapeExtendStatus::Done1),
                    1 => self.my_status = encode_status(ShapeExtendStatus::Done2),
                    2 => self.my_status = encode_status(ShapeExtendStatus::Done3),
                    3 => self.my_status = encode_status(ShapeExtendStatus::Done4),
                    _ => {}
                }
            }
            return true;
        }

        // OCCT L194.
        false
    }
}

/// OCCT `NCollection_UBTree<int, Bnd_Box>` +
/// `NCollection_UBTreeFiller<int, Bnd_Box>` (TKFoundation, untranslated):
/// the int/Bnd_Box tree instantiation used by the wire-connection flows
/// (FreeBounds.cxx L237-238).  Add/Fill keep the OCCT filler form; Select
/// walks the filled entries in insertion order invoking Reject/Accept and
/// counting the accepted ones — the documented reduction of the OCCT tree
/// traversal (see the module GAP note).
#[allow(dead_code)]
pub struct NCollectionUBTree {
    entries: Vec<(i32, BndBox)>,
}

#[allow(dead_code)]
impl NCollectionUBTree {
    /// OCCT NCollection_UBTree<int, Bnd_Box>().
    pub fn new() -> Self {
        NCollectionUBTree { entries: Vec::new() }
    }
}

/// OCCT NCollection_UBTreeFiller<int, Bnd_Box>(aBBTree) — Add queues the
/// (object, box) pair; Fill() commits them into the tree.
#[allow(dead_code)]
pub struct NCollectionUBTreeFiller {
    tree: NCollectionUBTree,
}

#[allow(dead_code)]
impl NCollectionUBTreeFiller {
    /// OCCT NCollection_UBTreeFiller(aBBTree).
    pub fn new(a_bb_tree: NCollectionUBTree) -> Self {
        NCollectionUBTreeFiller { tree: a_bb_tree }
    }

    /// OCCT Add(theObj, theBnd).
    pub fn add(&mut self, the_obj: i32, the_bnd: BndBox) {
        self.tree.entries.push((the_obj, the_bnd));
    }

    /// OCCT Fill() — commits the queued pairs (the carrier tree stores them
    /// in insertion order).
    pub fn fill(&mut self) {}

    /// OCCT aBBTree.Select(aSel): the tree traversal invoking Reject/Accept;
    /// returns the number of accepted objects.
    pub fn select(
        &mut self,
        brep: &mut BRep,
        a_sel: &mut ShapeAnalysisBoxBndTreeSelector,
        arrwires: &[Shape],
    ) -> i32 {
        let mut nsel = 0;
        for (the_obj, the_bnd) in &self.tree.entries {
            if a_sel.reject(the_bnd) {
                continue;
            }
            if a_sel.accept(brep, arrwires, *the_obj) {
                nsel += 1;
            }
        }
        nsel
    }
}
