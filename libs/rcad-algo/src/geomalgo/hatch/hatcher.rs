// OCCT Geom2dHatch_Hatcher (TKGeomAlgo/Geom2dHatch/Geom2dHatch_Hatcher.hxx
// L1-209 + .cxx L1-1919 + .lxx L1-257) — the hatching main engine over the
// landed leaf layer (Elements / Intersector / Hatching / Classifier).

use glam::{DVec2, DVec3};
use rcad_kernel::geom::Curve2d;
use rcad_kernel::precision::{ANGULAR, CONFUSION};
use rcad_kernel::topods::{Orientation, State};

use crate::geomalgo::geom2d_int::Curve2dAdaptor;
use crate::geomalgo::hatch::classifier::Classifier;
use crate::geomalgo::hatch::element::HatchElement;
use crate::geomalgo::hatch::elements::HatchElements;
use crate::geomalgo::hatch::hatch_gen::{
    Domain as HatchGenDomain, ErrorStatus, IntersectionType, PointOnElement, PointOnHatching,
};
use crate::geomalgo::hatch::hatching::Hatching;
use crate::geomalgo::hatch::intersector::HatchIntersector;
use crate::geomalgo::int_res2d::{Position, TypeTrans};
use crate::geomalgo::top_trans::CurveTransition;

/// OCCT NCollection_DataMap<int, Geom2dHatch_Hatching> (hxx L204) — the
/// order-preserving Vec<(usize, Hatching)> precedent of Geom2dHatch_Elements;
/// the DataMap Bind/UnBind/Find/ChangeFind semantics are kept.
#[derive(Debug)]
struct HatchingsMap {
    my_map: Vec<(usize, Hatching)>,
}

impl HatchingsMap {
    fn new() -> Self {
        HatchingsMap { my_map: Vec::new() }
    }

    /// OCCT NCollection_DataMap::Clear.
    fn clear(&mut self) {
        self.my_map.clear();
    }

    /// OCCT NCollection_DataMap::IsBound.
    fn is_bound(&self, k: usize) -> bool {
        self.my_map.iter().any(|(key, _)| *key == k)
    }

    /// OCCT NCollection_DataMap::UnBind.
    fn un_bind(&mut self, k: usize) -> bool {
        if let Some(pos) = self.my_map.iter().position(|(key, _)| *key == k) {
            self.my_map.remove(pos);
            true
        } else {
            false
        }
    }

    /// OCCT NCollection_DataMap::Bind — an already bound key gets its item
    /// overridden and returns false.
    fn bind(&mut self, k: usize, hatching: Hatching) -> bool {
        if let Some(pos) = self.my_map.iter().position(|(key, _)| *key == k) {
            self.my_map[pos].1 = hatching;
            return false;
        }
        self.my_map.push((k, hatching));
        true
    }

    /// OCCT NCollection_DataMap::Find — raises NoSuchObject when unbound.
    fn find(&self, k: usize) -> &Hatching {
        self.my_map
            .iter()
            .find(|(key, _)| *key == k)
            .map(|(_, h)| h)
            .expect("Geom2dHatch_Hatcher myHatchings Find - Standard_NoSuchObject")
    }

    /// OCCT NCollection_DataMap::ChangeFind — raises NoSuchObject when
    /// unbound.
    fn change_find(&mut self, k: usize) -> &mut Hatching {
        self.my_map
            .iter_mut()
            .find(|(key, _)| *key == k)
            .map(|(_, h)| h)
            .expect("Geom2dHatch_Hatcher myHatchings ChangeFind - Standard_NoSuchObject")
    }
}

/// OCCT Geom2dHatch_Hatcher — the hatching engine: it trims hatchings by the
/// boundary elements, computes the global transitions and the trimmed
/// domains.
pub struct Hatcher {
    /// OCCT Geom2dHatch_Intersector myIntersector (hxx L196) — a value
    /// member; OCCT copies it on construction/assignment, rcad transfers the
    /// ownership (no Clone is required on the GInter base).
    my_intersector: HatchIntersector,
    /// OCCT double myConfusion2d (hxx L197).
    my_confusion2d: f64,
    /// OCCT double myConfusion3d (hxx L198).
    my_confusion3d: f64,
    /// OCCT bool myKeepPoints (hxx L199).
    my_keep_points: bool,
    /// OCCT bool myKeepSegments (hxx L200).
    my_keep_segments: bool,
    /// OCCT int myNbElements (hxx L201).
    my_nb_elements: usize,
    /// OCCT Geom2dHatch_Elements myElements (hxx L202).
    my_elements: HatchElements,
    /// OCCT int myNbHatchings (hxx L203).
    my_nb_hatchings: usize,
    /// OCCT NCollection_DataMap<int, Geom2dHatch_Hatching> myHatchings (hxx
    /// L204).
    my_hatchings: HatchingsMap,
}

// ===========================================================================
//  OCCT Geom2dHatch_Hatcher.cxx L37-58 — Category : General use.
// ===========================================================================

impl Hatcher {
    /// OCCT Geom2dHatch_Hatcher::Geom2dHatch_Hatcher (cxx L45-58) — returns
    /// an empty hatcher.  The OCCT default arguments are KeepPnt = false,
    /// KeepSeg = false (hxx L48-49); rcad has no default arguments so they
    /// are passed explicitly.
    pub fn new(
        intersector: HatchIntersector,
        confusion2d: f64,
        confusion3d: f64,
        keep_pnt: bool,
        keep_seg: bool,
    ) -> Self {
        Hatcher {
            my_intersector: intersector,
            my_confusion2d: confusion2d,
            my_confusion3d: confusion3d,
            my_keep_points: keep_pnt,
            my_keep_segments: keep_seg,
            my_nb_elements: 0,
            my_elements: HatchElements::new(),
            my_nb_hatchings: 0,
            my_hatchings: HatchingsMap::new(),
        }
    }

    // -- Geom2dHatch_Hatcher.lxx inline methods ------------------------------

    /// OCCT Geom2dHatch_Hatcher::Intersector (cxx L65-76) — sets the
    /// associated intersector and clears the intersection points of all the
    /// hatchings.
    pub fn set_intersector(&mut self, intersector: HatchIntersector) {
        self.my_intersector = intersector;
        for indh in 1..=self.my_nb_hatchings {
            if self.my_hatchings.is_bound(indh) {
                let hatching = self.my_hatchings.change_find(indh);
                hatching.clr_points();
            }
        }
    }

    /// OCCT Geom2dHatch_Hatcher::Intersector (lxx L32-35) — returns the
    /// associated intersector.
    pub fn intersector(&self) -> &HatchIntersector {
        &self.my_intersector
    }

    /// OCCT Geom2dHatch_Hatcher::ChangeIntersector (lxx L42-45).
    pub fn change_intersector(&mut self) -> &mut HatchIntersector {
        &mut self.my_intersector
    }

    /// OCCT Geom2dHatch_Hatcher::Confusion2d (cxx L83-94) — sets the 2d
    /// confusion tolerance and clears the intersection points of all the
    /// hatchings.
    pub fn set_confusion2d(&mut self, confusion: f64) {
        self.my_confusion2d = confusion;
        for indh in 1..=self.my_nb_hatchings {
            if self.my_hatchings.is_bound(indh) {
                let hatching = self.my_hatchings.change_find(indh);
                hatching.clr_points();
            }
        }
    }

    /// OCCT Geom2dHatch_Hatcher::Confusion2d (lxx L52-55) — returns the 2d
    /// confusion tolerance.
    pub fn confusion2d(&self) -> f64 {
        self.my_confusion2d
    }

    /// OCCT Geom2dHatch_Hatcher::Confusion3d (cxx L101-112) — sets the 3d
    /// confusion tolerance and clears the intersection points of all the
    /// hatchings.
    pub fn set_confusion3d(&mut self, confusion: f64) {
        self.my_confusion3d = confusion;
        for indh in 1..=self.my_nb_hatchings {
            if self.my_hatchings.is_bound(indh) {
                let hatching = self.my_hatchings.change_find(indh);
                hatching.clr_points();
            }
        }
    }

    /// OCCT Geom2dHatch_Hatcher::Confusion3d (lxx L62-65) — returns the 3d
    /// confusion tolerance.
    pub fn confusion3d(&self) -> f64 {
        self.my_confusion3d
    }

    /// OCCT Geom2dHatch_Hatcher::KeepPoints (cxx L116-127) — sets the points
    /// consideration flag and clears the domains of all the hatchings.
    pub fn set_keep_points(&mut self, keep: bool) {
        self.my_keep_points = keep;
        for indh in 1..=self.my_nb_hatchings {
            if self.my_hatchings.is_bound(indh) {
                let hatching = self.my_hatchings.change_find(indh);
                hatching.clr_domains();
            }
        }
    }

    /// OCCT Geom2dHatch_Hatcher::KeepPoints (lxx L72-75).
    pub fn keep_points(&self) -> bool {
        self.my_keep_points
    }

    /// OCCT Geom2dHatch_Hatcher::KeepSegments (cxx L131-142) — sets the
    /// segments consideration flag and clears the domains of all the
    /// hatchings.
    pub fn set_keep_segments(&mut self, keep: bool) {
        self.my_keep_segments = keep;
        for indh in 1..=self.my_nb_hatchings {
            if self.my_hatchings.is_bound(indh) {
                let hatching = self.my_hatchings.change_find(indh);
                hatching.clr_domains();
            }
        }
    }

    /// OCCT Geom2dHatch_Hatcher::KeepSegments (lxx L82-85).
    pub fn keep_segments(&self) -> bool {
        self.my_keep_segments
    }

    /// OCCT Geom2dHatch_Hatcher::Clear (lxx L92-98) — removes all the
    /// hatchings and all the elements.
    pub fn clear(&mut self) {
        if self.my_nb_hatchings != 0 {
            self.clr_hatchings();
        }
        if self.my_nb_elements != 0 {
            self.clr_elements();
        }
    }

    /// OCCT Geom2dHatch_Hatcher::Element (lxx L102-109) — returns the IndE-th
    /// element (the OCCT protected accessor; Rust has no protected section).
    pub fn element(&mut self, inde: usize) -> &mut HatchElement {
        self.my_elements.change_find(inde)
    }

    /// OCCT Geom2dHatch_Hatcher::ElementCurve (lxx L116-123) — returns the
    /// curve associated to the IndE-th element.
    pub fn element_curve(&self, inde: usize) -> &Curve2d {
        let element = self.my_elements.find(inde);
        element.curve()
    }

    /// OCCT Geom2dHatch_Hatcher::Hatching (lxx L127-134) — returns the IndH-th
    /// hatching (the OCCT protected accessor).
    pub fn hatching(&mut self, indh: usize) -> &mut Hatching {
        self.my_hatchings.change_find(indh)
    }

    /// OCCT Geom2dHatch_Hatcher::HatchingCurve (lxx L141-148) — returns the
    /// curve associated to the IndH-th hatching.
    pub fn hatching_curve(&self, indh: usize) -> &Curve2d {
        let hatching = self.my_hatchings.find(indh);
        hatching.curve()
    }

    /// OCCT Geom2dHatch_Hatcher::NbPoints (lxx L156-163) — returns the number
    /// of intersection points of the IndH-th hatching.
    pub fn nb_points(&self, indh: usize) -> usize {
        let hatching = self.my_hatchings.find(indh);
        hatching.nb_points()
    }

    /// OCCT Geom2dHatch_Hatcher::Point (lxx L171-183) — returns the IndP-th
    /// intersection point of the IndH-th hatching.
    pub fn point(&self, indh: usize, indp: usize) -> &PointOnHatching {
        let hatching = self.my_hatchings.find(indh);
        hatching.point(indp)
    }

    /// OCCT Geom2dHatch_Hatcher::TrimDone (lxx L191-198).
    pub fn trim_done(&self, indh: usize) -> bool {
        let hatching = self.my_hatchings.find(indh);
        hatching.trim_done()
    }

    /// OCCT Geom2dHatch_Hatcher::TrimFailed (lxx L206-213).
    pub fn trim_failed(&self, indh: usize) -> bool {
        let hatching = self.my_hatchings.find(indh);
        hatching.trim_failed()
    }

    /// OCCT Geom2dHatch_Hatcher::IsDone (lxx L221-228) — returns the fact
    /// that the domains were computed for the IndH-th hatching.  (The no-arg
    /// IsDone() declared at hxx L160 has no definition anywhere in OCCT — a
    /// dead declaration; rcad cannot translate a missing body.)
    pub fn is_done(&self, indh: usize) -> bool {
        let hatching = self.my_hatchings.find(indh);
        hatching.is_done()
    }

    /// OCCT Geom2dHatch_Hatcher::Status (lxx L235-242) — returns the status
    /// about the IndH-th hatching.
    pub fn status(&self, indh: usize) -> ErrorStatus {
        let hatching = self.my_hatchings.find(indh);
        hatching.status()
    }

    /// OCCT Geom2dHatch_Hatcher::NbDomains (lxx L249-257) — returns the
    /// number of domains of the IndH-th hatching; raises StdFail_NotDone when
    /// the domains were not computed.
    pub fn nb_domains(&self, indh: usize) -> usize {
        let hatching = self.my_hatchings.find(indh);
        if !hatching.is_done() {
            panic!("Geom2dHatch_Hatcher::NbDomains");
        }
        hatching.nb_domains()
    }

    // ===========================================================================
    //  OCCT Geom2dHatch_Hatcher.cxx L144-179 — Category : Element.
    // ===========================================================================

    /// OCCT Geom2dHatch_Hatcher::AddElement (cxx L155-179) — adds an element
    /// to the hatcher and returns its index.
    pub fn add_element(&mut self, curve: &Curve2d, orientation: Orientation) -> usize {
        // for (IndE = 1; IndE <= myNbElements && myElements.IsBound(IndE); IndE++) { ; }
        let mut inde = 1usize;
        while inde <= self.my_nb_elements && self.my_elements.is_bound(inde) {
            inde += 1;
        }
        if inde > self.my_nb_elements {
            self.my_nb_elements += 1;
            inde = self.my_nb_elements;
        }
        let element = HatchElement::new(curve.clone(), orientation);
        self.my_elements.bind(inde, element);
        for indh in 1..=self.my_nb_hatchings {
            if self.my_hatchings.is_bound(indh) {
                let hatching = self.my_hatchings.change_find(indh);
                hatching.clr_points();
            }
        }
        inde
    }

    /// OCCT Geom2dHatch_Hatcher::AddElement (hxx L98-104) — the
    /// handle<Geom2d_Curve> overload: wraps the curve in an adaptor and
    /// delegates to the adaptor overload.  rcad `Curve2d` is its own adaptor,
    /// so the wrap is a copy.
    pub fn add_element_from_curve(&mut self, curve: &Curve2d, orientation: Orientation) -> usize {
        // Geom2dAdaptor_Curve aGAC(Curve);
        // return AddElement(aGAC, Orientation);
        let a_gac = curve.clone();
        self.add_element(&a_gac, orientation)
    }

    /// OCCT Geom2dHatch_Hatcher::RemElement (cxx L186-224) — removes the
    /// IndE-th element from the hatcher.
    pub fn rem_element(&mut self, inde: usize) {
        for indh in 1..=self.my_nb_hatchings {
            if self.my_hatchings.is_bound(indh) {
                let hatching = self.my_hatchings.change_find(indh);
                let mut domains_to_clear = false;
                // for (int IPntH = Hatching.NbPoints(); IPntH > 0; IPntH--)
                let nb_pnt_h = hatching.nb_points();
                for ipnt_h in (1..=nb_pnt_h).rev() {
                    // OCCT L199: `HatchGen_PointOnHatching PntH =
                    // Hatching.ChangePoint(IPntH);` — copy-initialization
                    // from the reference: the point is a COPY, so the
                    // RemPoint below edits the copy (OCCT form kept).
                    let mut pnt_h = hatching.change_point(ipnt_h).clone();
                    // for (int IPntE = PntH.NbPoints(); IPntE > 0; IPntE--)
                    let nb_pnt_e = pnt_h.nb_points();
                    for ipnt_e in (1..=nb_pnt_e).rev() {
                        if pnt_h.point(ipnt_e).index() == inde as i32 {
                            pnt_h.rem_point(ipnt_e);
                            domains_to_clear = true;
                        }
                    }
                    if pnt_h.nb_points() == 0 {
                        hatching.rem_point(ipnt_h);
                    }
                }
                if domains_to_clear {
                    hatching.clr_domains();
                }
            }
        }
        self.my_elements.un_bind(inde);
        if inde == self.my_nb_elements {
            self.my_nb_elements -= 1;
        }
    }

    /// OCCT Geom2dHatch_Hatcher::ClrElements (cxx L231-249) — removes all the
    /// elements from the hatcher.
    pub fn clr_elements(&mut self) {
        if self.my_nb_elements != 0 {
            if self.my_nb_hatchings != 0 {
                for indh in 1..=self.my_nb_hatchings {
                    if self.my_hatchings.is_bound(indh) {
                        let hatching = self.my_hatchings.change_find(indh);
                        hatching.clr_points();
                    }
                }
            }
            self.my_elements.clear();
            self.my_nb_elements = 0;
        }
    }

    // ===========================================================================
    //  OCCT Geom2dHatch_Hatcher.cxx L251-318 — Category : Hatching.
    // ===========================================================================

    /// OCCT Geom2dHatch_Hatcher::AddHatching (cxx L262-277) — adds a hatching
    /// to the hatcher and returns its index.
    pub fn add_hatching(&mut self, curve: &Curve2d) -> usize {
        // for (IndH = 1; IndH <= myNbHatchings && myHatchings.IsBound(IndH); IndH++) { ; }
        let mut indh = 1usize;
        while indh <= self.my_nb_hatchings && self.my_hatchings.is_bound(indh) {
            indh += 1;
        }
        if indh > self.my_nb_hatchings {
            self.my_nb_hatchings += 1;
            indh = self.my_nb_hatchings;
        }
        let hatching = Hatching::new(curve.clone());
        self.my_hatchings.bind(indh, hatching);
        indh
    }

    /// OCCT Geom2dHatch_Hatcher::RemHatching (cxx L284-296) — removes the
    /// IndH-th hatching from the hatcher.
    pub fn rem_hatching(&mut self, indh: usize) {
        self.my_hatchings.change_find(indh).clr_points();
        self.my_hatchings.un_bind(indh);
        if indh == self.my_nb_hatchings {
            self.my_nb_hatchings -= 1;
        }
    }

    /// OCCT Geom2dHatch_Hatcher::ClrHatchings (cxx L303-318) — removes all
    /// the hatchings from the hatcher.
    pub fn clr_hatchings(&mut self) {
        if self.my_nb_hatchings != 0 {
            for indh in 1..=self.my_nb_hatchings {
                if self.my_hatchings.is_bound(indh) {
                    let hatching = self.my_hatchings.change_find(indh);
                    hatching.clr_points();
                }
            }
            self.my_hatchings.clear();
            self.my_nb_hatchings = 0;
        }
    }

    // ===========================================================================
    //  OCCT Geom2dHatch_Hatcher.cxx L320-840 — Category : Computation -
    //  Trimming.
    // ===========================================================================

    /// OCCT Geom2dHatch_Hatcher::Trim (cxx L332-341) — trims all the
    /// hatchings of the hatcher by all the elements of the hatcher.
    pub fn trim(&mut self) {
        for indh in 1..=self.my_nb_hatchings {
            if self.my_hatchings.is_bound(indh) {
                self.trim_hatching(indh);
            }
        }
    }

    /// OCCT Geom2dHatch_Hatcher::Trim (cxx L349-354) — adds a hatching to the
    /// hatcher and trims it by the elements already given and returns its
    /// index.
    pub fn trim_curve(&mut self, curve: &Curve2d) -> usize {
        let indh = self.add_hatching(curve);
        self.trim_hatching(indh);
        indh
    }

    /// OCCT Geom2dHatch_Hatcher::Trim (cxx L361-395) — trims the IndH-th
    /// hatching by the elements already given.
    pub fn trim_hatching(&mut self, indh: usize) {
        self.my_hatchings.change_find(indh).clr_points();

        // bool OK, AllOK; AllOK = true;
        let mut all_ok = true;
        let mut ok;
        for inde in 1..=self.my_nb_elements {
            if self.my_elements.is_bound(inde) {
                ok = self.trim_by_element(indh, inde);
                all_ok = all_ok && ok;
            }
        }
        {
            let hatching = self.my_hatchings.change_find(indh);
            hatching.set_trim_done(true);
            hatching.set_trim_failed(!all_ok);
        }

        if all_ok {
            // The OCCT loop re-evaluates Hatching.NbPoints() at each step;
            // GlobalTransition never changes the point count.  Rust cannot
            // pass a reference into myHatchings to GlobalTransition while it
            // also reads self, so the point is cloned, processed, and written
            // back on success (GlobalTransition mutates the point only on its
            // success path).
            let mut ipnt = 1usize;
            while ipnt <= self.my_hatchings.find(indh).nb_points() {
                let mut pnt_h = self.my_hatchings.find(indh).point(ipnt).clone();
                let ok = self.global_transition(&mut pnt_h);
                all_ok = all_ok && ok;
                if ok {
                    *self.my_hatchings.change_find(indh).change_point(ipnt) = pnt_h;
                }
                ipnt += 1;
            }
            // Hatching.Status(AllOK ? HatchGen_NoProblem : HatchGen_TransitionFailure);
            let status = if all_ok {
                ErrorStatus::NoProblem
            } else {
                ErrorStatus::TransitionFailure
            };
            self.my_hatchings.change_find(indh).set_status(status);
        }
    }

    /// OCCT Geom2dHatch_Hatcher::Trim (cxx L522-840) — trims the IndH-th
    /// hatching of the hatcher by the IndE-th element.
    fn trim_by_element(&mut self, indh: usize, inde: usize) -> bool {
        let hatching = self.my_hatchings.change_find(indh);
        let element = self.my_elements.change_find(inde);

        // Geom2dAdaptor_Curve hatching = Hatching.ChangeCurve();
        // Geom2dAdaptor_Curve element  = Element.ChangeCurve();
        let hatching_curve = hatching.change_curve().clone();
        let element_curve = element.change_curve().clone();

        self.my_intersector.intersect(&hatching_curve, &element_curve);

        if !self.my_intersector.is_done() {
            print!(" Intersector -> Done = False ");
            return false;
        }

        if intersector_is_empty(&self.my_intersector) {
            return true;
        }

        //-----------------------------------------------------------------------
        // Traitement des points d intersection.
        //-----------------------------------------------------------------------

        for ipnt_i in 1..=self.my_intersector.nb_points() {
            let pnt_i = self.my_intersector.point(ipnt_i);

            let mut pnt_e = PointOnElement::new_intersection(pnt_i);
            pnt_e.set_index(inde as i32);

            let mut pnt_h = PointOnHatching::new_intersection(pnt_i);
            pnt_h.set_index(indh as i32);
            pnt_h.add_point(&pnt_e, self.my_confusion2d);

            hatching.add_point(&pnt_h, self.my_confusion2d);
        }

        //-----------------------------------------------------------------------
        // Traitement des segments d intersection.
        //-----------------------------------------------------------------------

        for iseg in 1..=self.my_intersector.nb_segments() {
            let seg = self.my_intersector.segment(iseg);

            let first_point = seg.has_first_point();
            let last_point = seg.has_last_point();

            //-----------------------------------------------------------------------
            // Les deux points peuvent etre confondus.
            //-----------------------------------------------------------------------

            if first_point && last_point {
                let pnt1 = seg.first_point();
                let pnt2 = seg.last_point();

                let trs_pnt1_h = pnt1.transition_of_first();
                let trs_pnt1_e = pnt1.transition_of_second();
                let trs_pnt2_h = pnt2.transition_of_first();
                let trs_pnt2_e = pnt2.transition_of_second();

                let type_pnt1_h = trs_pnt1_h.transition_type();
                let type_pnt1_e = trs_pnt1_e.transition_type();
                let type_pnt2_h = trs_pnt2_h.transition_type();
                let type_pnt2_e = trs_pnt2_e.transition_type();

                //-----------------------------------------------------------------------
                // Les deux points peuvent etre confondus au regard de la precision du
                // `hatcher'.
                //-----------------------------------------------------------------------

                let conf_2d = (pnt1.param_on_first() - pnt2.param_on_first()).abs()
                    <= self.my_confusion2d;

                //-----------------------------------------------------------------------
                // Les deux points peuvent etre `confondus' au regard des intersections.
                //-----------------------------------------------------------------------

                let mut conf_3d = false;

                if !conf_2d {
                    conf_3d = true;
                    if conf_3d {
                        conf_3d =
                            type_pnt1_h != TypeTrans::Touch && type_pnt1_h != TypeTrans::Undecided;
                    }
                    if conf_3d {
                        conf_3d =
                            type_pnt1_e != TypeTrans::Touch && type_pnt1_e != TypeTrans::Undecided;
                    }
                    if conf_3d {
                        conf_3d =
                            type_pnt2_h != TypeTrans::Touch && type_pnt2_h != TypeTrans::Undecided;
                    }
                    if conf_3d {
                        conf_3d =
                            type_pnt2_e != TypeTrans::Touch && type_pnt2_e != TypeTrans::Undecided;
                    }
                    if conf_3d {
                        conf_3d = type_pnt1_h == type_pnt2_h && type_pnt1_e == type_pnt2_e;
                    }
                    if conf_3d {
                        conf_3d = pnt1.value().distance(pnt2.value()) <= self.my_confusion3d;
                    }
                }

                if conf_2d || conf_3d {
                    let mut pnt_e = PointOnElement::new();
                    pnt_e.set_index(inde as i32);
                    pnt_e.set_parameter((pnt1.param_on_second() + pnt2.param_on_second()) / 2.);
                    match trs_pnt1_e.position_on_curve() {
                        Position::Head => {
                            pnt_e.set_position(Orientation::Forward);
                        }
                        Position::Middle => {
                            match trs_pnt2_e.position_on_curve() {
                                Position::Head => {
                                    pnt_e.set_position(Orientation::Forward);
                                }
                                Position::Middle => {
                                    pnt_e.set_position(Orientation::Internal);
                                }
                                Position::End => {
                                    pnt_e.set_position(Orientation::Reversed);
                                }
                            }
                        }
                        Position::End => {
                            pnt_e.set_position(Orientation::Reversed);
                        }
                    }
                    pnt_e.set_intersection_type(if pnt_e.position() == Orientation::Internal {
                        IntersectionType::True
                    } else {
                        IntersectionType::Touch
                    });
                    pnt_e.set_state_before(if type_pnt1_h == TypeTrans::In {
                        State::Out
                    } else {
                        State::In
                    });
                    pnt_e.set_state_after(if type_pnt2_h == TypeTrans::In {
                        State::Out
                    } else {
                        State::In
                    });

                    let mut pnt_h = PointOnHatching::new();
                    pnt_h.set_index(indh as i32);
                    pnt_h.set_parameter((pnt1.param_on_first() + pnt2.param_on_first()) / 2.);
                    match trs_pnt1_h.position_on_curve() {
                        Position::Head => {
                            pnt_h.set_position(Orientation::Forward);
                        }
                        Position::Middle => {
                            match trs_pnt2_h.position_on_curve() {
                                Position::Head => {
                                    pnt_h.set_position(Orientation::Forward);
                                }
                                Position::Middle => {
                                    pnt_h.set_position(Orientation::Internal);
                                }
                                Position::End => {
                                    pnt_h.set_position(Orientation::Reversed);
                                }
                            }
                        }
                        Position::End => {
                            pnt_h.set_position(Orientation::Reversed);
                        }
                    }

                    pnt_h.add_point(&pnt_e, self.my_confusion2d);
                    hatching.add_point(&pnt_h, self.my_confusion2d);

                    continue;
                }

                //-----------------------------------------------------------------------
                // Traitement du premier point du segment.
                //-----------------------------------------------------------------------

                if first_point {
                    let pnt_i = seg.first_point();

                    let mut pnt_e = PointOnElement::new_intersection(pnt_i);
                    pnt_e.set_index(inde as i32);
                    pnt_e.set_segment_beginning(true);
                    pnt_e.set_segment_end(false);

                    let mut pnt_h = PointOnHatching::new_intersection(pnt_i);
                    pnt_h.set_index(indh as i32);
                    pnt_h.add_point(&pnt_e, self.my_confusion2d);

                    hatching.add_point(&pnt_h, self.my_confusion2d);
                }

                //-----------------------------------------------------------------------
                // Traitement du deuxieme point du segment.
                //-----------------------------------------------------------------------

                if last_point {
                    let pnt_i = seg.last_point();

                    let mut pnt_e = PointOnElement::new_intersection(pnt_i);
                    pnt_e.set_index(inde as i32);
                    pnt_e.set_segment_beginning(false);
                    pnt_e.set_segment_end(true);

                    let mut pnt_h = PointOnHatching::new_intersection(pnt_i);
                    pnt_h.set_index(indh as i32);
                    pnt_h.add_point(&pnt_e, self.my_confusion2d);

                    hatching.add_point(&pnt_h, self.my_confusion2d);
                }
            }
        }
        true
    }

    // ===========================================================================
    //  OCCT Geom2dHatch_Hatcher.cxx L842-1121 — Category : Computation -
    //  Domains.
    // ===========================================================================

    /// OCCT Geom2dHatch_Hatcher::GlobalTransition (cxx L855-1121) — returns
    /// the before and after states of the complex transition of the point
    /// (the point is owned by the caller; GlobalTransition reads
    /// myIntersector / myElements / myHatchings immutably, which Rust only
    /// allows while the point is not borrowed from them).
    #[allow(unused_assignments)]
    fn global_transition(&self, point: &mut PointOnHatching) -> bool {
        // TopAbs_State StateBefore = TopAbs_UNKNOWN;
        let mut state_before = State::Unknown;
        // TopAbs_State StateAfter = TopAbs_UNKNOWN;
        let mut state_after = State::Unknown;
        // bool SegmentBegin = false;
        let mut segment_begin = false;
        // bool SegmentEnd = false;
        let mut segment_end = false;

        // gp_Dir2d Tangente2d, Normale2d — uninitialized in OCCT; rcad uses
        // neutral values so the code stays deterministic.
        let mut tangente2d = DVec2::X;
        let mut normale2d = DVec2::ZERO;
        // double Courbure — uninitialized in OCCT.
        let mut courbure = 0.0;

        // const Geom2dAdaptor_Curve& CurveH = HatchingCurve(Point.Index());
        let curve_h = self.hatching_curve(point.index() as usize);

        // myIntersector.LocalGeometry(CurveH.Curve(), Point.Parameter(),
        //                             Tangente2d, Normale2d, Courbure);
        // (rcad Curve2d folds the adaptor and its geometry: CurveH.Curve()
        // is CurveH itself.)
        (tangente2d, normale2d, courbure) =
            self.my_intersector
                .local_geometry(curve_h, point.parameter());

        // Tangente.SetCoord(Tangente2d.X(), Tangente2d.Y(), 0.0);
        let mut tangente = DVec3::new(tangente2d.x, tangente2d.y, 0.0);
        let mut normale;
        if courbure < CONFUSION {
            // Normale.SetCoord(-Tangente2d.Y(), Tangente2d.X(), 0.0);
            normale = DVec3::new(-tangente2d.y, tangente2d.x, 0.0);
        } else {
            // Normale.SetCoord(Normale2d.X(), Normale2d.Y(), 0.0);
            normale = DVec3::new(normale2d.x, normale2d.y, 0.0);
        }

        // TopTrans_CurveTransition ComplexTransition;
        // ComplexTransition.Reset(Tangente, Normale, Courbure);
        let mut complex_transition = CurveTransition::<DVec3>::new();
        complex_transition.reset_3d(tangente, normale, courbure);

        for ipnte in 1..=point.nb_points() {
            let pnte = point.point(ipnte);

            segment_begin = segment_begin || pnte.segment_beginning();
            segment_end = segment_end || pnte.segment_end();

            // const Geom2dHatch_Element& Element = myElements.Find(PntE.Index());
            let element = self.my_elements.find(pnte.index() as usize);
            // const Geom2dAdaptor_Curve& CurveE = Element.Curve();
            let curve_e = element.curve();

            // TopAbs_Orientation ElementOrientation = Element.Orientation();
            let element_orientation = element.orientation();
            // bool ToReverse = (ElementOrientation == TopAbs_REVERSED);
            let to_reverse = element_orientation == Orientation::Reversed;
            // double Param;
            let mut param;
            match pnte.position() {
                Orientation::Forward => {
                    param = if to_reverse {
                        curve_e.last_parameter()
                    } else {
                        curve_e.first_parameter()
                    };
                }
                Orientation::Internal => {
                    param = pnte.parameter();
                }
                Orientation::Reversed => {
                    param = if to_reverse {
                        curve_e.first_parameter()
                    } else {
                        curve_e.last_parameter()
                    };
                }
                Orientation::External => {}
            }

            //--
            // OCCT tree note (finding #25): this overwrites the Param
            // computed by the switch above, making the position-based logic
            // dead code; removing it caused blend regressions.  The switch
            // stays for the OCCT form.
            param = pnte.parameter();

            // myIntersector.LocalGeometry(CurveE.Curve(), Param,
            //                             Tangente2d, Normale2d, Courbure);
            (tangente2d, normale2d, courbure) =
                self.my_intersector.local_geometry(curve_e, param);

            //-----------------------------------------------------------------------
            // Calcul de la transition locale. On suppose les relations suivantes :
            //  - Si l orientation de l element est INTERNAL ==> INTERNAL
            //  - Si l orientation de l element est EXTERNAL ==> EXTERNAL
            //  - Si tangence, on a IN-IN  ou OUT-OUT ==> INTERNAL/EXTERNAL
            //  - Sinon,       on a IN-OUT ou OUT-IN  ==> REVERSED/FORWARD
            // Les deux dernieres conditions avec l element vu en FORWARD.
            //-----------------------------------------------------------------------
            let mut local_transition = Orientation::External;

            if element_orientation == Orientation::Internal {
                local_transition = Orientation::Internal;
            } else if element_orientation == Orientation::External {
                local_transition = Orientation::External;
            } else if pnte.intersection_type() == IntersectionType::Tangent {
                if pnte.position() == Orientation::Internal {
                    match pnte.state_before() {
                        State::In => {
                            local_transition = if to_reverse {
                                Orientation::External
                            } else {
                                Orientation::Internal
                            };
                        }
                        State::Out => {
                            local_transition = if to_reverse {
                                Orientation::Internal
                            } else {
                                Orientation::External
                            };
                        }
                        _ => {}
                    }
                } else {
                    match pnte.state_before() {
                        State::In => {
                            local_transition = if to_reverse {
                                Orientation::Forward
                            } else {
                                Orientation::Reversed
                            };
                        }
                        State::Out => {
                            local_transition = if to_reverse {
                                Orientation::Reversed
                            } else {
                                Orientation::Forward
                            };
                        }
                        _ => {}
                    }
                }
            } else {
                match pnte.state_before() {
                    State::In => {
                        local_transition = if to_reverse {
                            Orientation::Forward
                        } else {
                            Orientation::Reversed
                        };
                    }
                    State::Out => {
                        local_transition = if to_reverse {
                            Orientation::Reversed
                        } else {
                            Orientation::Forward
                        };
                    }
                    _ => {}
                }
            }

            //-----------------------------------------------------------------------
            // Orientation de la tangente au point d interference.
            //-----------------------------------------------------------------------
            let mut tangente_orientation = Orientation::Forward;
            match pnte.position() {
                Orientation::Forward => {
                    tangente_orientation = if to_reverse {
                        Orientation::Reversed
                    } else {
                        Orientation::Forward
                    };
                }
                Orientation::Internal => {
                    tangente_orientation = Orientation::Internal;
                }
                Orientation::Reversed => {
                    tangente_orientation = if to_reverse {
                        Orientation::Forward
                    } else {
                        Orientation::Reversed
                    };
                }
                Orientation::External => {}
            }

            //-----------------------------------------------------------------------
            // Proprietes geometriques.
            //-----------------------------------------------------------------------

            if to_reverse {
                tangente = DVec3::new(-tangente2d.x, -tangente2d.y, 0.0);
            } else {
                tangente = DVec3::new(tangente2d.x, tangente2d.y, 0.0);
            }
            normale = DVec3::new(normale2d.x, normale2d.y, 0.0);

            // ComplexTransition.Compare(Precision::Angular(), Tangente,
            //                           Normale, Courbure, LocalTransition,
            //                           TangenteOrientation);
            complex_transition.compare(
                ANGULAR,
                tangente,
                normale,
                courbure,
                local_transition,
                tangente_orientation,
            );
        }

        match complex_transition.state_before() {
            State::In => state_before = State::In,
            State::Out => state_before = State::Out,
            State::On => return false,
            State::Unknown => return false,
        }
        match complex_transition.state_after() {
            State::In => state_after = State::In,
            State::Out => state_after = State::Out,
            State::On => return false,
            State::Unknown => return false,
        }

        point.set_state_before(state_before);
        point.set_state_after(state_after);
        point.set_segment_beginning(segment_begin);
        point.set_segment_end(segment_end);
        true
    }

    /// OCCT Geom2dHatch_Hatcher::ComputeDomains (cxx L1128-1137) — computes
    /// the domains of all the hatchings.
    pub fn compute_domains(&mut self) {
        for indh in 1..=self.my_nb_hatchings {
            if self.my_hatchings.is_bound(indh) {
                self.compute_domains_hatching(indh);
            }
        }
    }

    /// OCCT Geom2dHatch_Hatcher::ComputeDomains (cxx L1144-1806) — computes
    /// the domains of the IndH-th hatching.
    fn compute_domains_hatching(&mut self, indh: usize) {
        {
            let hatching = self.my_hatchings.change_find(indh);
            hatching.clr_domains();
            hatching.set_is_done(false);
        }
        // The OCCT code holds the Hatching reference across Trim(IndH); Rust
        // forces the re-fetch around the &mut self call (same data flow).
        if !self.my_hatchings.find(indh).trim_done() {
            self.trim_hatching(indh);
        }
        if self.my_hatchings.find(indh).status() != ErrorStatus::NoProblem {
            return;
        }

        // bool Points = myKeepPoints; bool Segments = myKeepSegments;
        let points = self.my_keep_points;
        let segments = self.my_keep_segments;
        // int ISav = 0; bool SavPnt = false; int NbOpenedSegments = 0;
        let mut isav = 0usize;
        let mut sav_pnt = false;
        let mut nb_opened_segments = 0i32;
        // int NbPnt = Hatching.NbPoints(); int IPnt = 1;
        let nb_pnt = self.my_hatchings.find(indh).nb_points();

        if nb_pnt == 0 {
            //-- The hatching has to be classified.
            // Geom2dHatch_Classifier Classifier(myElements,
            //                                   Hatching.ClassificationPoint(),
            //                                   0.0000001);
            let classification_pnt = self.my_hatchings.find(indh).classification_point();
            let classifier = Classifier::with_perform(&mut self.my_elements, classification_pnt, 0.0000001);
            if classifier.state() == State::In {
                let domain = HatchGenDomain::new();
                self.my_hatchings.change_find(indh).add_domain(&domain);
            }

            self.my_hatchings.change_find(indh).set_is_done(true);
            return;
        }

        // The OCCT loop holds the Hatching reference throughout; the loop
        // body never calls &mut self methods, so the single borrow is kept.
        let hatching = self.my_hatchings.change_find(indh);
        // for (IPnt = 1; IPnt <= NbPnt; IPnt++)
        for ipnt in 1..=nb_pnt {
            let no_domain = hatching.nb_domains() == 0;
            let first_point = ipnt == 1;
            let last_point = ipnt == nb_pnt;

            // const HatchGen_PointOnHatching& CurPnt = Hatching.Point(IPnt);
            // (cloned: Rust cannot keep the reference across the mutable
            // Hatching calls below; CurPnt is read-only in OCCT.)
            let cur_pnt = hatching.point(ipnt).clone();

            //-----------------------------------------------------------------------
            // Calcul des domaines.
            //-----------------------------------------------------------------------

            let mut state_before = cur_pnt.state_before();
            let mut state_after = cur_pnt.state_after();
            let segment_begin = cur_pnt.segment_beginning();
            let segment_end = cur_pnt.segment_end();

            let mut domain = HatchGenDomain::new();

            //-----------------------------------------------------------------------
            // Initialisations dues au premier point.
            //-----------------------------------------------------------------------

            if first_point {
                sav_pnt = false;
                isav = 0;
                nb_opened_segments = 0;
                if segment_end && segment_begin {
                    if state_after == State::Unknown {
                        state_after = State::In;
                    }
                    if state_before == State::Unknown {
                        state_before = State::In;
                    }

                    if segments {
                        sav_pnt = true;
                        isav = 0;
                    }
                } else if segment_end {
                    if state_after == State::Unknown {
                        state_after = State::In;
                    }

                    if segments {
                        sav_pnt = true;
                        isav = 0;
                    }
                } else if segment_begin {
                    if state_before == State::Unknown {
                        state_before = State::In;
                    }
                    if state_before == State::In {
                        sav_pnt = true;
                        isav = 0;
                    }
                } else {
                    if state_before == State::In {
                        sav_pnt = true;
                        isav = 0;
                    }
                }
            }

            //-----------------------------------------------------------------------
            // Initialisations dues au dernier point.
            //-----------------------------------------------------------------------

            if last_point {
                if segment_end && segment_begin {
                    if state_after == State::Unknown {
                        state_after = State::In;
                    }

                    if state_before == State::Unknown {
                        state_before = State::In;
                    }
                } else if segment_end {
                    if state_after == State::Unknown {
                        state_after = State::In;
                    }
                } else if segment_begin {
                    if state_before == State::Unknown {
                        state_before = State::In;
                    }
                } else {
                }
            }

            //-----------------------------------------------------------------------
            // Cas general.
            //-----------------------------------------------------------------------

            let mut to_append = false;

            if segment_end && segment_begin {
                if state_before != State::In && state_after != State::In {
                    hatching.set_status(ErrorStatus::IncompatibleStates);
                    return;
                }

                if points {
                    if segments {
                        if !sav_pnt {
                            if no_domain {
                                hatching.set_status(ErrorStatus::IncoherentParity);
                            } else {
                                hatching.set_is_done(true);
                            }
                            return;
                        }

                        if isav != 0 {
                            domain.set_first_point(hatching.point(isav));
                        }

                        domain.set_second_point(&cur_pnt);
                        to_append = true;
                        sav_pnt = true;
                        isav = ipnt;
                    } else {
                        let is_in_in = state_before == State::In && state_after == State::In;
                        if sav_pnt && !is_in_in {
                            if no_domain {
                                hatching.set_status(ErrorStatus::IncoherentParity);
                            } else {
                                hatching.set_is_done(true);
                            }

                            return;
                        }

                        domain.set_points(&cur_pnt, &cur_pnt);
                        to_append = true;
                        sav_pnt = false;
                        isav = 0;
                    }
                }
            } else if segment_end {
                if segments {
                    if state_after == State::Out {
                        if !sav_pnt {
                            if no_domain {
                                hatching.set_status(ErrorStatus::IncoherentParity);
                            } else {
                                hatching.set_is_done(true);
                            }
                            return;
                        }

                        if isav != 0 {
                            domain.set_first_point(hatching.point(isav));
                        }

                        domain.set_second_point(&cur_pnt);
                        to_append = true;
                    } else {
                        if points {
                            if isav != 0 {
                                domain.set_first_point(hatching.point(isav));
                            }

                            domain.set_second_point(&cur_pnt);
                            to_append = true;
                            sav_pnt = true;
                            isav = ipnt;
                        }
                    }
                } else {
                    if state_after == State::In {
                        sav_pnt = true;
                        isav = ipnt;
                    }
                }

                nb_opened_segments -= 1;
            } else if segment_begin {
                if segments {
                    if state_before == State::Out {
                        sav_pnt = true;
                        isav = ipnt;
                    } else {
                        if points {
                            if !sav_pnt {
                                if no_domain {
                                    hatching.set_status(ErrorStatus::IncoherentParity);
                                } else {
                                    hatching.set_is_done(true);
                                }

                                return;
                            }

                            if isav != 0 {
                                domain.set_first_point(hatching.point(isav));
                            }

                            domain.set_second_point(&cur_pnt);
                            to_append = true;
                            sav_pnt = true;
                            isav = ipnt;
                        }
                    }
                } else {
                    if state_before == State::In {
                        if !sav_pnt {
                            if no_domain {
                                hatching.set_status(ErrorStatus::IncoherentParity);
                            } else {
                                hatching.set_is_done(true);
                            }

                            return;
                        }

                        if isav != 0 {
                            domain.set_first_point(hatching.point(isav));
                        }

                        domain.set_second_point(&cur_pnt);
                        to_append = true;

                        // Modified by Sergey KHROMOV - Fri Jan  5 12:05:30 2001
                        // SavPnt = false ;
                        // ISav = 0 ;

                        sav_pnt = true;
                        isav = ipnt;
                        // Modified by Sergey KHROMOV - Fri Jan  5 12:05:31 2001
                    }
                }

                nb_opened_segments += 1;
            } else {
                //-- Solution provisoire (lbr le 11 Aout 97)
                //-- si On a 2 points dont des points OUT OUT ou IN IN qui delimitent une isos
                //-- on transforme les transitions
                if state_before == State::Out && state_after == State::Out {
                    if nb_pnt == 2 {
                        if first_point {
                            state_after = State::In;
                        } else {
                            state_before = State::In;
                        }
                    }
                }
                if state_before == State::Out && state_after == State::Out {
                    if sav_pnt {
                        if no_domain {
                            hatching.set_status(ErrorStatus::IncoherentParity);
                        } else {
                            hatching.set_is_done(true);
                        }

                        return;
                    }

                    if points {
                        domain.set_points(&cur_pnt, &cur_pnt);
                        to_append = true;
                        sav_pnt = true;
                        isav = ipnt;
                    }
                } else if state_before == State::Out && state_after == State::In {
                    sav_pnt = true;
                    isav = ipnt;
                } else if state_before == State::In && state_after == State::Out {
                    if !sav_pnt {
                        if no_domain {
                            hatching.set_status(ErrorStatus::IncoherentParity);
                        } else {
                            hatching.set_is_done(true);
                        }

                        return;
                    }

                    if isav != 0 {
                        domain.set_first_point(hatching.point(isav));
                    }

                    domain.set_second_point(&cur_pnt);
                    to_append = true;
                    sav_pnt = false;
                    isav = 0;
                } else if state_before == State::In && state_after == State::In {
                    if points {
                        if nb_opened_segments == 0 {
                            if !sav_pnt {
                                if no_domain {
                                    hatching.set_status(ErrorStatus::IncoherentParity);
                                } else {
                                    hatching.set_is_done(true);
                                }

                                // return;
                                continue;
                            }

                            if isav != 0 {
                                domain.set_first_point(hatching.point(isav));
                            }

                            domain.set_second_point(&cur_pnt);
                            to_append = true;
                            sav_pnt = true;
                            isav = ipnt;
                        } else {
                            if segments {
                                if !sav_pnt {
                                    if no_domain {
                                        hatching.set_status(ErrorStatus::IncoherentParity);
                                    } else {
                                        hatching.set_is_done(true);
                                    }

                                    return;
                                }

                                if isav != 0 {
                                    domain.set_first_point(hatching.point(isav));
                                }

                                domain.set_second_point(&cur_pnt);
                                to_append = true;
                                sav_pnt = true;
                                isav = ipnt;
                            } else {
                                if sav_pnt {
                                    if no_domain {
                                        hatching.set_status(ErrorStatus::IncoherentParity);
                                    } else {
                                        hatching.set_is_done(true);
                                    }

                                    return;
                                }

                                domain.set_points(&cur_pnt, &cur_pnt);
                                to_append = true;
                                sav_pnt = false;
                                isav = 0;
                            }
                        }
                    }
                } else {
                    hatching.set_status(ErrorStatus::IncompatibleStates);
                    return;
                }
            }

            //-----------------------------------------------------------------------
            // Ajout du domaine.
            //-----------------------------------------------------------------------

            if to_append {
                hatching.add_domain(&domain);
            }

            //-----------------------------------------------------------------------
            // Traitement lie au dernier point.
            //-----------------------------------------------------------------------

            if last_point {
                domain.set_points_infinite();
                to_append = false;

                if segment_end && segment_begin {
                    if segments {
                        if !sav_pnt {
                            if no_domain {
                                hatching.set_status(ErrorStatus::IncoherentParity);
                            } else {
                                hatching.set_is_done(true);
                            }

                            return;
                        }

                        if isav != 0 {
                            domain.set_first_point(hatching.point(isav));
                        }

                        to_append = true;
                    }
                } else if segment_end {
                    if state_after == State::In {
                        if !sav_pnt {
                            if no_domain {
                                hatching.set_status(ErrorStatus::IncoherentParity);
                            } else {
                                hatching.set_is_done(true);
                            }

                            return;
                        }

                        if isav != 0 {
                            domain.set_first_point(hatching.point(isav));
                        }

                        to_append = true;
                    }
                } else if segment_begin {
                    if segments {
                        if !sav_pnt {
                            if no_domain {
                                hatching.set_status(ErrorStatus::IncoherentParity);
                            } else {
                                hatching.set_is_done(true);
                            }

                            return;
                        }

                        if isav != 0 {
                            domain.set_first_point(hatching.point(isav));
                        }

                        to_append = true;
                    }
                } else {
                    if state_after == State::In {
                        if !sav_pnt {
                            if no_domain {
                                hatching.set_status(ErrorStatus::IncoherentParity);
                            } else {
                                hatching.set_is_done(true);
                            }

                            return;
                        }

                        if isav != 0 {
                            domain.set_first_point(hatching.point(isav));
                        }

                        to_append = true;
                    }
                }

                if to_append {
                    hatching.add_domain(&domain);
                }
            }
        }

        hatching.set_is_done(true);
    }

    // ===========================================================================
    //  OCCT Geom2dHatch_Hatcher.cxx L1808-1831 — Category : Results.
    // ===========================================================================

    /// OCCT Geom2dHatch_Hatcher::Domain (cxx L1819-1831) — returns the IDom-th
    /// domain of the IndH-th hatching; raises StdFail_NotDone when the
    /// domains were not computed.
    pub fn domain(&self, indh: usize, idom: usize) -> &HatchGenDomain {
        let hatching = self.my_hatchings.find(indh);
        if !hatching.is_done() {
            panic!("Geom2dHatch_Hatcher::Domain");
        }
        hatching.domain(idom)
    }

    // ===========================================================================
    //  OCCT Geom2dHatch_Hatcher.cxx L1833-1919 — Category : Dump.
    // ===========================================================================

    /// OCCT Geom2dHatch_Hatcher::Dump (cxx L1841-1919).
    pub fn dump(&self) {
        println!();
        println!("========================================================");
        println!("=== Dump of the hatcher ================================");
        println!("========================================================");
        println!();

        println!(
            "The points   are {}considered.",
            if self.my_keep_points { "    " } else { "not " }
        );
        println!(
            "The segments are {}considered.",
            if self.my_keep_segments { "    " } else { "not " }
        );
        println!("2D Confusion tolerance : {}", self.my_confusion2d);
        println!("3D Confusion tolerance : {}", self.my_confusion3d);

        println!(
            "{} hatching{}",
            self.my_nb_hatchings,
            if self.my_nb_hatchings == 1 { "" } else { "s" }
        );
        println!(
            "{} element{}",
            self.my_nb_elements,
            if self.my_nb_elements == 1 { "" } else { "s" }
        );

        println!();
        println!("========================================================");
        println!("=== Hatchings ==========================================");
        println!("========================================================");
        println!();

        for indh in 1..=self.my_nb_hatchings {
            print!("Hatching # {}", indh);
            if !self.my_hatchings.is_bound(indh) {
                println!(" is not bound");
            } else {
                let hatching = self.my_hatchings.find(indh);
                let nb_pnt = hatching.nb_points();
                println!(" contains {} restriction points :", nb_pnt);
                for ipnt in 1..=nb_pnt {
                    let pnt_h = hatching.point(ipnt);
                    dump_point_on_hatching(pnt_h, ipnt as i32);
                }
                println!("----------------------------------------------");
            }
        }

        println!();
        println!("========================================================");
        println!("=== Elements ===========================================");
        println!("========================================================");
        println!();

        for inde in 1..=self.my_nb_elements {
            print!("Element # {}", inde);
            if !self.my_elements.is_bound(inde) {
                println!(" is not bound");
            } else {
                let element = self.my_elements.find(inde);
                match element.orientation() {
                    Orientation::Forward => println!(" is FORWARD"),
                    Orientation::Reversed => println!(" is REVERSED"),
                    Orientation::Internal => println!(" is INTERNAL"),
                    Orientation::External => println!(" is EXTERNAL"),
                }
            }
        }

        println!();
    }
}

/// OCCT HatchGen_PointOnHatching::Dump (HatchGen_PointOnHatching.cxx L152-237)
/// — mirrored here because the landed hatch_gen facade carries no Dump.
fn dump_point_on_hatching(pnt_h: &PointOnHatching, index: i32) {
    print!("--- Point on hatching ");
    if index > 0 {
        print!("# {:>3} ", index);
    } else {
        print!("------");
    }
    println!("------------------");

    println!("    Index of the hatching = {}", pnt_h.index());
    println!("    Parameter on hatching = {}", pnt_h.parameter());
    print!("    Position  on hatching = ");
    match pnt_h.position() {
        Orientation::Forward => println!("FORWARD  (i.e. BEGIN  )"),
        Orientation::Internal => println!("INTERNAL (i.e. MIDDLE )"),
        Orientation::Reversed => println!("REVERSED (i.e. END    )"),
        Orientation::External => println!("EXTERNAL (i.e. UNKNOWN)"),
    }
    print!("    State Before          = ");
    match pnt_h.state_before() {
        State::In => println!("IN"),
        State::Out => println!("OUT"),
        State::On => println!("ON"),
        State::Unknown => println!("UNKNOWN"),
    }
    print!("    State After           = ");
    match pnt_h.state_after() {
        State::In => println!("IN"),
        State::Out => println!("OUT"),
        State::On => println!("ON"),
        State::Unknown => println!("UNKNOWN"),
    }
    println!(
        "    Beginning of segment  = {}",
        if pnt_h.segment_beginning() { "TRUE" } else { "FALSE" }
    );
    println!(
        "    End       of segment  = {}",
        if pnt_h.segment_end() { "TRUE" } else { "FALSE" }
    );

    let nb_pnt = pnt_h.nb_points();
    if nb_pnt == 0 {
        println!("    No points on element");
    } else {
        println!("    Contains {} points on element", nb_pnt);
        for ipnt in 1..=nb_pnt {
            let point = pnt_h.point(ipnt);
            dump_point_on_element(point, ipnt as i32);
        }
    }

    println!("----------------------------------------------");
}

/// OCCT HatchGen_PointOnElement::Dump (HatchGen_PointOnElement.cxx L208-293).
fn dump_point_on_element(pnt_e: &PointOnElement, index: i32) {
    print!("    --- Point on element ");
    if index > 0 {
        print!("# {:>3} ", index);
    } else {
        print!("------");
    }
    println!("---------------");

    println!("        Index of the element = {}", pnt_e.index());
    println!("        Parameter on element = {}", pnt_e.parameter());
    print!("        Position  on element = ");
    match pnt_e.position() {
        Orientation::Forward => println!("FORWARD  (i.e. BEGIN  )"),
        Orientation::Internal => println!("INTERNAL (i.e. MIDDLE )"),
        Orientation::Reversed => println!("REVERSED (i.e. END    )"),
        Orientation::External => println!("EXTERNAL (i.e. UNKNOWN)"),
    }
    print!("        Intersection Type    = ");
    match pnt_e.intersection_type() {
        IntersectionType::True => println!("TRUE"),
        IntersectionType::Touch => println!("TOUCH"),
        IntersectionType::Tangent => println!("TANGENT"),
        IntersectionType::Undetermined => println!("UNDETERMINED"),
    }
    print!("        State Before         = ");
    match pnt_e.state_before() {
        State::In => println!("IN"),
        State::Out => println!("OUT"),
        State::On => println!("ON"),
        State::Unknown => println!("UNKNOWN"),
    }
    print!("        State After          = ");
    match pnt_e.state_after() {
        State::In => println!("IN"),
        State::Out => println!("OUT"),
        State::On => println!("ON"),
        State::Unknown => println!("UNKNOWN"),
    }
    println!(
        "        Beginning of segment = {}",
        if pnt_e.segment_beginning() { "TRUE" } else { "FALSE" }
    );
    println!(
        "        End       of segment = {}",
        if pnt_e.segment_end() { "TRUE" } else { "FALSE" }
    );

    println!("    ------------------------------------------");
}

/// OCCT IntRes2d_Intersection::IsEmpty (IntRes2d_Intersection.lxx L40-46) —
/// the accessor is not forwarded by the landed HatchIntersector facade, so
/// the body is mirrored here (raises StdFail_NotDone when not done).
fn intersector_is_empty(intersector: &HatchIntersector) -> bool {
    if !intersector.is_done() {
        panic!("IntRes2d_Intersection::IsEmpty - StdFail_NotDone");
    }
    intersector.nb_points() == 0 && intersector.nb_segments() == 0
}

// ===========================================================================
//  Anchor tests (analytically assertable).
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Line2d, TrimmedCurve2};

    /// A bounded element on the segment [p1, p2] (the rcad analogue of a
    /// Geom2dAdaptor_Curve over a trimmed Geom2d_Line).
    fn seg_curve(p1: DVec2, p2: DVec2) -> Curve2d {
        let dir = (p2 - p1).normalize();
        Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(Curve2d::Line(Line2d::new(p1, dir))),
            t_min: 0.0,
            t_max: 1.0,
        })
    }

    /// An unbounded straight hatching line.
    fn line_curve(origin: DVec2, dir: DVec2) -> Curve2d {
        Curve2d::Line(Line2d::new(origin, dir))
    }

    /// The unit square boundary, CCW, all edges FORWARD: element 1 = bottom,
    /// 2 = right, 3 = top, 4 = left.
    fn square_hatcher() -> Hatcher {
        let mut hatcher = Hatcher::new(
            HatchIntersector::with_tolerances(1.0e-7, 1.0e-7),
            1.0e-7,
            1.0e-7,
            false,
            false,
        );
        hatcher.add_element(&seg_curve(DVec2::new(0., 0.), DVec2::new(1., 0.)), Orientation::Forward);
        hatcher.add_element(&seg_curve(DVec2::new(1., 0.), DVec2::new(1., 1.)), Orientation::Forward);
        hatcher.add_element(&seg_curve(DVec2::new(1., 1.), DVec2::new(0., 1.)), Orientation::Forward);
        hatcher.add_element(&seg_curve(DVec2::new(0., 1.), DVec2::new(0., 0.)), Orientation::Forward);
        hatcher
    }

    /// The hatching line y = 0.5 crossing the square from left to right; it
    /// enters at line parameter 1.0 (left edge, element 4) and exits at 2.0
    /// (right edge, element 2).
    fn hatch_line() -> Curve2d {
        line_curve(DVec2::new(-1., 0.5), DVec2::X)
    }

    // -- 1. AddElement / AddHatching counters, accessors round-trip. --------

    #[test]
    fn hatcher_add_and_accessors_roundtrip() {
        let mut h = Hatcher::new(
            HatchIntersector::with_tolerances(1.0e-7, 1.0e-7),
            1.0e-7,
            3.0e-7,
            false,
            false,
        );
        assert_eq!(h.confusion2d(), 1.0e-7);
        assert_eq!(h.confusion3d(), 3.0e-7);
        assert!(!h.keep_points() && !h.keep_segments());
        assert_eq!(h.intersector().confusion_tolerance(), 1.0e-7);

        // Intersector setter and ChangeIntersector.
        h.set_intersector(HatchIntersector::with_tolerances(2.0e-7, 2.0e-7));
        assert_eq!(h.intersector().confusion_tolerance(), 2.0e-7);
        h.change_intersector().set_confusion_tolerance(4.0e-7);
        assert_eq!(h.intersector().confusion_tolerance(), 4.0e-7);

        // Elements: sequential indices, both AddElement entries.
        let bottom = seg_curve(DVec2::new(0., 0.), DVec2::new(1., 0.));
        let right = seg_curve(DVec2::new(1., 0.), DVec2::new(1., 1.));
        assert_eq!(h.add_element(&bottom, Orientation::Forward), 1);
        assert_eq!(h.add_element_from_curve(&right, Orientation::Reversed), 2);
        assert_eq!(h.element_curve(1).value(0.0), DVec2::new(0., 0.));
        assert_eq!(h.element_curve(2).value(1.0), DVec2::new(1., 1.));
        assert_eq!(h.element(2).orientation(), Orientation::Reversed);
        h.element(2).set_orientation(Orientation::Forward);
        assert_eq!(h.element(2).orientation(), Orientation::Forward);

        // Hatching: index, curve, empty point count.
        let hl = line_curve(DVec2::new(-1., 0.5), DVec2::X);
        let indh = h.add_hatching(&hl);
        assert_eq!(indh, 1);
        assert_eq!(h.hatching_curve(indh).value(0.5), DVec2::new(-0.5, 0.5));
        assert_eq!(h.nb_points(indh), 0);
        h.hatching(indh).set_is_done(true);
        assert!(h.is_done(indh));
        assert_eq!(h.status(indh), ErrorStatus::NoProblem);
        assert!(!h.trim_done(indh) && !h.trim_failed(indh));

        // Confusion setters round-trip.
        h.set_confusion2d(2.0e-7);
        assert_eq!(h.confusion2d(), 2.0e-7);
        h.set_confusion3d(5.0e-7);
        assert_eq!(h.confusion3d(), 5.0e-7);

        // KeepPoints / KeepSegments setters -> ClrDomains -> IsDone false.
        h.set_keep_points(true);
        assert!(h.keep_points());
        assert!(!h.is_done(indh));
        h.set_keep_segments(true);
        assert!(h.keep_segments());

        // Clear removes everything; the counters restart at 1.
        h.clear();
        assert_eq!(h.add_element(&bottom, Orientation::Forward), 1);
        assert_eq!(h.add_hatching(&hl), 1);
    }

    // -- 2. Trim: straight hatching across the square. -----------------------

    #[test]
    fn hatcher_trim_line_through_square() {
        let mut h = square_hatcher();
        let indh = h.add_hatching(&hatch_line());
        h.trim_hatching(indh);
        assert!(h.trim_done(indh));
        assert!(!h.trim_failed(indh));
        assert_eq!(h.status(indh), ErrorStatus::NoProblem);
        assert_eq!(h.nb_points(indh), 2);

        // Entering point: line parameter 1.0 on the left edge (element 4).
        let p1 = h.point(indh, 1);
        assert!((p1.parameter() - 1.0).abs() < 1.0e-9);
        assert_eq!(p1.nb_points(), 1);
        let pnte1 = p1.point(1);
        assert_eq!(pnte1.index(), 4);
        assert!((pnte1.parameter() - 0.5).abs() < 1.0e-9);
        assert_eq!(pnte1.position(), Orientation::Internal);
        assert_eq!(pnte1.intersection_type(), IntersectionType::True);
        assert_eq!(pnte1.state_before(), State::Out);
        assert_eq!(pnte1.state_after(), State::In);

        // Exiting point: line parameter 2.0 on the right edge (element 2).
        let p2 = h.point(indh, 2);
        assert!((p2.parameter() - 2.0).abs() < 1.0e-9);
        assert_eq!(p2.nb_points(), 1);
        let pnte2 = p2.point(1);
        assert_eq!(pnte2.index(), 2);
        assert!((pnte2.parameter() - 0.5).abs() < 1.0e-9);
        assert_eq!(pnte2.position(), Orientation::Internal);
        assert_eq!(pnte2.intersection_type(), IntersectionType::True);
        assert_eq!(pnte2.state_before(), State::In);
        assert_eq!(pnte2.state_after(), State::Out);

        // GlobalTransition wrote the IN/OUT states on the hatching points.
        assert_eq!(p1.state_before(), State::Out);
        assert_eq!(p1.state_after(), State::In);
        assert_eq!(p2.state_before(), State::In);
        assert_eq!(p2.state_after(), State::Out);

        // Trim(Curve) overload: adds and trims in one step (hatching #2).
        let indh2 = h.trim_curve(&hatch_line());
        assert_eq!(indh2, 2);
        assert_eq!(h.nb_points(indh2), 2);
    }

    // -- 3. ComputeDomains / Domain. ------------------------------------------

    #[test]
    fn hatcher_compute_domains_square() {
        let mut h = square_hatcher();
        let indh = h.add_hatching(&hatch_line());
        h.compute_domains_hatching(indh);
        assert!(h.is_done(indh));
        assert_eq!(h.nb_domains(indh), 1);
        let d = h.domain(indh, 1);
        assert!(d.has_first_point() && d.has_second_point());
        assert!((d.first_point().parameter() - 1.0).abs() < 1.0e-9);
        assert!((d.second_point().parameter() - 2.0).abs() < 1.0e-9);

        // A hatching fully outside the contour: no intersection points, so
        // the classifier branch runs, the state is OUT and no domain appears.
        let mut h2 = square_hatcher();
        let indh2 = h2.add_hatching(&line_curve(DVec2::new(-1., 2.0), DVec2::X));
        h2.compute_domains();
        assert!(h2.is_done(indh2));
        assert_eq!(h2.nb_points(indh2), 0);
        assert_eq!(h2.nb_domains(indh2), 0);
    }

    // -- 4. KeepPoints / KeepSegments two-state difference. -------------------

    #[test]
    fn hatcher_keep_flags_invalidate_domains() {
        let mut h = square_hatcher();
        let indh = h.add_hatching(&hatch_line());
        h.compute_domains();
        assert!(h.is_done(indh));
        assert_eq!(h.nb_domains(indh), 1);

        // KeepSegments(true) clears the computed domains (ClrDomains).
        h.set_keep_segments(true);
        assert!(h.keep_segments());
        assert!(!h.is_done(indh));

        // NbDomains raises StdFail_NotDone while the domains are not
        // computed.
        let panicked = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = h.nb_domains(indh);
        }));
        assert!(panicked.is_err());

        // Recomputing with KeepSegments=true yields the same analytic domain
        // for this crossing case (no tangent segment on the hatching).
        h.compute_domains_hatching(indh);
        assert!(h.is_done(indh));
        assert_eq!(h.nb_domains(indh), 1);
        assert!((h.domain(indh, 1).first_point().parameter() - 1.0).abs() < 1.0e-9);
        assert!((h.domain(indh, 1).second_point().parameter() - 2.0).abs() < 1.0e-9);
    }

    // -- 5. RemElement / re-trim / removal bookkeeping. -----------------------

    #[test]
    fn hatcher_rem_element_and_retrim() {
        let mut h = square_hatcher();
        let indh = h.add_hatching(&hatch_line());
        h.trim_hatching(indh);
        assert_eq!(h.nb_points(indh), 2);

        // Removing the left edge (4) prunes the entering point (its only
        // point on element references element 4).
        h.rem_element(4);
        assert_eq!(h.nb_points(indh), 1);
        assert!((h.point(indh, 1).parameter() - 2.0).abs() < 1.0e-9);

        // The removed last slot is reclaimed: myNbElements dropped 4 -> 3 and
        // the next AddElement reuses the index 4.
        let left = seg_curve(DVec2::new(0., 1.), DVec2::new(0., 0.));
        assert_eq!(h.add_element(&left, Orientation::Forward), 4);

        // Re-trim of the same hatching by the elements already given.
        h.trim_hatching(indh);
        assert_eq!(h.nb_points(indh), 2);
        assert!((h.point(indh, 1).parameter() - 1.0).abs() < 1.0e-9);

        // Removing an inner element (3, top edge, no intersections) keeps the
        // points and the high-water count; the freed slot 3 is handed out by
        // the next AddElement.
        h.rem_element(3);
        assert_eq!(h.nb_points(indh), 2);
        let top = seg_curve(DVec2::new(1., 1.), DVec2::new(0., 1.));
        assert_eq!(h.add_element(&top, Orientation::Forward), 3);

        // ClrElements resets the counter and clears the hatching points.
        h.clr_elements();
        assert_eq!(h.nb_points(indh), 0);
        assert!(!h.trim_done(indh) && !h.trim_failed(indh));
        assert_eq!(h.add_element(&left, Orientation::Forward), 1);

        // RemHatching / ClrHatchings bookkeeping (the removed last slot is
        // reclaimed as well).
        let indh2 = h.add_hatching(&hatch_line());
        assert_eq!(indh2, 2);
        h.rem_hatching(indh2);
        assert_eq!(h.add_hatching(&hatch_line()), 2);
        h.clr_hatchings();
        assert_eq!(h.add_hatching(&hatch_line()), 1);
    }
}
