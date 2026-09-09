//! ChFiDS gap-fill — methods of `ChFiDS_Spine`, `ChFiDS_FilSpine` and
//! `ChFiDS_ElSpine` that were missing from `chfi_ds.rs`, translated 1:1
//! from OCCT TKFillet/ChFiDS.  Kept in a dedicated file because
//! `chfi_ds.rs` exceeds the 2000-line guideline.

use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

use glam::{DVec2, DVec3};
use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::topo::topods::{Orientation, Shape};
use rcad_kernel::geom::CurveEval as _;

use crate::geomalgo::law::{
    LawComposite, LawConstant, LawFunction, LawFunctionHandle, LawInterpol, LawS,
};

use super::chfi_ds::{ChFiDSElSpine, ChFiDSFilSpine, ChFiDSSpine, elclib_in_period};

// =========================================================================
// OCCT ChFiDS_Spine — missing method translations
// =========================================================================

impl ChFiDSSpine {
    /// OCCT ChFiDS_Spine.lxx L117-120 — OffsetEdges(I) (1-based).
    pub fn offset_edges(&self, i: usize) -> &Shape {
        &self.offsetspine[i - 1]
    }

    /// OCCT ChFiDS_Spine.cxx L89-93 — AppendElSpine(Els).
    pub fn append_el_spine(&mut self, els: ChFiDSElSpine) {
        self.elspines.push(els);
    }

    /// OCCT ChFiDS_Spine.cxx L96-99 — AppendOffsetElSpine(Els).
    pub fn append_offset_el_spine(&mut self, els: ChFiDSElSpine) {
        self.offset_elspines.push(els);
    }

    /// OCCT ChFiDS_Spine.cxx L103-106 — ElSpine(const TopoDS_Edge& E).
    pub fn el_spine_of_edge(&self, e: &Shape) -> Option<&ChFiDSElSpine> {
        self.el_spine_of_index(self.index_of_edge(e))
    }

    /// OCCT ChFiDS_Spine.cxx L108-116 — ElSpine(const int IE) (1-based).
    pub fn el_spine_of_index(&self, ie: usize) -> Option<&ChFiDSElSpine> {
        let mut wmil = 0.5 * (self.first_parameter_of(ie) + self.last_parameter_of(ie));
        if self.is_periodic() {
            wmil = elclib_in_period(wmil, self.first_parameter(), self.last_parameter());
        }
        self.el_spine_of_param(wmil)
    }

    /// OCCT ChFiDS_Spine.cxx L118-139 — ElSpine(const double W).
    pub fn el_spine_of_param(&self, w: f64) -> Option<&ChFiDSElSpine> {
        if self.elspines.len() == 1 {
            return self.elspines.first();
        }
        for cur in &self.elspines {
            let uf = cur.first_parameter();
            let ul = cur.last_parameter();
            if uf <= w && w <= ul {
                return Some(cur);
            }
        }
        None
    }

    /// OCCT ChFiDS_Spine.cxx L143-146 — ChangeElSpines().
    pub fn change_el_spines(&mut self) -> &mut Vec<ChFiDSElSpine> {
        &mut self.elspines
    }

    /// OCCT ChFiDS_Spine.cxx L150-153 — ChangeOffsetElSpines().
    pub fn change_offset_el_spines(&mut self) -> &mut Vec<ChFiDSElSpine> {
        &mut self.offset_elspines
    }

    /// OCCT ChFiDS_Spine.cxx L345-348 — Resolution(R3d).
    pub fn resolution(&self, r3d: f64) -> f64 {
        r3d
    }

    /// OCCT ChFiDS_Spine.cxx L396-399 — HasFirstTgt().
    pub fn has_first_tgt(&self) -> bool {
        self.hasfirsttgt
    }

    /// OCCT ChFiDS_Spine.cxx L403-406 — HasLastTgt().
    pub fn has_last_tgt(&self) -> bool {
        self.haslasttgt
    }

    /// OCCT ChFiDS_Spine.cxx L496-499 — UnsetReference().
    pub fn unset_reference(&mut self) {
        self.hasref = false;
    }

    /// OCCT ChFiDS_Spine.cxx L743-746 — D0(AbsC, P).
    pub fn d0(&mut self, absc: f64) -> DVec3 {
        self.value_at(absc)
    }

    /// OCCT ChFiDS_Spine.cxx L798-851 — D2(AbsC, P, V1, V2): point, tangent
    /// and second derivative (normalized, orientation-adjusted) on the
    /// composite spine.
    pub fn d2(&mut self, absc: f64) -> (DVec3, DVec3, DVec3) {
        let mut l = absc;
        let index = self.prepare(&mut l);

        if index == -1 {
            let p = self.firstori + self.firsttgt * l;
            return (p, self.firsttgt, DVec3::ZERO);
        } else if index as i32 == self.abscissa.as_ref().map_or(0, |a| a.len()) as i32 + 1 {
            let p = self.lastori + self.lasttgt * l;
            return (p, self.lasttgt, DVec3::ZERO);
        }
        let index = index as usize;
        self.indexofcurve = index as i32;
        let e = &self.spine[index - 1];
        let ed = e.as_edge().expect("not an edge");
        let curve = ed.curve.as_ref().expect("edge curve").clone();
        let t = l / self.length_of(index);
        let (cf, cl) = (ed.range[0], ed.range[1]);
        let uapp = (1.0 - t) * cf + t * cl;
        let u = rcad_kernel::base::gcpnts::abscissa_point::abscissa_point_parameter(
            &curve, cf, cl, l, uapp,
        );
        let p = curve.point_at(u);
        let mut v1 = curve.derivative_at(u);
        let mut v2 = curve.derivative2_at(u);
        // OCCT: N1 = V1.SquareMagnitude()
        let mut n1 = v1.dot(v1);
        // OCCT: D2 = -(V1.Dot(V2)) * (1./N1) * (1./N1)
        let d2c = -(v1.dot(v2)) * (1.0 / n1) * (1.0 / n1);
        v2 = v2 * (1.0 / n1);
        n1 = n1.sqrt();
        let va = v1 * d2c;
        v2 = v2 + va;
        let mut d1 = 1.0 / n1;
        if e.orientation == Orientation::Reversed {
            d1 = -d1;
        }
        v1 = v1 * d1;
        (p, v1, v2)
    }

    /// OCCT ChFiDS_Spine.cxx L855-862 — SetCurrent(Index).
    /// The OCCT myCurve.Initialize(TopoDS::Edge(...)) re-initialization is a
    /// pending boundary: rcad tracks the current elementary spine by
    /// `indexofcurve` alone (BRepAdaptor_Curve pending translation).
    pub fn set_current(&mut self, index: i32) {
        if index != self.indexofcurve {
            self.indexofcurve = index;
        }
    }

    /// OCCT ChFiDS_Spine.lxx L152-155 — CurrentIndexOfElementarySpine().
    pub fn current_index_of_elementary_spine(&self) -> i32 {
        self.indexofcurve
    }

    /// OCCT ChFiDS_Spine.lxx L159-163 — Mode().
    pub fn mode(&self) -> super::chfi_ds::ChFiDS_ChamfMode {
        self.my_mode
    }

    /// OCCT ChFiDS_Spine.lxx L167-170 — GetTolesp().
    pub fn get_tolesp(&self) -> f64 {
        self.tolesp
    }
}

// =========================================================================
// OCCT ChFiDS_FilSpine — missing method translations
// =========================================================================

impl ChFiDSFilSpine {
    /// OCCT ChFiDS_FilSpine.cxx L96-119 — UnSetRadius(const TopoDS_Edge& E)
    /// (1-based IE).
    pub fn unset_radius_on_edge(&mut self, e: &Shape) {
        self.base.splitdone = false;
        let ie = self.base.index_of_edge(e);

        let uf = self.base.first_parameter_of(ie);
        let ul = self.base.last_parameter_of(ie);
        let mut ifirst = 0usize;
        let mut ilast = 0usize;
        for i in 1..=self.parandrad.len() {
            if (self.parandrad[i - 1].x - uf).abs() <= f64::MIN_POSITIVE {
                ifirst = i;
            }
            if (self.parandrad[i - 1].x - ul).abs() <= f64::MIN_POSITIVE {
                ilast = i;
            }
        }
        if ifirst != 0 && ilast != 0 {
            // OCCT: parandrad.Remove(ifirst, ilast) — the inclusive range.
            self.parandrad.drain((ifirst - 1)..ilast);
        }
    }

    /// OCCT ChFiDS_FilSpine.cxx L220-231 — UnSetRadius(const TopoDS_Vertex& V).
    pub fn unset_radius_at_vertex(&mut self, v: &Shape) {
        let npar = self.base.absc_of_vertex(v);
        for i in 1..=self.parandrad.len() {
            if self.parandrad[i - 1].x == npar {
                self.parandrad.remove(i - 1);
                break;
            }
        }
    }
}

// =========================================================================
// OCCT ChFiDS_ElSpine — missing accessor translations.
// (SetFirstPointAndTgt / SetLastPointAndTgt / the default constructor live
// in chfi3d.rs — pre-existing.)
// =========================================================================

impl ChFiDSElSpine {
    /// OCCT ChFiDS_ElSpine.cxx L76-79 — FirstParameter().
    pub fn first_parameter(&self) -> f64 {
        self.firstparam
    }

    /// OCCT ChFiDS_ElSpine.cxx L83-86 — LastParameter().
    pub fn last_parameter(&self) -> f64 {
        self.lastparam
    }

    /// OCCT ChFiDS_ElSpine.cxx L269-273 — FirstPointAndTgt(P, T).
    pub fn first_point_and_tgt(&self) -> (DVec3, DVec3) {
        (self.firstpnt, self.firsttgt)
    }

    /// OCCT ChFiDS_ElSpine.cxx L277-281 — LastPointAndTgt(P, T).
    pub fn last_point_and_tgt(&self) -> (DVec3, DVec3) {
        (self.lastpnt, self.lasttgt)
    }

    /// OCCT ChFiDS_ElSpine.cxx L90-93 — GetSavedFirstParameter().
    pub fn get_saved_first_parameter(&self) -> f64 {
        self.pfirstsav
    }

    /// OCCT ChFiDS_ElSpine.cxx L97-100 — GetSavedLastParameter().
    pub fn get_saved_last_parameter(&self) -> f64 {
        self.plastsav
    }

    /// OCCT ChFiDS_ElSpine.cxx L155-159 — SetPeriodic(I).
    pub fn set_periodic(&mut self, i: bool) {
        self.periodic = i;
        self.period = self.lastparam - self.firstparam;
    }

    /// OCCT ChFiDS_ElSpine.cxx L163-170 — Period().
    pub fn period(&self) -> f64 {
        if !self.periodic {
            panic!("Standard_Failure: ElSpine non periodique");
        }
        self.period
    }

    /// OCCT ChFiDS_ElSpine.cxx L216-219 — SaveFirstParameter().
    pub fn save_first_parameter(&mut self) {
        self.pfirstsav = self.firstparam;
    }

    /// OCCT ChFiDS_ElSpine.cxx L223-226 — SaveLastParameter().
    pub fn save_last_parameter(&mut self) {
        self.plastsav = self.lastparam;
    }

    /// OCCT ChFiDS_ElSpine.cxx L230-242 — SetOrigin(O): re-origin of the
    /// underlying periodic BSpline (downcast to Geom_BSplineCurve, then
    /// bs->SetOrigin(O, Precision::PConfusion()) + curve.Load(bs)).  The
    /// Geom_BSplineCurve::SetOrigin(U, Tol) body lives with the other
    /// BSpline primitives in chfi3d_perform_elspine (performed by
    /// ChFi3d_PerformElSpine).  The null-curve state cannot exist in OCCT;
    /// the no-op keeps the pre-PerformElSpine consumer behavior.
    pub fn set_origin(&mut self, o: f64) {
        if let Some(rcad_kernel::geom::Curve3::BSpline(bs)) = self.curve.as_mut() {
            super::chfi3d_perform_elspine::bspline_set_origin_u_tol(bs, o, CONFUSION);
        }
    }

    /// OCCT ChFiDS_ElSpine.cxx L134-137 — Resolution(R3d) (the
    /// GeomAdaptor_Curve resolution of the loaded curve).  Null-curve
    /// fallback: the pre-PerformElSpine consumers keep the input tolerance
    /// (the stub semantics this method replaces).
    pub fn resolution(&self, r3d: f64) -> f64 {
        match &self.curve {
            Some(c) => c.resolution(r3d),
            None => r3d,
        }
    }

    /// OCCT ChFiDS_ElSpine.cxx L299-302 — SetCurve(C) (curve.Load(C)).
    pub fn set_curve(&mut self, c: rcad_kernel::geom::Curve3) {
        self.curve = Some(c);
    }

    /// OCCT ChFiDS_ElSpine — the loaded curve as the Adaptor3d_Curve basis
    /// (Adaptor3d_Curve::Value/D1/... route to the GeomAdaptor_Curve
    /// member).  None encodes the not-yet-performed state (OCCT has no null
    /// curve).
    pub fn adaptor_curve(&self) -> Option<rcad_kernel::geom::Curve3> {
        self.curve.clone()
    }

    /// OCCT ChFiDS_ElSpine (Adaptor3d_Curve::Value → curve member Value).
    pub fn value(&self, u: f64) -> DVec3 {
        match &self.curve {
            Some(c) => c.point_at(u),
            None => self.firstpnt + self.firsttgt * (u - self.firstparam),
        }
    }

    /// OCCT ChFiDS_ElSpine (Adaptor3d_Curve::D1 → curve member D1).
    pub fn d1(&self, u: f64) -> (DVec3, DVec3) {
        match &self.curve {
            Some(c) => (c.point_at(u), c.derivative_at(u)),
            None => (DVec3::ZERO, DVec3::ZERO),
        }
    }

    /// OCCT ChFiDS_ElSpine (Adaptor3d_Curve::D2 → curve member D2).
    pub fn d2(&self, u: f64) -> (DVec3, DVec3, DVec3) {
        match &self.curve {
            Some(c) => (c.point_at(u), c.derivative_at(u), c.derivative2_at(u)),
            None => (DVec3::ZERO, DVec3::ZERO, DVec3::ZERO),
        }
    }

    /// OCCT ChFiDS_ElSpine.cxx L262-265 — AddVertexWithTangent(anAx1).
    pub fn add_vertex_with_tangent(&mut self, an_ax1: rcad_kernel::math::gp::Ax1) {
        self.vertices_with_tangents.push(an_ax1);
    }

    /// OCCT ChFiDS_ElSpine.cxx L285-288 — NbVertices().
    pub fn nb_vertices(&self) -> usize {
        self.vertices_with_tangents.len()
    }

    /// OCCT ChFiDS_ElSpine.cxx L292-295 — VertexWithTangent(Index) (1-based).
    pub fn vertex_with_tangent(&self, index: usize) -> rcad_kernel::math::gp::Ax1 {
        self.vertices_with_tangents[index - 1]
    }

    // OCCT ChFiDS_ElSpine.cxx L339-384 — Line()/Circle()/Ellipse()/
    // Hyperbola()/Parabola()/Bezier()/BSpline() (the typed curve queries
    // behind GetType()).  Pending consumer: no translated ChFi3d code calls
    // them yet; each routes to the matching Curve3 variant of `curve` when
    // a consumer appears.
}

// =========================================================================
// OCCT ChFiDS_FilSpine — the Law mechanism (ChFiDS_FilSpine.cxx
// L143-216 replay tail, L235-242, L369-381, L383-526 static mklaw,
// L530-854 ComputeLaw, L858-870 Law, L874-889 ChangeLaw,
// L893-929 MaxRadFromSeqAndLaws).
//
// Architecture deviation for the OCCT field
// `NCollection_List<occ::handle<Law_Function>> myLaws`
// (ChFiDS_FilSpine.hxx): the rcad ChFiDSFilSpine struct lives in
// chfi_ds.rs (frozen for this batch) whose placeholder
// `laws: Vec<LawFunction>` field cannot change type, so the real law
// list is stored in the thread-local side table below, keyed by the
// spine's first edge identity (Shape::ptr_id — the same ptr-based
// identity the chfi3d_ds side tables use).  The key is value-stable
// across the by-value ChFiDSSpineHandle clones that the rcad handle
// semantics produce, so the laws follow the spine through
// take/write-back cycles (st.my_spine = Some(spine)) exactly like the
// other ChFiDS_FilSpine fields.
// =========================================================================

thread_local! {
    static FILSPINE_LAWS: RefCell<HashMap<u64, Vec<Rc<RefCell<LawComposite>>>>> =
        RefCell::new(HashMap::new());
}

impl ChFiDSFilSpine {
    /// The side-table key (OCCT: the spine object identity).
    fn filspine_laws_key(&self) -> u64 {
        self.base.spine[0].ptr_id()
    }

    /// Access to this spine's law list (OCCT: the myLaws member — stored as
    /// Law_Composite handles, the concrete type ComputeLaw produces).
    fn filspine_laws<R>(&self, f: impl FnOnce(&mut Vec<Rc<RefCell<LawComposite>>>) -> R) -> R {
        FILSPINE_LAWS.with(|regs| {
            let mut regs = regs.borrow_mut();
            f(regs.entry(self.filspine_laws_key()).or_default())
        })
    }

    /// OCCT ChFiDS_FilSpine.cxx L235-242 — SetRadius(handle<Law_Function>&,
    /// IinC).  Literal OCCT body: the Law_Composite is a function-local
    /// object that is never stored (only parandrad.Clear() has an effect —
    /// the laws themselves are (re)built by AppendLaw/ComputeLaw from
    /// parandrad).
    pub fn set_radius_law_fn(&mut self, c: LawFunctionHandle, _iinc: usize) {
        self.base.splitdone = false;
        // OCCT: prout = new Law_Composite(); lst = prout->ChangeLaws();
        //       lst.Append(C);
        let mut prout = LawComposite::new();
        prout.change_laws().push(c);
        // OCCT: parandrad.Clear();
        self.parandrad.clear();
    }

    /// OCCT ChFiDS_FilSpine.cxx L50-53 (the laws.Clear() line of Reset) —
    /// exposed for the reset paths that live in chfi_ds.rs (frozen file;
    /// callers will be switched when it reopens).
    pub fn clear_laws(&mut self) {
        self.filspine_laws(|laws| laws.clear());
    }

    /// OCCT ChFiDS_FilSpine.cxx L193-215 — the splitdone tail of
    /// SetRadius(UandR, IinC): "si le split est done il faut rejouer la law
    /// correspondant au parametre W".  Standalone here because the
    /// set_radius_uandr body that calls it lives in chfi_ds.rs (frozen
    /// file, tail marked pending there); the replay is wired by its
    /// callers.
    pub fn perform_law_replay(&mut self, w: f64) {
        self.filspine_laws(|laws| {
            // OCCT: NCollection_List<...>::Iterator It(elspines);
            //       NCollection_List<...>::Iterator Itl(laws);
            //       Els = It.Value();
            let els_periodic = self.base.elspines.first().map(|e| e.periodic);
            if els_periodic == Some(true) {
                // OCCT: if (Els->IsPeriodic()) Itl.ChangeValue() = ComputeLaw(Els);
                if let Some(els) = self.base.elspines.first() {
                    let newlaw = self.compute_law(els);
                    if let Some(slot) = laws.first_mut() {
                        *slot = newlaw;
                    }
                }
            } else {
                // OCCT: for (; It.More(); It.Next(), Itl.Next()) { ... }
                for (i, els) in self.base.elspines.iter().enumerate() {
                    let uf = els.first_parameter();
                    let ul = els.last_parameter();
                    if uf <= w && w <= ul {
                        // OCCT: Itl.ChangeValue() = ComputeLaw(Els);
                        let newlaw = self.compute_law(els);
                        if let Some(slot) = laws.get_mut(i) {
                            *slot = newlaw;
                        }
                    }
                }
            }
        });
    }

    /// OCCT ChFiDS_FilSpine.cxx L369-373 — AppendElSpine(Els) (the
    /// ChFiDS_FilSpine override: base append + AppendLaw).
    pub fn append_el_spine_fil(&mut self, els: &ChFiDSElSpine) {
        // OCCT: ChFiDS_Spine::AppendElSpine(Els);
        self.base.elspines.push(els.clone());
        // OCCT: AppendLaw(Els);
        self.append_law(els);
    }

    /// OCCT ChFiDS_FilSpine.cxx L377-381 — AppendLaw(Els).
    pub fn append_law(&mut self, els: &ChFiDSElSpine) {
        // OCCT: occ::handle<Law_Composite> l = ComputeLaw(Els); laws.Append(l);
        let l = self.compute_law(els);
        self.filspine_laws(|laws| laws.push(l));
    }

    /// OCCT ChFiDS_FilSpine.cxx L858-870 — Law(Els): the law attached to
    /// the given elspine, as a Law_Composite handle (down_cast<Law_Composite>
    /// in OCCT; every stored law is a ComputeLaw composite).  The elspine
    /// match uses pointer identity with the (first, last) parameter range
    /// as identity substitute (architecture note: OCCT compares
    /// occ::handle identity).
    pub fn law_of(&self, els: &ChFiDSElSpine) -> Option<Rc<RefCell<LawComposite>>> {
        let laws = self.filspine_laws(|l| l.clone());
        for (i, sp) in self.base.elspines.iter().enumerate() {
            let same = std::ptr::eq(sp, els)
                || (sp.first_parameter() == els.first_parameter()
                    && sp.last_parameter() == els.last_parameter());
            if same {
                return laws.get(i).cloned();
            }
        }
        None
    }

    /// OCCT ChFiDS_FilSpine.cxx L530-854 — ComputeLaw(Els) — returns the
    /// new composite as a Law_Composite handle (handle<Law_Composite> in
    /// OCCT).
    pub fn compute_law(&self, els: &ChFiDSElSpine) -> Rc<RefCell<LawComposite>> {
        let tol3d = CONFUSION;
        let mut deb;
        let mut fin;
        let mut curdeb;
        let mut curfin;
        curdeb = els.first_parameter();
        deb = curdeb;
        curfin = els.last_parameter();
        fin = curfin;
        let ideb = self.base.index_of_param(deb, true);
        let ifin = self.base.index_of_param(fin, false);
        let len = self.base.nb_edges();
        // if the spine is periodic, attention to the index and parameters
        let mut spinedeb = self.base.first_parameter();
        let mut spinefin = self.base.last_parameter();

        let mut nbed = ifin - ideb + 1;
        let mut bidfin = ifin;

        let mut loi = LawComposite::new();
        // OCCT: NCollection_List<occ::handle<Law_Function>>& list =
        //           loi->ChangeLaws();
        // (built as a local list and installed into the composite at return
        // — the Rust mutable borrow of loi.change_laws() cannot be held
        // across the SetPeriodic call below.)
        let mut funclist: Vec<LawFunctionHandle> = Vec::new();
        let mut rdeb = 0.0f64;
        let mut rfin = 0.0f64;
        let mut rcur;
        let mut icur = 1usize;
        let mut lastloi: Option<LawFunctionHandle> = None;
        let mut lawencours = false;
        let mut loi_periodic = false;

        if self.base.is_periodic() {
            if deb < 0.0 && ideb > ifin {
                bidfin += len;
            } else if fin > self.base.last_parameter_of(len) && ideb > ifin {
                bidfin += len;
            }
            nbed = bidfin - ideb + 1;
        }
        // OCCT: NCollection_Array1<int> ind(1, nbed); int j = 1;
        let mut ind = vec![0usize; nbed];
        let mut j = 1usize;
        for i in ideb..=bidfin {
            ind[j - 1] = ((i - 1) % len) + 1;
            j += 1;
        }

        if els.periodic {
            // A pereodic composite is created at range, which is eventually
            // offset relatively to the elspine, to avoid a single point at
            // origin.
            // OCCT: loi->SetPeriodic();
            loi_periodic = true;
            // Is there a constant edge?
            let mut k = 1usize;
            while k <= len {
                if self.is_constant_on(k) {
                    // yes  !
                    curdeb = self.base.first_parameter_of(k);
                    deb = curdeb;
                    spinedeb = curdeb;
                    fin = deb + self.base.period();
                    spinefin = fin;
                    for l in 1..=len {
                        ind[l - 1] = ((k + l - 2) % len) + 1;
                    }
                    rdeb = self.radius_on(k);
                    rfin = rdeb;
                    icur += 1;
                    if len == 1 {
                        // because InPeriod will make 0.!!!
                        curfin = self.base.last_parameter_of(k);
                    } else {
                        curfin = elclib_in_period(
                            self.base.last_parameter_of(k),
                            spinedeb,
                            spinefin,
                        );
                    }
                    let mut curloi = LawConstant::new();
                    curloi.set(rdeb, curdeb, curfin);
                    funclist.push(Rc::new(RefCell::new(curloi)));
                    curdeb = curfin;
                    break;
                }
                k += 1;
            }
            if k > len {
                // no !
                if self.parandrad.is_empty() {
                    panic!("Standard_DomainError: Radius not defined");
                }
                let mut nbp = self.parandrad.len();
                if nbp > 1 {
                    deb = self.parandrad[0].x;
                    fin = deb + self.base.period();
                    if self.parandrad[self.parandrad.len() - 1].x - fin < -tol3d {
                        nbp += 1;
                    }
                } else {
                    nbp += 1;
                }
                // OCCT: NCollection_Array1<gp_Pnt2d> pr(1, nbp);
                let mut pr = vec![DVec2::ZERO; nbp];
                for l in 1..nbp {
                    pr[l - 1] = self.parandrad[l - 1];
                }
                pr[nbp - 1] = DVec2::new(fin, pr[0].y);
                let mut curloi = LawInterpol::new();
                // OCCT: curloi->Set(pr, true);
                curloi.set(&pr, true);
                funclist.push(Rc::new(RefCell::new(curloi)));
                // OCCT: return loi;
                *loi.change_laws() = funclist;
                loi.set_periodic();
                return Rc::new(RefCell::new(loi));
            }
        } else if self.base.is_periodic() {
            // start radius.
            if self.is_constant_on(ind[0]) {
                rdeb = self.radius_on(ind[0]);
                curfin = self.base.last_parameter_of(ind[0]);
                curfin = elclib_in_period(curfin, spinedeb + tol3d, spinefin + tol3d);
                curfin = fin.min(curfin);
                let mut curloi = LawConstant::new();
                curloi.set(rdeb, curdeb, curfin);
                funclist.push(Rc::new(RefCell::new(curloi)));
                curdeb = curfin;
                icur += 1;
            } else {
                // There is inevitably kpart right before!
                let mut iprec = ind[0] - 1;
                if iprec == 0 {
                    iprec = len;
                }
                if self.is_constant_on(iprec) {
                    rdeb = self.radius_on(iprec);
                } else {
                    panic!("Standard_DomainError: AppendLaw : previous constant is missing!");
                }
                lawencours = true;
            }
            // the raduis at end.
            if self.is_constant_on(ind[nbed - 1]) {
                rfin = self.radius_on(ind[nbed - 1]);
            } else {
                // There is inevitably kpart right after!
                let mut isuiv = ind[nbed - 1] + 1;
                if isuiv == len + 1 {
                    isuiv = 1;
                }
                if self.is_constant_on(isuiv) {
                    rfin = self.radius_on(isuiv);
                } else {
                    panic!("Standard_DomainError: AppendLaw : next constant is missing!");
                }
            }
        } else {
            // the radius at start.
            if self.is_constant_on(ind[0]) {
                rdeb = self.radius_on(ind[0]);
                curfin = fin.min(self.base.last_parameter_of(ind[0]));
                let mut curloi = LawConstant::new();
                curloi.set(rdeb, curdeb, curfin);
                funclist.push(Rc::new(RefCell::new(curloi)));
                curdeb = curfin;
                icur += 1;
            } else {
                if ind[0] > 1 {
                    if self.is_constant_on(ind[0] - 1) {
                        rdeb = self.radius_on(ind[0] - 1);
                    } else {
                        panic!(
                            "Standard_DomainError: AppendLaw : previous constant is missing"
                        );
                    }
                } else if self.parandrad.is_empty() {
                    panic!("Standard_DomainError: AppendLaw : no radius on vertex");
                } else {
                    rdeb = -1.0;
                }
                lawencours = true;
            }
            // the radius at end.
            if self.is_constant_on(ind[nbed - 1]) {
                rfin = self.radius_on(ind[nbed - 1]);
            } else {
                if ind[nbed - 1] < len {
                    if self.is_constant_on(ind[nbed - 1] + 1) {
                        rfin = self.radius_on(ind[nbed - 1] + 1);
                    } else {
                        panic!(
                            "Standard_DomainError: AppendLaw : next constant is missing"
                        );
                    }
                } else if self.parandrad.is_empty() {
                    panic!("Standard_DomainError: AppendLaw : no radius on vertex");
                } else {
                    rfin = -1.0;
                }
            }
        }

        // There are infos on the extremities of the elspine,
        // all edges are parsed
        while icur <= nbed {
            if self.is_constant_on(ind[icur - 1]) {
                rcur = self.radius_on(ind[icur - 1]);
                if lawencours {
                    let mut temp: Vec<LawFunctionHandle> = Vec::new();
                    mklaw(
                        &mut temp,
                        &self.parandrad,
                        curdeb,
                        curfin,
                        rdeb,
                        rcur,
                        self.base.is_periodic(),
                        spinedeb,
                        spinefin,
                        tol3d,
                    );
                    // OCCT: list.Append(temp); — list append of the whole
                    // temporary list.
                    funclist.extend(temp);
                    lawencours = false;
                    curdeb = curfin;
                }
                curfin = self.base.last_parameter_of(ind[icur - 1]);
                if self.base.is_periodic() {
                    curfin = elclib_in_period(curfin, spinedeb + tol3d, spinefin + tol3d);
                    if ind[icur - 1] == ind[nbed - 1] {
                        // Attention the curfin can be wrong if the last edge
                        // passes above the  origin periodic spline.
                        let mut biddeb = self.base.first_parameter_of(ind[icur - 1]);
                        biddeb =
                            elclib_in_period(biddeb, spinedeb + tol3d, spinefin + tol3d);
                        if biddeb >= curfin {
                            curfin = fin;
                        } else {
                            curfin = fin.min(curfin);
                        }
                    } else {
                        curfin = fin.min(curfin);
                    }
                }
                if (curfin - curdeb) > tol3d {
                    rdeb = rcur;
                    let mut curloi = LawConstant::new();
                    curloi.set(rdeb, curdeb, curfin);
                    funclist.push(Rc::new(RefCell::new(curloi)));
                    curdeb = curfin;
                }
            } else {
                curfin = self.base.last_parameter_of(ind[icur - 1]);
                if self.base.is_periodic() {
                    curfin = elclib_in_period(curfin, spinedeb + tol3d, spinefin + tol3d);
                }
                curfin = fin.min(curfin);
                lawencours = true;
                if ind[icur - 1] == ind[nbed - 1] {
                    // Attention the curfin can be wrong if the last edge
                    // passes above the  origin periodic spline.
                    if self.base.is_periodic() {
                        let mut biddeb = self.base.first_parameter_of(ind[icur - 1]);
                        curfin = self.base.last_parameter_of(ind[icur - 1]);
                        biddeb =
                            elclib_in_period(biddeb, spinedeb + tol3d, spinefin + tol3d);
                        curfin = elclib_in_period(curfin, spinedeb + tol3d, spinefin + tol3d);
                        if biddeb >= curfin {
                            curfin = fin;
                        } else {
                            curfin = fin.min(curfin);
                        }
                    }
                    // or if it is the end of spine with extension.
                    else if ind[icur - 1] == len {
                        curfin = fin;
                    }
                    let mut temp: Vec<LawFunctionHandle> = Vec::new();
                    mklaw(
                        &mut temp,
                        &self.parandrad,
                        curdeb,
                        curfin,
                        rdeb,
                        rfin,
                        self.base.is_periodic(),
                        spinedeb,
                        spinefin,
                        tol3d,
                    );
                    // OCCT: list.Append(temp);
                    funclist.extend(temp);
                }
            }
            icur += 1;
        }
        // OCCT: if (!lastloi.IsNull()) list.Append(lastloi); — lastloi is
        // never assigned in the OCCT body; the branch is kept literally.
        if let Some(ll) = lastloi.take() {
            funclist.push(ll);
        }
        *loi.change_laws() = funclist;
        if loi_periodic {
            loi.set_periodic();
        }
        Rc::new(RefCell::new(loi))
    }

    /// OCCT ChFiDS_FilSpine.cxx L874-889 — ChangeLaw(E).
    pub fn change_law(&mut self, e: &Shape) -> Option<LawFunctionHandle> {
        if !self.base.splitdone {
            panic!(
                "Standard_DomainError: ChFiDS_FilSpine::ChangeLaw : the limits are not up-to-date"
            );
        }
        let ie = self.base.index_of_edge(e);
        if self.is_constant_on(ie) {
            panic!(
                "Standard_DomainError: ChFiDS_FilSpine::ChangeLaw : no law on constant edges"
            );
        }
        // OCCT: occ::handle<ChFiDS_ElSpine> hsp = ElSpine(IE);
        let hsp = self
            .base
            .el_spine_of_index(ie)
            .expect("ElSpine(IE) (ChFiDS_Spine::ElSpine)");
        let w = 0.5 * (self.base.first_parameter_of(ie) + self.base.last_parameter_of(ie));
        // OCCT: occ::handle<Law_Composite> lc = Law(hsp);
        //       return lc->ChangeElementaryLaw(w);
        let lc = self.law_of(hsp);
        lc.and_then(|c| c.borrow_mut().change_elementary_law(w).cloned())
    }

    /// OCCT ChFiDS_FilSpine.cxx L893-929 — MaxRadFromSeqAndLaws().
    pub fn max_rad_from_seq_and_laws(&self) -> f64 {
        let mut max_rad = 0.0f64;

        for i in 1..=self.parandrad.len() {
            if self.parandrad[i - 1].y > max_rad {
                max_rad = self.parandrad[i - 1].y;
            }
        }

        let laws = self.filspine_laws(|l| l.clone());
        for law in &laws {
            let (mut fpar, mut lpar) = (0.0f64, 0.0f64);
            // OCCT: law->Bounds(fpar, lpar);
            law.borrow().bounds(&mut fpar, &mut lpar);
            let delta = (lpar - fpar) * 0.2;
            for i in 0..=4usize {
                let par = fpar + i as f64 * delta;
                // OCCT: rad = law->Value(par);
                let rad = law.borrow_mut().value(par);
                if rad > max_rad {
                    max_rad = rad;
                }
            }
            let rad = law.borrow_mut().value(lpar);
            if rad > max_rad {
                max_rad = rad;
            }
        }

        max_rad
    }
}

/// OCCT ChFiDS_FilSpine.cxx L383-526 — static mklaw(res, pr, curdeb, curfin,
/// Rdeb, Rfin, recadre, deb, fin, tol3d).
#[allow(clippy::too_many_arguments)]
fn mklaw(
    res: &mut Vec<LawFunctionHandle>,
    pr: &[DVec2],
    curdeb: f64,
    curfin: f64,
    rdeb: f64,
    rfin: f64,
    recadre: bool,
    deb: f64,
    fin: f64,
    tol3d: f64,
) {
    let mut npr: Vec<DVec2> = Vec::new();
    let mut rad = rdeb;
    let mut raf = rfin;
    let mut yaunpointsurledeb = false;
    let mut yaunpointsurlefin = false;
    if !pr.is_empty() {
        for i in 1..=pr.len() {
            let cur = pr[i - 1];
            let mut wcur = cur.x;
            if recadre {
                wcur = elclib_in_period(wcur, deb, fin);
            }
            if curdeb - tol3d <= wcur && wcur <= curfin + tol3d {
                if wcur - curdeb < tol3d {
                    yaunpointsurledeb = true;
                    let mut ncur = cur;
                    if rdeb < 0.0 {
                        rad = cur.y;
                    }
                    ncur.x = curdeb;
                    ncur.y = rad;
                    npr.push(ncur);
                } else if curfin - wcur < tol3d {
                    yaunpointsurlefin = true;
                    let mut ncur = cur;
                    if rfin < 0.0 {
                        raf = cur.y;
                    }
                    ncur.x = curfin;
                    ncur.y = raf;
                    npr.push(ncur);
                } else {
                    npr.push(DVec2::new(wcur, cur.y));
                }
            }
        }
    }

    if npr.is_empty() {
        if rdeb < 0.0 && rfin < 0.0 {
            panic!("Standard_DomainError: Impossible to create the law");
        } else if rdeb < 0.0 || rfin < 0.0 {
            let r = if rfin < 0.0 { rdeb } else { rfin };
            let mut loi = LawConstant::new();
            loi.set(r, curdeb, curfin);
            res.push(Rc::new(RefCell::new(loi)));
        } else {
            let mut loi = LawS::new();
            loi.set(curdeb, rdeb, curfin, rfin);
            res.push(Rc::new(RefCell::new(loi)));
        }
    } else {
        if !yaunpointsurledeb && rdeb >= 0.0 {
            npr.push(DVec2::new(curdeb, rdeb));
        }
        if !yaunpointsurlefin && rfin >= 0.0 {
            npr.push(DVec2::new(curfin, rfin));
        }
        let mut nbp = npr.len();
        // OCCT L473-484: the bubble sort on X.
        for i in 1..nbp {
            for j in (i + 1)..=nbp {
                if npr[i - 1].x > npr[j - 1].x {
                    npr.swap(i - 1, j - 1);
                }
            }
        }
        // Duplicates are removed. (OCCT L486-500)
        let mut fini = nbp <= 1;
        let mut i = 1usize;
        while !fini {
            if (npr[i - 1].x - npr[i].x).abs() < tol3d {
                npr.remove(i - 1);
                nbp -= 1;
            } else {
                i += 1;
            }
            fini = i >= nbp;
        }

        if rad < 0.0 {
            let mut loi = LawConstant::new();
            loi.set(npr[0].y, curdeb, npr[0].x);
            res.push(Rc::new(RefCell::new(loi)));
        }
        if nbp > 1 {
            // OCCT L510-517: NCollection_Array1<gp_Pnt2d> tpr(1, nbp);
            // curloi->Set(tpr, 0., 0., false);
            let mut curloi = LawInterpol::new();
            curloi.set_with_tangents(&npr[..nbp], 0.0, 0.0, false);
            res.push(Rc::new(RefCell::new(curloi)));
        }
        if raf < 0.0 {
            let mut loi = LawConstant::new();
            loi.set(npr[nbp - 1].y, npr[nbp - 1].x, curfin);
            res.push(Rc::new(RefCell::new(loi)));
        }
    }
}
