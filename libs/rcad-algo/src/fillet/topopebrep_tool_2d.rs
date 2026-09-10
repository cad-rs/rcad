//! OCCT TopOpeBRepTool_2d.cxx (TKBool/TopOpeBRepTool) — the FC2D
//! pcurve-on-face family, landed as the E3-N item (1) component of the
//! TopOpeBRepBuild_HBuilder carrier (the SplitFace1 re-landing
//! prerequisite; sanctioned by the D6 ruling that TKBool code translates
//! as a component of the consuming carrier).
//!
//! Translated here (all FC2D functions of TopOpeBRepTool_2d.cxx):
//! FC2D_FancestorE, FC2D_Prepare, FC2D_HasC3D, FC2D_HasCurveOnSurface,
//! FC2D_HasOldCurveOnSurface (both overloads), FC2D_PNewCurveOnSurface,
//! FC2D_HasNewCurveOnSurface (both overloads), FC2D_AddNewCurveOnSurface,
//! FC2D_make2d (both overloads), FC2D_MakeCurveOnSurface,
//! FC2D_CurveOnSurface (both overloads), FC2D_EditableCurveOnSurface,
//! FC2D_translate.
//!
//! Architecture differences (Rust <-> C++):
//! - Every function carries the owning `brep: &BRep` first parameter:
//!   OCCT BRep_Tool reads TShape data through global shape pointers,
//!   rcad reads through the BRep owner (the pcurve lookup composes the
//!   face/edge locations the way `bb_update_edge_pcurve` stores them).
//! - OCCT `occ::handle<Geom2d_Curve>` nullability -> `Option<Curve2d>`.
//! - OCCT out-parameters (`double&`) -> `&mut f64` / `&mut Option<...>`;
//!   OCCT leaves them uninitialized on the not-found path, Rust
//!   initializes the locals to 0.0 at the call sites.
//! - The OCCT C++ overloads (3-arg vs 6-arg HasOld/HasNew, 6-arg vs 7-arg
//!   CurveOnSurface, 5-arg vs 6-arg make2d) map to distinct Rust names
//!   (`_fltol` for the f2d/l2d/tol out-parameter form, `_ef` for the
//!   EF-neighbor form); the `f` out-parameter of the (E, F, f, l, tol,
//!   trim3d) overloads is named `f_out` (collides with the face `F`).
//! - OCCT function-static thread-locals -> `thread_local!` statics with
//!   the OCCT names.
//! - `FC2D_PNewCurveOnSurface` returns a mutable pointer into the static
//!   list in OCCT; Rust cannot hand out that reference across the
//!   thread-local borrow, so it returns a snapshot copy and the single
//!   mutating consumer (`FC2D_CurveOnSurface` 7-arg) writes back through
//!   a re-lookup of the stored element.
//! - The OCCT `#ifdef OCCT_DEBUG` trace blocks are compiled out and are
//!   not carried.

#![allow(dead_code)]
#![allow(non_upper_case_globals)]
// The OCCT bodies declare `bool hasold = false;` before an unconditional
// reassignment (TopOpeBRepTool_2d.cxx L164 / L386 / L396) — the initial
// value is never read there either.
#![allow(unused_assignments)]

use std::cell::{Cell, RefCell};
use std::collections::HashMap;

use glam::{DAffine3, DVec2};

use rcad_kernel::geom::{
    translate_curve2d, transform_curve, Curve2d, Curve2dEval as _, SurfaceEval as _,
};
use rcad_kernel::topods::{BRep, BRepTool as _, Orientation, Shape};

use super::topopebrep_tool_curve_tool::{
    brep_adaptor_surface, curve_tool_make_pcurve, curve_tool_make_pcurve_on_face,
    BRepAdaptorCurveOnFace,
};
use rcad_kernel::base::proj_lib::proj_lib_projected_curve_b::{GeomCurveHandle, GeomSurfaceHandle, ProjLibProjectedCurve};

// =========================================================================
// OCCT TopOpeBRepTool_2d.cxx L43-56 — the file-static thread-locals.
// =========================================================================

// OCCT L44-46: structure e -> C2D/F.  NCollection_DataMap<TopoDS_Shape,
// NCollection_List<TopOpeBRepTool_C2DF>, TopTools_ShapeMapHasher> keyed by
// the ShapeMapHasher identity (TShape + Location, orientation ignored) ->
// HashMap<(u64, u32), Vec<C2df>> (the hbuilder shape-key convention).
// The OCCT `*_ptr == nullptr` state -> `None`.
thread_local! {
    static GLOBAL_pmosloc2df: RefCell<Option<HashMap<(u64, u32), Vec<C2df>>>> =
        const { RefCell::new(None) };
    // OCCT L47: static thread_local int GLOBAL_C2D_i = 0; // DEB
    static GLOBAL_C2D_i: Cell<i32> = const { Cell::new(0) }; // DEB

    // OCCT L50-53: structure ancetre.  The IndexedDataMap insertion order
    // is never consumed (Contains / FindFromKey / Extent only) -> HashMap.
    static GLOBAL_pidmoslosc2df: RefCell<Option<HashMap<(u64, u32), Vec<Shape>>>> =
        const { RefCell::new(None) };
    // OCCT L54-56: GLOBAL_pFc2df / GLOBAL_pS1c2df / GLOBAL_pS2c2df.  The
    // OCCT nullified-face state (TopoDS_Face::Nullify) -> None.
    static GLOBAL_pFc2df: RefCell<Option<Shape>> = const { RefCell::new(None) };
    static GLOBAL_pS1c2df: RefCell<Option<Shape>> = const { RefCell::new(None) };
    static GLOBAL_pS2c2df: RefCell<Option<Shape>> = const { RefCell::new(None) };
}

/// OCCT TopTools_ShapeMapHasher identity: TShape + Location (orientation
/// ignored) — the rcad `(ptr_id, location)` key (classify.rs convention).
fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT TopoDS_Shape::IsEqual (TopoDS_Shape.hxx L276-280): TShape +
/// Location + Orientation.
fn shapes_equal(a: &Shape, b: &Shape) -> bool {
    shape_key(a) == shape_key(b) && a.orientation == b.orientation
}

/// OCCT TopoDS_Shape::Oriented(theOr) — a copy carrying the orientation.
fn shape_oriented(s: &Shape, the_or: Orientation) -> Shape {
    let mut c = s.clone();
    c.orientation = the_or;
    c
}

// =========================================================================
// GAP carrier of TopOpeBRepTool_C2DF (TKBool/TopOpeBRepTool untranslated;
// TopOpeBRepTool_C2DF.hxx L27-61, TopOpeBRepTool_C2DF.cxx L37-107).
//
// GAP: TopOpeBRepTool_C2DF (TKBool/TopOpeBRepTool untranslated) — the
// value type of the FC2D global pcurve map (TopOpeBRepTool_2d.cxx
// L44-46).  The FC2D family bodies need the type to exist to stay 1:1,
// so this carrier holds the OCCT private members myPC/myf2d/myl2d/mytol/
// myFace verbatim and its accessors reproduce the OCCT one-line bodies;
// the independent translation of TopOpeBRepTool_C2DF.cxx remains a
// pending dependency (see the D6 judgment note in the delivery report).
// =========================================================================

/// GAP carrier of `TopOpeBRepTool_C2DF` (TKBool/TopOpeBRepTool
/// untranslated; TopOpeBRepTool_C2DF.hxx L55-60 members).
#[derive(Debug, Clone)]
pub(crate) struct C2df {
    /// OCCT myPC: handle<Geom2d_Curve> (null -> None).
    pub(crate) my_pc: Option<Curve2d>,
    /// OCCT myf2d.
    pub(crate) myf2d: f64,
    /// OCCT myl2d.
    pub(crate) myl2d: f64,
    /// OCCT mytol.
    pub(crate) mytol: f64,
    /// OCCT myFace.
    pub(crate) my_face: Shape,
}

impl C2df {
    /// GAP carrier of the OCCT constructor
    /// TopOpeBRepTool_C2DF(PC, f2d, l2d, tol, F) (TopOpeBRepTool_C2DF.cxx
    /// L41-52): member assignments only.
    pub(crate) fn new(
        my_pc: Option<Curve2d>,
        myf2d: f64,
        myl2d: f64,
        mytol: f64,
        my_face: &Shape,
    ) -> Self {
        C2df {
            my_pc,
            myf2d,
            myl2d,
            mytol,
            my_face: my_face.clone(),
        }
    }

    /// GAP carrier of TopOpeBRepTool_C2DF::SetPC (TopOpeBRepTool_C2DF.cxx
    /// L56-65).
    pub(crate) fn set_pc(&mut self, my_pc: Option<Curve2d>, myf2d: f64, myl2d: f64, mytol: f64) {
        self.my_pc = my_pc;
        self.myf2d = myf2d;
        self.myl2d = myl2d;
        self.mytol = mytol;
    }

    /// GAP carrier of TopOpeBRepTool_C2DF::PC (TopOpeBRepTool_C2DF.cxx
    /// L76-84): the out-parameters receive the stored range/tolerance and
    /// the pcurve handle is returned (None when null).
    pub(crate) fn pc(&self, f2d: &mut f64, l2d: &mut f64, tol: &mut f64) -> Option<Curve2d> {
        *f2d = self.myf2d;
        *l2d = self.myl2d;
        *tol = self.mytol;
        self.my_pc.clone()
    }

    /// GAP carrier of TopOpeBRepTool_C2DF::IsFace (TopOpeBRepTool_C2DF.cxx
    /// L103-107): F.IsEqual(myFace).
    pub(crate) fn is_face(&self, f: &Shape) -> bool {
        let b = shapes_equal(f, &self.my_face);
        b
    }
}

// =========================================================================
// GAP carriers of the projection machinery (external untranslated
// dependencies of the FC2D projection arm).
// =========================================================================

/// GAP: TopExp::MapShapesAndAncestors (TKBRep/TopExp.cxx L80-120 — no
/// shared rcad rehost; the existing per-file copies are private).  The
/// carrier leaves the ancestor map empty, which preserves the OCCT
/// outcome for an edge without recorded ancestors: FC2D_FancestorE
/// returns the null face (TopOpeBRepTool_2d.cxx L82-90) and the no-3D-
/// curve projection arm returns a null pcurve.
fn gap_map_shapes_and_ancestors(
    _brep: &BRep,
    _s1: Option<&Shape>,
    _s2: Option<&Shape>,
    _map: &mut HashMap<(u64, u32), Vec<Shape>>,
) {
    // GAP body: no rcad equivalent of TopExp::MapShapesAndAncestors; the
    // OCCT ancestor map stays empty (the preserved OCCT failure path is
    // the FC2D_FancestorE null-face return).
}

/// GAP: TopOpeBRepTool_TOOL::UVISO (TKBool/TopOpeBRepTool_TOOL.cxx
/// L1191-1221 untranslated).  OCCT failure path preserved: a null PC
/// (and the non-line / non-iso outcomes) returns false with isoU = isoV
/// = false (OCCT L1194-1196); the carrier always takes that path, so
/// FC2D_translate never shifts a periodic-surface pcurve (the OCCT
/// behavior for every false outcome of UVISO).
fn gap_tool_uviso(
    pc: Option<&Curve2d>,
    iso_u: &mut bool,
    iso_v: &mut bool,
    d2d: &mut DVec2,
    o2d: &mut DVec2,
) -> bool {
    *iso_u = false;
    *iso_v = false;
    let _ = (pc, d2d, o2d);
    false
}

/// GAP: FTOL_FaceTolerances3d (TKBool/TopOpeBRepTool_tol.cxx L221-260
/// untranslated; it reads the FBOX_GetHBoxTool face boxes and reduces
/// through FTOL_FaceTolerances).  OCCT failure path preserved: the carrier
/// returns the 0.0 initialization of the OCCT call-site local `double
/// tolin;` (TopOpeBRepTool_2d.cxx L320-321 / L516-517); the
/// ProjLib_ProjectedCurve
/// tolerance constructor clamps with max(Tol, Precision::Confusion())
/// (ProjLib_ProjectedCurve.cxx L317), so the projection runs at the
/// Confusion tolerance — no FTOL box-based reduction is applied.
fn gap_ftol_face_tolerances3d(_brep: &BRep, _f: &Shape, _fe: &Shape) -> f64 {
    // GAP body: the box-based tolerance reduction is untranslated; the
    // value feeds the ProjLib_ProjectedCurve tolerance constructor.
    0.0
}

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L61-93
fn fc2d_fancestor_e(brep: &BRep, e: &Shape) -> Shape {
    // OCCT L63-68: lazy allocation of GLOBAL_pmosloc2df.
    if GLOBAL_pmosloc2df.with(|m| m.borrow().is_none()) {
        GLOBAL_pmosloc2df.with(|m| *m.borrow_mut() = Some(HashMap::new()));
    }
    // OCCT L69: int ancemp = (*GLOBAL_pidmoslosc2df).Extent();
    let ancemp = GLOBAL_pidmoslosc2df.with(|m| m.borrow().as_ref().map_or(0, |x| x.len()));
    if ancemp == 0 {
        // OCCT L72-79: TopExp::MapShapesAndAncestors(S1/S2, TopAbs_EDGE,
        // TopAbs_FACE, map) — GAP carrier (the map stays empty).
        let s1 = GLOBAL_pS1c2df.with(|g| g.borrow().clone());
        let s2 = GLOBAL_pS2c2df.with(|g| g.borrow().clone());
        GLOBAL_pidmoslosc2df.with(|m| {
            let mut m = m.borrow_mut();
            let map = m.get_or_insert_with(HashMap::new);
            gap_map_shapes_and_ancestors(brep, s1.as_ref(), s2.as_ref(), map);
        });
    }
    // OCCT L81-85: bool Eb = map.Contains(E); if (!Eb) return
    // *GLOBAL_pFc2df (the nullified face).
    let eb = GLOBAL_pidmoslosc2df
        .with(|m| m.borrow().as_ref().map_or(false, |x| x.contains_key(&shape_key(e))));
    if !eb {
        return GLOBAL_pFc2df.with(|g| g.borrow().clone().unwrap_or_else(Shape::null));
    }
    // OCCT L86-90: const NCollection_List<TopoDS_Shape>& lf =
    // map.FindFromKey(E); if (lf.IsEmpty()) return *GLOBAL_pFc2df;
    let lf = GLOBAL_pidmoslosc2df.with(|m| {
        m.borrow()
            .as_ref()
            .and_then(|x| x.get(&shape_key(e)).cloned())
    });
    let Some(lf) = lf else {
        return GLOBAL_pFc2df.with(|g| g.borrow().clone().unwrap_or_else(Shape::null));
    };
    if lf.is_empty() {
        return GLOBAL_pFc2df.with(|g| g.borrow().clone().unwrap_or_else(Shape::null));
    }
    // OCCT L91-92: const TopoDS_Face& F = TopoDS::Face(lf.First());
    // return F;  (rcad Shape is untyped until its data is read — the
    // cast is the type assertion, architecture difference.)
    let f = lf[0].clone();
    f
}

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L96-134
pub(crate) fn fc2d_prepare(s1: &Shape, s2: &Shape) -> i32 {
    // OCCT L98-104: allocate then Clear GLOBAL_pmosloc2df.
    GLOBAL_pmosloc2df.with(|m| {
        let mut m = m.borrow_mut();
        if m.is_none() {
            *m = Some(HashMap::new());
        }
        m.as_mut().unwrap().clear();
    });
    // OCCT L105: GLOBAL_C2D_i = 0;
    GLOBAL_C2D_i.with(|i| i.set(0));
    // OCCT L107-113: allocate then Clear GLOBAL_pidmoslosc2df.
    GLOBAL_pidmoslosc2df.with(|m| {
        let mut m = m.borrow_mut();
        if m.is_none() {
            *m = Some(HashMap::new());
        }
        m.as_mut().unwrap().clear();
    });
    // OCCT L115-119: allocate GLOBAL_pFc2df then Nullify().
    GLOBAL_pFc2df.with(|g| *g.borrow_mut() = None);
    // OCCT L121-125: allocate GLOBAL_pS1c2df then *GLOBAL_pS1c2df = S1;
    GLOBAL_pS1c2df.with(|g| *g.borrow_mut() = Some(s1.clone()));
    // OCCT L127-131: allocate GLOBAL_pS2c2df then *GLOBAL_pS2c2df = S2;
    GLOBAL_pS2c2df.with(|g| *g.borrow_mut() = Some(s2.clone()));
    // OCCT L133: return 0;
    0
}

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L137-144
pub(crate) fn fc2d_has_c3d(brep: &BRep, e: &Shape) -> bool {
    // OCCT L139-141: TopLoc_Location loc; double f3d, l3d;
    // C3D = BRep_Tool::Curve(E, loc, f3d, l3d);
    let c3d = brep.edge_curve_data(e);
    // OCCT L142-143: bool b = (!C3D.IsNull()); return b;
    let b = c3d.is_some();
    b
}

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L147-154
pub(crate) fn fc2d_has_curve_on_surface(brep: &BRep, e: &Shape, f: &Shape) -> bool {
    // OCCT L149: occ::handle<Geom2d_Curve> C2D;
    let mut c2d: Option<Curve2d> = None;
    // OCCT L150-151.
    let hasold = fc2d_has_old_curve_on_surface(brep, e, f, &mut c2d);
    let hasnew = fc2d_has_new_curve_on_surface(e, f, &mut c2d);
    // OCCT L152-153: bool b = hasold || hasnew; return b;
    let b = hasold || hasnew;
    b
}

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L157-169 — the 6-argument overload (C2D,
// f2d, l2d, tol out-parameters; Rust name adds `_fltol`).
fn fc2d_has_old_curve_on_surface_fltol(
    brep: &BRep,
    e: &Shape,
    f: &Shape,
    c2d: &mut Option<Curve2d>,
    f2d: &mut f64,
    l2d: &mut f64,
    tol: &mut f64,
) -> bool {
    // OCCT L164: bool hasold = false;
    let mut hasold = false;
    // OCCT L165: tol = BRep_Tool::Tolerance(E);
    *tol = brep.tolerance(e);
    // OCCT L166: C2D = BRep_Tool::CurveOnSurface(E, F, f2d, l2d);
    match brep.curve_on_surface(e, f) {
        Some((curv, first, last)) => {
            *c2d = Some(curv);
            *f2d = first;
            *l2d = last;
        }
        None => {
            *c2d = None;
        }
    }
    // OCCT L167-168: hasold = (!C2D.IsNull()); return hasold;
    hasold = c2d.is_some();
    hasold
}

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L171-178 — the 3-argument overload.
fn fc2d_has_old_curve_on_surface(
    brep: &BRep,
    e: &Shape,
    f: &Shape,
    c2d: &mut Option<Curve2d>,
) -> bool {
    // OCCT L175: double f2d, l2d, tol;
    let mut f2d: f64 = 0.0; // OCCT leaves the locals uninitialized
    let mut l2d: f64 = 0.0;
    let mut tol: f64 = 0.0;
    // OCCT L176-177.
    let hasold = fc2d_has_old_curve_on_surface_fltol(brep, e, f, c2d, &mut f2d, &mut l2d, &mut tol);
    hasold
}

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L181-205
//
// Architecture difference: OCCT returns `TopOpeBRepTool_C2DF*` — a
// mutable pointer INTO the thread-local list; Rust cannot hand that
// reference across the thread-local borrow, so this returns a snapshot
// copy.  The only OCCT consumers read PC() (a copy) or SetPC() — the
// mutating consumer (FC2D_CurveOnSurface 7-arg) writes back through a
// re-lookup of the stored element.
fn fc2d_p_new_curve_on_surface(e: &Shape, f: &Shape) -> Option<C2df> {
    // OCCT L183: TopOpeBRepTool_C2DF* pc2df = nullptr;
    let mut pc2df: Option<C2df> = None;
    // OCCT L184-187: if (GLOBAL_pmosloc2df == nullptr) return nullptr;
    let map_none = GLOBAL_pmosloc2df.with(|m| m.borrow().is_none());
    if map_none {
        return None;
    }
    // OCCT L188-192: bool Eisb = GLOBAL_pmosloc2df->IsBound(E);
    // if (!Eisb) return nullptr;
    let eisb = GLOBAL_pmosloc2df
        .with(|m| m.borrow().as_ref().map_or(false, |x| x.contains_key(&shape_key(e))));
    if !eisb {
        return None;
    }
    // OCCT L193-204: iterate the bound list; the c2df with
    // IsFace(F) -> pc2df = &c2df; break;
    GLOBAL_pmosloc2df.with(|m| {
        let m = m.borrow();
        let lc2df = m.as_ref().unwrap().get(&shape_key(e)).unwrap();
        // OCCT: NCollection_List<TopOpeBRepTool_C2DF>::Iterator it(...).
        for c2df in lc2df {
            let isf = c2df.is_face(f);
            if isf {
                pc2df = Some(c2df.clone());
                break;
            }
        }
    });
    // OCCT L204: return pc2df;
    pc2df
}

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L207-221 — the 6-argument overload.
fn fc2d_has_new_curve_on_surface_fltol(
    e: &Shape,
    f: &Shape,
    c2d: &mut Option<Curve2d>,
    f2d: &mut f64,
    l2d: &mut f64,
    tol: &mut f64,
) -> bool {
    // OCCT L214: const TopOpeBRepTool_C2DF* pc2df = FC2D_PNewCurveOnSurface(E, F);
    let pc2df = fc2d_p_new_curve_on_surface(e, f);
    // OCCT L215: bool hasnew = (pc2df != nullptr);
    let hasnew = pc2df.is_some();
    // OCCT L216-219: if (hasnew) { C2D = pc2df->PC(f2d, l2d, tol); }
    if let Some(c2df) = &pc2df {
        *c2d = c2df.pc(f2d, l2d, tol);
    }
    // OCCT L220: return hasnew;
    hasnew
}

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L223-230 — the 3-argument overload.
fn fc2d_has_new_curve_on_surface(e: &Shape, f: &Shape, c2d: &mut Option<Curve2d>) -> bool {
    // OCCT L227: double f2d, l2d, tol;
    let mut f2d: f64 = 0.0; // OCCT leaves the locals uninitialized
    let mut l2d: f64 = 0.0;
    let mut tol: f64 = 0.0;
    // OCCT L228-229.
    let b = fc2d_has_new_curve_on_surface_fltol(e, f, c2d, &mut f2d, &mut l2d, &mut tol);
    b
}

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L233-254
fn fc2d_add_new_curve_on_surface(
    c2d: &Option<Curve2d>,
    e: &Shape,
    f: &Shape,
    f2d: f64,
    l2d: f64,
    tol: f64,
) -> i32 {
    // OCCT L240-243: if (C2D.IsNull()) return 1;
    if c2d.is_none() {
        return 1;
    }
    // OCCT L244: TopOpeBRepTool_C2DF c2df(C2D, f2d, l2d, tol, F); — the
    // GAP-carrier constructor (C2df::new).
    let c2df = C2df::new(c2d.clone(), f2d, l2d, tol, f);
    // OCCT L245-248: if (GLOBAL_pmosloc2df == nullptr) return 1;
    let map_none = GLOBAL_pmosloc2df.with(|m| m.borrow().is_none());
    if map_none {
        return 1;
    }
    // OCCT L249-252: NCollection_List<TopOpeBRepTool_C2DF> thelist;
    // GLOBAL_pmosloc2df->Bind(E, thelist);  (Bind rebinds — any previous
    // binding of E is replaced.)  NCollection_List<TopOpeBRepTool_C2DF>&
    // lc2df = GLOBAL_pmosloc2df->ChangeFind(E); lc2df.Append(c2df);
    GLOBAL_pmosloc2df.with(|m| {
        let mut m = m.borrow_mut();
        let map = m.as_mut().unwrap();
        let key = shape_key(e);
        map.insert(key, Vec::new());
        let lc2df = map.get_mut(&key).unwrap();
        lc2df.push(c2df);
    });
    // OCCT L253: return 0;
    0
}

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L257-338 (forward declaration L257-262 is
// the C++ prototype; Rust needs none).  The `trim3d` default argument
// (false) is explicit at the call sites.
fn fc2d_make2d(
    brep: &BRep,
    e: &Shape,
    f: &Shape,
    f2d: &mut f64,
    l2d: &mut f64,
    tol: &mut f64,
    trim3d: bool,
) -> Option<Curve2d> {
    // OCCT L271: C2D = BRep_Tool::CurveOnSurface(E, F, f2d, l2d);
    // OCCT L272-275: if (!C2D.IsNull()) return C2D;
    if let Some((c, first, last)) = brep.curve_on_surface(e, f) {
        *f2d = first;
        *l2d = last;
        return Some(c);
    }
    // OCCT L270: occ::handle<Geom2d_Curve> C2D;
    let mut c2d: Option<Curve2d> = None;

    // pas de 2D
    // OCCT L278-281: double f3d, l3d; TopLoc_Location eloc;
    // C1 = BRep_Tool::Curve(E, eloc, f3d, l3d); bool hasC3D = (!C1.IsNull());
    let c1 = brep.edge_curve_data(e);
    let eloc: DAffine3 = brep.get_location(e.location);
    let (f3d, l3d) = {
        let range = brep.edge_range(e);
        (range[0], range[1])
    };
    let has_c3d = c1.is_some();

    if has_c3d {
        // OCCT L285: bool elocid = eloc.IsIdentity();
        let elocid = eloc == DAffine3::IDENTITY;
        // OCCT L286-294: C2 = elocid ? C1 : down_cast<Geom_Curve>(C1->Transformed(eloc.Transformation()));
        let c1 = c1.unwrap();
        let c2 = if elocid {
            c1
        } else {
            transform_curve(&c1, &eloc)
        };
        // OCCT L295-300: double f = 0., l = 0.; if (trim3d) { f = f3d; l = l3d; }
        // (the OCCT locals f/l are renamed f_par/l_par — they collide with
        // the face parameter `f`.)
        let mut f_par: f64 = 0.0;
        let mut l_par: f64 = 0.0;
        if trim3d {
            f_par = f3d;
            l_par = l3d;
        }
        // OCCT L301: C2D = TopOpeBRepTool_CurveTool::MakePCurveOnFace(F, C2, tol, f, l);
        c2d = curve_tool_make_pcurve_on_face(brep, f, &c2, tol, f_par, l_par);
        // OCCT L302-303: f2d = f3d; l2d = l3d;
        *f2d = f3d;
        *l2d = l3d;
        // OCCT L304: return C2D;
        return c2d;
    } else {
        // E sans courbe 2d sur F, E sans courbe 3d
        // une face accedant a E : FE
        // OCCT L310: const TopoDS_Face& FE = FC2D_FancestorE(E);
        let fe = fc2d_fancestor_e(brep, e);
        // OCCT L311-314: if (FE.IsNull()) return C2D;
        if fe.is_null() {
            return c2d;
        }
        // OCCT L315: bool compminmaxUV = false;
        let compminmax_uv = false;
        // OCCT L316: BRepAdaptor_Surface BAS(F, compminmaxUV);
        let bas = brep_adaptor_surface(brep, f, compminmax_uv);
        // OCCT L317: occ::handle<BRepAdaptor_Surface> BAHS = new BRepAdaptor_Surface(BAS);
        let bahs: GeomSurfaceHandle = std::sync::Arc::new(bas.clone());
        // OCCT L318: BRepAdaptor_Curve AC(E, FE);
        let ac = BRepAdaptorCurveOnFace::new(brep, e, &fe);
        // OCCT L319: occ::handle<BRepAdaptor_Curve> AHC = new BRepAdaptor_Curve(AC);
        let ahc: GeomCurveHandle = std::sync::Arc::new(ac);
        // OCCT L320-321: double tolin; FTOL_FaceTolerances3d(F, FE, tolin);
        let tolin = gap_ftol_face_tolerances3d(brep, f, &fe);
        // OCCT L322: ProjLib_ProjectedCurve projcurv(BAHS, AHC, tolin);
        let projcurv = ProjLibProjectedCurve::with_surface_curve_tol(bahs, ahc, tolin);
        // OCCT L323: C2D = MakePCurve(projcurv);
        c2d = curve_tool_make_pcurve(&projcurv);
        // OCCT L324-327: double f, l; BRep_Tool::Range(E, f, l); f2d = f; l2d = l;
        let (f_par, l_par) = {
            let range = brep.edge_range(e);
            (range[0], range[1])
        };
        *f2d = f_par;
        *l2d = l_par;
    }

    // OCCT L330-335: the #ifdef OCCT_DEBUG trace block is compiled out.
    // OCCT L337: return C2D;
    c2d
} // make2d1

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L342-353
// modified by NIZHNY-MZV  Mon Oct  4 10:37:36 1999
pub(crate) fn fc2d_make_curve_on_surface(
    brep: &BRep,
    e: &Shape,
    f: &Shape,
    f_out: &mut f64,
    l_out: &mut f64,
    tol: &mut f64,
    trim3d: bool,
) -> Option<Curve2d> {
    // OCCT L350: C2D = FC2D_make2d(E, F, f, l, tol, trim3d);
    let c2d = fc2d_make2d(brep, e, f, f_out, l_out, tol, trim3d);
    // OCCT L351: FC2D_AddNewCurveOnSurface(C2D, E, F, f, l, tol);
    fc2d_add_new_curve_on_surface(&c2d, e, f, *f_out, *l_out, *tol);
    // OCCT L352: return C2D;
    c2d
}

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L356-376 — the 6-argument overload (no EF;
// `trim3d` defaults to false in OCCT and is explicit here).
pub(crate) fn fc2d_curve_on_surface(
    brep: &BRep,
    e: &Shape,
    f: &Shape,
    f_out: &mut f64,
    l_out: &mut f64,
    tol: &mut f64,
    trim3d: bool,
) -> Option<Curve2d> {
    // OCCT L363: occ::handle<Geom2d_Curve> C2D;
    let mut c2d: Option<Curve2d> = None;
    // OCCT L364-368.
    let hasold = fc2d_has_old_curve_on_surface_fltol(brep, e, f, &mut c2d, f_out, l_out, tol);
    if hasold {
        return c2d;
    }
    // OCCT L369-373.
    let hasnew = fc2d_has_new_curve_on_surface_fltol(e, f, &mut c2d, f_out, l_out, tol);
    if hasnew {
        return c2d;
    }
    // OCCT L374-375: C2D = FC2D_MakeCurveOnSurface(E, F, f, l, tol, trim3d); return C2D;
    c2d = fc2d_make_curve_on_surface(brep, e, f, f_out, l_out, tol, trim3d);
    c2d
}

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L379-407
pub(crate) fn fc2d_editable_curve_on_surface(
    brep: &BRep,
    e: &Shape,
    f: &Shape,
    f_out: &mut f64,
    l_out: &mut f64,
    tol: &mut f64,
    trim3d: bool,
) -> Option<Curve2d> {
    // OCCT L386: bool hasold = false;
    let mut hasold = false;
    {
        // OCCT L388: occ::handle<Geom2d_Curve> C2D;
        let mut c2d: Option<Curve2d> = None;
        // OCCT L389-394: hasold = FC2D_HasOldCurveOnSurface(E, F, C2D, f, l, tol);
        // if (hasold) { copC2D = down_cast<Geom2d_Curve>(C2D->Copy()); return copC2D; }
        hasold = fc2d_has_old_curve_on_surface_fltol(brep, e, f, &mut c2d, f_out, l_out, tol);
        if hasold {
            // OCCT L392: C2D->Copy() — rcad Curve2d is a value type; the
            // clone is the copy (architecture difference).
            let cop_c2d = c2d.clone();
            return cop_c2d;
        }
    }
    // OCCT L396: bool hasnew = false;
    let mut hasnew = false;
    {
        // OCCT L398: occ::handle<Geom2d_Curve> newC2D;
        let mut new_c2d: Option<Curve2d> = None;
        // OCCT L399-403.
        hasnew = fc2d_has_new_curve_on_surface_fltol(e, f, &mut new_c2d, f_out, l_out, tol);
        if hasnew {
            return new_c2d;
        }
    }
    // OCCT L405-406: makC2D = FC2D_MakeCurveOnSurface(E, F, f, l, tol, trim3d); return makC2D;
    let mak_c2d = fc2d_make_curve_on_surface(brep, e, f, f_out, l_out, tol, trim3d);
    mak_c2d
}

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L410-447 — the edge parameter is unnamed in
// the OCCT definition (L412: `const TopoDS_Edge&,`).
fn fc2d_translate(brep: &BRep, c2d: &mut Option<Curve2d>, _e: &Shape, f: &Shape, ef: &Shape) {
    // OCCT L416-418: TopLoc_Location sloc;
    // const handle<Geom_Surface>& S1 = BRep_Tool::Surface(F, sloc);
    // bool isperio = S1->IsUPeriodic() || S1->IsVPeriodic();
    // (a face without surface would dereference null in OCCT; rcad reads
    // through the Option and keeps OCCT's false outcome — architecture
    // difference.)
    let isperio = brep
        .face_surface(f)
        .map(|s1| s1.is_u_periodic() || s1.is_v_periodic())
        .unwrap_or(false);
    // OCCT L419-421: gp_Dir2d d2d; gp_Pnt2d O2d; bool isuiso, isviso;
    let mut d2d: DVec2 = DVec2::ZERO;
    let mut o2d: DVec2 = DVec2::ZERO;
    let mut isuiso = false;
    let mut isviso = false;
    // OCCT L422: bool uviso = TopOpeBRepTool_TOOL::UVISO(C2D, isuiso, isviso, d2d, O2d);
    // — GAP carrier (TopOpeBRepTool_TOOL untranslated).
    let uviso = gap_tool_uviso(c2d.as_ref(), &mut isuiso, &mut isviso, &mut d2d, &mut o2d);
    // OCCT L423: bool EFnull = EF.IsNull();
    let efnull = ef.is_null();

    // OCCT L425: if (isperio && uviso && !EFnull)
    if isperio && uviso && !efnull {
        // C2D prend comme origine dans F l'origine de la pcurve de EF dans F
        // OCCT L428-429: TopoDS_Face FFOR = F; FFOR.Orientation(TopAbs_FORWARD);
        let ffor = shape_oriented(f, Orientation::Forward);
        // OCCT L431: BRep_Tool::UVPoints(EF, FFOR, p1, p2) (BRep_Tool.cxx
        // L1090-1108) — the pcurve endpoint read (PFirst/PLast).  The read
        // goes through the same composed-location pcurve lookup the
        // hasold arm uses (the brep_algo::tool rehost keys the pcurve by
        // the raw face key — architecture difference).
        let (p1, _p2) = match brep.curve_on_surface(ef, &ffor) {
            Some((pc, first, last)) => (pc.point_at(first), pc.point_at(last)),
            None => (DVec2::ZERO, DVec2::ZERO),
        };
        // OCCT L432-435.
        let p_ef = if isuiso { p1.x } else { p1.y };
        let p_c2d = if isuiso { o2d.x } else { o2d.y };
        let factor = p_ef - p_c2d;
        let b = factor.abs() > 1.0e-6;
        // OCCT L436-445.
        if b {
            // OCCT L438-442: gp_Vec2d transl(1., 0.); if (isviso) transl = gp_Vec2d(0., 1.);
            let mut transl = DVec2::new(1.0, 0.0);
            if isviso {
                transl = DVec2::new(0.0, 1.0);
            }
            // OCCT L443: transl.Multiply(factor);
            transl = transl * factor;
            // OCCT L444: C2D->Translate(transl);
            if let Some(curv) = c2d.as_mut() {
                *curv = translate_curve2d(curv, transl);
            }
        }
    }
}

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L450-535 (forward declaration L450-456 is
// the C++ prototype; Rust needs none).
fn fc2d_make2d_ef(
    brep: &BRep,
    e: &Shape,
    f: &Shape,
    ef: &Shape,
    f2d: &mut f64,
    l2d: &mut f64,
    tol: &mut f64,
    trim3d: bool,
) -> Option<Curve2d> {
    // OCCT L466-469: C2D = BRep_Tool::CurveOnSurface(E, F, f2d, l2d);
    // if (!C2D.IsNull()) return C2D;
    if let Some((c, first, last)) = brep.curve_on_surface(e, f) {
        *f2d = first;
        *l2d = last;
        return Some(c);
    }
    // OCCT L465: occ::handle<Geom2d_Curve> C2D;
    let mut c2d: Option<Curve2d> = None;

    // pas de 2D
    // OCCT L473-476.
    let c1 = brep.edge_curve_data(e);
    let eloc: DAffine3 = brep.get_location(e.location);
    let (f3d, l3d) = {
        let range = brep.edge_range(e);
        (range[0], range[1])
    };
    let has_c3d = c1.is_some();

    if has_c3d {
        // OCCT L480: bool elocid = eloc.IsIdentity();
        let elocid = eloc == DAffine3::IDENTITY;
        // OCCT L481-489.
        let c1 = c1.unwrap();
        let c2 = if elocid {
            c1
        } else {
            transform_curve(&c1, &eloc)
        };
        // OCCT L490-495: double f = 0., l = 0.; if (trim3d) { f = f3d; l = l3d; }
        // (the OCCT locals f/l are renamed f_par/l_par — they collide with
        // the face parameter `f`.)
        let mut f_par: f64 = 0.0;
        let mut l_par: f64 = 0.0;
        if trim3d {
            f_par = f3d;
            l_par = l3d;
        }
        // OCCT L496: C2D = TopOpeBRepTool_CurveTool::MakePCurveOnFace(F, C2, tol, f, l);
        c2d = curve_tool_make_pcurve_on_face(brep, f, &c2, tol, f_par, l_par);
        // OCCT L497-498: f2d = f3d; l2d = l3d;
        *f2d = f3d;
        *l2d = l3d;
        // OCCT L499: FC2D_translate(C2D, E, F, EF);
        fc2d_translate(brep, &mut c2d, e, f, ef);
        // OCCT L500: return C2D;
        return c2d;
    } else {
        // E sans courbe 2d sur F, E sans courbe 3d
        // une face accedant a E : FE
        // OCCT L506: const TopoDS_Face& FE = FC2D_FancestorE(E);
        let fe = fc2d_fancestor_e(brep, e);
        // OCCT L507-510.
        if fe.is_null() {
            return c2d;
        }
        // OCCT L511: bool compminmaxUV = false;
        let compminmax_uv = false;
        // OCCT L512: BRepAdaptor_Surface BAS(F, compminmaxUV);
        let bas = brep_adaptor_surface(brep, f, compminmax_uv);
        // OCCT L513: occ::handle<BRepAdaptor_Surface> BAHS = new BRepAdaptor_Surface(BAS);
        let bahs: GeomSurfaceHandle = std::sync::Arc::new(bas.clone());
        // OCCT L514: BRepAdaptor_Curve AC(E, FE);
        let ac = BRepAdaptorCurveOnFace::new(brep, e, &fe);
        // OCCT L515: occ::handle<BRepAdaptor_Curve> AHC = new BRepAdaptor_Curve(AC);
        let ahc: GeomCurveHandle = std::sync::Arc::new(ac);
        // OCCT L516-517: double tolin; FTOL_FaceTolerances3d(F, FE, tolin);
        let tolin = gap_ftol_face_tolerances3d(brep, f, &fe);
        // OCCT L518: ProjLib_ProjectedCurve projcurv(BAHS, AHC, tolin);
        let projcurv = ProjLibProjectedCurve::with_surface_curve_tol(bahs, ahc, tolin);
        // OCCT L519: C2D = MakePCurve(projcurv);
        c2d = curve_tool_make_pcurve(&projcurv);
        // OCCT L520-523: double f, l; BRep_Tool::Range(E, f, l); f2d = f; l2d = l;
        let (f_par, l_par) = {
            let range = brep.edge_range(e);
            (range[0], range[1])
        };
        *f2d = f_par;
        *l2d = l_par;
        // OCCT L524: FC2D_translate(C2D, E, F, EF);
        fc2d_translate(brep, &mut c2d, e, f, ef);
    }

    // OCCT L527-532: the #ifdef OCCT_DEBUG trace block is compiled out.
    // OCCT L534: return C2D;
    c2d
} // make2d2

// ------------------------------------------------------------------------------------
// OCCT TopOpeBRepTool_2d.cxx L538-566 — the 7-argument overload (with the
// EF neighbor edge).
pub(crate) fn fc2d_curve_on_surface_ef(
    brep: &BRep,
    e: &Shape,
    f: &Shape,
    ef: &Shape,
    f2d: &mut f64,
    l2d: &mut f64,
    tol: &mut f64,
    trim3d: bool,
) -> Option<Curve2d> {
    // OCCT L546: occ::handle<Geom2d_Curve> C2D;
    let mut c2d: Option<Curve2d> = None;

    // OCCT L548-552.
    let hasold = fc2d_has_old_curve_on_surface_fltol(brep, e, f, &mut c2d, f2d, l2d, tol);
    if hasold {
        return c2d;
    }

    // OCCT L554: TopOpeBRepTool_C2DF* pc2df = FC2D_PNewCurveOnSurface(E, F);
    // (the snapshot-copy architecture difference is documented on
    // fc2d_p_new_curve_on_surface).
    let pc2df = fc2d_p_new_curve_on_surface(e, f);
    // OCCT L555-561: if (pc2df != nullptr) { C2D = pc2df->PC(f2d, l2d, tol);
    // FC2D_translate(C2D, E, F, EF); pc2df->SetPC(C2D, f2d, l2d, tol); return C2D; }
    if let Some(c2df) = &pc2df {
        c2d = c2df.pc(f2d, l2d, tol);
        fc2d_translate(brep, &mut c2d, e, f, ef);
        // OCCT L559: pc2df->SetPC(C2D, f2d, l2d, tol) — the write-back
        // through the stored element (OCCT holds the raw pointer; rcad
        // re-finds the list element — architecture difference).
        GLOBAL_pmosloc2df.with(|m| {
            let mut m = m.borrow_mut();
            let lc2df = m.as_mut().unwrap().get_mut(&shape_key(e)).unwrap();
            for stored in lc2df.iter_mut() {
                if stored.is_face(f) {
                    stored.set_pc(c2d.clone(), *f2d, *l2d, *tol);
                    break;
                }
            }
        });
        return c2d;
    }

    // OCCT L563-565: C2D = FC2D_make2d(E, F, EF, f2d, l2d, tol, trim3d);
    // FC2D_AddNewCurveOnSurface(C2D, E, F, f2d, l2d, tol); return C2D;
    c2d = fc2d_make2d_ef(brep, e, f, ef, f2d, l2d, tol, trim3d);
    fc2d_add_new_curve_on_surface(&c2d, e, f, *f2d, *l2d, *tol);
    c2d
}
