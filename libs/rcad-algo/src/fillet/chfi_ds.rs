//! ChFiDS data structures — OCCT TKFillet/ChFiDS 1:1 translation.
//!
//! Sources:
//!   - ChFiDS_State.hxx, ChFiDS_ErrorStatus.hxx, ChFiDS_ChamfMethod.hxx,
//!     ChFiDS_ChamfMode.hxx, ChFiDS_TypeOfConcavity.hxx (enums)
//!   - ChFiDS_Spine.cxx/.hxx/.lxx (Spine)
//!   - ChFiDS_FilSpine.cxx/.hxx (FilSpine)
//!   - ChFiDS_ChamfSpine.cxx/.hxx (ChamfSpine)
//!   - ChFiDS_Stripe.hxx (Stripe)
//!   - ChFiDS_StripeMap.hxx (StripeMap)
//!   - ChFiDS_CommonPoint.hxx (CommonPoint)
//!
//! OCCT sequences are 1-based; the Vec translations keep the OCCT order and
//! the translated code adjusts the index arithmetic (OCCT `Value(i)` /
//! `Length()` becomes `self[i-1]` / `self.len()`).

use glam::{DVec2, DVec3};
use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::topo::topods::{Orientation, Shape, TEdgeData, TShape};
use rcad_kernel::topods;
use std::sync::Arc;

// =========================================================================
// OCCT ChFiDS_State.hxx — enum ChFiDS_State
// =========================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ChFiDS_State {
    OnSame,
    OnDiff,
    AllSame,
    BreakPoint,
    FreeBoundary,
    Closed,
    Tangent,
}

// =========================================================================
// OCCT ChFiDS_ErrorStatus.hxx — enum ChFiDS_ErrorStatus
// =========================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChFiDS_ErrorStatus {
    Ok,
    Error,
    WalkingFailure,
    StartsolFailure,
    TwistedSurface,
}

// =========================================================================
// OCCT ChFi3d_FilletShape.hxx — enum ChFi3d_FilletShape
// =========================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChFi3dFilletShape {
    Rational,
    QuasiAngular,
    Polynomial,
}

// =========================================================================
// OCCT ChFiDS_ChamfMethod.hxx — enum ChFiDS_ChamfMethod
// =========================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChFiDS_ChamfMethod {
    Sym,
    TwoDist,
    DistAngle,
}

// =========================================================================
// OCCT ChFiDS_ChamfMode.hxx — enum ChFiDS_ChamfMode
// =========================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChFiDS_ChamfMode {
    /// chamfer with constant distance from spine to one of the two surfaces
    ClassicChamfer,
    /// symmetric chamfer with constant throat
    ConstThroatChamfer,
    /// chamfer with constant throat with penetration
    ConstThroatWithPenetrationChamfer,
}

// =========================================================================
// OCCT ChFiDS_TypeOfConcavity.hxx — enum ChFiDS_TypeOfConcavity
// =========================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChFiDS_TypeOfConcavity {
    Concave,
    Convex,
    Tangential,
    FreeBound,
    Other,
    Mixed,
}

// =========================================================================
// OCCT Law_Function / Law_Composite — pending TKMath Law package.
// The ChFiDS_FilSpine law list and SetRadius(Law) overloads reference the
// law objects opaquely until the Law package is translated.
// =========================================================================

#[derive(Debug, Clone)]
pub struct LawFunction;

// =========================================================================
// OCCT ChFiDS_ElSpine — elementary spine (pending full field translation,
// ChFiDS_ElSpine.hxx).  Referenced by ChFiDS_Spine::elspines.
// =========================================================================

#[derive(Debug, Clone)]
pub struct ChFiDSElSpine;

// =========================================================================
// OCCT ChFiDS_SurfData + NCollection_HSequence — pending ChFiDS_SurfData
// translation.  The Stripe holds the sequence; the skeleton keeps it empty
// until PerformSetOfSurf (the numerical core) is translated.
// =========================================================================

#[derive(Debug, Clone)]
pub struct ChFiDSSurfData;

/// OCCT NCollection_HArray1<ChFiDS_CircSection> — pending ChFiDS_CircSection
/// translation; ChFi3d_FilBuilder::Sect returns this opaquely.
pub type ChFiDSCircSectionArray = Vec<ChFiDSCircSection>;

#[derive(Debug, Clone)]
pub struct ChFiDSCircSection;

// =========================================================================
// OCCT ChFiDS_CommonPoint.hxx L139-151 — private fields.
// Accessor bodies pending full ChFiDS_CommonPoint.cxx translation.
// =========================================================================

#[derive(Debug, Clone)]
pub struct ChFiDS_CommonPoint {
    pub arc: Shape,
    pub vtx: Shape,
    pub point: DVec3,
    pub vector: DVec3,
    pub tol: f64,
    pub prmarc: f64,
    pub prmtg: f64,
    pub traarc: Orientation,
    pub isonarc: bool,
    pub isvtx: bool,
    pub hasvector: bool,
}

impl Default for ChFiDS_CommonPoint {
    fn default() -> Self {
        ChFiDS_CommonPoint {
            arc: Shape::null(),
            vtx: Shape::null(),
            point: DVec3::ZERO,
            vector: DVec3::ZERO,
            tol: 0.0,
            prmarc: 0.0,
            prmtg: 0.0,
            traarc: Orientation::Forward,
            isonarc: false,
            isvtx: false,
            hasvector: false,
        }
    }
}

impl ChFiDS_CommonPoint {
    /// OCCT ChFiDS_CommonPoint::Point() — the 3D point of the common point.
    pub fn point(&self) -> DVec3 {
        self.point
    }

    /// OCCT ChFiDS_CommonPoint::IsEqual(other, Tol) — tolerance compare.
    pub fn is_equal(&self, other: &ChFiDS_CommonPoint, tol: f64) -> bool {
        self.point.distance(other.point) <= tol
    }
}

// =========================================================================
// OCCT ChFiDS_Spine — ChFiDS_Spine.hxx L?  fields (private:) +
// ChFiDS_Spine.cxx L37-67 (constructors), L171-184 (Reset),
// L186-264 (parameters), L266-302 (periodicity/vertices), L482-492
// (Index(E)), L503-540 (Load), L919-932 (error status),
// ChFiDS_Spine.lxx L124-127 (SetEdges) and the status inlines.
//
// OCCT inheritance (ChFiDS_FilSpine / ChFiDS_ChamfSpine derive from this
// class) is modeled by composition: the derived structs embed `ChFiDSSpine`.
// =========================================================================

#[derive(Debug, Clone)]
pub struct ChFiDSSpine {
    // OCCT: BRepAdaptor_Curve myCurve / myOffsetCurve — pending BRepAdaptor
    // translation; the current edge is tracked by `indexofcurve` alone.
    /// OCCT: int indexofcurve
    pub indexofcurve: i32,
    /// OCCT: ChFiDS_TypeOfConcavity myTypeOfConcavity
    pub my_type_of_concavity: ChFiDS_TypeOfConcavity,
    /// OCCT: ChFiDS_State firstState
    pub first_state: ChFiDS_State,
    /// OCCT: ChFiDS_State lastState
    pub last_state: ChFiDS_State,
    /// OCCT: NCollection_Sequence<TopoDS_Shape> spine
    pub spine: Vec<Shape>,
    /// OCCT: NCollection_Sequence<TopoDS_Shape> offsetspine
    pub offsetspine: Vec<Shape>,
    /// OCCT: occ::handle<NCollection_HArray1<double>> abscissa (1-based
    /// array of length = NbEdges; index i-1 here).
    pub abscissa: Option<Vec<f64>>,
    /// OCCT: occ::handle<NCollection_HArray1<double>> offset_abscissa
    pub offset_abscissa: Option<Vec<f64>>,
    /// OCCT: double tolesp
    pub tolesp: f64,
    /// OCCT: double firstparam
    pub firstparam: f64,
    /// OCCT: double lastparam
    pub lastparam: f64,
    /// OCCT: bool firstprolon
    pub firstprolon: bool,
    /// OCCT: bool lastprolon
    pub lastprolon: bool,
    /// OCCT: bool firstistgt
    pub firstistgt: bool,
    /// OCCT: bool lastistgt
    pub lastistgt: bool,
    /// OCCT: double firsttgtpar
    pub firsttgtpar: f64,
    /// OCCT: double lasttgtpar
    pub lasttgtpar: f64,
    /// OCCT: bool hasfirsttgt
    pub hasfirsttgt: bool,
    /// OCCT: bool haslasttgt
    pub haslasttgt: bool,
    /// OCCT: gp_Pnt firstori
    pub firstori: DVec3,
    /// OCCT: gp_Pnt lastori
    pub lastori: DVec3,
    /// OCCT: gp_Vec firsttgt
    pub firsttgt: DVec3,
    /// OCCT: gp_Vec lasttgt
    pub lasttgt: DVec3,
    /// OCCT: double valref
    pub valref: f64,
    /// OCCT: bool hasref
    pub hasref: bool,
    /// OCCT: ChFiDS_ErrorStatus errorstate
    pub errorstate: ChFiDS_ErrorStatus,
    /// OCCT: bool splitdone (ChFiDS_Spine.hxx public section)
    pub splitdone: bool,
    /// OCCT: ChFiDS_ChamfMode myMode
    pub my_mode: ChFiDS_ChamfMode,
    /// OCCT: NCollection_List<occ::handle<ChFiDS_ElSpine>> elspines
    pub elspines: Vec<ChFiDSElSpine>,
    /// OCCT: NCollection_List<occ::handle<ChFiDS_ElSpine>> offset_elspines
    pub offset_elspines: Vec<ChFiDSElSpine>,
}

impl Default for ChFiDSSpine {
    fn default() -> Self {
        Self::new()
    }
}

impl ChFiDSSpine {
    /// OCCT ChFiDS_Spine.cxx L37-61.
    pub fn new() -> Self {
        ChFiDSSpine {
            indexofcurve: 0,
            my_type_of_concavity: ChFiDS_TypeOfConcavity::Other,
            first_state: ChFiDS_State::OnSame,
            last_state: ChFiDS_State::OnSame,
            spine: Vec::new(),
            offsetspine: Vec::new(),
            abscissa: None,
            offset_abscissa: None,
            tolesp: CONFUSION, // Precision::Confusion()
            firstparam: 0.0,
            lastparam: 0.0,
            firstprolon: false,
            lastprolon: false,
            firstistgt: false,
            lastistgt: false,
            firsttgtpar: 0.0,
            lasttgtpar: 0.0,
            hasfirsttgt: false,
            haslasttgt: false,
            firstori: DVec3::ZERO,
            lastori: DVec3::ZERO,
            firsttgt: DVec3::ZERO,
            lasttgt: DVec3::ZERO,
            valref: 0.0,
            hasref: false,
            errorstate: ChFiDS_ErrorStatus::Ok,
            splitdone: false,
            my_mode: ChFiDS_ChamfMode::ClassicChamfer,
            elspines: Vec::new(),
            offset_elspines: Vec::new(),
        }
    }

    /// OCCT ChFiDS_Spine.cxx L63-67.
    pub fn with_tol(tol: f64) -> Self {
        let mut sp = Self::new();
        sp.tolesp = tol;
        sp
    }

    /// OCCT ChFiDS_Spine.lxx L124-127.
    pub fn set_edges(&mut self, e: Shape) {
        self.spine.push(e);
    }

    /// OCCT ChFiDS_Spine.lxx — NbEdges().
    pub fn nb_edges(&self) -> usize {
        self.spine.len()
    }

    /// OCCT ChFiDS_Spine.lxx — Edges(I) (1-based).
    pub fn edges(&self, i: usize) -> &Shape {
        &self.spine[i - 1]
    }

    /// OCCT ChFiDS_Spine.cxx L171-184.
    pub fn reset(&mut self, all_data: bool) {
        self.splitdone = false;
        self.elspines.clear();
        if all_data {
            self.firstparam = 0.0;
            self.lastparam = self.abscissa.as_ref().map_or(0.0, |a| a[a.len() - 1]);
            self.firstprolon = false;
            self.lastprolon = false;
        }
    }

    /// OCCT ChFiDS_Spine.cxx L186-195.
    pub fn first_parameter(&self) -> f64 {
        if self.firstprolon {
            return self.firstparam;
        }
        0.0
    }

    /// OCCT ChFiDS_Spine.cxx L197-206.
    pub fn last_parameter(&self) -> f64 {
        if self.lastprolon {
            return self.lastparam;
        }
        self.abscissa.as_ref().map_or(0.0, |a| a[a.len() - 1])
    }

    /// OCCT ChFiDS_Spine.cxx L208-220.
    pub fn set_first_parameter(&mut self, par: f64) {
        self.firstprolon = true;
        self.firstparam = par;
    }

    /// OCCT ChFiDS_Spine.cxx L222-235.
    pub fn set_last_parameter(&mut self, par: f64) {
        self.lastprolon = true;
        self.lastparam = par;
    }

    /// OCCT ChFiDS_Spine.cxx L237-246.
    pub fn first_parameter_of(&self, index_spine: usize) -> f64 {
        if index_spine == 1 {
            return 0.0;
        }
        self.abscissa.as_ref().map_or(0.0, |a| a[index_spine - 2])
    }

    /// OCCT ChFiDS_Spine.cxx L248-253.
    pub fn last_parameter_of(&self, index_spine: usize) -> f64 {
        self.abscissa.as_ref().map_or(0.0, |a| a[index_spine - 1])
    }

    /// OCCT ChFiDS_Spine.cxx L255-264.
    pub fn length_of(&self, index_spine: usize) -> f64 {
        if index_spine == 1 {
            return self.abscissa.as_ref().map_or(0.0, |a| a[index_spine - 1]);
        }
        self.abscissa.as_ref().map_or(0.0, |a| {
            a[index_spine - 1] - a[index_spine - 2]
        })
    }

    /// OCCT ChFiDS_Spine.cxx L266-269.
    pub fn is_periodic(&self) -> bool {
        self.first_state == ChFiDS_State::Closed
    }

    /// OCCT ChFiDS_Spine.cxx L271-274 — IsClosed via FirstVertex/LastVertex.
    pub fn is_closed(&self) -> bool {
        self.first_vertex().is_same(&self.last_vertex())
    }

    /// OCCT ChFiDS_Spine.cxx L280-290.
    pub fn first_vertex(&self) -> Shape {
        let e = &self.spine[0];
        if e.orientation == Orientation::Forward {
            edge_first_vertex(e)
        } else {
            edge_last_vertex(e)
        }
    }

    /// OCCT ChFiDS_Spine.cxx L292-302.
    pub fn last_vertex(&self) -> Shape {
        let e = &self.spine[self.spine.len() - 1];
        if e.orientation == Orientation::Forward {
            edge_last_vertex(e)
        } else {
            edge_first_vertex(e)
        }
    }

    /// OCCT ChFiDS_Spine.cxx L482-492.
    pub fn index_of_edge(&self, e: &Shape) -> usize {
        for (ie, s) in self.spine.iter().enumerate() {
            if e.is_same(s) {
                return ie + 1;
            }
        }
        0
    }

    /// OCCT ChFiDS_Spine.cxx L503-540.
    pub fn load(&mut self) {
        let len = self.spine.len();
        let mut abscissa = vec![0.0; len];
        let mut a1 = 0.0;
        for (i, s) in self.spine.iter().enumerate() {
            // OCCT: myCurve.Initialize(TopoDS::Edge(spine.Value(i)));
            // a1 += GCPnts_AbscissaPoint::Length(myCurve);
            let ed: &TEdgeData = s.as_edge().expect("ChFiDS_Spine::Load: not an edge");
            if let Some(c) = &ed.curve {
                a1 += rcad_kernel::base::gcpnts::abscissa_point::arc_length(
                    c, ed.range[0], ed.range[1],
                );
            }
            abscissa[i] = a1;
        }
        self.abscissa = Some(abscissa);
        self.indexofcurve = 1;

        // Here, we should update tolesp according to curve parameter range
        // if tolesp candidate less than default initial value.
        let umin = self.first_parameter();
        let umax = self.last_parameter();

        let new_tolesp = 5.0e-5 * (umax - umin);
        if self.tolesp > new_tolesp {
            self.tolesp = new_tolesp;
        }
    }

    /// OCCT ChFiDS_Spine.lxx — SetFirstStatus.
    pub fn set_first_status(&mut self, s: ChFiDS_State) {
        self.first_state = s;
    }

    /// OCCT ChFiDS_Spine.lxx — SetLastStatus.
    pub fn set_last_status(&mut self, s: ChFiDS_State) {
        self.last_state = s;
    }

    /// OCCT ChFiDS_Spine.lxx — FirstStatus.
    pub fn first_status(&self) -> ChFiDS_State {
        self.first_state
    }

    /// OCCT ChFiDS_Spine.lxx — LastStatus.
    pub fn last_status(&self) -> ChFiDS_State {
        self.last_state
    }

    /// OCCT ChFiDS_Spine.lxx — SetStatus(S, IsFirst).
    pub fn set_status(&mut self, s: ChFiDS_State, is_first: bool) {
        if is_first {
            self.first_state = s;
        } else {
            self.last_state = s;
        }
    }

    /// OCCT ChFiDS_Spine.lxx — Status(IsFirst).
    pub fn status(&self, is_first: bool) -> ChFiDS_State {
        if is_first {
            self.first_state
        } else {
            self.last_state
        }
    }

    /// OCCT ChFiDS_Spine.cxx L919-922.
    pub fn set_error_status(&mut self, state: ChFiDS_ErrorStatus) {
        self.errorstate = state;
    }

    /// OCCT ChFiDS_Spine.cxx L924-932.
    pub fn error_status(&self) -> ChFiDS_ErrorStatus {
        self.errorstate
    }

    /// OCCT ChFiDS_Spine.cxx Absc(const TopoDS_Vertex&) — the abscissa of a
    /// vertex on the composite spine.  The body runs through Prepare() and
    /// the BRepAdaptor_Curve composite traversal (pending translation); the
    /// call sites (SetRadius at a vertex) are marked pending boundaries.
    pub fn absc_of_vertex(&self, _v: &Shape) -> f64 {
        // OCCT: double npar = Absc(V); — pending BRepAdaptor_Curve.
        0.0
    }
}

/// OCCT TopExp::FirstVertex(E) — the edge TShape stores its start vertex.
fn edge_first_vertex(e: &Shape) -> Shape {
    let ed = e.as_edge().expect("not an edge");
    ed.first.clone()
}

/// OCCT TopExp::LastVertex(E).
fn edge_last_vertex(e: &Shape) -> Shape {
    let ed = e.as_edge().expect("not an edge");
    ed.last.clone()
}

// =========================================================================
// OCCT ChFiDS_FilSpine — ChFiDS_FilSpine.cxx L41-45 (ctors), L50-83 (Reset),
// L85-94 / L96-121 / L123-130 / L132-141 / L143-218 (SetRadius overloads),
// L246-266 / L268-307 (IsConstant), L309-315 / L317-356 / L358-366 (Radius).
// Fields: ChFiDS_FilSpine.hxx (parandrad, laws).
// =========================================================================

#[derive(Debug, Clone)]
pub struct ChFiDSFilSpine {
    /// OCCT: ChFiDS_Spine base class subobject.
    pub base: ChFiDSSpine,
    /// OCCT: NCollection_Sequence<gp_XY> parandrad
    pub parandrad: Vec<DVec2>,
    /// OCCT: NCollection_List<occ::handle<Law_Function>> laws (Law package
    /// pending; stored opaquely).
    pub laws: Vec<LawFunction>,
}

impl Default for ChFiDSFilSpine {
    fn default() -> Self {
        Self::new()
    }
}

impl ChFiDSFilSpine {
    /// OCCT ChFiDS_FilSpine.cxx L41-43.
    pub fn new() -> Self {
        ChFiDSFilSpine {
            base: ChFiDSSpine::new(),
            parandrad: Vec::new(),
            laws: Vec::new(),
        }
    }

    /// OCCT ChFiDS_FilSpine.cxx L43-45.
    pub fn with_tol(tol: f64) -> Self {
        ChFiDSFilSpine {
            base: ChFiDSSpine::with_tol(tol),
            parandrad: Vec::new(),
            laws: Vec::new(),
        }
    }

    /// OCCT ChFiDS_FilSpine.cxx L50-83.
    pub fn reset(&mut self, all_data: bool) {
        self.base.reset(all_data);
        self.laws.clear();
        if all_data {
            self.parandrad.clear();
        } else {
            // Complete parandrad
            let spinedeb = self.base.first_parameter();
            let spinefin = self.base.last_parameter();

            if let (Some(first), Some(last)) = (self.parandrad.first(), self.parandrad.last()) {
                let mut first_uandr = *first;
                let mut last_uandr = *last;
                let mut prepend = None;
                let mut append = None;
                if (spinedeb - first_uandr.x).abs() > f64::MIN_POSITIVE {
                    first_uandr.x = spinedeb;
                    prepend = Some(first_uandr);
                }
                if (spinefin - last_uandr.x).abs() > f64::MIN_POSITIVE {
                    last_uandr.x = spinefin;
                    append = Some(last_uandr);
                }
                if let Some(v) = prepend {
                    self.parandrad.insert(0, v);
                }
                if let Some(v) = append {
                    self.parandrad.push(v);
                }
            }

            if self.base.is_periodic() {
                let n = self.parandrad.len();
                if n >= 1 {
                    let y = self.parandrad[0].y;
                    self.parandrad[n - 1].y = y;
                }
            }
        }
    }

    /// OCCT ChFiDS_FilSpine.cxx L85-94.
    pub fn set_radius_on_edge(&mut self, radius: f64, e: &Shape) {
        self.base.splitdone = false;
        let ie = self.base.index_of_edge(e);
        let first_uandr = DVec2::new(0.0, radius);
        let last_uandr = DVec2::new(1.0, radius);
        self.set_radius_uandr(first_uandr, ie);
        self.set_radius_uandr(last_uandr, ie);
    }

    /// OCCT ChFiDS_FilSpine.cxx L132-141.
    pub fn set_radius(&mut self, radius: f64) {
        self.parandrad.clear();
        let first_uandr = DVec2::new(self.base.first_parameter(), radius);
        let last_uandr = DVec2::new(self.base.last_parameter(), radius);
        self.set_radius_uandr(first_uandr, 0);
        self.set_radius_uandr(last_uandr, 0);
    }

    /// OCCT ChFiDS_FilSpine.cxx L143-218 (the splitdone law-replay tail
    /// depends on Law_Composite/ChFiDS_ElSpine internals — pending).
    pub fn set_radius_uandr(&mut self, uandr: DVec2, iinc: usize) {
        let w;
        if iinc == 0 {
            w = uandr.x;
        } else {
            let uf = self.base.first_parameter_of(iinc);
            let ul = self.base.last_parameter_of(iinc);
            w = uf + uandr.x * (ul - uf);
        }

        let pr = DVec2::new(w, uandr.y);
        let mut i = 0usize; // OCCT 1-based loop variable i
        while i < self.parandrad.len() {
            let xi = self.parandrad[i].x;
            if xi == w {
                self.parandrad[i].y = uandr.y;
                if !self.base.splitdone {
                    return;
                } else {
                    i += 1;
                    break;
                }
            } else if xi > w {
                self.parandrad.insert(i, pr);
                if !self.base.splitdone {
                    return;
                } else {
                    i += 1;
                    break;
                }
            }
            i += 1;
        }
        if i == self.parandrad.len() {
            self.parandrad.push(pr);
        }
        // si le split est done il faut rejouer la law correspondant au
        // parametre W — pending Law_Composite translation.
    }

    /// OCCT ChFiDS_FilSpine.cxx L246-266.
    pub fn is_constant(&self) -> bool {
        if self.parandrad.is_empty() {
            return false;
        }
        let radius = self.parandrad[0].y;
        for i in 1..self.parandrad.len() {
            if (radius - self.parandrad[i].y).abs() > CONFUSION {
                return false;
            }
        }
        true
    }

    /// OCCT ChFiDS_FilSpine.cxx L268-307 (1-based IE).
    pub fn is_constant_on(&self, ie: usize) -> bool {
        let uf = self.base.first_parameter_of(ie);
        let ul = self.base.last_parameter_of(ie);

        let mut start_rad = 0.0f64;
        let mut i = 0usize;
        while i + 1 < self.parandrad.len() {
            let par = self.parandrad[i].x;
            let rad = self.parandrad[i].y;
            let nextpar = self.parandrad[i + 1].x;
            if (uf - par).abs() <= f64::MIN_POSITIVE
                || (par < uf && uf < nextpar && nextpar - uf > f64::MIN_POSITIVE)
            {
                start_rad = rad;
                break;
            }
            i += 1;
        }
        let mut i = i + 1;
        while i < self.parandrad.len() {
            let par = self.parandrad[i].x;
            let rad = self.parandrad[i].y;
            if (rad - start_rad).abs() > CONFUSION {
                return false;
            }
            if (ul - par).abs() <= f64::MIN_POSITIVE {
                return true;
            }
            if par > ul {
                return true;
            }
            i += 1;
        }
        true
    }

    /// OCCT ChFiDS_FilSpine.cxx L317-356 (1-based IE).
    pub fn radius_on(&self, ie: usize) -> f64 {
        let uf = self.base.first_parameter_of(ie);
        let ul = self.base.last_parameter_of(ie);

        let mut start_rad = 0.0f64;
        let mut i = 0usize;
        while i + 1 < self.parandrad.len() {
            let par = self.parandrad[i].x;
            let rad = self.parandrad[i].y;
            let nextpar = self.parandrad[i + 1].x;
            if (uf - par).abs() <= f64::MIN_POSITIVE
                || (par < uf && uf < nextpar && nextpar - uf > f64::MIN_POSITIVE)
            {
                start_rad = rad;
                break;
            }
            i += 1;
        }
        let mut i = i + 1;
        while i < self.parandrad.len() {
            let par = self.parandrad[i].x;
            let rad = self.parandrad[i].y;
            if (rad - start_rad).abs() > CONFUSION {
                panic!("Standard_DomainError: Edge is not constant");
            }
            if (ul - par).abs() <= f64::MIN_POSITIVE {
                return start_rad;
            }
            if par > ul {
                return start_rad;
            }
            i += 1;
        }
        start_rad
    }

    /// OCCT ChFiDS_FilSpine.cxx L358-366.
    pub fn radius(&self) -> f64 {
        if !self.is_constant() {
            panic!("Standard_DomainError: Spine is not constant");
        }
        self.parandrad[0].y
    }

    /// OCCT ChFiDS_FilSpine.cxx L309-315.
    pub fn radius_on_edge(&self, e: &Shape) -> f64 {
        let ie = self.base.index_of_edge(e);
        self.radius_on(ie)
    }

    /// OCCT ChFiDS_FilSpine.cxx L123-130.
    ///
    /// The `Absc(V)` spine query (abscissa of a vertex on the composite
    /// spine curve) depends on the BRepAdaptor_Curve machinery
    /// (ChFiDS_Spine.cxx Prepare/Absc) that is pending translation, so the
    /// parameter resolution point is a marked boundary.
    pub fn set_radius_at_vertex(&mut self, radius: f64, v: &Shape) {
        let npar = self.base.absc_of_vertex(v);
        let uandr = DVec2::new(npar, radius);
        self.set_radius_uandr(uandr, 0);
    }
}

// =========================================================================
// OCCT ChFiDS_ChamfSpine — ChFiDS_ChamfSpine.cxx L24-45 (ctors) and
// L47-126 (dist/angle accessors).  Fields: ChFiDS_ChamfSpine.hxx L58-64.
// =========================================================================

#[derive(Debug, Clone)]
pub struct ChFiDSChamfSpine {
    /// OCCT: ChFiDS_Spine base class subobject.
    pub base: ChFiDSSpine,
    /// OCCT: double d1
    pub d1: f64,
    /// OCCT: double d2
    pub d2: f64,
    /// OCCT: double angle
    pub angle: f64,
    /// OCCT: ChFiDS_ChamfMethod mChamf
    pub mchamf: ChFiDS_ChamfMethod,
}

impl Default for ChFiDSChamfSpine {
    fn default() -> Self {
        Self::new()
    }
}

impl ChFiDSChamfSpine {
    /// OCCT ChFiDS_ChamfSpine.cxx L24-36.
    pub fn new() -> Self {
        let mut sp = ChFiDSChamfSpine {
            base: ChFiDSSpine::new(),
            d1: 0.0,
            d2: 0.0,
            angle: 0.0,
            mchamf: ChFiDS_ChamfMethod::Sym,
        };
        sp.base.my_mode = ChFiDS_ChamfMode::ClassicChamfer;
        sp
    }

    /// OCCT ChFiDS_ChamfSpine.cxx L30-36.
    pub fn with_tol(tol: f64) -> Self {
        let mut sp = ChFiDSChamfSpine {
            base: ChFiDSSpine::with_tol(tol),
            d1: 0.0,
            d2: 0.0,
            angle: 0.0,
            mchamf: ChFiDS_ChamfMethod::Sym,
        };
        sp.base.my_mode = ChFiDS_ChamfMode::ClassicChamfer;
        sp
    }

    /// OCCT ChFiDS_ChamfSpine.cxx L38-45.
    pub fn get_dist(&self) -> f64 {
        if self.mchamf != ChFiDS_ChamfMethod::Sym {
            panic!("Standard_Failure: Chamfer is not symmetric");
        }
        self.d1
    }

    /// OCCT ChFiDS_ChamfSpine.cxx L47-54.
    pub fn set_dist(&mut self, dis: f64) {
        self.mchamf = ChFiDS_ChamfMethod::Sym;
        self.d1 = dis;
    }

    /// OCCT ChFiDS_ChamfSpine.cxx L56-67.
    pub fn dists(&self) -> (f64, f64) {
        if self.mchamf != ChFiDS_ChamfMethod::TwoDist {
            panic!("Standard_Failure: Chamfer is not a Two Dists Chamfer");
        }
        (self.d1, self.d2)
    }

    /// OCCT ChFiDS_ChamfSpine.cxx L69-77.
    pub fn set_dists(&mut self, dis1: f64, dis2: f64) {
        self.mchamf = ChFiDS_ChamfMethod::TwoDist;
        self.d1 = dis1;
        self.d2 = dis2;
    }

    /// OCCT ChFiDS_ChamfSpine.cxx L79-95.
    pub fn get_dist_angle(&self) -> (f64, f64) {
        if self.mchamf != ChFiDS_ChamfMethod::DistAngle {
            panic!("Standard_Failure: Chamfer is not a Two Dists Chamfer");
        }
        (self.d1, self.angle)
    }

    /// OCCT ChFiDS_ChamfSpine.cxx L97-107.
    pub fn set_dist_angle(&mut self, dis: f64, angle: f64) {
        self.mchamf = ChFiDS_ChamfMethod::DistAngle;
        self.d1 = dis;
        self.angle = angle;
    }

    /// OCCT ChFiDS_ChamfSpine.cxx (SetMode — sets the base myMode).
    pub fn set_mode(&mut self, the_mode: ChFiDS_ChamfMode) {
        self.base.my_mode = the_mode;
    }

    /// OCCT ChFiDS_ChamfSpine.cxx L116-126.
    pub fn is_chamfer(&self) -> ChFiDS_ChamfMethod {
        self.mchamf
    }
}

// =========================================================================
// OCCT ChFiDS_Stripe — ChFiDS_Stripe.hxx L146-168 (fields).  Accessor
// bodies live in ChFiDS_Stripe.lxx (inline); only the accessors used by the
// translated chain are provided here.
// =========================================================================

pub type SharedSurfData = Arc<ChFiDSSurfData>;

#[derive(Debug, Clone)]
pub struct ChFiDSStripe {
    /// OCCT: double pardeb1
    pub pardeb1: f64,
    /// OCCT: double parfin1
    pub parfin1: f64,
    /// OCCT: double pardeb2
    pub pardeb2: f64,
    /// OCCT: double parfin2
    pub parfin2: f64,
    /// OCCT: occ::handle<ChFiDS_Spine> mySpine — polymorphic Fil/Chamf spine.
    pub my_spine: Option<ChFiDSSpineHandle>,
    /// OCCT: occ::handle<NCollection_HSequence<occ::handle<ChFiDS_SurfData>>>
    /// myHdata
    pub my_hdata: Vec<SharedSurfData>,
    /// OCCT: int myChoix
    pub my_choix: i32,
    /// OCCT: int indexOfSolid
    pub index_of_solid: i32,
    /// OCCT: int indexOfcurve1
    pub index_ofcurve1: i32,
    /// OCCT: int indexOfcurve2
    pub index_ofcurve2: i32,
    /// OCCT: int indexfirstPOnS1
    pub indexfirst_pon_s1: i32,
    /// OCCT: int indexlastPOnS1
    pub indexlast_pon_s1: i32,
    /// OCCT: int indexfirstPOnS2
    pub indexfirst_pon_s2: i32,
    /// OCCT: int indexlastPOnS2
    pub indexlast_pon_s2: i32,
    /// OCCT: int begfilled
    pub begfilled: i32,
    /// OCCT: int endfilled
    pub endfilled: i32,
    /// OCCT: TopAbs_Orientation myOr1
    pub my_or1: Orientation,
    /// OCCT: TopAbs_Orientation myOr2
    pub my_or2: Orientation,
    /// OCCT: TopAbs_Orientation orcurv1
    pub orcurv1: Orientation,
    /// OCCT: TopAbs_Orientation orcurv2
    pub orcurv2: Orientation,
}

impl Default for ChFiDSStripe {
    fn default() -> Self {
        ChFiDSStripe {
            pardeb1: 0.0,
            parfin1: 0.0,
            pardeb2: 0.0,
            parfin2: 0.0,
            my_spine: None,
            my_hdata: Vec::new(),
            my_choix: 0,
            index_of_solid: 0,
            index_ofcurve1: 0,
            index_ofcurve2: 0,
            indexfirst_pon_s1: 0,
            indexlast_pon_s1: 0,
            indexfirst_pon_s2: 0,
            indexlast_pon_s2: 0,
            begfilled: 0,
            endfilled: 0,
            my_or1: Orientation::Forward,
            my_or2: Orientation::Forward,
            orcurv1: Orientation::Forward,
            orcurv2: Orientation::Forward,
        }
    }
}

impl ChFiDSStripe {
    /// OCCT ChFiDS_Stripe.lxx — Spine().
    pub fn spine(&self) -> Option<&ChFiDSSpineHandle> {
        self.my_spine.as_ref()
    }

    /// OCCT ChFiDS_Stripe.lxx — ChangeSpine() (assignable handle slot).
    pub fn change_spine(&mut self, sp: ChFiDSSpineHandle) {
        self.my_spine = Some(sp);
    }

    /// OCCT ChFiDS_Stripe.lxx — Reset().
    pub fn reset(&mut self) {
        if let Some(sp) = &mut self.my_spine {
            sp.base_mut().reset(false);
        }
    }
}

/// OCCT handle<ChFiDS_Spine> polymorphic slot: down_cast chooses the
/// concrete derived spine, matching ChFi3d_FilBuilder / ChFi3d_ChBuilder.
#[derive(Debug, Clone)]
pub enum ChFiDSSpineHandle {
    Fil(ChFiDSFilSpine),
    Chamf(ChFiDSChamfSpine),
}

impl ChFiDSSpineHandle {
    /// Common base access (OCCT: the handle is used as ChFiDS_Spine&).
    pub fn base(&self) -> &ChFiDSSpine {
        match self {
            ChFiDSSpineHandle::Fil(sp) => &sp.base,
            ChFiDSSpineHandle::Chamf(sp) => &sp.base,
        }
    }

    pub fn base_mut(&mut self) -> &mut ChFiDSSpine {
        match self {
            ChFiDSSpineHandle::Fil(sp) => &mut sp.base,
            ChFiDSSpineHandle::Chamf(sp) => &mut sp.base,
        }
    }

    /// OCCT occ::down_cast<ChFiDS_FilSpine>.
    pub fn down_cast_fil(&self) -> Option<&ChFiDSFilSpine> {
        match self {
            ChFiDSSpineHandle::Fil(sp) => Some(sp),
            _ => None,
        }
    }

    /// OCCT occ::down_cast<ChFiDS_ChamfSpine>.
    pub fn down_cast_chamf(&self) -> Option<&ChFiDSChamfSpine> {
        match self {
            ChFiDSSpineHandle::Chamf(sp) => Some(sp),
            _ => None,
        }
    }

    pub fn down_cast_fil_mut(&mut self) -> Option<&mut ChFiDSFilSpine> {
        match self {
            ChFiDSSpineHandle::Fil(sp) => Some(sp),
            _ => None,
        }
    }

    pub fn down_cast_chamf_mut(&mut self) -> Option<&mut ChFiDSChamfSpine> {
        match self {
            ChFiDSSpineHandle::Chamf(sp) => Some(sp),
            _ => None,
        }
    }
}

pub type SharedStripe = Arc<std::sync::RwLock<ChFiDSStripe>>;

// =========================================================================
// OCCT ChFiDS_StripeMap — ChFiDS_StripeMap.hxx L65-69
// (NCollection_IndexedDataMap<TopoDS_Vertex, List<Stripe>>).
// Keyed by TShape pointer identity (IsSame semantics).
// =========================================================================

#[derive(Debug, Clone, Default)]
pub struct ChFiDSStripeMap {
    /// Insertion-ordered keys (FindKey(i) is 1-based in OCCT).
    my_keys: Vec<Shape>,
    /// FindFromIndex values keyed by TShape pointer.
    my_map: std::collections::HashMap<u64, Vec<SharedStripe>>,
}

impl ChFiDSStripeMap {
    pub fn new() -> Self {
        ChFiDSStripeMap {
            my_keys: Vec::new(),
            my_map: std::collections::HashMap::new(),
        }
    }

    /// OCCT ChFiDS_StripeMap.cxx — Add(V, Stripe).
    pub fn add(&mut self, v: &Shape, stripe: SharedStripe) {
        let key = v.ptr_id();
        if !self.my_map.contains_key(&key) {
            self.my_keys.push(v.clone());
            self.my_map.insert(key, Vec::new());
        }
        self.my_map.get_mut(&key).unwrap().push(stripe);
    }

    /// OCCT — Extent().
    pub fn extent(&self) -> usize {
        self.my_keys.len()
    }

    /// OCCT — FindKey(i) (1-based).
    pub fn find_key(&self, i: usize) -> &Shape {
        &self.my_keys[i - 1]
    }

    /// OCCT — FindFromIndex(i) (1-based).
    pub fn find_from_index(&self, i: usize) -> &Vec<SharedStripe> {
        self.my_map.get(&self.my_keys[i - 1].ptr_id()).expect("stripe map index")
    }

    /// OCCT — Clear().
    pub fn clear(&mut self) {
        self.my_keys.clear();
        self.my_map.clear();
    }
}

// =========================================================================
// OCCT ChFiDS_Map (TKFillet/ChFiDS/ChFiDS_Map.hxx) — ancestor map filled by
// TopExp::MapShapesAndAncestors(S, TOR, TOS).  Keyed by TShape pointer
// identity (IsSame semantics), insertion-ordered for indexed access.
// =========================================================================

#[derive(Debug, Clone, Default)]
pub struct ChFiDSMap {
    my_keys: Vec<Shape>,
    my_map: std::collections::HashMap<u64, Vec<Shape>>,
}

impl ChFiDSMap {
    pub fn new() -> Self {
        ChFiDSMap {
            my_keys: Vec::new(),
            my_map: std::collections::HashMap::new(),
        }
    }

    /// OCCT ChFiDS_Map::Fill(S, TOR, TOS) — TopExp::MapShapesAndAncestors:
    /// every shape of type `tor` (the key) is mapped to the shapes of type
    /// `tos` that contain it.
    pub fn fill(
        &mut self,
        brep: &topods::BRep,
        tor: topods::ShapeType,
        tos: topods::ShapeType,
    ) {
        self.my_keys.clear();
        self.my_map.clear();

        // Emit (key, ancestor) pairs per TopExp::MapShapesAndAncestors.
        let mut pairs: Vec<(Shape, Shape)> = Vec::new();
        for (fi, ts) in brep.tshapes.iter().enumerate() {
            let anc = Shape::from_parts(ts.clone(), fi, 0, Orientation::Forward);
            let child_shape =
                |idx: usize| -> Option<Shape> {
                    brep.tshapes.get(idx).map(|t| {
                        Shape::from_parts(t.clone(), idx, 0, Orientation::Forward)
                    })
                };
            match ts.as_ref() {
                TShape::Shell(sd) => {
                    if tos == topods::ShapeType::Face && tor == topods::ShapeType::Edge {
                        for fs in &sd.faces {
                            if let Some(fts) = brep.tshapes.get(fs.index) {
                                if let TShape::Face(fd) = fts.as_ref() {
                                    if let Some(wt) = brep.tshapes.get(fd.outer_wire.index) {
                                        if let TShape::Wire(wd) = wt.as_ref() {
                                            for we in &wd.edges {
                                                if let Some(es) = child_shape(we.index) {
                                                    pairs.push((es, anc.clone()));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                TShape::Solid(sd) => {
                    if tos == topods::ShapeType::Shell && tor == topods::ShapeType::Edge {
                        for shs in &sd.shells {
                            if let Some(shts) = brep.tshapes.get(shs.index) {
                                if let TShape::Shell(shd) = shts.as_ref() {
                                    for fs in &shd.faces {
                                        if let Some(fts) = brep.tshapes.get(fs.index) {
                                            if let TShape::Face(fd) = fts.as_ref() {
                                                if let Some(wt) =
                                                    brep.tshapes.get(fd.outer_wire.index)
                                                {
                                                    if let TShape::Wire(wd) = wt.as_ref() {
                                                        for we in &wd.edges {
                                                            if let Some(es) = child_shape(we.index)
                                                            {
                                                                pairs.push((es, anc.clone()));
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                TShape::Face(fd) => {
                    // Face ancestors: faces over edges (EDGE,FACE) and over
                    // vertices (VERTEX,FACE).
                    let want = matches!(
                        (tor, tos),
                        (topods::ShapeType::Edge, topods::ShapeType::Face)
                            | (topods::ShapeType::Vertex, topods::ShapeType::Face)
                    );
                    if !want {
                        continue;
                    }
                    if let Some(wt) = brep.tshapes.get(fd.outer_wire.index) {
                        if let TShape::Wire(wd) = wt.as_ref() {
                            for we in &wd.edges {
                                if let Some(es) = child_shape(we.index) {
                                    if tor == topods::ShapeType::Edge {
                                        pairs.push((es.clone(), anc.clone()));
                                    }
                                    if tor == topods::ShapeType::Vertex {
                                        if let Some(ed) = es.as_edge() {
                                            for v in [&ed.first, &ed.last] {
                                                pairs.push((v.clone(), anc.clone()));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                TShape::Edge(ed) => {
                    // Edge ancestors: edges over vertices (VERTEX,EDGE).
                    if tos == topods::ShapeType::Edge && tor == topods::ShapeType::Vertex {
                        for v in [&ed.first, &ed.last] {
                            pairs.push((v.clone(), anc.clone()));
                        }
                    }
                }
                _ => {}
            }
        }

        for (key, anc) in pairs {
            let k = key.ptr_id();
            if !self.my_map.contains_key(&k) {
                self.my_keys.push(key.clone());
                self.my_map.insert(k, Vec::new());
            }
            let list = self.my_map.get_mut(&k).unwrap();
            // MapShapesAndAncestors does not append the ancestor twice for
            // the same key entry.
            if !list.iter().any(|a| a.ptr_id() == anc.ptr_id()) {
                list.push(anc);
            }
        }
    }

    /// OCCT ChFiDS_Map::Contains(key).
    pub fn contains(&self, key: &Shape) -> bool {
        self.my_map.contains_key(&key.ptr_id())
    }

    /// OCCT ChFiDS_Map::operator()(key) — ancestor list of key.
    pub fn find(&self, key: &Shape) -> &Vec<Shape> {
        self.my_map.get(&key.ptr_id()).expect("ChFiDSMap: key not bound")
    }
}
