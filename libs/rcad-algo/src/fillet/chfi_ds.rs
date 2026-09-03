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
pub struct ChFiDSElSpine {
    /// OCCT ChFiDS_ElSpine.hxx: double firstparam / lastparam
    pub firstparam: f64,
    pub lastparam: f64,
    /// OCCT: gp_Pnt firstPnt / lastPnt, gp_Vec firstTgt / lastTgt
    pub firstpnt: DVec3,
    pub firsttgt: DVec3,
    pub lastpnt: DVec3,
    pub lasttgt: DVec3,
    /// OCCT: bool periodic
    pub periodic: bool,
    /// OCCT: ChangeNext() / ChangePrevious() — SurfData links
    pub next: Option<SharedSurfData>,
    pub previous: Option<SharedSurfData>,
}

// =========================================================================
// OCCT ChFiDS_FaceInterference.hxx L72-77 — fields; bodies from
// ChFiDS_FaceInterference.cxx / .lxx.
// =========================================================================

/// OCCT ChFiDS_FaceInterference — the interference data of the fillet
/// surface with one support face (3D tangency line + the two pcurves).
#[derive(Debug, Clone)]
pub struct ChFiDS_FaceInterference {
    /// OCCT: double firstParam
    pub first_param: f64,
    /// OCCT: double lastParam
    pub last_param: f64,
    /// OCCT: occ::handle<Geom2d_Curve> pCurveOnFace
    pub pcurve_on_face: Option<rcad_kernel::geom::Curve2d>,
    /// OCCT: occ::handle<Geom2d_Curve> pCurveOnSurf
    pub pcurve_on_surf: Option<rcad_kernel::geom::Curve2d>,
    /// OCCT: int lineindex (DS index of the 3D tangency line)
    pub lineindex: i32,
    /// OCCT: TopAbs_Orientation LineTransition
    pub line_transition: Orientation,
}

impl Default for ChFiDS_FaceInterference {
    fn default() -> Self {
        ChFiDS_FaceInterference {
            first_param: 0.0,
            last_param: 0.0,
            pcurve_on_face: None,
            pcurve_on_surf: None,
            lineindex: 0,
            line_transition: Orientation::Forward,
        }
    }
}

impl ChFiDS_FaceInterference {
    /// OCCT ChFiDS_FaceInterference.cxx SetInterference(Index, Or, PC1, PC2).
    pub fn set_interference(
        &mut self,
        index: i32,
        or: Orientation,
        pcurve_on_face: Option<rcad_kernel::geom::Curve2d>,
        pcurve_on_surf: Option<rcad_kernel::geom::Curve2d>,
    ) {
        self.lineindex = index;
        self.line_transition = or;
        self.pcurve_on_face = pcurve_on_face;
        self.pcurve_on_surf = pcurve_on_surf;
    }

    /// OCCT ChFiDS_FaceInterference.cxx Parameter(IsFirst).
    pub fn parameter(&self, is_first: bool) -> f64 {
        if is_first {
            self.first_param
        } else {
            self.last_param
        }
    }

    /// OCCT .lxx — PCurveOnFace().
    pub fn pcurve_on_face(&self) -> Option<&rcad_kernel::geom::Curve2d> {
        self.pcurve_on_face.as_ref()
    }

    /// OCCT .lxx — PCurveOnSurf().
    pub fn pcurve_on_surf(&self) -> Option<&rcad_kernel::geom::Curve2d> {
        self.pcurve_on_surf.as_ref()
    }

    /// OCCT .lxx — first/last parameter setter.
    pub fn set_parameter(&mut self, is_first: bool, par: f64) {
        if is_first {
            self.first_param = par;
        } else {
            self.last_param = par;
        }
    }

    /// OCCT .lxx — FirstParameter().
    pub fn parameter_first(&self) -> f64 {
        self.first_param
    }

    /// OCCT .lxx — LastParameter().
    pub fn parameter_last(&self) -> f64 {
        self.last_param
    }

    /// OCCT .lxx — LineIndex().
    pub fn line_index(&self) -> i32 {
        self.lineindex
    }

    /// OCCT .lxx — Transition() (the line transition).
    pub fn transition(&self) -> Orientation {
        self.line_transition
    }

    /// OCCT .lxx — ChangePCurveOnFace() (assignable slot).
    pub fn change_pcurve_on_face(&mut self) -> &mut Option<rcad_kernel::geom::Curve2d> {
        &mut self.pcurve_on_face
    }

    /// OCCT .lxx — ChangePCurveOnSurf() (assignable slot).
    pub fn change_pcurve_on_surf(&mut self) -> &mut Option<rcad_kernel::geom::Curve2d> {
        &mut self.pcurve_on_surf
    }
}

// =========================================================================
// OCCT ChFiDS_SurfData.hxx L154-179 — fields; accessors from
// ChFiDS_SurfData.cxx / .lxx.  The Simul handle (simul) is pending the
// ChFiDS_CircSection translation.
// =========================================================================

/// OCCT ChFiDS_SurfData — all the data of one fillet surface patch along
/// the spine.
#[derive(Debug, Clone)]
pub struct ChFiDSSurfData {
    /// OCCT: ChFiDS_CommonPoint pfirstOnS1
    pub pfirst_on_s1: ChFiDS_CommonPoint,
    /// OCCT: ChFiDS_CommonPoint plastOnS1
    pub plast_on_s1: ChFiDS_CommonPoint,
    /// OCCT: ChFiDS_CommonPoint pfirstOnS2
    pub pfirst_on_s2: ChFiDS_CommonPoint,
    /// OCCT: ChFiDS_CommonPoint plastOnS2
    pub plast_on_s2: ChFiDS_CommonPoint,
    /// OCCT: ChFiDS_FaceInterference intf1
    pub intf1: ChFiDS_FaceInterference,
    /// OCCT: ChFiDS_FaceInterference intf2
    pub intf2: ChFiDS_FaceInterference,
    /// OCCT: gp_Pnt2d p2df1
    pub p2df1: DVec2,
    /// OCCT: gp_Pnt2d p2dl1
    pub p2dl1: DVec2,
    /// OCCT: gp_Pnt2d p2df2
    pub p2df2: DVec2,
    /// OCCT: gp_Pnt2d p2dl2
    pub p2dl2: DVec2,
    /// OCCT: double ufspine
    pub ufspine: f64,
    /// OCCT: double ulspine
    pub ulspine: f64,
    /// OCCT: double myfirstextend
    pub myfirstextend: f64,
    /// OCCT: double mylastextend
    pub mylastextend: f64,
    /// OCCT: int indexOfS1
    pub index_of_s1: i32,
    /// OCCT: int indexOfC1
    pub index_of_c1: i32,
    /// OCCT: int indexOfS2
    pub index_of_s2: i32,
    /// OCCT: int indexOfC2
    pub index_of_c2: i32,
    /// OCCT: int indexOfConge
    pub index_of_conge: i32,
    /// OCCT: bool isoncurv1
    pub isoncurv1: bool,
    /// OCCT: bool isoncurv2
    pub isoncurv2: bool,
    /// OCCT: bool twistons1
    pub twistons1: bool,
    /// OCCT: bool twistons2
    pub twistons2: bool,
    /// OCCT: TopAbs_Orientation orientation
    pub orientation: Orientation,
    /// rcad extension: the fillet surface's DS index (OCCT stores the
    /// surface in the TopOpeBRepDS; the DS placeholder mirrors it).
    pub surf_index: i32,
}

impl Default for ChFiDSSurfData {
    fn default() -> Self {
        ChFiDSSurfData {
            pfirst_on_s1: ChFiDS_CommonPoint::default(),
            plast_on_s1: ChFiDS_CommonPoint::default(),
            pfirst_on_s2: ChFiDS_CommonPoint::default(),
            plast_on_s2: ChFiDS_CommonPoint::default(),
            intf1: ChFiDS_FaceInterference::default(),
            intf2: ChFiDS_FaceInterference::default(),
            p2df1: DVec2::ZERO,
            p2dl1: DVec2::ZERO,
            p2df2: DVec2::ZERO,
            p2dl2: DVec2::ZERO,
            ufspine: 0.0,
            ulspine: 0.0,
            myfirstextend: 0.0,
            mylastextend: 0.0,
            index_of_s1: 0,
            index_of_c1: 0,
            index_of_s2: 0,
            index_of_c2: 0,
            index_of_conge: 0,
            isoncurv1: false,
            isoncurv2: false,
            twistons1: false,
            twistons2: false,
            orientation: Orientation::Forward,
            surf_index: 0,
        }
    }
}

impl ChFiDSSurfData {
    /// OCCT ChFiDS_SurfData.cxx ChangeSurf(Index).
    pub fn change_surf(&mut self, index: i32) {
        self.surf_index = index;
    }

    /// OCCT .lxx — Surf().
    pub fn surf(&self) -> i32 {
        self.surf_index
    }

    /// OCCT .lxx — ChangeOrientation().
    pub fn change_orientation(&mut self) -> &mut Orientation {
        &mut self.orientation
    }

    /// OCCT .lxx — InterferenceOnS1().
    pub fn interference_on_s1(&self) -> &ChFiDS_FaceInterference {
        &self.intf1
    }

    /// OCCT .lxx — ChangeInterferenceOnS1().
    pub fn change_interference_on_s1(&mut self) -> &mut ChFiDS_FaceInterference {
        &mut self.intf1
    }

    /// OCCT .lxx — InterferenceOnS2().
    pub fn interference_on_s2(&self) -> &ChFiDS_FaceInterference {
        &self.intf2
    }

    /// OCCT .lxx — ChangeInterferenceOnS2().
    pub fn change_interference_on_s2(&mut self) -> &mut ChFiDS_FaceInterference {
        &mut self.intf2
    }

    /// OCCT .cxx FirstSpineParam() (ufspine).
    pub fn first_spine_param(&self) -> f64 {
        self.ufspine
    }

    /// OCCT .cxx LastSpineParam() (ulspine).
    pub fn last_spine_param(&self) -> f64 {
        self.ulspine
    }

    /// OCCT .cxx FirstSpineParam(Par).
    pub fn set_first_spine_param(&mut self, par: f64) {
        self.ufspine = par;
    }

    /// OCCT .cxx LastSpineParam(Par).
    pub fn set_last_spine_param(&mut self, par: f64) {
        self.ulspine = par;
    }

    /// OCCT .cxx Vertex(First, OnS).
    pub fn vertex(&self, first: bool, on_s: i32) -> &ChFiDS_CommonPoint {
        match (first, on_s) {
            (true, 1) => &self.pfirst_on_s1,
            (true, _) => &self.pfirst_on_s2,
            (false, 1) => &self.plast_on_s1,
            (false, _) => &self.plast_on_s2,
        }
    }

    /// OCCT .cxx ChangeVertex(First, OnS).
    pub fn change_vertex(&mut self, first: bool, on_s: i32) -> &mut ChFiDS_CommonPoint {
        match (first, on_s) {
            (true, 1) => &mut self.pfirst_on_s1,
            (true, _) => &mut self.pfirst_on_s2,
            (false, 1) => &mut self.plast_on_s1,
            (false, _) => &mut self.plast_on_s2,
        }
    }

    /// OCCT .cxx ChangeIndexOfS1(Index).
    pub fn change_index_of_s1(&mut self, index: i32) {
        self.index_of_s1 = index;
    }

    /// OCCT .cxx ChangeIndexOfS2(Index).
    pub fn change_index_of_s2(&mut self, index: i32) {
        self.index_of_s2 = index;
    }

    /// OCCT .cxx Index(OfS).
    pub fn index_of(&self, of_s: i32) -> i32 {
        if of_s == 1 {
            self.index_of_s1
        } else {
            self.index_of_s2
        }
    }

    /// OCCT .cxx Interference(OnS).
    pub fn interference(&self, on_s: i32) -> &ChFiDS_FaceInterference {
        if on_s == 1 {
            &self.intf1
        } else {
            &self.intf2
        }
    }

    /// OCCT .cxx ChangeInterference(OnS).
    pub fn change_interference(&mut self, on_s: i32) -> &mut ChFiDS_FaceInterference {
        if on_s == 1 {
            &mut self.intf1
        } else {
            &mut self.intf2
        }
    }

    /// OCCT .lxx VertexFirstOnS1().
    pub fn vertex_first_on_s1(&self) -> &ChFiDS_CommonPoint {
        &self.pfirst_on_s1
    }

    /// OCCT .lxx VertexLastOnS1().
    pub fn vertex_last_on_s1(&self) -> &ChFiDS_CommonPoint {
        &self.plast_on_s1
    }

    /// OCCT .lxx VertexFirstOnS2().
    pub fn vertex_first_on_s2(&self) -> &ChFiDS_CommonPoint {
        &self.pfirst_on_s2
    }

    /// OCCT .lxx VertexLastOnS2().
    pub fn vertex_last_on_s2(&self) -> &ChFiDS_CommonPoint {
        &self.plast_on_s2
    }

    /// OCCT .lxx IsOnCurve1().
    pub fn is_on_curve1(&self) -> bool {
        self.isoncurv1
    }

    /// OCCT .lxx IsOnCurve2().
    pub fn is_on_curve2(&self) -> bool {
        self.isoncurv2
    }

    /// OCCT .lxx Orientation().
    pub fn orientation(&self) -> Orientation {
        self.orientation
    }
}
// =========================================================================

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
    /// OCCT ChFiDS_CommonPoint.cxx L29-36 — Reset().
    pub fn reset(&mut self) {
        self.tol = 0.0;
        self.isvtx = false;
        self.isonarc = false;
        self.hasvector = false;
    }

    /// OCCT ChFiDS_CommonPoint.hxx — SetVertex(theVertex).
    pub fn set_vertex(&mut self, v: Shape) {
        self.vtx = v;
        self.isvtx = true;
    }

    /// OCCT ChFiDS_CommonPoint.cxx L52-64 — SetArc(Tol, A, Param, TArc);
    /// the tolerance is not crushed (only grown).
    pub fn set_arc(&mut self, tol: f64, arc: Shape, param: f64, tarc: Orientation) {
        self.isonarc = true;
        if tol > self.tol {
            self.tol = tol;
        }
        self.arc = arc;
        self.prmarc = param;
        self.traarc = tarc;
    }

    /// OCCT ChFiDS_CommonPoint.cxx L69-72 — SetParameter(Param) (the
    /// parameter in the tangency line).
    pub fn set_parameter(&mut self, param: f64) {
        self.prmtg = param;
    }

    /// OCCT ChFiDS_CommonPoint.hxx — SetPoint(thePoint).
    pub fn set_point(&mut self, p: DVec3) {
        self.point = p;
    }

    /// OCCT ChFiDS_CommonPoint.hxx — SetTolerance(Tol) (fuzziness, only
    /// grows).
    pub fn set_tolerance(&mut self, tol: f64) {
        if tol > self.tol {
            self.tol = tol;
        }
    }

    /// OCCT ChFiDS_CommonPoint.hxx — SetVector(theVector).
    pub fn set_vector(&mut self, v: DVec3) {
        self.hasvector = true;
        self.vector = v;
    }

    /// OCCT ChFiDS_CommonPoint.cxx L101-104 — Parameter() (in the
    /// tangency line).
    pub fn parameter(&self) -> f64 {
        self.prmtg
    }

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

    /// OCCT ChFiDS_Spine.lxx — SetTypeOfConcavity.
    pub fn set_type_of_concavity(&mut self, the_type: ChFiDS_TypeOfConcavity) {
        self.my_type_of_concavity = the_type;
    }

    /// OCCT ChFiDS_Spine.lxx — GetTypeOfConcavity.
    pub fn type_of_concavity(&self) -> ChFiDS_TypeOfConcavity {
        self.my_type_of_concavity
    }

    /// OCCT ChFiDS_Spine.lxx L131-134.
    pub fn set_offset_edges(&mut self, e: Shape) {
        self.offsetspine.push(e);
    }

    /// OCCT ChFiDS_Spine.lxx L138-142 — spine.InsertBefore(1, E).
    pub fn put_in_first(&mut self, e: Shape) {
        self.spine.insert(0, e);
    }

    /// OCCT ChFiDS_Spine.lxx L145-148 — offsetspine.InsertBefore(1, E).
    pub fn put_in_first_offset(&mut self, e: Shape) {
        self.offsetspine.insert(0, e);
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

    /// OCCT ChFiDS_Spine.lxx — SetTangencyExtremity(IsTangency, IsFirst).
    pub fn set_tangency_extremity(&mut self, is_tangency: bool, is_first: bool) {
        if is_first {
            self.firstistgt = is_tangency;
        } else {
            self.lastistgt = is_tangency;
        }
    }

    /// OCCT ChFiDS_Spine.lxx — IsTangencyExtremity(IsFirst).
    pub fn is_tangency_extremity(&self, is_first: bool) -> bool {
        if is_first {
            self.firstistgt
        } else {
            self.lastistgt
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

pub type SharedSurfData = Arc<std::sync::RwLock<ChFiDSSurfData>>;

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
    /// OCCT: occ::handle<Geom2d_Curve> pcrv1
    pub pcrv1: Option<rcad_kernel::geom::Curve2d>,
    /// OCCT: occ::handle<Geom2d_Curve> pcrv2
    pub pcrv2: Option<rcad_kernel::geom::Curve2d>,
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
            pcrv1: None,
            pcrv2: None,
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

    /// OCCT ChFiDS_Stripe.cxx L268-282 — InDS(First, Nb).
    pub fn in_ds(&mut self, first: bool, nb: i32) {
        if first {
            self.begfilled = nb;
        } else {
            self.endfilled = nb;
        }
    }

    /// OCCT ChFiDS_Stripe.cxx L285-292 — IsInDS(First).
    pub fn is_in_ds(&self, first: bool) -> i32 {
        if first {
            self.begfilled
        } else {
            self.endfilled
        }
    }

    /// OCCT ChFiDS_Stripe.lxx — FirstParameters(Pdeb, Pfin).
    pub fn first_parameters(&self) -> (f64, f64) {
        (self.pardeb1, self.parfin1)
    }

    /// OCCT ChFiDS_Stripe.lxx — LastParameters(Pdeb, Pfin).
    pub fn last_parameters(&self) -> (f64, f64) {
        (self.pardeb2, self.parfin2)
    }

    /// OCCT ChFiDS_Stripe.lxx — FirstCurve().
    pub fn first_curve(&self) -> i32 {
        self.index_ofcurve1
    }

    /// OCCT ChFiDS_Stripe.lxx — LastCurve().
    pub fn last_curve(&self) -> i32 {
        self.index_ofcurve2
    }

    /// OCCT ChFiDS_Stripe.lxx — FirstPCurve().
    pub fn first_pcurve(&self) -> Option<&rcad_kernel::geom::Curve2d> {
        self.pcrv1.as_ref()
    }

    /// OCCT ChFiDS_Stripe.lxx — LastPCurve().
    pub fn last_pcurve(&self) -> Option<&rcad_kernel::geom::Curve2d> {
        self.pcrv2.as_ref()
    }

    /// OCCT ChFiDS_Stripe.cxx — IndexPoint(First, OnS).
    pub fn index_point(&self, first: bool, on_s: i32) -> i32 {
        if first {
            if on_s == 1 {
                self.indexfirst_pon_s1
            } else {
                self.indexfirst_pon_s2
            }
        } else if on_s == 1 {
            self.indexlast_pon_s1
        } else {
            self.indexlast_pon_s2
        }
    }

    /// OCCT ChFiDS_Stripe.cxx — Orientation(OnS) (the face-side
    /// orientations myOr1/myOr2 set by StripeOrientations).
    pub fn orientation_on_s(&self, on_s: i32) -> Orientation {
        if on_s == 1 {
            self.my_or1
        } else {
            self.my_or2
        }
    }

    /// OCCT ChFiDS_Stripe.cxx — SetOrientation(Or, OnS).
    pub fn set_orientation_on_s(&mut self, or: Orientation, on_s: i32) {
        if on_s == 1 {
            self.my_or1 = or;
        } else {
            self.my_or2 = or;
        }
    }

    /// OCCT ChFiDS_Stripe.cxx — Orientation(First) (the pcurve
    /// orientation of the end curve, orcurv1/orcurv2).
    pub fn orientation(&self, first: bool) -> Orientation {
        if first {
            self.orcurv1
        } else {
            self.orcurv2
        }
    }

    /// OCCT ChFiDS_Stripe.cxx — SetOrientation(Or, First).
    pub fn set_orientation(&mut self, or: Orientation, first: bool) {
        if first {
            self.orcurv1 = or;
        } else {
            self.orcurv2 = or;
        }
    }

    /// OCCT ChFiDS_Stripe.lxx — FirstPCurveOrientation().
    pub fn first_pcurve_orientation(&self) -> Orientation {
        self.orcurv1
    }

    /// OCCT ChFiDS_Stripe.lxx — LastPCurveOrientation().
    pub fn last_pcurve_orientation(&self) -> Orientation {
        self.orcurv2
    }

    /// OCCT ChFiDS_Stripe.lxx — SolidIndex().
    pub fn solid_index(&self) -> i32 {
        self.index_of_solid
    }
}

// =========================================================================
// OCCT ChFiDS_Regul (ChFiDS_Regul.hxx + .cxx) — storage of a curve and its
// 2 faces or surfaces of support.  A negative S index encodes a surface
// (IsSurfaceN).
// =========================================================================
#[derive(Debug, Clone, Copy)]
pub struct ChFiDSRegul {
    icurv: i32,
    is1: i32,
    is2: i32,
}

impl Default for ChFiDSRegul {
    fn default() -> Self {
        ChFiDSRegul::new()
    }
}

impl ChFiDSRegul {
    /// OCCT ChFiDS_Regul.cxx L21-27.
    pub fn new() -> Self {
        ChFiDSRegul { icurv: 0, is1: 0, is2: 0 }
    }

    /// OCCT ChFiDS_Regul.cxx L30-33 (icurv = |IC|).
    pub fn set_curve(&mut self, ic: i32) {
        self.icurv = ic.abs();
    }

    /// OCCT ChFiDS_Regul.cxx L36-45 (face keeps the sign, surface negates).
    pub fn set_s1(&mut self, is1: i32, is_face: bool) {
        self.is1 = if is_face { is1.abs() } else { -is1.abs() };
    }

    /// OCCT ChFiDS_Regul.cxx L48-57.
    pub fn set_s2(&mut self, is2: i32, is_face: bool) {
        self.is2 = if is_face { is2.abs() } else { -is2.abs() };
    }

    /// OCCT ChFiDS_Regul.cxx L60-63.
    pub fn is_surface1(&self) -> bool {
        self.is1 < 0
    }

    /// OCCT ChFiDS_Regul.cxx L66-69.
    pub fn is_surface2(&self) -> bool {
        self.is2 < 0
    }

    /// OCCT ChFiDS_Regul.cxx L72-75.
    pub fn curve(&self) -> i32 {
        self.icurv
    }

    /// OCCT ChFiDS_Regul.cxx L78-81.
    pub fn s1(&self) -> i32 {
        self.is1.abs()
    }

    /// OCCT ChFiDS_Regul.cxx L84-87.
    pub fn s2(&self) -> i32 {
        self.is2.abs()
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

    // OCCT NCollection_IndexedDataMap — the key list (1-based FindKey).
    pub fn keys(&self) -> &Vec<Shape> {
        &self.my_keys
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
    // OCCT NCollection_IndexedDataMap — the insertion-ordered key list.
    pub fn keys(&self) -> &Vec<Shape> {
        &self.my_keys
    }

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
                    if tos == topods::ShapeType::Shell && tor == topods::ShapeType::Edge {
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
                    if tos == topods::ShapeType::Solid && tor == topods::ShapeType::Edge {
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


impl ChFiDSSpine {
    /// OCCT ChFiDS_Spine.cxx L352-373 — SetFirstTgt(W).
    pub fn set_first_tgt(&mut self, w: f64) {
        if self.is_periodic() {
            panic!("Standard_Failure: No extension by tangent on periodic contours");
        }
        // The flag is suspended if already positioned to avoid stopping d1.
        self.hasfirsttgt = false;
        let (p, t) = self.d1(w);
        self.firstori = p;
        self.firsttgt = t;
        // and it is reset.
        self.hasfirsttgt = true;
        self.firsttgtpar = w;
    }

    /// OCCT ChFiDS_Spine.cxx L375-395 — SetLastTgt(W).
    pub fn set_last_tgt(&mut self, w: f64) {
        if self.is_periodic() {
            panic!("Standard_Failure: No extension by tangent periodic contours");
        }
        self.haslasttgt = false;
        let (p, t) = self.d1(w);
        self.lastori = p;
        self.lasttgt = t;
        self.haslasttgt = true;
        self.lasttgtpar = w;
    }

    /// OCCT ChFiDS_Spine.cxx L410-424 — SetReference(W).
    pub fn set_reference(&mut self, w: f64) {
        self.hasref = true;
        let lll = self.abscissa.as_ref().map_or(0.0, |a| a[a.len() - 1]);
        if self.is_periodic() {
            self.valref = elclib_in_period(w, 0.0, lll);
        } else {
            self.valref = w;
        }
    }

    /// OCCT ChFiDS_Spine.cxx L426-439 — SetReference(I).
    pub fn set_reference_index(&mut self, i: usize) {
        self.hasref = true;
        let a = self.abscissa.as_ref();
        if i == 1 {
            self.valref = a.map_or(0.0, |a| a[0]) * 0.5;
        } else {
            self.valref = a.map_or(0.0, |a| (a[i - 1] + a[i - 2]) * 0.5);
        }
    }

    /// OCCT ChFiDS_Spine.cxx L441-480 — Index(W, Forward).
    pub fn index_of_param(&self, w: f64, forward: bool) -> usize {
        let t = self.tolesp.max(CONFUSION);
        let last = self.abscissa.as_ref().map_or(0.0, |a| a[a.len() - 1]);
        let len = self.abscissa.as_ref().map_or(0, |a| a.len());
        let mut par = w;
        if self.is_periodic() && par.abs() >= t && (par - last).abs() >= t {
            par = elclib_in_period(par, 0.0, last);
        }
        let mut f = 0.0;
        let mut l = 0.0;
        let mut ind = 1usize;
        while ind <= len {
            f = l;
            l = self.abscissa.as_ref().map_or(0.0, |a| a[ind - 1]);
            if par < l || ind == len {
                break;
            }
            ind += 1;
        }
        if forward && ind < len && (par - l).abs() < t {
            ind += 1;
        } else if !forward && ind > 1 && (par - f).abs() < t {
            ind -= 1;
        } else if forward && self.is_periodic() && ind == len && (par - l).abs() < t {
            ind = 1;
        } else if !forward && self.is_periodic() && ind == 1 && (par - f).abs() < t {
            ind = len;
        }
        ind
    }

    /// OCCT ChFiDS_Spine.cxx L619-698 — Prepare(L, Ind): resolves the
    /// elementary-spine index for the abscissa L, adjusting L for periodic
    /// contours and tangent extensions.  Returns the index; -1 encodes the
    /// first tangent extension, len+1 the last one (OCCT conventions).
    fn prepare(&self, l: &mut f64) -> i32 {
        let tol = self.tolesp.max(CONFUSION);
        let last = self.abscissa.as_ref().map_or(0.0, |a| a[a.len() - 1]);
        let len = self.abscissa.as_ref().map_or(0, |a| a.len()) as i32;
        let mut lval = *l;
        if self.is_periodic() && lval.abs() >= tol && (lval - last).abs() >= tol {
            lval = elclib_in_period(lval, 0.0, last);
        }

        let ind: i32;
        if self.hasfirsttgt && lval <= self.firsttgtpar {
            if self.hasref && self.valref >= lval && (lval - self.firsttgtpar).abs() <= tol {
                ind = self.index_of_param(lval, true) as i32;
            } else {
                ind = -1;
                lval -= self.firsttgtpar;
            }
        } else if lval <= 0.0 {
            ind = 1;
        } else if self.haslasttgt && lval >= self.lasttgtpar {
            if self.hasref && self.valref <= lval && (lval - self.lasttgtpar).abs() <= tol {
                ind = self.index_of_param(lval, true) as i32;
            } else {
                ind = len + 1;
                lval -= self.lasttgtpar;
            }
        } else if lval >= last {
            ind = len;
        } else {
            ind = self.index_of_param(lval, true) as i32;
        }
        *l = lval;
        ind
    }

    /// OCCT ChFiDS_Spine.cxx L537-560 — Absc(U, I): abscissa of the
    /// parameter U on the elementary spine I.
    pub fn absc_of_param(&mut self, u: f64, i: usize) -> f64 {
        self.indexofcurve = i as i32;
        let e = &self.spine[i - 1];
        let ed = e.as_edge().expect("not an edge");
        let curve = ed.curve.as_ref().expect("edge curve");
        let l0 = self.first_parameter_of(i);
        if e.orientation == Orientation::Reversed {
            l0
                + rcad_kernel::base::gcpnts::abscissa_point::arc_length(curve, u, ed.range[1])
        } else {
            l0
                + rcad_kernel::base::gcpnts::abscissa_point::arc_length(curve, ed.range[0], u)
        }
    }

    /// OCCT ChFiDS_Spine.cxx L562-601 — Parameter(AbsC, U, Oriented).
    pub fn parameter_at(&mut self, absc: f64, oriented: bool) -> f64 {
        let len = self.abscissa.as_ref().map_or(0, |a| a.len());
        let mut index = 1usize;
        while index < len {
            if absc < self.abscissa.as_ref().unwrap()[index - 1] {
                break;
            }
            index += 1;
        }
        self.parameter_on(absc, index, oriented)
    }

    /// OCCT ChFiDS_Spine.cxx L603-617 — Parameter(Index, AbsC, U, Oriented).
    pub fn parameter_on(&mut self, absc: f64, index: usize, oriented: bool) -> f64 {
        self.indexofcurve = index as i32;
        let e = &self.spine[index - 1];
        let ed = e.as_edge().expect("not an edge");
        let curve = ed.curve.as_ref().expect("edge curve").clone();
        let or = e.orientation;
        let l;
        if or == Orientation::Reversed {
            l = self.abscissa.as_ref().unwrap()[index - 1] - absc;
        } else if index == 1 {
            l = absc;
        } else {
            l = absc - self.abscissa.as_ref().unwrap()[index - 2];
        }
        let t = l / self.length_of(index);
        let (cf, cl) = (ed.range[0], ed.range[1]);
        let uapp = (1.0 - t) * cf + t * cl;
        let mut u = rcad_kernel::base::gcpnts::abscissa_point::abscissa_point_parameter(
            &curve, cf, cl, l, uapp,
        );
        if or == Orientation::Reversed && oriented {
            u = (cl + cf) - u;
        }
        u
    }

    /// OCCT ChFiDS_Spine.cxx L703-748 — Value(AbsC).
    pub fn value_at(&mut self, absc: f64) -> DVec3 {
        use rcad_kernel::geom::CurveEval as _;
        let mut l = absc;
        let index = self.prepare(&mut l);

        if index == -1 {
            let vp = self.firsttgt * l;
            return self.firstori + vp;
        } else if index as i32 == self.abscissa.as_ref().map_or(0, |a| a.len()) as i32 + 1 {
            let vp = self.lasttgt * l;
            return self.lastori + vp;
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
        curve.point_at(u)
    }

    /// OCCT ChFiDS_Spine.cxx L750-798 — D1(AbsC, P, V1): point and tangent
    /// (normalized, orientation-adjusted) on the composite spine.
    pub fn d1(&mut self, absc: f64) -> (DVec3, DVec3) {
        use rcad_kernel::geom::CurveEval as _;
        let mut l = absc;
        let index = self.prepare(&mut l);

        if index == -1 {
            let p = self.firstori + self.firsttgt * l;
            return (p, self.firsttgt);
        } else if index as i32 == self.abscissa.as_ref().map_or(0, |a| a.len()) as i32 + 1 {
            let p = self.lastori + self.lasttgt * l;
            return (p, self.lasttgt);
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
        let d1scale = 1.0 / v1.length();
        if e.orientation == Orientation::Reversed {
            v1 = -(v1 * d1scale);
        } else {
            v1 = v1 * d1scale;
        }
        (p, v1)
    }

    /// OCCT ChFiDS_Spine.hxx — SplitDone(B).
    pub fn set_split_done(&mut self, b: bool) {
        self.splitdone = b;
    }

    /// OCCT ChFiDS_Spine.hxx — SplitDone().
    pub fn split_done(&self) -> bool {
        self.splitdone
    }
}

/// OCCT ElCLib::InPeriod(U, UFirst, ULast).
pub fn elclib_in_period(u: f64, ufirst: f64, ulast: f64) -> f64 {
    let period = ulast - ufirst;
    if period <= 0.0 {
        return u;
    }
    let mut r = u - ((u - ufirst) / period).floor() * period;
    if r < ufirst {
        r += period;
    }
    if r >= ulast {
        r -= period;
    }
    r
}

impl ChFiDSSpine {
    /// OCCT ChFiDS_Spine.cxx Period() — the spine period (last abscissa).
    pub fn period(&self) -> f64 {
        self.abscissa.as_ref().map_or(0.0, |a| a[a.len() - 1])
    }
}
impl ChFiDSStripe {
    /// OCCT ChFiDS_Stripe.lxx — SetCurve(Icurv, IsFirst).
    pub fn set_curve(&mut self, icurv: i32, is_first: bool) {
        if is_first {
            self.index_ofcurve1 = icurv;
        } else {
            self.index_ofcurve2 = icurv;
        }
    }

    /// OCCT ChFiDS_Stripe.lxx — SetParameters(IsFirst, Pardeb, Parfin).
    pub fn set_parameters(&mut self, is_first: bool, pardeb: f64, parfin: f64) {
        if is_first {
            self.pardeb1 = pardeb;
            self.parfin1 = parfin;
        } else {
            self.pardeb2 = pardeb;
            self.parfin2 = parfin;
        }
    }

    /// OCCT ChFiDS_Stripe.lxx — ChangePCurve(IsFirst) (assignable slot).
    pub fn change_pcurve(&mut self, is_first: bool, pc: rcad_kernel::geom::Curve2d) {
        if is_first {
            self.pcrv1 = Some(pc);
        } else {
            self.pcrv2 = Some(pc);
        }
    }

    /// OCCT ChFiDS_Stripe.lxx — SetIndexPoint(Index, IsFirst, OnS).
    pub fn set_index_point(&mut self, index: i32, is_first: bool, on_s: i32) {
        if on_s == 1 {
            if is_first {
                self.indexfirst_pon_s1 = index;
            } else {
                self.indexlast_pon_s1 = index;
            }
        } else if is_first {
            self.indexfirst_pon_s2 = index;
        } else {
            self.indexlast_pon_s2 = index;
        }
    }

    /// OCCT ChFiDS_Stripe.lxx — SetSolidIndex(I).
    pub fn set_solid_index(&mut self, i: i32) {
        self.index_of_solid = i;
    }
}
impl ChFiDS_CommonPoint {
    /// OCCT ChFiDS_CommonPoint.lxx — IsOnArc().
    pub fn is_on_arc(&self) -> bool {
        self.isonarc
    }

    /// OCCT ChFiDS_CommonPoint.lxx — Arc().
    pub fn arc(&self) -> &Shape {
        &self.arc
    }

    /// OCCT ChFiDS_CommonPoint.lxx — ParameterOnArc().
    pub fn parameter_on_arc(&self) -> f64 {
        self.prmarc
    }

    /// OCCT ChFiDS_CommonPoint.lxx — TransitionOnArc().
    pub fn transition_on_arc(&self) -> Orientation {
        self.traarc
    }

    /// OCCT ChFiDS_CommonPoint.lxx — IsVertex().
    pub fn is_vertex(&self) -> bool {
        self.isvtx
    }

    /// OCCT ChFiDS_CommonPoint.lxx — Vertex().
    pub fn vertex(&self) -> &Shape {
        &self.vtx
    }

    /// OCCT ChFiDS_CommonPoint.lxx — Tolerance().
    pub fn tolerance(&self) -> f64 {
        self.tol
    }
}
