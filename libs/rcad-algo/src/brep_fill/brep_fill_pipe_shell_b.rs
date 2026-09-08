//! OCCT BRepFill_PipeShell — GAP placeholder types and split helpers
//! (TKBool/BRepFill sweep stack + TKGeomAlgo/TKMath dependencies not
//! translated yet, plan §0.6).
//!
//! The main file ([`super::brep_fill_pipe_shell`]) keeps the OCCT control
//! flow 1:1; the classes below carry the OCCT constructor/method surface so
//! the call sites stay line-by-line, and panic until the owning translation
//! units land:
//! - BRepFill_Section / BRepFill_SectionPlacement / BRepFill_Sweep /
//!   BRepFill_LocationLaw / BRepFill_SectionLaw / BRepFill_Edge3DLaw /
//!   BRepFill_EdgeOnSurfLaw / BRepFill_ACRLaw / BRepFill_ShapeLaw /
//!   BRepFill_NSections (TKBool/BRepFill — separate units);
//! - GeomFill_CurveAndTrihedron / GeomFill_Fixed / GeomFill_ConstantBiNormal /
//!   GeomFill_DiscreteTrihedron / GeomFill_GuideTrihedronAC /
//!   GeomFill_GuideTrihedronPlan / GeomFill_LocationGuide (TKGeomAlgo sweep
//!   machinery — the geomfill module currently covers the Filling/BSpline
//!   and trihedron-law units);
//! - Law_Interpol (TKMath/Law);
//! - BRepAdaptor_CompCurve (TKBRep adaptor);
//! - BRepBuilderAPI_Transform / BRepBuilderAPI_Copy (TKTopAlgo);
//! - IntCurveSurface_HInter (TKG3d/TKGeomAlgo);
//! - BRepGProp::LinearProperties (same gap as CompatibleWires);
//! - BRepLib_MakeFace(Wire, OnlyPlane) — reused from [`super::brep_fill_axe`].
//!
//! Real translations in this split unit: the plain OCCT enums
//! GeomFill_Trihedron / BRepFill_TransitionStyle / BRepFill_TypeOfContact.
//!
//! Note: OCCT GeomFill_PipeError maps to the rcad
//! `geomfill::trihedron_law::PipeError` (variant set PipeOk / PipeNotOk /
//! PlaneNotIntersectGuide / ImpossibleContact vs the OCCT PipeOk / PipeNotOk /
//! PipeNoSolution / PipeNotPlan).

use std::collections::HashMap;

use glam::DVec3;

use rcad_kernel::geom::Surface3;
use rcad_kernel::math::gp::Trsf;
use rcad_kernel::topo::topods::{BRep, GeomAbsShape, Shape};

use crate::brep_fill::generator::ShapeKey;

// ===========================================================================
// Real OCCT enums
// ===========================================================================

/// OCCT GeomFill_Trihedron (GeomFill_Trihedron.hxx) — the trihedron mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomFillTrihedron {
    /// OCCT GeomFill_IsFrenet.
    IsFrenet,
    /// OCCT GeomFill_IsCorrectedFrenet.
    IsCorrectedFrenet,
    /// OCCT GeomFill_IsDarboux.
    IsDarboux,
    /// OCCT GeomFill_IsDiscreteTrihedron.
    IsDiscreteTrihedron,
    /// OCCT GeomFill_IsFixed.
    IsFixed,
    /// OCCT GeomFill_IsConstantNormal.
    IsConstantNormal,
    /// OCCT GeomFill_IsGuideAC.
    IsGuideAC,
    /// OCCT GeomFill_IsGuideACWithContact.
    IsGuideACWithContact,
    /// OCCT GeomFill_IsGuidePlan.
    IsGuidePlan,
    /// OCCT GeomFill_IsGuidePlanWithContact.
    IsGuidePlanWithContact,
}

/// OCCT BRepFill_TransitionStyle (BRepFill_TransitionStyle.hxx).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BRepFillTransitionStyle {
    /// OCCT BRepFill_Modified.
    Modified,
    /// OCCT BRepFill_Round.
    Round,
    /// OCCT BRepFill_RightCorner.
    RightCorner,
    /// OCCT BRepFill_RoundCorner.
    RoundCorner,
}

/// OCCT BRepFill_TypeOfContact (BRepFill_TypeOfContact.hxx).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BRepFillTypeOfContact {
    /// OCCT BRepFill_NoContact.
    NoContact,
    /// OCCT BRepFill_Contact.
    Contact,
    /// OCCT BRepFill_ContactOnBorder.
    ContactOnBorder,
}

/// OCCT `NCollection_HArray2<TopoDS_Shape>` mapping (the sweep result arrays).
pub type ShapeHArray2 = Vec<Vec<Shape>>;

/// OCCT `NCollection_Map<TopoDS_Shape>` mapping.
pub type ShapeSet = std::collections::HashSet<ShapeKey>;

/// OCCT `NCollection_DataMap<TopoDS_Shape, handle<HArray2<TopoDS_Shape>>>`
/// mapping.
pub type ShapeToArray2Map = HashMap<ShapeKey, ShapeHArray2>;

// ===========================================================================
// Law_Interpol / Law_Function (TKMath/Law)
// ===========================================================================

/// OCCT Law_Interpol — GAP: TKMath/Law is not translated (plan §0.6).
pub struct LawInterpol;

impl LawInterpol {
    /// OCCT Law_Interpol().
    pub fn new() -> Self {
        panic!("GAP: Law_Interpol (TKMath/Law) is not translated — see file header")
    }

    /// OCCT Law_Interpol::Set(ParAndRad, IsPeriodic).
    pub fn set(&self, _par_and_rad: &[glam::DVec2], _is_periodic: bool) {
        panic!("GAP: Law_Interpol::Set (TKMath/Law) is not translated — see file header")
    }
}

impl Default for LawInterpol {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT handle<Law_Function> myLaw — the base-class handle; the only concrete
/// law used by PipeShell is Law_Interpol.
#[derive(Clone)]
pub struct LawFunction;

impl LawFunction {
    /// OCCT myLaw = new Law_Interpol().
    pub fn new_interpol() -> Self {
        LawFunction
    }

    /// OCCT occ::down_cast<Law_Interpol>(myLaw).
    pub fn downcast_interpol(&self) -> &LawInterpol {
        panic!(
            "GAP: Law_Interpol (TKMath/Law) is not translated — see file header"
        )
    }

    /// OCCT `myLaw = L` (the SetLaw argument).
    pub fn from_handle(_l: &LawFunction) -> Self {
        LawFunction
    }
}

// ===========================================================================
// BRepFill_Section (TKBool/BRepFill — separate translation unit)
// ===========================================================================

/// OCCT BRepFill_Section — GAP: own file in TKBool/BRepFill, not part of this
/// batch (plan §0.6).
#[derive(Clone)]
pub struct BRepFillSection;

impl BRepFillSection {
    /// OCCT BRepFill_Section(Profile, Location, WithContact, WithCorrection).
    pub fn new(
        _profile: &Shape,
        _location: &Shape,
        _with_contact: bool,
        _with_correction: bool,
    ) -> Self {
        panic!(
            "GAP: BRepFill_Section (TKBool/BRepFill) is not translated — see file header"
        )
    }

    /// OCCT BRepFill_Section::Set(IsLaw).
    pub fn set(&mut self, _is_law: bool) {
        panic!(
            "GAP: BRepFill_Section (TKBool/BRepFill) is not translated — see file header"
        )
    }

    /// OCCT BRepFill_Section::OriginalShape().
    pub fn original_shape(&self) -> Shape {
        panic!(
            "GAP: BRepFill_Section (TKBool/BRepFill) is not translated — see file header"
        )
    }

    /// OCCT BRepFill_Section::Wire().
    pub fn wire(&self) -> Shape {
        panic!(
            "GAP: BRepFill_Section (TKBool/BRepFill) is not translated — see file header"
        )
    }

    /// OCCT BRepFill_Section::Vertex().
    pub fn vertex(&self) -> Shape {
        panic!(
            "GAP: BRepFill_Section (TKBool/BRepFill) is not translated — see file header"
        )
    }

    /// OCCT BRepFill_Section::WithContact().
    pub fn with_contact(&self) -> bool {
        panic!(
            "GAP: BRepFill_Section (TKBool/BRepFill) is not translated — see file header"
        )
    }

    /// OCCT BRepFill_Section::WithCorrection().
    pub fn with_correction(&self) -> bool {
        panic!(
            "GAP: BRepFill_Section (TKBool/BRepFill) is not translated — see file header"
        )
    }

    /// OCCT BRepFill_Section::IsLaw().
    pub fn is_law(&self) -> bool {
        panic!(
            "GAP: BRepFill_Section (TKBool/BRepFill) is not translated — see file header"
        )
    }

    /// OCCT BRepFill_Section::IsPunctual().
    pub fn is_punctual(&self) -> bool {
        panic!(
            "GAP: BRepFill_Section (TKBool/BRepFill) is not translated — see file header"
        )
    }

    /// OCCT BRepFill_Section::ModifiedShape(TheShape).
    pub fn modified_shape(&self, _the_shape: &Shape) -> Shape {
        panic!(
            "GAP: BRepFill_Section (TKBool/BRepFill) is not translated — see file header"
        )
    }
}

// ===========================================================================
// BRepFill_LocationLaw / BRepFill_SectionLaw and derived laws
// (TKBool/BRepFill — separate translation units)
// ===========================================================================

/// OCCT handle<BRepFill_LocationLaw> — GAP (plan §0.6).
pub struct BRepFillLocationLaw;

impl BRepFillLocationLaw {
    /// OCCT DeleteTransform().
    pub fn delete_transform(&mut self) {
        panic!("GAP: BRepFill_LocationLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT TransformInG0Law().
    pub fn transform_in_g0_law(&mut self) {
        panic!("GAP: BRepFill_LocationLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT TransformInCompatibleLaw(Angmin).
    pub fn transform_in_compatible_law(&mut self, _angmin: f64) {
        panic!("GAP: BRepFill_LocationLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT NbLaw().
    pub fn nb_law(&self) -> usize {
        panic!("GAP: BRepFill_LocationLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT CurvilinearBounds(Index, First, Last).
    pub fn curvilinear_bounds(&self, _index: usize, _first: &mut f64, _last: &mut f64) {
        panic!("GAP: BRepFill_LocationLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT D0(U, W).
    pub fn d0(&self, _u: f64, _w: &mut Shape) {
        panic!("GAP: BRepFill_LocationLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT IsClosed().
    pub fn is_closed(&self) -> bool {
        panic!("GAP: BRepFill_LocationLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT IsG1(Index) — the continuity rank (>= 0 means G1).
    pub fn is_g1(&self, _index: f64) -> i32 {
        panic!("GAP: BRepFill_LocationLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT GetStatus().
    pub fn get_status(&self) -> crate::geomalgo::geomfill::trihedron_law::PipeError {
        panic!("GAP: BRepFill_LocationLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT Law(Index) — the polymorphic GeomFill location-law handle.
    pub fn law(&self, _index: usize) -> GeomFillLocationLawHandle {
        panic!("GAP: BRepFill_LocationLaw (TKBool/BRepFill) is not translated — see file header")
    }
}

/// OCCT handle<BRepFill_SectionLaw> — GAP (plan §0.6).
pub struct BRepFillSectionLaw;

impl BRepFillSectionLaw {
    /// OCCT Law(1)->GetDomain(First, Last) (the single-section path).
    pub fn get_domain(&self, _first: &mut f64, _last: &mut f64) {
        panic!("GAP: BRepFill_SectionLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT D0(U, W).
    pub fn d0(&self, _u: f64, _w: &mut Shape) {
        panic!("GAP: BRepFill_SectionLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT IsVClosed().
    pub fn is_vclosed(&self) -> bool {
        panic!("GAP: BRepFill_SectionLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT IsUClosed().
    pub fn is_uclosed(&self) -> bool {
        panic!("GAP: BRepFill_SectionLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT NbLaw().
    pub fn nb_law(&self) -> usize {
        panic!("GAP: BRepFill_SectionLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT IndexOfEdge(E) — signed index (0 when the edge is not a law edge).
    pub fn index_of_edge(&self, _e: &Shape) -> i32 {
        panic!("GAP: BRepFill_SectionLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT ConcatenedLaw().
    pub fn concatened_law(&self) -> GeomFillSectionLawHandle {
        panic!("GAP: BRepFill_SectionLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        panic!("GAP: BRepFill_SectionLaw (TKBool/BRepFill) is not translated — see file header")
    }
}

/// OCCT new BRepFill_Edge3DLaw(Spine, Law) — GAP (plan §0.6).  The OCCT
/// parameter is the base handle<GeomFill_LocationLaw>; the rcad enum
/// models the two concrete subclasses reaching this constructor
/// (PipeShell.cxx L271/287/301/315 CurveAndTrihedron, L429/1020
/// LocationGuide).
pub fn edge3d_law_new(
    _spine: &Shape,
    _loc: GeomFillLocationLawHandle,
) -> BRepFillLocationLaw {
    panic!("GAP: BRepFill_Edge3DLaw (TKBool/BRepFill) is not translated — see file header")
}

/// OCCT new BRepFill_ACRLaw(Spine, Loc) — GAP (plan §0.6).
pub fn acr_law_new(_spine: &Shape, _loc: GeomFillLocationGuide) -> BRepFillLocationLaw {
    panic!("GAP: BRepFill_ACRLaw (TKBool/BRepFill) is not translated — see file header")
}

/// OCCT new BRepFill_EdgeOnSurfLaw(Spine, Support) — GAP (plan §0.6).
pub fn edge_on_surf_law_new(_spine: &Shape, _support: &Shape) -> BRepFillEdgeOnSurfLaw {
    panic!("GAP: BRepFill_EdgeOnSurfLaw (TKBool/BRepFill) is not translated — see file header")
}

/// OCCT handle<BRepFill_EdgeOnSurfLaw> — GAP (plan §0.6).
pub struct BRepFillEdgeOnSurfLaw;

impl BRepFillEdgeOnSurfLaw {
    /// OCCT HasResult().
    pub fn has_result(&self) -> bool {
        panic!("GAP: BRepFill_EdgeOnSurfLaw (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT upcast to the location-law handle (`myLocation = loc`).
    pub fn into_location_law(self) -> BRepFillLocationLaw {
        panic!("GAP: BRepFill_EdgeOnSurfLaw (TKBool/BRepFill) is not translated — see file header")
    }
}

/// OCCT new BRepFill_ShapeLaw(Wire) — GAP (plan §0.6).
pub fn shape_law_new(_wire: &Shape) -> BRepFillSectionLaw {
    panic!("GAP: BRepFill_ShapeLaw (TKBool/BRepFill) is not translated — see file header")
}

/// OCCT new BRepFill_ShapeLaw(Wire, Law) — GAP (plan §0.6).
pub fn shape_law_new_with_law(_wire: &Shape, _law: &LawFunction) -> BRepFillSectionLaw {
    panic!("GAP: BRepFill_ShapeLaw (TKBool/BRepFill) is not translated — see file header")
}

/// OCCT new BRepFill_NSections(Sections, Trsfs, Params, First, Last) — GAP
/// (plan §0.6).
pub fn nsections_law_new(
    _sections: &[Shape],
    _trsfs: &[Trsf],
    _params: &[f64],
    _first: f64,
    _last: f64,
) -> BRepFillSectionLaw {
    panic!("GAP: BRepFill_NSections (TKBool/BRepFill) is not translated — see file header")
}

// ===========================================================================
// GeomFill sweep machinery (TKGeomAlgo)
// ===========================================================================

/// OCCT handle<GeomFill_CurveAndTrihedron> — GAP (plan §0.6).
pub struct GeomFillCurveAndTrihedron;

impl GeomFillCurveAndTrihedron {
    /// OCCT new GeomFill_CurveAndTrihedron(TLaw).
    pub fn new(_t_law: GeomFillTrihedronLawHandle) -> Self {
        panic!(
            "GAP: GeomFill_CurveAndTrihedron (TKGeomAlgo) is not translated — see file header"
        )
    }
}

/// OCCT handle<GeomFill_Fixed> — GAP (plan §0.6).
pub struct GeomFillFixed;

impl GeomFillFixed {
    /// OCCT new GeomFill_Fixed(V1, V2).
    pub fn new(_v1: DVec3, _v2: DVec3) -> Self {
        panic!("GAP: GeomFill_Fixed (TKGeomAlgo) is not translated — see file header")
    }

    /// OCCT upcast to GeomFill_TrihedronLaw.
    pub fn into_trihedron_law(self) -> GeomFillTrihedronLawHandle {
        panic!("GAP: GeomFill_Fixed (TKGeomAlgo) is not translated — see file header")
    }
}

/// OCCT handle<GeomFill_ConstantBiNormal> — GAP (plan §0.6).
pub struct GeomFillConstantBiNormal;

impl GeomFillConstantBiNormal {
    /// OCCT new GeomFill_ConstantBiNormal(BiNormal).
    pub fn new(_binormal: DVec3) -> Self {
        panic!("GAP: GeomFill_ConstantBiNormal (TKGeomAlgo) is not translated — see file header")
    }

    /// OCCT upcast to GeomFill_TrihedronLaw.
    pub fn into_trihedron_law(self) -> GeomFillTrihedronLawHandle {
        panic!("GAP: GeomFill_ConstantBiNormal (TKGeomAlgo) is not translated — see file header")
    }
}

/// OCCT handle<GeomFill_TrihedronLaw> — GAP carrier over the trihedron kinds
/// used by PipeShell (Frenet / CorrectedFrenet exist in rcad geomfill, but the
/// CurveAndTrihedron consumption path is not ported).
pub enum GeomFillTrihedronLawHandle {
    /// OCCT new GeomFill_Frenet().
    Frenet,
    /// OCCT new GeomFill_CorrectedFrenet().
    CorrectedFrenet,
    /// OCCT new GeomFill_DiscreteTrihedron().
    Discrete,
    /// OCCT GeomFill_Fixed.
    Fixed(GeomFillFixed),
    /// OCCT GeomFill_ConstantBiNormal.
    ConstantBiNormal(GeomFillConstantBiNormal),
    /// OCCT GeomFill_GuideTrihedronAC.
    GuideAC(GeomFillGuideTrihedronAC),
}

/// OCCT handle<GeomFill_DiscreteTrihedron> — GAP (plan §0.6).
pub struct GeomFillDiscreteTrihedron;

impl GeomFillDiscreteTrihedron {
    /// OCCT new GeomFill_DiscreteTrihedron().
    pub fn new() -> Self {
        panic!("GAP: GeomFill_DiscreteTrihedron (TKGeomAlgo) is not translated — see file header")
    }

    /// OCCT upcast to GeomFill_TrihedronLaw.
    pub fn into_trihedron_law(self) -> GeomFillTrihedronLawHandle {
        panic!("GAP: GeomFill_DiscreteTrihedron (TKGeomAlgo) is not translated — see file header")
    }
}

impl Default for GeomFillDiscreteTrihedron {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT handle<GeomFill_GuideTrihedronAC> — GAP (plan §0.6).
pub struct GeomFillGuideTrihedronAC;

impl GeomFillGuideTrihedronAC {
    /// OCCT new GeomFill_GuideTrihedronAC(Guide).
    pub fn new(_guide: &BRepAdaptorCompCurve) -> Self {
        panic!("GAP: GeomFill_GuideTrihedronAC (TKGeomAlgo) is not translated — see file header")
    }
}

/// OCCT handle<GeomFill_GuideTrihedronPlan> — GAP (plan §0.6).
pub struct GeomFillGuideTrihedronPlan;

impl GeomFillGuideTrihedronPlan {
    /// OCCT new GeomFill_GuideTrihedronPlan(Guide).
    pub fn new(_guide: &BRepAdaptorCompCurve) -> Self {
        panic!("GAP: GeomFill_GuideTrihedronPlan (TKGeomAlgo) is not translated — see file header")
    }
}

/// OCCT handle<GeomFill_LocationGuide> — GAP (plan §0.6).
pub struct GeomFillLocationGuide;

impl GeomFillLocationGuide {
    /// OCCT new GeomFill_LocationGuide(TLaw).
    pub fn new(_t_law: &GeomFillGuideTrihedronAC) -> Self {
        panic!("GAP: GeomFill_LocationGuide (TKGeomAlgo) is not translated — see file header")
    }

    /// OCCT new GeomFill_LocationGuide(TLaw) — the plan-guide overload form.
    pub fn new_plan(_t_law: &GeomFillGuideTrihedronPlan) -> Self {
        panic!("GAP: GeomFill_LocationGuide (TKGeomAlgo) is not translated — see file header")
    }

    /// OCCT ComputeAutomaticLaw(ParAndRad).
    pub fn compute_automatic_law(&self) -> Vec<glam::DVec2> {
        panic!("GAP: GeomFill_LocationGuide (TKGeomAlgo) is not translated — see file header")
    }

    /// OCCT GetCurve() — the spine curve of the guide.
    pub fn get_curve(&self) -> IntCurveSurfaceHCurve {
        panic!("GAP: GeomFill_LocationGuide (TKGeomAlgo) is not translated — see file header")
    }

    /// OCCT Guide() — the guide curve.
    pub fn guide(&self) -> IntCurveSurfaceHCurve {
        panic!("GAP: GeomFill_LocationGuide (TKGeomAlgo) is not translated — see file header")
    }

    /// OCCT Set(Sec, Bool, First, Last, OldAngle, Angle).
    #[allow(clippy::too_many_arguments)]
    pub fn set_law(
        &self,
        _sec: &GeomFillSectionLawHandle,
        _boolean: bool,
        _first: f64,
        _last: f64,
        _old_angle: f64,
        _angle: &mut f64,
    ) {
        panic!("GAP: GeomFill_LocationGuide (TKGeomAlgo) is not translated — see file header")
    }

    /// OCCT EraseRotation().
    pub fn erase_rotation(&self) {
        panic!("GAP: GeomFill_LocationGuide (TKGeomAlgo) is not translated — see file header")
    }
}

/// The polymorphic `myLocation->Law(i)` handle (GeomFill_LocationGuide /
/// GeomFill_CurveAndTrihedron downcasts).
pub enum GeomFillLocationLawHandle {
    /// OCCT occ::down_cast<GeomFill_LocationGuide>(...).
    LocationGuide(GeomFillLocationGuide),
    /// OCCT the CurveAndTrihedron form.
    CurveAndTrihedron(GeomFillCurveAndTrihedron),
}

/// OCCT handle<GeomFill_SectionLaw> (ConcatenedLaw result) — GAP (plan §0.6).
pub struct GeomFillSectionLawHandle;

impl GeomFillSectionLawHandle {
    /// OCCT GetDomain(First, Last).
    pub fn get_domain(&self, _first: &mut f64, _last: &mut f64) {
        panic!("GAP: GeomFill_SectionLaw (TKGeomAlgo) is not translated — see file header")
    }
}

// ===========================================================================
// BRepAdaptor_CompCurve / BRepBuilderAPI_Transform / Copy / IntCurveSurface
// ===========================================================================

/// OCCT handle<BRepAdaptor_CompCurve> — GAP: the TKBRep adaptor is not
/// translated (plan §0.6).
pub struct BRepAdaptorCompCurve;

impl BRepAdaptorCompCurve {
    /// OCCT new BRepAdaptor_CompCurve(W).
    pub fn new(_w: &Shape) -> Self {
        panic!("GAP: BRepAdaptor_CompCurve (TKBRep) is not translated — see file header")
    }

    /// OCCT D1(U, P, V).
    pub fn d1(&self, _u: f64, _p: &mut DVec3, _v: &mut DVec3) {
        panic!("GAP: BRepAdaptor_CompCurve (TKBRep) is not translated — see file header")
    }
}

/// OCCT handle<Adaptor3d_Curve> as consumed by IntCurveSurface_HInter — GAP.
pub struct IntCurveSurfaceHCurve;

/// OCCT BRepBuilderAPI_Transform — GAP (TKTopAlgo; the BRepBuilderAPI facade
/// unit is separate).
pub struct BRepBuilderAPITransform;

impl BRepBuilderAPITransform {
    /// OCCT BRepBuilderAPI_Transform(S, T, Copy).
    pub fn new(_s: &Shape, _t: &Trsf, _copy: bool) -> Self {
        panic!("GAP: BRepBuilderAPI_Transform (TKTopAlgo) is not translated — see file header")
    }

    /// OCCT Shape() — the transformed shape.
    pub fn shape(&self) -> Shape {
        panic!("GAP: BRepBuilderAPI_Transform (TKTopAlgo) is not translated — see file header")
    }
}

/// OCCT BRepBuilderAPI_Copy — GAP (TKTopAlgo).
pub struct BRepBuilderAPICopy;

impl BRepBuilderAPICopy {
    /// OCCT BRepBuilderAPI_Copy(S).
    pub fn new(_s: &Shape) -> Self {
        panic!("GAP: BRepBuilderAPI_Copy (TKTopAlgo) is not translated — see file header")
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        panic!("GAP: BRepBuilderAPI_Copy (TKTopAlgo) is not translated — see file header")
    }

    /// OCCT Shape().
    pub fn shape(&self) -> Shape {
        panic!("GAP: BRepBuilderAPI_Copy (TKTopAlgo) is not translated — see file header")
    }
}

/// OCCT IntCurveSurface_HInter — GAP (plan §0.6).
pub struct IntCurveSurfaceHInter;

impl IntCurveSurfaceHInter {
    /// OCCT Perform(HCurve, SurfAdaptor).
    pub fn perform(&mut self, _curve: &IntCurveSurfaceHCurve, _plane: &Surface3) {
        panic!("GAP: IntCurveSurface_HInter (TKG3d) is not translated — see file header")
    }

    /// OCCT NbPoints().
    pub fn nb_points(&self) -> usize {
        panic!("GAP: IntCurveSurface_HInter (TKG3d) is not translated — see file header")
    }

    /// OCCT Point(j).Pnt().
    pub fn point_pnt(&self, _j: usize) -> DVec3 {
        panic!("GAP: IntCurveSurface_HInter (TKG3d) is not translated — see file header")
    }
}

/// OCCT BRepGProp::LinearProperties — GAP (same as CompatibleWires).
pub fn brep_gprop_linear_properties_centre(_profile: &Shape) -> DVec3 {
    panic!(
        "GAP: BRepGProp::LinearProperties (TKTopAlgo) is not translated — see file header"
    )
}

// ===========================================================================
// BRepFill_SectionPlacement / BRepFill_Sweep (TKBool/BRepFill — separate units)
// ===========================================================================

/// OCCT BRepFill_SectionPlacement — GAP (plan §0.6).
pub struct BRepFillSectionPlacement;

impl BRepFillSectionPlacement {
    /// OCCT BRepFill_SectionPlacement(Location, Section, Vertex, WithContact,
    /// WithCorrection).
    pub fn new(
        _location: &BRepFillLocationLaw,
        _section: &Shape,
        _vertex: &Shape,
        _with_contact: bool,
        _with_correction: bool,
    ) -> Self {
        panic!(
            "GAP: BRepFill_SectionPlacement (TKBool/BRepFill) is not translated — see file header"
        )
    }

    /// OCCT Transformation().
    pub fn transformation(&self) -> Trsf {
        panic!(
            "GAP: BRepFill_SectionPlacement (TKBool/BRepFill) is not translated — see file header"
        )
    }

    /// OCCT AbscissaOnPath().
    pub fn abscissa_on_path(&self) -> f64 {
        panic!(
            "GAP: BRepFill_SectionPlacement (TKBool/BRepFill) is not translated — see file header"
        )
    }
}

/// OCCT BRepFill_Sweep — GAP (plan §0.6).
pub struct BRepFillSweep;

impl BRepFillSweep {
    /// OCCT BRepFill_Sweep(Section, Location, Boolean).
    pub fn new(
        _section: &BRepFillSectionLaw,
        _location: &BRepFillLocationLaw,
        _boolean: bool,
    ) -> Self {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT SetTolerance(Tol3d, BoundTol, TolCurve, TolAngular).
    pub fn set_tolerance(&mut self, _tol3d: f64, _bound_tol: f64, _tol_curve: f64, _tol_angular: f64) {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT SetAngularControl(Angmin, Angmax).
    pub fn set_angular_control(&mut self, _angmin: f64, _angmax: f64) {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT SetForceApproxC1(ForceApproxC1).
    pub fn set_force_approx_c1(&mut self, _force_approx_c1: bool) {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT SetBounds(BoundFirst, BoundLast).
    pub fn set_bounds(&mut self, _bound_first: &Shape, _bound_last: &Shape) {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT Build(WDone, WInter, WSubS, Transition, Continuity, Whatdegree,
    /// Degmax, Segmax).
    #[allow(clippy::too_many_arguments)]
    pub fn build(
        &mut self,
        _w_done: &mut ShapeSet,
        _w_inter: &mut ShapeToArray2Map,
        _w_sub_s: &mut ShapeToArray2Map,
        _transition: BRepFillTransitionStyle,
        _continuity: GeomAbsShape,
        _what_degree: i32,
        _degmax: i32,
        _segmax: i32,
    ) {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT Shape().
    pub fn shape(&self) -> Shape {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT ErrorOnSurface().
    pub fn error_on_surface(&self) -> f64 {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT InterFaces() — the u-edges array.
    pub fn inter_faces(&self) -> ShapeHArray2 {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT SubShape() — the spine faces array.
    pub fn sub_shape(&self) -> ShapeHArray2 {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT Sections() — the v-edges array.
    pub fn sections(&self) -> ShapeHArray2 {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill) is not translated — see file header")
    }

    /// OCCT Tape(I) — the tape (shell) generated by the I-th section edge.
    pub fn tape(&self, _i: usize) -> Shape {
        panic!("GAP: BRepFill_Sweep (TKBool/BRepFill) is not translated — see file header")
    }
}

/// OCCT GeomFill_Location constant for BRepFill_Sweep::Build.
pub const GEOM_FILL_LOCATION: i32 = 1;

/// OCCT static BRepFill::SearchOrigin(W, P, Dir, Tol) (BRepFill.cxx L885+) —
/// GAP: the BRepFill.cxx statics beyond Axe are a separate unit (the
/// CompatibleWires port carries its own private SearchOrigin method).
pub fn brep_fill_search_origin(
    _brep: &mut BRep,
    _w: &mut Shape,
    _p: DVec3,
    _dir: DVec3,
    _tol: f64,
) {
    panic!(
        "GAP: BRepFill::SearchOrigin (TKBool/BRepFill, BRepFill.cxx) is not translated — \
         see file header"
    )
}
