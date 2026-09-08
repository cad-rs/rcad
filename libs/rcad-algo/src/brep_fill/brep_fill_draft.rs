//! OCCT BRepFill_Draft — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepFill/
//!         BRepFill_Draft.hxx (L14-102) +
//!         BRepFill_Draft.cxx (L1-975).
//!
//! OCCT inheritance chain: none — BRepFill_Draft is a standalone class.
//!
//! Architecture differences (referenced from the affected functions):
//! 1. The sweep engine stack (BRepFill_Sweep, BRepFill_ShapeLaw /
//!    BRepFill_SectionLaw, BRepFill_DraftLaw / BRepFill_LocationLaw,
//!    GeomFill_LocationDraft, BndLib_Add3dCurve / BndLib_AddSurface,
//!    BRepAdaptor_Surface, GeomLProp_SLProps, BRepLib_MakeFace,
//!    BRepBuilderAPI_Sewing, BRepExtrema_DistShapeShape, BRepAlgoAPI_Section
//!    + the BOPAlgo PaveFiller / Builder public API, BRepTools::UVBounds)
//!    is not translated yet — the OCCT-class-named GAP carriers below keep
//!    the call form (plan D3).  The sweep / laws / arrays are shared with
//!    the brep_fill_pipe.rs translation.
//! 2. NCollection_List / NCollection_IndexedDataMap ->
//!    Vec<Shape> / IndexMap keyed by the (TShape ptr, location) identity.
//! 3. TopoDS_Iterator -> `topods_iterator` (the brep_fill_pipe.rs shared
//!    helper); TopExp_Explorer -> feat::brep_feat_builder::explorer.
//! 4. gp_Trsf / gp_Ax3 / gp_Mat / gp_Dir -> glam DAffine3 / DMat3 / DVec3;
//!    the pure-math helpers are re-hosted below (the loc_ope_pipe.rs
//!    precedent).
//! 5. Bnd_Box -> rcad_kernel::math::bnd::BndBox.  The closed-flag mutation
//!    of myWire (cxx L271) requires the owning BRep pool in rcad; the flag
//!    update goes through a local pool (the mutation is carried by the flag
//!    of the shape the pool created — the caller pool is untouched).
//!
//! first consumer: BRepOffsetAPI_MakeDraft (the 2e translation; the GAP
//! carrier there is rewired to this module).

use std::collections::HashMap;

use glam::{DAffine3, DMat3, DVec3};
use indexmap::IndexMap;

use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::geom::{Curve2dEval, Curve3, CurveEval, Line3, Plane, Surface3, TrimmedSurface};
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, ShapeType, TShape, tshape_flags};

use std::cell::RefCell;
use std::rc::Rc;

use crate::brep_fill::brep_fill_pipe::{
    add_to_shell, add_to_solid, set_closed_flag, shape_reversed, topods_iterator,
    BRepFillLocationLaw, BRepFillSectionLaw, BRepFillSweep, GeomAbsShape,
    GeomFillApproxStyle, ShapeArray2, BRepFillTransitionStyle, TOPABS_IN,
};
use crate::brep_fill::brep_fill_section_law::BRepFillSectionLawOps;
// The real OCCT BRepFill_ShapeLaw (TKBool/BRepFill) — the section law of
// BRepFill_Draft::Init (cxx L459).  The location-law chain (the real
// BRepFill_DraftLaw / BRepFill_Sweep) stays on the placeholder carriers
// below: the real laws consume the owning `BRep` pool while BRepFill_Draft
// works on detached shapes (architecture difference #5), and the
// GeomFill_LocationDraft translation has not landed.
use crate::brep_fill::brep_fill_shape_law::BRepFillShapeLaw;
use crate::brep_fill::generator::ShapeKey;

/// OCCT Precision::Confusion().
const TOL_CONFUSION: f64 = CONFUSION;

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
fn shape_key(s: &Shape) -> ShapeKey {
    ShapeKey(s.ptr_id())
}

// ---------------------------------------------------------------------------
// GAP carriers — the draft sweep stack (architecture difference #1)
// ---------------------------------------------------------------------------

/// GAP: GeomFill_LocationDraft (TKGeomAlgo/GeomFill, parallel batch) — the
/// draft location law (D0 returns the frame + draft direction).
#[derive(Debug, Clone)]
pub struct GeomFillLocationDraft {
    #[allow(dead_code)]
    my_dir: DVec3,
    #[allow(dead_code)]
    my_angle: f64,
}

impl GeomFillLocationDraft {
    /// OCCT new GeomFill_LocationDraft(Dir, Angle).
    pub fn new(the_dir: DVec3, the_angle: f64) -> Self {
        GeomFillLocationDraft {
            my_dir: the_dir,
            my_angle: the_angle,
        }
    }

    /// OCCT GeomFill_LocationDraft::SetAngle(Angle).
    pub fn set_angle(&mut self, _the_angle: f64) {
        panic!("GAP: GeomFill_LocationDraft (GeomFill batch) — see file header")
    }

    /// OCCT GeomFill_LocationLaw::GetDomain(f, l).
    pub fn get_domain(&self, _first: &mut f64, _last: &mut f64) {
        panic!("GAP: GeomFill_LocationDraft (GeomFill batch) — see file header")
    }

    /// OCCT GeomFill_LocationLaw::D0(param, M, V).
    pub fn d0(&self, _param: f64, _m: &mut DMat3, _v: &mut DVec3) {
        panic!("GAP: GeomFill_LocationDraft (GeomFill batch) — see file header")
    }
}

/// GAP: Adaptor3d_Curve handle — the curve handle consumed by
/// GoodOrientation (Law(Ind)->GetCurve()).
pub struct Adaptor3dCurve;

impl Adaptor3dCurve {
    /// OCCT Adaptor3d_Curve::D0(t, P).
    pub fn d0(&self, _t: f64, _p: &mut DVec3) {
        panic!("GAP: Adaptor3d_Curve (Adaptor3d not translated) — see file header")
    }
}

/// GAP: BRepFill_DraftLaw (TKBool/BRepFill) — the draft location law over
/// the spine (BRepFill_LocationLaw subclass); not translated (plan D3).
pub struct BRepFillDraftLaw;

impl BRepFillDraftLaw {
    /// OCCT new BRepFill_DraftLaw(Spine, Law).
    pub fn new(_the_spine: &Shape, _the_law: &GeomFillLocationDraft) -> Self {
        panic!("GAP: BRepFill_DraftLaw (TKBool/BRepFill not translated) — see file header")
    }

    /// OCCT BRepFill_LocationLaw::NbLaw().
    pub fn nb_law(&self) -> usize {
        panic!("GAP: BRepFill_DraftLaw (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_LocationLaw::Edge(ii).
    pub fn edge(&self, _ii: usize) -> Shape {
        panic!("GAP: BRepFill_DraftLaw (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_LocationLaw::Vertex(ii).
    pub fn vertex(&self, _ii: usize) -> Shape {
        panic!("GAP: BRepFill_DraftLaw (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_LocationLaw::Law(ii) — the location law frame access.
    pub fn law(&self, _ii: usize) -> &GeomFillLocationDraft {
        panic!("GAP: BRepFill_DraftLaw (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_LocationLaw::CleanLaw(Tol) — clean small
    /// discontinuities.
    pub fn clean_law(&self, _tol: f64) {
        panic!("GAP: BRepFill_DraftLaw (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_LocationLaw::CurvilinearBounds(Index, f, l).
    pub fn curvilinear_bounds(&self, _index: usize, _first: &mut f64, _last: &mut f64) {
        panic!("GAP: BRepFill_DraftLaw (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_LocationLaw::Parameter(abspars, Ind, par).
    pub fn parameter(&self, _abs_pars: f64, _ind: &mut i32, _par: &mut f64) {
        panic!("GAP: BRepFill_DraftLaw (TKBool/BRepFill not translated) — see file header")
    }
    /// OCCT BRepFill_LocationLaw::Law(Ind)->GetCurve() — the adaptor curve.
    pub fn get_curve(&self, _ind: usize) -> Adaptor3dCurve {
        panic!("GAP: BRepFill_DraftLaw (TKBool/BRepFill not translated) — see file header")
    }
}

/// OCCT handle up-cast Handle(BRepFill_DraftLaw) -> Handle(BRepFill_
/// LocationLaw) — GAP (the carriers are distinct types in rcad).
fn upcast_law(_l: &BRepFillDraftLaw) -> &BRepFillLocationLaw {
    panic!("GAP: handle up-cast BRepFill_DraftLaw -> BRepFill_LocationLaw — see file header")
}

/// GAP: BndLib_Add3dCurve (TKMath/BndLib).
pub struct BndLibAdd3dCurve;

impl BndLibAdd3dCurve {
    /// OCCT BndLib_Add3dCurve::Add(C, Tol, B).
    pub fn add(_c: &Adaptor3dCurve, _tol: f64, _b: &mut BndBox) {
        panic!("GAP: BndLib_Add3dCurve (TKMath/BndLib not translated) — see file header")
    }
}

/// GAP: BndLib_AddSurface (TKMath/BndLib).
pub struct BndLibAddSurface;

impl BndLibAddSurface {
    /// OCCT BndLib_AddSurface::Add(S, Tol, B).
    pub fn add(_s: &GeomAdaptorSurface, _tol: f64, _b: &mut BndBox) {
        panic!("GAP: BndLib_AddSurface (TKMath/BndLib not translated) — see file header")
    }
    /// OCCT BndLib_AddSurface::Add(S, UMin, UMax, VMin, VMax, Tol, B).
    #[allow(clippy::too_many_arguments)]
    pub fn add_with_bounds(
        _s: &GeomAdaptorSurface,
        _u_min: f64,
        _u_max: f64,
        _v_min: f64,
        _v_max: f64,
        _tol: f64,
        _b: &mut BndBox,
    ) {
        panic!("GAP: BndLib_AddSurface (TKMath/BndLib not translated) — see file header")
    }
}

/// GAP: GeomAdaptor_Surface (TKGeomBase/GeomAdaptor).
pub struct GeomAdaptorSurface;

impl GeomAdaptorSurface {
    /// OCCT GeomAdaptor_Surface(S).
    pub fn new(_s: &Surface3) -> Self {
        panic!("GAP: GeomAdaptor_Surface (TKGeomBase not translated) — see file header")
    }
}

/// GAP: BRepAdaptor_Surface (TKTopAlgo/BRepAdaptor) — the face surface
/// adaptor consumed by the BuildShell orientation control.
pub struct BRepAdaptorSurface;

impl BRepAdaptorSurface {
    /// OCCT BRepAdaptor_Surface(F).
    pub fn new(_f: &Shape) -> Self {
        panic!("GAP: BRepAdaptor_Surface (TKTopAlgo not translated) — see file header")
    }
    /// OCCT BRepAdaptor_Surface::FirstUParameter().
    pub fn first_u_parameter(&self) -> f64 {
        panic!("GAP: BRepAdaptor_Surface (TKTopAlgo not translated) — see file header")
    }
    /// OCCT BRepAdaptor_Surface::FirstVParameter().
    pub fn first_v_parameter(&self) -> f64 {
        panic!("GAP: BRepAdaptor_Surface (TKTopAlgo not translated) — see file header")
    }
    /// OCCT BRepAdaptor_Surface::D1(U, V, P, D1U, D1V).
    pub fn d1(&self, _u: f64, _v: f64, _p: &mut DVec3, _d1u: &mut DVec3, _d1v: &mut DVec3) {
        panic!("GAP: BRepAdaptor_Surface (TKTopAlgo not translated) — see file header")
    }
}

/// GAP: GeomLProp_SLProps (TKGeomBase/GeomLProp).
pub struct GeomLPropSLProps;

impl GeomLPropSLProps {
    /// OCCT GeomLProp_SLProps(S, U, V, N, Tol).
    pub fn new(_s: &Surface3, _u: f64, _v: f64, _n: i32, _tol: f64) -> Self {
        panic!("GAP: GeomLProp_SLProps (TKGeomBase not translated) — see file header")
    }
    /// OCCT GeomLProp_SLProps::SetParameters(U, V).
    pub fn set_parameters(&mut self, _u: f64, _v: f64) {
        panic!("GAP: GeomLProp_SLProps (TKGeomBase not translated) — see file header")
    }
    /// OCCT GeomLProp_SLProps::IsNormalDefined().
    pub fn is_normal_defined(&self) -> bool {
        panic!("GAP: GeomLProp_SLProps (TKGeomBase not translated) — see file header")
    }
    /// OCCT GeomLProp_SLProps::Normal().
    pub fn normal(&self) -> DVec3 {
        panic!("GAP: GeomLProp_SLProps (TKGeomBase not translated) — see file header")
    }
}

/// GAP: BRepLib_MakeFace (TKBRep/BRepLib) — the natural-bound face of
/// BuildShell.
pub struct BRepLibMakeFace;

impl BRepLibMakeFace {
    /// OCCT BRepLib_MakeFace::Init(S, Bound, TolDegen).
    pub fn init(_surf: &Surface3, _bound: bool, _tol_degen: f64) -> Self {
        panic!("GAP: BRepLib_MakeFace (TKBRep/BRepLib not translated) — see file header")
    }
    /// OCCT BRepLib_MakeFace::Face().
    pub fn face(&self) -> Shape {
        panic!("GAP: BRepLib_MakeFace (TKBRep/BRepLib not translated) — see file header")
    }
}

/// GAP: BRepBuilderAPI_Sewing (TKShHealing/BRepBuilderAPI).
pub struct BRepBuilderAPISewing;

impl BRepBuilderAPISewing {
    /// OCCT BRepBuilderAPI_Sewing(Tol, Option1, Option2, Option3).
    pub fn new(_tol: f64, _option1: bool, _option2: bool, _option3: bool) -> Self {
        panic!("GAP: BRepBuilderAPI_Sewing (TKShHealing not translated) — see file header")
    }
    /// OCCT BRepBuilderAPI_Sewing::Add(S).
    pub fn add(&mut self, _s: &Shape) {
        panic!("GAP: BRepBuilderAPI_Sewing (TKShHealing not translated) — see file header")
    }
    /// OCCT BRepBuilderAPI_Sewing::Perform().
    pub fn perform(&mut self) {
        panic!("GAP: BRepBuilderAPI_Sewing (TKShHealing not translated) — see file header")
    }
    /// OCCT BRepBuilderAPI_Sewing::NbContigousEdges().
    pub fn nb_contigous_edges(&self) -> usize {
        panic!("GAP: BRepBuilderAPI_Sewing (TKShHealing not translated) — see file header")
    }
    /// OCCT BRepBuilderAPI_Sewing::SewedShape().
    pub fn sewed_shape(&self) -> Shape {
        panic!("GAP: BRepBuilderAPI_Sewing (TKShHealing not translated) — see file header")
    }
    /// OCCT BRepBuilderAPI_Sewing::IsModified(S).
    pub fn is_modified(&self, _s: &Shape) -> bool {
        panic!("GAP: BRepBuilderAPI_Sewing (TKShHealing not translated) — see file header")
    }
    /// OCCT BRepBuilderAPI_Sewing::Modified(S).
    pub fn modified(&self, _s: &Shape) -> Shape {
        panic!("GAP: BRepBuilderAPI_Sewing (TKShHealing not translated) — see file header")
    }
}

/// GAP: BRepExtrema_DistShapeShape (TKTopAlgo/BRepExtrema).
pub struct BRepExtremaDistShapeShape;

impl BRepExtremaDistShapeShape {
    /// OCCT BRepExtrema_DistShapeShape().
    pub fn new() -> Self {
        panic!("GAP: BRepExtrema_DistShapeShape (TKTopAlgo not translated) — see file header")
    }
    /// OCCT BRepExtrema_DistShapeShape::LoadS1(S).
    pub fn load_s1(&mut self, _s: &Shape) {
        panic!("GAP: BRepExtrema_DistShapeShape (TKTopAlgo not translated) — see file header")
    }
    /// OCCT BRepExtrema_DistShapeShape::LoadS2(S).
    pub fn load_s2(&mut self, _s: &Shape) {
        panic!("GAP: BRepExtrema_DistShapeShape (TKTopAlgo not translated) — see file header")
    }
    /// OCCT BRepExtrema_DistShapeShape::Perform().
    pub fn perform(&mut self) {
        panic!("GAP: BRepExtrema_DistShapeShape (TKTopAlgo not translated) — see file header")
    }
    /// OCCT BRepExtrema_DistShapeShape::IsDone().
    pub fn is_done(&self) -> bool {
        panic!("GAP: BRepExtrema_DistShapeShape (TKTopAlgo not translated) — see file header")
    }
    /// OCCT BRepExtrema_DistShapeShape::Value().
    pub fn value(&self) -> f64 {
        panic!("GAP: BRepExtrema_DistShapeShape (TKTopAlgo not translated) — see file header")
    }
}

impl Default for BRepExtremaDistShapeShape {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT BOPAlgo_Operation (BOPAlgo_Operation.hxx) — the CUT code consumed
/// by BuildBOP.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BOPAlgoOperation {
    BOPAlgo_CUT,
}

/// OCCT BOPAlgo_GlueEnum (BOPAlgo_GlueEnum.hxx).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BOPAlgoGlue {
    BOPAlgo_GlueShift,
}

/// GAP: BOPAlgo_PaveFiller (the public-API form; the low-level rcad
/// bop::algo::pave_filler::PaveFiller exists but has no SetArguments /
/// Perform public facade yet).
pub struct BOPAlgoPaveFiller;

impl BOPAlgoPaveFiller {
    /// OCCT BOPAlgo_PaveFiller::SetArguments(Args).
    pub fn set_arguments(&mut self, _args: &[Shape]) {
        panic!("GAP: BOPAlgo_PaveFiller public facade — see file header")
    }
    /// OCCT BOPAlgo_PaveFiller::Perform().
    pub fn perform(&mut self) {
        panic!("GAP: BOPAlgo_PaveFiller public facade — see file header")
    }
    /// OCCT BOPAlgo_Options::HasErrors().
    pub fn has_errors(&self) -> bool {
        panic!("GAP: BOPAlgo_PaveFiller public facade — see file header")
    }
}

/// GAP: BRepAlgoAPI_Section (TKBoolean/BRepAlgoAPI) — the (S1, S2, Filler)
/// constructor form.
pub struct BRepAlgoAPISection;

impl BRepAlgoAPISection {
    /// OCCT BRepAlgoAPI_Section(S1, S2, Filler).
    pub fn new(_s1: &Shape, _s2: &Shape, _filler: &BOPAlgoPaveFiller) -> Self {
        panic!("GAP: BRepAlgoAPI_Section (TKBoolean not translated) — see file header")
    }
    /// OCCT BRepAlgoAPI_BooleanOperation::Shape().
    pub fn shape(&self) -> Shape {
        panic!("GAP: BRepAlgoAPI_Section (TKBoolean not translated) — see file header")
    }
}

/// GAP: BOPAlgo_Builder (the public-API form: AddArgument /
/// PerformWithFiller / BuildBOP / History).
pub struct BOPAlgoBuilder;

impl BOPAlgoBuilder {
    /// OCCT BOPAlgo_Builder::AddArgument(S).
    pub fn add_argument(&mut self, _s: &Shape) {
        panic!("GAP: BOPAlgo_Builder public facade — see file header")
    }
    /// OCCT BOPAlgo_Builder::SetGlue(Glue).
    pub fn set_glue(&mut self, _glue: BOPAlgoGlue) {
        panic!("GAP: BOPAlgo_Builder public facade — see file header")
    }
    /// OCCT BOPAlgo_Builder::Perform().
    pub fn perform(&mut self) {
        panic!("GAP: BOPAlgo_Builder public facade — see file header")
    }
    /// OCCT BOPAlgo_Builder::PerformWithFiller(PF).
    pub fn perform_with_filler(&mut self, _pf: &BOPAlgoPaveFiller) {
        panic!("GAP: BOPAlgo_Builder public facade — see file header")
    }
    /// OCCT BOPAlgo_BOP::BuildBOP(ListObjects, ListTools, Op, Range) — the
    /// plain-operation form.
    pub fn build_bop(&mut self, _lo: &[Shape], _lt: &[Shape], _op: BOPAlgoOperation) {
        panic!("GAP: BOPAlgo_Builder::BuildBOP public facade — see file header")
    }
    /// OCCT BOPAlgo_BOP::BuildBOP(ListObjects, StateObjects, ListTools,
    /// StateTools, Range) — the state form.
    pub fn build_bop_states(
        &mut self,
        _lo: &[Shape],
        _state_objects: u8,
        _lt: &[Shape],
        _state_tools: u8,
    ) {
        panic!("GAP: BOPAlgo_Builder::BuildBOP public facade — see file header")
    }
    /// OCCT BOPAlgo_Options::HasErrors().
    pub fn has_errors(&self) -> bool {
        panic!("GAP: BOPAlgo_Builder public facade — see file header")
    }
    /// OCCT BOPAlgo_Builder::Shape().
    pub fn shape(&self) -> Shape {
        panic!("GAP: BOPAlgo_Builder public facade — see file header")
    }
    /// OCCT BOPAlgo_Builder::History().
    pub fn history(&self) -> BRepToolsHistory {
        panic!("GAP: BOPAlgo_Builder public facade — see file header")
    }
}

/// GAP: BRepTools_History (TKTopAlgo/BRepTools — the rcad
/// bop::history::BRepToolsHistory lacks Merge).
#[derive(Debug, Default)]
pub struct BRepToolsHistory;

impl BRepToolsHistory {
    /// OCCT BRepTools_History::Merge(List, Builder).
    pub fn merge_with_builder(&mut self, _the_list: &[Shape], _the_builder: &BOPAlgoBuilder) {
        panic!("GAP: BRepTools_History::Merge (the rcad bop history lacks Merge) — see file header")
    }
    /// OCCT BRepTools_History::Merge(History).
    pub fn merge(&mut self, _the_history: &BRepToolsHistory) {
        panic!("GAP: BRepTools_History::Merge (the rcad bop history lacks Merge) — see file header")
    }
    /// OCCT BRepTools_History::Clear().
    pub fn clear(&mut self) {}
    /// OCCT BRepTools_History::Modified(S).
    pub fn modified(&self, _s: &Shape) -> Vec<Shape> {
        panic!("GAP: BRepTools_History::Modified — see file header")
    }
}

// ---------------------------------------------------------------------------
// BRep_Tool / TopExp re-hosts
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Pnt(V).
fn brep_tool_pnt(v: &Shape) -> DVec3 {
    match v.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => DVec3::ZERO,
    }
}

/// OCCT BRep_Tool::Degenerated(E).
fn brep_tool_degenerated(e: &Shape) -> bool {
    match e.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT BRep_Tool::IsClosed(S) / TopoDS_Shape::Closed() — the closed flag
/// of the shape TShape.
fn brep_tool_is_closed(s: &Shape) -> bool {
    match s.data.as_ref() {
        TShape::Wire(wd) => wd.flags & tshape_flags::CLOSED != 0,
        TShape::Shell(sd) => sd.flags & tshape_flags::CLOSED != 0,
        TShape::Edge(ed) => ed.flags & tshape_flags::CLOSED != 0,
        _ => false,
    }
}

/// OCCT TopExp::MapShapesAndAncestors(S, EDGE, FACE, edgemap) — the
/// loc_ope_glued_shape reduction.
fn map_shapes_and_ancestors_edge_face(s: &Shape) -> IndexMap<(u64, u32), (Shape, Vec<Shape>)> {
    let mut the_map = IndexMap::new();
    crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
        s,
        ShapeType::Edge,
        ShapeType::Face,
        &mut the_map,
    );
    the_map
}

/// OCCT TopExp_Explorer(S, TopAbs_EDGE / FACE / SHELL / SOLID).
fn explored(s: &Shape, t: ShapeType) -> Vec<Shape> {
    crate::feat::brep_feat_builder::explorer(s, t, ShapeType::Shape)
}

// ---------------------------------------------------------------------------
// Static helpers (BRepFill_Draft.cxx)
// ---------------------------------------------------------------------------

/// OCCT gp_Dir::Angle(theOther) — the [0, PI] angle between directions
/// (pure-math helper).
fn dir_angle(d1: DVec3, d2: DVec3) -> f64 {
    let dot = d1.dot(d2).clamp(-1.0, 1.0);
    d1.cross(d2).length().atan2(dot)
}

/// OCCT gp_Vec::AngleWithRef(R, Normal) — the signed angle about the
/// reference (pure-math helper).
fn angle_with_ref(v1: &DVec3, v2: &DVec3, reference: &DVec3) -> f64 {
    let n1 = v1.length();
    let n2 = v2.length();
    if n1 <= 1.0e-15 || n2 <= 1.0e-15 {
        return 0.0;
    }
    let cross = v1.cross(*v2);
    let sign = if reference.dot(cross) < 0.0 { -1.0 } else { 1.0 };
    sign * cross.length().atan2(v1.dot(*v2))
}

/// OCCT gp_Ax3(P, V) + gp_Trsf::SetTransformation(Ax3) — the world-to-local
/// frame transformation (pure-math helper; the loc_ope_pipe.rs precedent).
fn trsf_set_transformation(origin: DVec3, dir: DVec3) -> DAffine3 {
    // gp_Ax3(P, V): a right-handed frame with the given main direction.
    let axis = dir.normalize_or_zero();
    let ref_dir = if axis.x.abs() > 1.0 - 1e-12 {
        DVec3::Z
    } else {
        DVec3::X
    };
    let x_axis = ref_dir - axis * ref_dir.dot(axis);
    let x_axis = if x_axis.length_squared() <= 1e-24 {
        rcad_kernel::geom::any_perpendicular(axis)
    } else {
        x_axis.normalize()
    };
    let y_axis = axis.cross(x_axis).normalize();
    // The frame affine (local -> world); SetTransformation gives the
    // world -> local form.
    let mut frame = DAffine3::from_mat3(DMat3::from_cols(x_axis, y_axis, axis));
    frame.translation = origin;
    frame.inverse()
}

/// OCCT BRepTools_WireExplorer(W) — the wire edges (the
/// brep_fill_pipe.rs reduction).
fn brep_tools_wire_edges(w: &Shape) -> Vec<Shape> {
    match w.data.as_ref() {
        TShape::Wire(wd) => wd.edges.clone(),
        _ => Vec::new(),
    }
}

/// OCCT BRepAdaptor_Curve() + Initialize(E) — the edge curve adaptor (the
/// ComputeTrsf / GoodOrientation GAP stand-in).
fn brep_adaptor_curve_of(_e: &Shape) -> Adaptor3dCurve {
    Adaptor3dCurve
}

/// OCCT static ComputeTrsf (L68-104): the approximate barycenter and the
/// box of the transformed wire.
fn compute_trsf(w: &Shape, d: DVec3, box_out: &mut BndBox, tf: &mut DAffine3) {
    // OCCT L71-83: the barycenter of the wire vertices.
    let exp = brep_tools_wire_edges(w);
    let mut bary = DVec3::ZERO;
    let mut nb: i32 = 0;
    for e in &exp {
        // OCCT L77-82: Bary += BRep_Tool::Pnt(Exp.CurrentVertex()).XYZ().
        let (v_f, _v_l) = super::brep_fill_pipe::top_exp_vertices(e);
        bary += brep_tool_pnt(&v_f);
        nb += 1;
    }
    bary /= nb as f64;

    // OCCT L86-87: gp_Ax3 N(Bary, D); Tf.SetTransformation(N).
    *tf = trsf_set_transformation(bary, d);

    // OCCT L92-94: the transformed wire (the location pool is BRep-bound;
    // the transformed evaluation is carried by the box walk below).
    let _ = tf;

    // OCCT L97-103: the box of the transformed wire edges (BndLib GAP).
    box_out.set_void();
    for e in &exp {
        // OCCT L100: AC.Initialize(Exp.Current()).
        let ac = brep_adaptor_curve_of(e);
        // OCCT L102: BndLib_Add3dCurve::Add(AC, 0.1, Box).
        BndLibAdd3dCurve::add(&ac, 0.1, box_out);
    }
}

/// OCCT static Longueur (L108-134): the maximum length of the rule from
/// the wire box to the stop-surface box along the direction.
fn longueur(w_box: &BndBox, s_box: &BndBox, d: &mut DVec3, p: &mut DVec3) -> f64 {
    // OCCT L115-117: WBox.Get(...) — WZmin / WZmax.
    let (_, _, wz_min, _, _, wz_max) = w_box.get().expect("WBox.Get");
    // OCCT L119-120: SBox.Get + P.SetCoord.
    let (x_min, y_min, z_min, x_max, y_max, z_max) = s_box.get().expect("SBox.Get");
    *p = DVec3::new((x_min + x_max) / 2.0, (y_min + y_max) / 2.0, z_max);

    // OCCT L122-133.
    if z_max < wz_min {
        // Skin in the wrong direction. Invert...
        *d = -*d;
        p.z = z_min;
        wz_max - z_min
    } else {
        z_max - wz_min
    }
}

/// OCCT static GoodOrientation (L140-190): check if the law is oriented to
/// have an exterior skin.
fn good_orientation(b: &BndBox, law: &BRepFillDraftLaw, d: DVec3) -> bool {
    // OCCT L147-149: the box diagonal vector.
    let (a_xmin, a_ymin, a_zmin, a_xmax, a_ymax, a_zmax) = b.get().expect("B.Get");
    let p1 = DVec3::new(a_xmin, a_ymin, a_zmin);
    let p2 = DVec3::new(a_xmax, a_ymax, a_zmax);
    let v = p2 - p1;

    // OCCT L151-152.
    let mut f = 0.0;
    let mut l = 0.0;
    let nb_law = law.nb_law();
    law.curvilinear_bounds(nb_law, &mut f, &mut l);
    let mut r = v.length() / l;

    // OCCT L156-162.
    let mut nb = (4.0 + (10.0 * r)) as i32;
    r = l / nb as f64;
    nb += 1; // Number of points

    // OCCT L164-174: the sampled points.
    let mut pnts: Vec<DVec3> = Vec::new();
    let mut bary = DVec3::ZERO;
    for ii in 1..=nb {
        let mut ind: i32 = 0;
        let mut t = 0.0;
        law.parameter((ii - 1) as f64 * r, &mut ind, &mut t);
        let ac = law.get_curve(ind as usize);
        let mut p = DVec3::ZERO;
        ac.d0(t, &mut p);
        pnts.push(p);
        bary += p;
    }

    // OCCT L176-189: the winding angle around the normal.
    bary /= nb as f64;
    let centre = bary;
    let normal = d;
    let mut angle = 0.0;
    let mut reference = pnts[0] - centre;

    for ii in 1..pnts.len() {
        let rr = pnts[ii] - centre;
        angle += angle_with_ref(&reference, &rr, &normal);
        reference = rr;
    }

    angle >= 0.0
}

/// OCCT BRepLib_MakeEdge(TC) + BRepLib_MakeWire(EG) re-host — the section
// generating line (the add_tedge / add_twire kernel forms; the vertices
// stand at the trimmed line extremities).
fn brep_lib_make_edge_wire(pool: &mut BRep, curve: Curve3, p1: DVec3, p2: DVec3, f: f64, l: f64) -> Shape {
    let v1 = pool.add_tvertex_unique(p1);
    let v2 = pool.add_tvertex_unique(p2);
    let edge = pool.add_tedge(
        Some(curve),
        crate::brep_fill::generator::shape_oriented(&v1, rcad_kernel::topods::Orientation::Forward),
        crate::brep_fill::generator::shape_oriented(&v2, rcad_kernel::topods::Orientation::Reversed),
        [f, l],
    );
    pool.add_twire(vec![edge])
}

// ---------------------------------------------------------------------------
// BRepFill_Draft (hxx L37-100)
// ---------------------------------------------------------------------------

/// OCCT BRepFill_Draft (hxx L37-100) — the draft surface engine.
pub struct BRepFillDraft {
    my_dir: DVec3,                      // OCCT: myDir (gp_Dir)
    my_angle: f64,                      // OCCT: myAngle
    angmin: f64,                        // OCCT: angmin
    angmax: f64,                        // OCCT: angmax
    my_tol: f64,                        // OCCT: myTol
    my_loc: Option<BRepFillDraftLaw>,   // OCCT: myLoc (None = null handle)
    my_sec: Option<Rc<RefCell<dyn BRepFillSectionLawOps>>>, // OCCT: mySec (None = null handle)
    my_sections: Option<ShapeArray2>,   // OCCT: mySections
    my_faces: Option<ShapeArray2>,      // OCCT: myFaces
    my_generated: Vec<Shape>,           // OCCT: myGenerated
    my_shape: Shape,                    // OCCT: myShape
    my_top: Shape,                      // OCCT: myTop
    my_shell: Shape,                    // OCCT: myShell (TopoDS_Shell)
    my_wire: Shape,                     // OCCT: myWire (TopoDS_Wire)
    my_cont: GeomAbsShape,              // OCCT: myCont
    my_style: BRepFillTransitionStyle,  // OCCT: myStyle
    is_internal: bool,                  // OCCT: IsInternal
    my_done: bool,                      // OCCT: myDone
}

impl BRepFillDraft {
    /// OCCT BRepFill_Draft::BRepFill_Draft(S, Dir, Angle) (cxx L194-283).
    pub fn new(s: &Shape, dir: DVec3, angle: f64) -> Self {
        // OCCT L196-199: myLoc/mySec/myFaces/mySections.Nullify().
        let mut r = BRepFillDraft {
            my_dir: dir,
            my_angle: 0.0,
            angmin: 0.0,
            angmax: 0.0,
            my_tol: 0.0,
            my_loc: None,
            my_sec: None,
            my_sections: None,
            my_faces: None,
            my_generated: Vec::new(),
            my_shape: Shape::null(),
            my_top: Shape::null(),
            my_shell: Shape::null(),
            my_wire: Shape::null(),
            my_cont: GeomAbsShape::AbsC1,
            my_style: BRepFillTransitionStyle::TransitionStyle_Right,
            is_internal: false,
            my_done: false,
        };

        // OCCT L201-262.
        match s.shape_type() {
            ShapeType::Wire => {
                // OCCT L204: myWire = TopoDS::Wire(S).
                r.my_wire = s.clone();
            }
            ShapeType::Face => {
                // OCCT L208-210: myWire = TopoDS::Wire(Exp.Value()).
                let exp = topods_iterator(s);
                r.my_wire = exp[0].clone();
            }
            ShapeType::Shell => {
                // OCCT L213-232: the free-border edges.
                let mut list: Vec<Shape> = Vec::new();
                let edgemap = map_shapes_and_ancestors_edge_face(s);
                for (_key, (the_edge, faces)) in edgemap.iter() {
                    let the_edge = the_edge.clone();
                    let nbf = faces.len();
                    // OCCT L224: skip degenerated edges.
                    if !brep_tool_degenerated(&the_edge) && nbf == 1 {
                        // OCCT L229: List.Append(theEdge).
                        list.push(the_edge);
                    }
                }

                // OCCT L234-258.
                if !list.is_empty() {
                    // OCCT L236-242: BRepLib_MakeWire MW; MW.Add(List);
                    // Err == BRepLib_WireDone -> myWire = MW.Wire().
                    let mut pool = BRep::new();
                    r.my_wire = pool.add_twire(list);
                } else {
                    // OCCT L256: throw Standard_ConstructionError.
                    panic!("Standard_ConstructionError: BRepFill_Draft");
                }
            }
            _ => {
                // OCCT L261: throw Standard_ConstructionError.
                panic!("Standard_ConstructionError: BRepFill_Draft");
            }
        }

        // OCCT L264-273: attention to closed non declared wires!
        if !brep_tool_is_closed(&r.my_wire) {
            let (v_f, v_l) = super::brep_fill_pipe::top_exp_vertices(&r.my_wire);
            if v_f.is_same(&v_l) {
                // OCCT L271: myWire.Closed(true) — the flag mutation is
                // pool-bound in rcad (architecture difference #5).
                let mut pool = BRep::new();
                set_closed_flag(&mut pool, &r.my_wire.clone(), true);
            }
        }

        // OCCT L275-280.
        r.my_angle = angle.abs();
        r.my_dir = dir;
        r.my_top = s.clone();
        r.my_done = false;
        r.my_tol = 1.0e-4;
        r.my_cont = GeomAbsShape::AbsC1;
        // OCCT L281-282: SetOptions(); SetDraft() — the hxx defaults.
        r.set_options(BRepFillTransitionStyle::TransitionStyle_Right, 0.01, 3.0);
        r.set_draft(false);
        r
    }

    /// OCCT BRepFill_Draft::SetOptions(Style, Min, Max) (cxx L287-294).
    pub fn set_options(&mut self, style: BRepFillTransitionStyle, min: f64, max: f64) {
        self.my_style = style;
        self.angmin = min;
        self.angmax = max;
    }

    /// OCCT BRepFill_Draft::SetDraft(Internal) (cxx L298-301).
    pub fn set_draft(&mut self, internal: bool) {
        self.is_internal = internal;
    }

    /// OCCT BRepFill_Draft::Perform(LengthMax) (cxx L307-318).
    pub fn perform(&mut self, length_max: f64) {
        // OCCT L309-312: S.Nullify(); Bnd_Box WBox; gp_Trsf Trsf.
        let mut w_box = BndBox::new();
        let mut trsf = DAffine3::IDENTITY;

        compute_trsf(&self.my_wire.clone(), self.my_dir, &mut w_box, &mut trsf);
        // OCCT L315-317: Init(S, LengthMax, WBox); BuildShell(S); Sewing().
        self.init(&None, length_max, &w_box);
        self.build_shell(&None, false);
        self.sewing();
    }

    /// OCCT BRepFill_Draft::Perform(Surface, KeepInsideSurface)
    /// (cxx L324-349).
    pub fn perform_with_surface(&mut self, surface: &Surface3, keep_inside_surface: bool) {
        let mut w_box = BndBox::new();
        let mut s_box = BndBox::new();
        let mut trsf = DAffine3::IDENTITY;
        let mut pt = DVec3::ZERO;

        compute_trsf(&self.my_wire.clone(), self.my_dir, &mut w_box, &mut trsf);

        // OCCT L334-339: the box of the transformed stop surface.
        let surf = rcad_kernel::geom::transform_surface(surface, &trsf);
        let s1 = GeomAdaptorSurface::new(&surf);
        BndLibAddSurface::add(&s1, 0.1, &mut s_box);

        // OCCT L342-343: the maximum length of the rule.
        let mut l = longueur(&w_box, &s_box, &mut self.my_dir, &mut pt);
        l /= self.my_angle.cos().abs();

        // OCCT L346-348: construction.
        self.init(&Some(surface.clone()), l, &w_box);
        self.build_shell(&Some(surface.clone()), !keep_inside_surface);
        self.sewing();
    }

    /// OCCT BRepFill_Draft::Perform(StopShape, KeepOutSide) (cxx L355-407).
    pub fn perform_with_shape(&mut self, stop_shape: &Shape, keep_out_side: bool) {
        let mut w_box = BndBox::new();
        let mut s_box = BndBox::new();
        let mut trsf = DAffine3::IDENTITY;
        let mut pt = DVec3::ZERO;

        compute_trsf(&self.my_wire.clone(), self.my_dir, &mut w_box, &mut trsf);

        // OCCT L365-389: the bounding box of the stop shape faces.
        let mut b_surf = BndBox::new();
        let mut u_min = 0.0;
        let mut u_max = 0.0;
        let mut v_min = 0.0;
        let mut v_max = 0.0;

        let ex = explored(stop_shape, ShapeType::Face);
        s_box.set_void();
        for current in &ex {
            // OCCT L379: BRepTools::UVBounds(...) — GAP (plan D3).
            brep_tools_uv_bounds(current, &mut u_min, &mut u_max, &mut v_min, &mut v_max);
            // OCCT L380-382: the transformed face surface.
            let f_surf = brep_tool_surface(current).expect("Surface");
            let surf = rcad_kernel::geom::transform_surface(&f_surf, &trsf);
            // OCCT L383-386: the box of the current face.
            let s1 = GeomAdaptorSurface::new(&surf);
            BndLibAddSurface::add_with_bounds(&s1, u_min, u_max, v_min, v_max, 0.1, &mut b_surf);
            // OCCT L387: SBox.Add(BSurf) — group boxes.
            s_box.add_box(&b_surf);
        }

        // OCCT L392-393: the maximum length of the rule.
        let mut l = longueur(&w_box, &s_box, &mut self.my_dir, &mut pt);
        l /= self.my_angle.cos().abs();

        // OCCT L396-400: the stop surface (a trimmed plane).
        let inv = trsf.inverse();
        let pt_inv = inv.transform_point3(pt);
        let plan = Plane::new(pt_inv, self.my_dir);
        let surf = Surface3::Trimmed(TrimmedSurface::new(
            Surface3::Plane(plan),
            -l,
            l,
            -l,
            l,
        ));

        // OCCT L403-406: sweeping and restriction.
        self.init(&Some(Surface3::Plane(plan)), l * 1.01, &w_box);
        self.build_shell(&Some(surf), false);
        self.fuse(stop_shape, keep_out_side);
        self.sewing();
    }

    /// OCCT BRepFill_Draft::Init(Surf, Length, Box) (cxx L411-460).
    fn init(&mut self, surf: &Option<Surface3>, length: f64, box_in: &BndBox) {
        // OCCT L416-417: the law of positioning.
        let mut loc = GeomFillLocationDraft::new(self.my_dir, self.my_angle);
        self.my_loc = Some(BRepFillDraftLaw::new(&self.my_wire.clone(), &loc));

        // OCCT L419: B = GoodOrientation(Box, myLoc, myDir).
        let b = good_orientation(box_in, self.my_loc.as_ref().expect("myLoc"), self.my_dir);

        // OCCT L421-426: if (IsInternal ^ (!B)) — invert the draft angle.
        if self.is_internal ^ !b {
            self.my_angle = -self.my_angle;
            loc.set_angle(self.my_angle);
            self.my_loc = Some(BRepFillDraftLaw::new(&self.my_wire.clone(), &loc));
        }

        // OCCT L428: myLoc->CleanLaw(angmin) — clean small discontinuities.
        self.my_loc.as_ref().expect("myLoc").clean_law(self.angmin);

        // OCCT L430-449: the law of section — the generating line is
        // straight and parallel to the binormal.
        let p = DVec3::ZERO;
        let mut d = DVec3::new(0.0, 1.0, 0.0);

        // OCCT L436-444: control of the orientation.
        let mut f = 0.0;
        let mut l = 0.0;
        let my_loc = self.my_loc.as_ref().expect("myLoc");
        my_loc.law(1).get_domain(&mut f, &mut l);
        let mut m = DMat3::IDENTITY;
        let mut bid = DVec3::ZERO;
        my_loc.law(1).d0((f + l) / 2.0, &mut m, &mut bid);
        let bn = m.col(2);

        let ang = dir_angle(self.my_dir, bn);
        if ang > std::f64::consts::PI / 2.0 {
            d = -d;
        }

        // OCCT L449-457: Geom_Line(P, D); Geom_TrimmedCurve(L, 0, Length);
        // BRepLib_MakeEdge(TC); BRepLib_MakeWire(EG).
        let line = Line3::new(p, d);
        let tc = Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3::new(
            Curve3::Line(line),
            0.0,
            length,
        ));
        let p_first = tc.point_at(0.0);
        let p_last = tc.point_at(length);
        let mut pool = BRep::new();
        let g = brep_lib_make_edge_wire(&mut pool, tc, p_first, p_last, 0.0, length);

        // OCCT L459: mySec = new BRepFill_ShapeLaw(G, true).
        self.my_sec = Some(Rc::new(RefCell::new(BRepFillShapeLaw::new(&pool, &g, true))));

        let _ = surf;
    }

    /// OCCT BRepFill_Draft::BuildShell(Surf, KeepOutSide) (cxx L466-531).
    fn build_shell(&mut self, surf: &Option<Surface3>, keep_out_side: bool) {
        // OCCT L469-471: BRepFill_Sweep Sweep(mySec, myLoc, true).  The
        // mySec member carries the real BRepFill_ShapeLaw slot; the sweep
        // construction stays on the placeholder carriers of brep_fill_pipe
        // until the owning-pool refactor + the GeomFill_LocationDraft
        // translation land (see the import note).
        let mut dummy: HashMap<ShapeKey, Shape> = HashMap::new();
        let mut dummy2: HashMap<ShapeKey, ShapeArray2> = HashMap::new();
        let mut dummy3: HashMap<ShapeKey, ShapeArray2> = HashMap::new();
        let my_sec = crate::brep_fill::brep_fill_pipe::BRepFillSectionLaw;
        let my_loc_draft = BRepFillDraftLaw::new(&Shape::null(), &GeomFillLocationDraft::new(DVec3::Y, 0.0));
        let my_loc: &BRepFillLocationLaw = upcast_law(&my_loc_draft);
        let mut sweep = BRepFillSweep::new(my_sec, my_loc, true);
        sweep.set_tolerance(self.my_tol);
        sweep.set_angular_control(self.angmin, self.angmax);
        sweep.build(
            &mut dummy,
            &mut dummy2,
            &mut dummy3,
            self.my_style,
            self.my_cont,
            GeomFillApproxStyle::GeomFill_Location,
            11,
            30,
        );
        // OCCT L482-520: the done branch — the orientation control.
        if sweep.is_done() {
            self.my_shape = sweep.shape();
            self.my_shell = self.my_shape.clone();
            self.my_faces = sweep.sub_shape();
            self.my_sections = sweep.sections();
            self.my_done = true;

            // OCCT L490-509: control of the orientation.
            let mut out = true;
            let ex = explored(&self.my_shell.clone(), ShapeType::Face);
            let f = ex[0].clone();
            let sf = BRepAdaptorSurface::new(&f);
            let u = sf.first_u_parameter();
            let v = sf.first_v_parameter();
            let mut p = DVec3::ZERO;
            let mut v1 = DVec3::ZERO;
            let mut v2 = DVec3::ZERO;
            sf.d1(u, v, &mut p, &mut v1, &mut v2);
            let mut v_n = v1.cross(v2);
            if f.orientation == rcad_kernel::topods::Orientation::Reversed {
                v_n = -v_n;
            }
            if v_n.length() > 1.0e-10 {
                out = dir_angle(self.my_dir, v_n) > std::f64::consts::PI / 2.0;
            }
            if out == self.is_internal {
                self.my_shell = shape_reversed(&self.my_shell);
                self.my_shape = shape_reversed(&self.my_shape);
            }
        } else {
            // OCCT L517-519.
            self.my_done = false;
            return;
        }

        // OCCT L522-530: add the face at end (the Fuse).
        if let Some(surf) = surf {
            let mk_f = BRepLibMakeFace::init(surf, true, TOL_CONFUSION);
            let mk_f_face = mk_f.face();
            self.fuse(&mk_f_face, keep_out_side);
        }
    }

    /// OCCT BRepFill_Draft::Fuse(S, KeepOutSide) (cxx L538-823).
    fn fuse(&mut self, stop_shape: &Shape, keep_out_side: bool) -> bool {
        let mut pool = BRep::new();
        let b = BRepBuilder::new();
        let mut issolid = false;
        let mut state1: u8 = 1; // TopAbs_OUT
        let mut state2: u8 = 1; // TopAbs_OUT

        // OCCT L545-554: Sol1 from myShape.
        let mut sol1 = Shape::null();
        if self.my_shape.shape_type() == ShapeType::Solid {
            sol1 = self.my_shape.clone();
            issolid = true;
        } else {
            // OCCT L552-553: shell => solid (for fusion).
            sol1 = pool.add_tsolid(Vec::new());
            add_to_solid(&mut pool, &b, &mut sol1, &self.my_shape.clone());
        }

        // OCCT L556-585: Sol2 from StopShape.
        let mut sol2 = Shape::null();
        match stop_shape.shape_type() {
            ShapeType::Compound => {
                // OCCT L558-561: recurse on the first part.
                let it = topods_iterator(stop_shape);
                return self.fuse(&it[0].clone(), keep_out_side);
            }
            ShapeType::Solid => {
                // OCCT L562-564.
                sol2 = stop_shape.clone();
            }
            ShapeType::Shell => {
                // OCCT L566-569: shell => solid (for fusion).
                sol2 = pool.add_tsolid(Vec::new());
                add_to_solid(&mut pool, &b, &mut sol2, &stop_shape.clone());
            }
            ShapeType::Face => {
                // OCCT L572-579.
                let mut s = pool.add_tshell(Vec::new());
                add_to_shell(&mut pool, &b, &mut s, &stop_shape.clone());
                set_closed_flag(&mut pool, &s, brep_tool_is_closed(&s));
                sol2 = pool.add_tsolid(Vec::new());
                add_to_solid(&mut pool, &b, &mut sol2, &s);
            }
            _ => {
                // OCCT L583: impossible to do.
                return false;
            }
        }

        // OCCT L588-597: the PaveFiller.
        let mut a_pf = BOPAlgoPaveFiller;
        let an_args = vec![sol1.clone(), sol2.clone()];
        a_pf.set_arguments(&an_args);
        a_pf.perform();
        if a_pf.has_errors() {
            return false;
        }

        // OCCT L599-607: the section.
        let a_sec = BRepAlgoAPISection::new(&sol1, &sol2, &a_pf);
        let a_section = a_sec.shape();

        let exp = explored(&a_section, ShapeType::Edge);
        if exp.is_empty() {
            // OCCT L604-606: no section edges produced.
            return false;
        }

        // OCCT L609-676: the state by the geometry.
        if stop_shape.shape_type() != ShapeType::Solid {
            // OCCT L614-615: the section edge closest to myWire.
            let mut a_se_min = Shape::null();
            let mut d_min = f64::MAX;
            let mut dist_tool = BRepExtremaDistShapeShape::new();
            dist_tool.load_s1(&self.my_wire.clone());

            for a_se in &exp {
                dist_tool.load_s2(a_se);
                dist_tool.perform();
                if dist_tool.is_done() {
                    let d = dist_tool.value();
                    if d < d_min {
                        d_min = d;
                        a_se_min = a_se.clone();
                        if d_min < TOL_CONFUSION {
                            break;
                        }
                    }
                }
            }

            if !a_se_min.is_null() {
                // OCCT L642-647: get geometry of StopShape — the second
                // pcurve on the section edge (the surface carrier is GAP;
                // the stand-in plane only feeds the SLProps GAP below).
                let Some((c2d, f_par, l_par)) =
                    brep_tool_curve_on_surface_index(&a_se_min, 2)
                else {
                    return false;
                };

                // OCCT L649-661: find a normal.
                let mut p2d = c2d.point_at((f_par + l_par) / 2.0);
                let mut sp = GeomLPropSLProps::new(
                    &Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Y)),
                    p2d.x,
                    p2d.y,
                    1,
                    1.0e-12,
                );
                if !sp.is_normal_defined() {
                    p2d = c2d.point_at((3.0 * f_par + l_par) / 4.0);
                    sp.set_parameters(p2d.x, p2d.y);
                    if !sp.is_normal_defined() {
                        p2d = c2d.point_at((f_par + 3.0 * l_par) / 4.0);
                        sp.set_parameters(p2d.x, p2d.y);
                    }
                }

                // OCCT L663-674: subtract State1.
                if sp.is_normal_defined() {
                    if dir_angle(self.my_dir, sp.normal()) < std::f64::consts::PI / 2.0 {
                        state1 = TOPABS_IN;
                    } else {
                        state1 = 1; // TopAbs_OUT
                    }
                }
            }
        }

        // OCCT L678-688: invert State2.
        if !keep_out_side {
            if state2 == TOPABS_IN {
                state2 = 1; // TopAbs_OUT
            } else {
                state2 = TOPABS_IN;
            }
        }

        // OCCT L691-698: the boolean operation builder.
        let mut a_builder = BOPAlgoBuilder;
        a_builder.add_argument(&sol1);
        a_builder.add_argument(&sol2);
        a_builder.perform_with_filler(&a_pf);
        if a_builder.has_errors() {
            return false;
        }

        let mut result = Shape::null();
        let mut a_history = BRepToolsHistory::default();

        let mut is_single_op_needed = true;
        // OCCT L703-769: to get rid of the unnecessary parts of the first
        // solid make the cutting first.
        if state1 == 1 {
            // TopAbs_OUT
            let a_lo = vec![sol1.clone()];
            let a_lt = vec![sol2.clone()];
            a_builder.build_bop(&a_lo, &a_lt, BOPAlgoOperation::BOPAlgo_CUT);
            if !a_builder.has_errors() {
                // OCCT L713-741: the closest cut solid.
                let mut a_cut_min = Shape::null();
                let an_exp_s = explored(&a_builder.shape(), ShapeType::Solid);
                let mut an_exp_s_iter = an_exp_s.into_iter();
                if let Some(first_solid) = an_exp_s_iter.next() {
                    a_cut_min = first_solid.clone();
                    let rest: Vec<Shape> = an_exp_s_iter.collect();
                    if !rest.is_empty() {
                        let mut a_d_min = f64::MAX;
                        let mut dist_tool = BRepExtremaDistShapeShape::new();
                        dist_tool.load_s1(&self.my_wire.clone());

                        for a_cut in rest {
                            dist_tool.load_s2(&a_cut);
                            dist_tool.perform();
                            if dist_tool.is_done() {
                                let d = dist_tool.value();
                                if d < a_d_min {
                                    a_d_min = d;
                                    a_cut_min = a_cut.clone();
                                }
                            }
                        }
                    }
                }

                if !a_cut_min.is_null() {
                    // OCCT L745-746: save history for the first argument
                    // only.
                    a_history.merge_with_builder(&a_lo, &a_builder);

                    // OCCT L748-757: the glue builder.
                    let mut a_gluer = BOPAlgoBuilder;
                    a_gluer.add_argument(&a_cut_min);
                    a_gluer.add_argument(&sol2);
                    a_gluer.set_glue(BOPAlgoGlue::BOPAlgo_GlueShift);
                    a_gluer.perform();

                    let a_lo = vec![a_cut_min.clone()];
                    a_gluer.build_bop_states(&a_lo, state1, &a_lt, state2);

                    if !a_gluer.has_errors() {
                        // OCCT L760-765.
                        let gluer_history = a_gluer.history();
                        a_history.merge(&gluer_history);

                        result = a_gluer.shape();
                        let solids = explored(&result, ShapeType::Solid);
                        is_single_op_needed = solids.is_empty();
                    }
                }
            }
        }

        // OCCT L771-787: the single operation.
        if is_single_op_needed {
            a_history.clear();

            let a_lo = vec![sol1.clone()];
            let a_lt = vec![sol2.clone()];

            a_builder.build_bop_states(&a_lo, state1, &a_lt, state2);
            if a_builder.has_errors() {
                return false;
            }

            let builder_history = a_builder.history();
            a_history.merge(&builder_history);
            result = a_builder.shape();
        }

        // OCCT L789-801.
        if issolid {
            self.my_shape = result;
        } else {
            let shells = explored(&result, ShapeType::Shell);
            if let Some(first_shell) = shells.first() {
                self.my_shape = first_shell.clone();
            }
        }

        // OCCT L803-820: update the history.
        let nb_law = self.my_loc.as_ref().expect("myLoc").nb_law();
        for ii in 1..=nb_law {
            let f = self.my_faces.as_ref().expect("myFaces").value(1, ii);
            let l = a_history.modified(&f);
            if !l.is_empty() {
                self.my_faces
                    .as_mut()
                    .expect("myFaces")
                    .set_value(1, ii, l[0].clone());
            }
        }
        for ii in 1..=(nb_law + 1) {
            let s = self.my_sections.as_ref().expect("mySections").value(1, ii);
            let l = a_history.modified(&s);
            if !l.is_empty() {
                self.my_sections
                    .as_mut()
                    .expect("mySections")
                    .set_value(1, ii, l[0].clone());
            }
        }

        true
    }

    /// OCCT BRepFill_Draft::Sewing() (cxx L829-916).
    fn sewing(&mut self) -> bool {
        let mut to_ass = self.my_top.shape_type() != ShapeType::Wire;
        let ok: bool;

        // OCCT L835-838.
        if !to_ass || !self.my_done {
            return false;
        }

        // OCCT L841-844: assembly make a shell from the faces of the shape
        // + the input shape.
        let mut ass = BRepBuilderAPISewing::new(5.0 * self.my_tol, true, true, false);
        ass.add(&self.my_shape.clone());
        ass.add(&self.my_top.clone());
        to_ass = true;
        let _ = to_ass;

        // OCCT L848-872: check if the assembly is real.
        ass.perform();
        let nb_ce = ass.nb_contigous_edges();

        let mut ok_now = false;
        if nb_ce > 0 {
            let mut res = ass.sewed_shape();
            if res.shape_type() == ShapeType::Shell || res.shape_type() == ShapeType::Solid {
                self.my_shape = res.clone();
                ok_now = true;
            } else if res.shape_type() == ShapeType::Compound {
                let it = topods_iterator(&res);
                if it.len() == 1 {
                    // OCCT L866-869: only one part => this is correct.
                    res = it[0].clone();
                    self.my_shape = res;
                    ok_now = true;
                }
            }
        }
        ok = ok_now;

        if ok {
            // OCCT L877-891: update the history.
            let nb_law = self.my_loc.as_ref().expect("myLoc").nb_law();
            for ii in 1..=nb_law {
                let f = self.my_faces.as_ref().expect("myFaces").value(1, ii);
                if ass.is_modified(&f) {
                    let m = ass.modified(&f);
                    self.my_faces
                        .as_mut()
                        .expect("myFaces")
                        .set_value(1, ii, m);
                }
            }
            for ii in 1..=(nb_law + 1) {
                let s = self.my_sections.as_ref().expect("mySections").value(1, ii);
                if ass.is_modified(&s) {
                    let m = ass.modified(&s);
                    self.my_sections
                        .as_mut()
                        .expect("mySections")
                        .set_value(1, ii, m);
                }
            }

            // OCCT L893-909: make a solid when the shape is closed.
            if brep_tool_is_closed(&self.my_shape) {
                let mut pool = BRep::new();
                let b = BRepBuilder::new();
                let mut solid = pool.add_tsolid(Vec::new());
                add_to_solid(&mut pool, &b, &mut solid, &self.my_shape.clone());

                let mut sc =
                    crate::topalgo::brep_class3d::solid_classifier::SolidClassifier::new();
                sc.load(&solid);
                sc.perform_infinite_point(TOL_CONFUSION);
                if sc.state() == TOPABS_IN {
                    solid = pool.add_tsolid(Vec::new());
                    self.my_shape = shape_reversed(&self.my_shape);
                    add_to_solid(&mut pool, &b, &mut solid, &self.my_shape.clone());
                }
                self.my_shape = solid;
            }
        }

        ok
    }

    /// OCCT BRepFill_Draft::Generated(S) (cxx L922-952).
    pub fn generated(&mut self, s: &Shape) -> Vec<Shape> {
        self.my_generated.clear();
        // OCCT L925-928: E = TopoDS::Edge(S); if (E.IsNull()) — a failed
        // cast yields the null edge (the non-edge branch walks the law
        // vertices, OCCT source as written).
        let e_is_null = s.shape_type() != ShapeType::Edge;

        if e_is_null {
            // OCCT L930-937.
            let nb_law = self.my_loc.as_ref().expect("myLoc").nb_law();
            for ii in 0..=nb_law {
                let v = self.my_loc.as_ref().expect("myLoc").vertex(ii);
                if s.is_same(&v) {
                    let sec = self.my_sections.as_ref().expect("mySections").value(1, ii + 1);
                    self.my_generated.push(sec);
                    break;
                }
            }
        } else {
            // OCCT L941-948.
            let nb_law = self.my_loc.as_ref().expect("myLoc").nb_law();
            for ii in 1..=nb_law {
                let e = self.my_loc.as_ref().expect("myLoc").edge(ii);
                if s.is_same(&e) {
                    let f = self.my_faces.as_ref().expect("myFaces").value(1, ii);
                    self.my_generated.push(f);
                    break;
                }
            }
        }

        self.my_generated.clone()
    }

    /// OCCT BRepFill_Draft::Shape() (cxx L956-959).
    pub fn shape(&self) -> Shape {
        self.my_shape.clone()
    }

    /// OCCT BRepFill_Draft::Shell() (cxx L965-968).
    pub fn shell(&self) -> Shape {
        self.my_shell.clone()
    }

    /// OCCT BRepFill_Draft::IsDone() (cxx L972-975).
    pub fn is_done(&self) -> bool {
        self.my_done
    }
}

// ---------------------------------------------------------------------------
// Small local helpers
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::CurveOnSurface(E, C, S, L, f, l, Index) — the indexed
/// form (Index = 2 selects the second pcurve representation; the surface
/// of the representation is GAP — the offset/draft stop faces carry the
/// surface in the owning pool).
fn brep_tool_curve_on_surface_index(
    edg: &Shape,
    index: usize,
) -> Option<(rcad_kernel::geom::Curve2d, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => {
            let reps: Vec<(&rcad_kernel::geom::Curve2d, f64, f64)> = ed
                .pcurves
                .values()
                .map(|(c, a, b)| (c, *a, *b))
                .collect();
            if index == 0 {
                return reps.first().map(|(c, a, b)| ((*c).clone(), *a, *b));
            }
            reps.get(index - 1).map(|(c, a, b)| ((*c).clone(), *a, *b))
        }
        _ => None,
    }
}

/// OCCT BRep_Tool::Surface(F) — the face surface.
fn brep_tool_surface(f: &Shape) -> Option<Surface3> {
    match f.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT BRepTools::UVBounds(F, UMin, UMax, VMin, VMax) — GAP (plan D3: the
/// pcurve-bounds walk is not translated).
fn brep_tools_uv_bounds(
    _f: &Shape,
    _u_min: &mut f64,
    _u_max: &mut f64,
    _v_min: &mut f64,
    _v_max: &mut f64,
) {
    panic!("GAP: BRepTools::UVBounds — see file header")
}
