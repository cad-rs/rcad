//! OCCT BRepFill_PipeShell — split helpers and the residual GAP carriers
//! (TKBool/BRepFill sweep stack + TKGeomAlgo/TKMath dependencies not
//! translated yet, plan §0.6).
//!
//! The main file ([`super::brep_fill_pipe_shell`]) keeps the OCCT control
//! flow 1:1.  The GeomFill sweep machinery and the BRepFill law/placement
//! family are the landed translations — the consumers now carry the real
//! classes directly:
//! - GeomFill_Fixed / ConstantBiNormal / CurveAndTrihedron /
//!   GuideTrihedronAC / GuideTrihedronPlan / LocationGuide
//!   ([`crate::geomalgo::geomfill`]);
//! - BRepFill_LocationLaw / SectionLaw / Edge3DLaw / EdgeOnSurfLaw /
//!   ACRLaw / ShapeLaw / NSections / SectionPlacement / Sweep(part A)
//!   ([`crate::brep_fill`]);
//! - Law_Interpol ([`crate::geomalgo::law::law_interpol`]).
//!
//! Residual GAP carriers kept here (no real counterpart yet):
//! - BRepFill_Section (TKBool/BRepFill — separate unit);
//! - GeomFill_DiscreteTrihedron (TKGeomAlgo sweep machinery);
//! - BRepAdaptor_CompCurve (TKBRep adaptor — re-exported from the
//!   SectionPlacement unit, the same class carrier);
//! - BRepBuilderAPI_Transform / BRepBuilderAPI_Copy (TKTopAlgo);
//! - IntCurveSurface_HInter (the geomfill GAP carrier re-export);
//! - BRepGProp::LinearProperties (same gap as CompatibleWires);
//! - BRepFill_Sweep::Build (the part-B backlog of BRepFill_Sweep.cxx);
//! - BRepFill::SearchOrigin (BRepFill.cxx statics beyond Axe).
//!
//! Real translations in this split unit: the plain OCCT enums
//! GeomFill_Trihedron / BRepFill_TransitionStyle / BRepFill_TypeOfContact.
//!
//! Note: OCCT GeomFill_PipeError maps to the rcad
//! `geomfill::trihedron_law::PipeError` (variant set PipeOk / PipeNotOk /
//! PlaneNotIntersectGuide / ImpossibleContact vs the OCCT PipeOk / PipeNotOk /
//! PipeNoSolution / PipeNotPlan).

use std::collections::HashMap;

use glam::{DVec3, DAffine3};

use rcad_kernel::topo::topods::{BRep, GeomAbsShape, Shape};

use crate::brep_fill::generator::ShapeKey;
use crate::geomalgo::geomfill::trihedron_law::TrihedronLaw;

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
// GeomFill sweep machinery — the residual GAPs (TKGeomAlgo)
// ===========================================================================

/// OCCT handle<GeomFill_DiscreteTrihedron> — GAP (plan §0.6).
pub struct GeomFillDiscreteTrihedron;

impl GeomFillDiscreteTrihedron {
    /// OCCT new GeomFill_DiscreteTrihedron().
    pub fn new() -> Self {
        panic!("GAP: GeomFill_DiscreteTrihedron (TKGeomAlgo) is not translated — see file header")
    }

    /// OCCT upcast to GeomFill_TrihedronLaw.
    pub fn into_trihedron_law(self) -> Box<dyn TrihedronLaw> {
        panic!("GAP: GeomFill_DiscreteTrihedron (TKGeomAlgo) is not translated — see file header")
    }
}

impl Default for GeomFillDiscreteTrihedron {
    fn default() -> Self {
        Self::new()
    }
}

// ===========================================================================
// BRepAdaptor_CompCurve / BRepBuilderAPI_Transform / Copy / IntCurveSurface
// ===========================================================================

/// OCCT handle<BRepAdaptor_CompCurve> — GAP: the TKBRep adaptor is not
/// translated (plan §0.6); the SectionPlacement unit carries the same class
/// carrier (BRepFill_SectionPlacement.cxx L168 construction site).
pub use crate::brep_fill::brep_fill_section_placement::BRepAdaptorCompCurve;

/// OCCT handle<Adaptor3d_Curve> as consumed by IntCurveSurface_HInter —
/// the rcad guide curve is the plain [`rcad_kernel::geom::Curve3`].
pub use rcad_kernel::geom::Curve3 as IntCurveSurfaceHCurve;

/// OCCT IntCurveSurface_HInter — the geomfill GAP carrier (the rcad
/// IntCurveSurface engine exists, the production HCurveTool/HSurfaceTool
/// host markers over Curve3/Surface3 do not — see the carrier file header).
pub use crate::geomalgo::geomfill::int_curve_surface_h_inter::IntCurveSurfaceHInter;

/// OCCT BRepBuilderAPI_Transform — GAP (TKTopAlgo; the BRepBuilderAPI facade
/// unit is separate).
pub struct BRepBuilderAPITransform;

impl BRepBuilderAPITransform {
    /// OCCT BRepBuilderAPI_Transform(S, T, Copy).
    pub fn new(_s: &Shape, _t: &DAffine3, _copy: bool) -> Self {
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

/// OCCT BRepGProp::LinearProperties — GAP (same as CompatibleWires).
pub fn brep_gprop_linear_properties_centre(_profile: &Shape) -> DVec3 {
    panic!(
        "GAP: BRepGProp::LinearProperties (TKTopAlgo) is not translated — see file header"
    )
}

// ===========================================================================
// BRepFill_Sweep::Build — the part-B backlog of BRepFill_Sweep.cxx
// (NumberOfPoles / BuildFace / BuildEdge / BuildShell / Build /
//  PerformCorner / ... are the remaining translation units; the real part-A
//  class lives in brep_fill_sweep.rs).
// ===========================================================================

/// OCCT BRepFill_Sweep::Build(WDone, WInter, WSubS, Transition, Continuity,
/// Whatdegree, Degmax, Segmax) — GAP: the part-B backlog of
/// BRepFill_Sweep.cxx is not translated (see the brep_fill_sweep.rs header);
/// the call site keeps the OCCT failure path.
#[allow(clippy::too_many_arguments)]
impl crate::brep_fill::brep_fill_sweep::BRepFillSweep {
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
        panic!(
            "GAP: BRepFill_Sweep::Build (TKBool/BRepFill, part-B backlog) is not translated — \
             see file header"
        )
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
