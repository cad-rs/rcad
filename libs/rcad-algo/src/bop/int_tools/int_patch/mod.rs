//! IntPatch (TKGeomAlgo) — intersection of surface patches.
//!
//! Rust translation of OCCT IntPatch_IType, IntPatch_Point, IntPatch_Line,
//! IntAna_QuadQuadGeo, IntPatch_ImpImpIntersection, IntPatch_Intersection.
//!
//! This module is a 1:1 translation of the OCCT analytic surface intersection
//! chain used by IntTools_FaceFace for the analytic (quadric) surface pairs:
//!   IntPatch_Intersection -> IntPatch_ImpImpIntersection -> IntAna_QuadQuadGeo
//!
//! rcad data-model notes:
//! - OCCT `Handle(IntPatch_Line)` (GLine/ALine/WLine) maps to `IntPatchLine`,
//!   a curve-carrying struct. The abstract-base transitions are omitted: they
//!   are not consumed by the rcad pipeline.
//! - OCCT `IntPatch_Point` maps to `IntPatchPoint { p1, p2, u1, v1, u2, v2 }`.

pub mod imp_imp_intersection;
pub mod int_quad_quad;
pub mod intersection;
pub mod point_line;
pub mod quad_quad_geo;
pub mod special_points;
pub mod a_line_to_w_line;
pub mod w_line_tool;

pub use imp_imp_intersection::ImpImpIntersection;
pub use intersection::IntPatchIntersection;
pub use quad_quad_geo::{AnaResultType, QuadQuadGeo};

use glam::DVec3;
use rcad_kernel::geom::{Curve2d, Curve3, Surface3};

/// OCCT GeomAbs_SurfaceType.hxx — surface type classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomAbsSurfaceType {
    Plane,
    Cylinder,
    Sphere,
    Cone,
    Torus,
    BezierSurface,
    BSplineSurface,
    SurfaceOfRevolution,
    SurfaceOfExtrusion,
    OffsetSurface,
    OtherSurface,
}

/// Map rcad QuadricType (topalgo::int_surf) to OCCT GeomAbs_SurfaceType.
pub fn geom_abs_of_quadric(typ: crate::topalgo::int_surf::quadric::QuadricType) -> GeomAbsSurfaceType {
    match typ {
        crate::topalgo::int_surf::quadric::QuadricType::Plane => GeomAbsSurfaceType::Plane,
        crate::topalgo::int_surf::quadric::QuadricType::Cylinder => GeomAbsSurfaceType::Cylinder,
        crate::topalgo::int_surf::quadric::QuadricType::Sphere => GeomAbsSurfaceType::Sphere,
        crate::topalgo::int_surf::quadric::QuadricType::Cone => GeomAbsSurfaceType::Cone,
        crate::topalgo::int_surf::quadric::QuadricType::Torus => GeomAbsSurfaceType::Torus,
        crate::topalgo::int_surf::quadric::QuadricType::Other => GeomAbsSurfaceType::OtherSurface,
    }
}

/// Convert rcad Surface3 to OCCT GeomAbs_SurfaceType.
pub fn classify_surface_type(surf: &Surface3) -> GeomAbsSurfaceType {
    match surf {
        Surface3::Plane(_) => GeomAbsSurfaceType::Plane,
        Surface3::Cylinder(_) => GeomAbsSurfaceType::Cylinder,
        Surface3::Sphere(_) => GeomAbsSurfaceType::Sphere,
        Surface3::Cone(_) => GeomAbsSurfaceType::Cone,
        Surface3::Torus(_) => GeomAbsSurfaceType::Torus,
        Surface3::BSpline(_) | Surface3::Bezier(_) => GeomAbsSurfaceType::BSplineSurface,
        _ => GeomAbsSurfaceType::OtherSurface,
    }
}

/// OCCT IntPatch_IType.hxx — type of intersection line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntPatchIType {
    Unknown,
    Line,
    Circle,
    Ellipse,
    Parabola,
    Hyperbola,
    /// OCCT IntPatch_ALine: an analytic line carrying an IntAna_Curve (the
    /// result of IntAna_IntQuadQuad for cylinder/cone pairs).  Converted to a
    /// WLine by IntPatch_ALineToWLine::MakeWLine inside GeomGeomPerfom.
    Analytic,
    Walking,
    Restriction,
}

/// IntPatch_Point — vertex on an IntPatch_Line marking a boundary
/// intersection. Stores its parameter on the line, the 3D position, and
/// the UV coordinates on both surfaces.
#[derive(Debug, Clone, Copy)]
pub struct IntPatchVertex {
    pub param_on_line: f64,
    pub p3d: DVec3,
    pub u1: f64,
    pub v1: f64,
    pub u2: f64,
    pub v2: f64,
}

/// OCCT IntPatch_Point.hxx — a point of intersection of two patches.
#[derive(Debug, Clone)]
pub struct IntPatchPoint {
    pub p1: DVec3,
    pub p2: DVec3,
    pub u1: f64,
    pub v1: f64,
    pub u2: f64,
    pub v2: f64,
    pub tolerance: f64,
}

/// OCCT IntSurf_PntOn2S.hxx — point on two surfaces.
#[derive(Debug, Clone)]
pub struct PntOn2S {
    pub p1: DVec3,
    pub p2: DVec3,
    pub u1: f64,
    pub v1: f64,
    pub u2: f64,
    pub v2: f64,
}

/// A sampled point of a walking line (IntPatch_WLine point).
#[derive(Debug, Clone, Copy)]
pub struct WLinePnt {
    pub p3d: DVec3,
    pub u1: f64,
    pub v1: f64,
    pub u2: f64,
    pub v2: f64,
}

/// Walking line kind — which sub-algorithm produced it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WLineType {
    Unknown,
    ImpImp,
    ImpPrm,
    PrmPrm,
}

/// OCCT IntPatch_Line — the abstract intersection line of two patches,
/// carrying the 3D curve and (when computed) the 2D pcurves.
#[derive(Debug, Clone)]
pub struct IntPatchLine {
    pub line_type: IntPatchIType,
    pub curve: Curve3,
    pub t_range: [f64; 2],
    pub pcurve1: Option<Curve2d>,
    pub pcurve2: Option<Curve2d>,
    pub tolerance: f64,
    pub tang_tolerance: f64,
    pub wline_pnts: Vec<WLinePnt>,
    pub is_purging_allowed: bool,
    pub wl_type: WLineType,
    /// Vertices on this line (IntPatch_Point / boundary crossings).
    pub vertices: Vec<IntPatchVertex>,
    /// OCCT IntPatch_ALine: the analytic IntAna_Curve carried by an analytic
    /// line (the result of IntAna_IntQuadQuad for cylinder/cone pairs).  It is
    /// converted to a WLine by IntPatch_ALineToWLine::MakeWLine.
    pub a_curve: Option<crate::bop::int_tools::int_patch::int_quad_quad::IntAnaCurve>,
}

impl IntPatchLine {
    pub fn analytic(lt: IntPatchIType, curve: Curve3, tr: [f64; 2]) -> Self {
        Self {
            line_type: lt,
            curve,
            t_range: tr,
            pcurve1: None,
            pcurve2: None,
            tolerance: 1e-7,
            tang_tolerance: 1e-7,
            wline_pnts: Vec::new(),
            is_purging_allowed: false,
            wl_type: WLineType::Unknown,
            vertices: Vec::new(),
            a_curve: None,
        }
    }
    pub fn walking(pnts: Vec<WLinePnt>, wt: WLineType) -> Self {
        let line = rcad_kernel::geom::Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        Self {
            line_type: IntPatchIType::Walking,
            curve: Curve3::Line(line),
            t_range: [0.0, 1.0],
            pcurve1: None,
            pcurve2: None,
            tolerance: 1e-7,
            tang_tolerance: 1e-7,
            wline_pnts: pnts,
            is_purging_allowed: true,
            wl_type: wt,
            vertices: Vec::new(),
            a_curve: None,
        }
    }
    pub fn is_wline(&self) -> bool {
        !self.wline_pnts.is_empty()
    }
    pub fn nb_points(&self) -> usize {
        self.wline_pnts.len()
    }
    pub fn point(&self, i: usize) -> &WLinePnt {
        &self.wline_pnts[i]
    }
}
