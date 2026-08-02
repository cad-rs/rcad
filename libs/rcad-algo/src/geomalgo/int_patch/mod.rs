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
pub mod transitions;
pub mod int_cs;
pub mod so_on_bounds;
pub mod restriction;
pub mod elclib;
pub mod int_xx;
pub mod int_cycy;
pub mod cycy_common;
pub mod cycy_coeffs;
pub mod cycy_boundaries;
pub mod cycy_walking;
pub mod curve_surface;

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

/// Map rcad QuadricType (geomalgo::int_surf) to OCCT GeomAbs_SurfaceType.
pub fn geom_abs_of_quadric(typ: crate::geomalgo::int_surf::quadric::QuadricType) -> GeomAbsSurfaceType {
    match typ {
        crate::geomalgo::int_surf::quadric::QuadricType::Plane => GeomAbsSurfaceType::Plane,
        crate::geomalgo::int_surf::quadric::QuadricType::Cylinder => GeomAbsSurfaceType::Cylinder,
        crate::geomalgo::int_surf::quadric::QuadricType::Sphere => GeomAbsSurfaceType::Sphere,
        crate::geomalgo::int_surf::quadric::QuadricType::Cone => GeomAbsSurfaceType::Cone,
        crate::geomalgo::int_surf::quadric::QuadricType::Torus => GeomAbsSurfaceType::Torus,
        crate::geomalgo::int_surf::quadric::QuadricType::Other => GeomAbsSurfaceType::OtherSurface,
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
#[derive(Debug, Clone)]
pub struct IntPatchVertex {
    pub param_on_line: f64,
    pub p3d: DVec3,
    pub u1: f64,
    pub v1: f64,
    pub u2: f64,
    pub v2: f64,
    // ---- OCCT IntPatch_Point fields (used by PutPointsOnLine / ProcessSegments /
    // ProcessRLine in the SOnBounds post-processing) ----
    pub tolerance: f64,
    pub multiple: bool,
    pub on_dom_s1: bool,
    pub on_dom_s2: bool,
    /// 2D arc on surface 1 / surface 2 this point lies on (IntPatch_Point ArcOnS1/S2).
    pub arc_on_s1: Option<Curve2d>,
    pub arc_on_s2: Option<Curve2d>,
    pub param_on_arc1: f64,
    pub param_on_arc2: f64,
}

impl Default for IntPatchVertex {
    fn default() -> Self {
        IntPatchVertex {
            param_on_line: 0.0,
            p3d: DVec3::ZERO,
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
            tolerance: 1e-7,
            multiple: false,
            on_dom_s1: false,
            on_dom_s2: false,
            arc_on_s1: None,
            arc_on_s2: None,
            param_on_arc1: 0.0,
            param_on_arc2: 0.0,
        }
    }
}

impl IntPatchVertex {
    /// OCCT IntPatch_Point::SetValue(Pt, Tol, Tangent).
    pub fn set_value(&mut self, pt: DVec3, tol: f64, _tangent: bool) {
        self.p3d = pt;
        self.tolerance = tol;
    }
    /// OCCT IntPatch_Point::SetParameters(U1, V1, U2, V2).
    pub fn set_parameters(&mut self, u1: f64, v1: f64, u2: f64, v2: f64) {
        self.u1 = u1;
        self.v1 = v1;
        self.u2 = u2;
        self.v2 = v2;
    }
    /// OCCT IntPatch_Point::SetTolerance(Tol).
    pub fn set_tolerance(&mut self, tol: f64) {
        self.tolerance = tol;
    }
    /// OCCT IntPatch_Point::SetParameter(Para).
    pub fn set_parameter(&mut self, para: f64) {
        self.param_on_line = para;
    }
    /// OCCT IntPatch_Point::SetMultiple(IsMult).
    pub fn set_multiple(&mut self, is_mult: bool) {
        self.multiple = is_mult;
    }
    /// OCCT IntPatch_Point::IsMultiple().
    pub fn is_multiple(&self) -> bool {
        self.multiple
    }
    /// OCCT IntPatch_Point::IsOnDomS1().
    pub fn is_on_dom_s1(&self) -> bool {
        self.on_dom_s1
    }
    /// OCCT IntPatch_Point::IsOnDomS2().
    pub fn is_on_dom_s2(&self) -> bool {
        self.on_dom_s2
    }
    /// OCCT IntPatch_Point::SetArc(OnFirst, A, Param, TLine, TArc).
    /// Transitions are not stored on the rcad vertex (they are only used
    /// transiently in the restriction processing).
    pub fn set_arc(&mut self, on_first: bool, arc: Curve2d, param: f64) {
        if on_first {
            self.arc_on_s1 = Some(arc);
            self.param_on_arc1 = param;
            self.on_dom_s1 = true;
        } else {
            self.arc_on_s2 = Some(arc);
            self.param_on_arc2 = param;
            self.on_dom_s2 = true;
        }
    }
    /// OCCT IntPatch_Point::ArcOnS1().
    pub fn arc_on_s1(&self) -> Option<&Curve2d> {
        self.arc_on_s1.as_ref()
    }
    /// OCCT IntPatch_Point::ArcOnS2().
    pub fn arc_on_s2(&self) -> Option<&Curve2d> {
        self.arc_on_s2.as_ref()
    }
    /// OCCT IntPatch_Point::ParameterOnArc1().
    pub fn parameter_on_arc1(&self) -> f64 {
        self.param_on_arc1
    }
    /// OCCT IntPatch_Point::ParameterOnArc2().
    pub fn parameter_on_arc2(&self) -> f64 {
        self.param_on_arc2
    }
    /// OCCT IntPatch_Point::ParametersOnS1(U1, V1).
    pub fn parameters_on_s1(&self) -> (f64, f64) {
        (self.u1, self.v1)
    }
    /// OCCT IntPatch_Point::ParametersOnS2(U2, V2).
    pub fn parameters_on_s2(&self) -> (f64, f64) {
        (self.u2, self.v2)
    }
    /// OCCT IntPatch_Point::ParameterOnLine().
    pub fn parameter_on_line(&self) -> f64 {
        self.param_on_line
    }
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
    pub a_curve: Option<crate::geomalgo::int_patch::int_quad_quad::IntAnaCurve>,
    // ---- OCCT IntPatch_RLine fields (restriction lines) ----
    /// The 2D restriction arc on surface 1 / surface 2.
    pub arc_on_s1: Option<Curve2d>,
    pub arc_on_s2: Option<Curve2d>,
    /// Transitions of the line on surface 1 / surface 2 (IntPatch_Line
    /// TransitionOnS1/S2).  OCCT GLine/ALine store a TypeTrans (In/Out/
    /// Undecided) or Situation (Touch) transition per surface; both map to the
    /// IntSurf_Transition representation.
    pub trans1: Option<transitions::Transition>,
    pub trans2: Option<transitions::Transition>,
    /// Indices into `vertices` of the first/last point of the line.
    pub first_point: Option<usize>,
    pub last_point: Option<usize>,
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
            arc_on_s1: None,
            arc_on_s2: None,
            trans1: None,
            trans2: None,
            first_point: None,
            last_point: None,
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
            arc_on_s1: None,
            arc_on_s2: None,
            trans1: None,
            trans2: None,
            first_point: None,
            last_point: None,
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

    // ---- OCCT IntPatch_GLine / IntPatch_ALine / IntPatch_RLine vertex ops ----

    /// OCCT GLine/ALine::AddVertex(IntPatch_Point).
    pub fn add_vertex(&mut self, v: IntPatchVertex) {
        self.vertices.push(v);
    }
    /// OCCT GLine/ALine::Replace(Index, IntPatch_Point) — 1-based index.
    pub fn replace_vertex(&mut self, index: usize, v: IntPatchVertex) {
        self.vertices[index - 1] = v;
    }
    /// OCCT GLine/ALine::NbVertex().
    pub fn nb_vertex(&self) -> usize {
        self.vertices.len()
    }
    /// OCCT GLine/ALine::Vertex(Index) — 1-based.
    pub fn vertex(&self, index: usize) -> &IntPatchVertex {
        &self.vertices[index - 1]
    }

    // ---- OCCT IntPatch_RLine fields ----

    /// OCCT RLine::SetArcOnS1(A).
    pub fn set_arc_on_s1(&mut self, a: Curve2d) {
        self.arc_on_s1 = Some(a);
    }
    /// OCCT RLine::SetArcOnS2(A).
    pub fn set_arc_on_s2(&mut self, a: Curve2d) {
        self.arc_on_s2 = Some(a);
    }
    /// OCCT RLine::IsArcOnS1().
    pub fn is_arc_on_s1(&self) -> bool {
        self.arc_on_s1.is_some()
    }
    /// OCCT RLine::IsArcOnS2().
    pub fn is_arc_on_s2(&self) -> bool {
        self.arc_on_s2.is_some()
    }
    /// OCCT RLine::ArcOnS1().
    pub fn arc_on_s1(&self) -> Option<&Curve2d> {
        self.arc_on_s1.as_ref()
    }
    /// OCCT RLine::ArcOnS2().
    pub fn arc_on_s2(&self) -> Option<&Curve2d> {
        self.arc_on_s2.as_ref()
    }
    /// OCCT RLine::SetFirstPoint(Index) — 1-based index into vertices.
    pub fn set_first_point(&mut self, index: usize) {
        self.first_point = Some(index);
    }
    /// OCCT RLine::SetLastPoint(Index).
    pub fn set_last_point(&mut self, index: usize) {
        self.last_point = Some(index);
    }
    /// OCCT RLine::HasFirstPoint().
    pub fn has_first_point(&self) -> bool {
        self.first_point.is_some()
    }
    /// OCCT RLine::HasLastPoint().
    pub fn has_last_point(&self) -> bool {
        self.last_point.is_some()
    }
    /// OCCT RLine::FirstPoint().
    pub fn first_point(&self) -> &IntPatchVertex {
        &self.vertices[self.first_point.unwrap() - 1]
    }
    /// OCCT RLine::LastPoint().
    pub fn last_point(&self) -> &IntPatchVertex {
        &self.vertices[self.last_point.unwrap() - 1]
    }
}
