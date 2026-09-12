//! OCCT Draft_Modification_1.cxx — the file-local statics and the pure-math
//! re-hosts / reduced-kernel carriers consumed by the Draft package
//! (Draft_Modification_1.cxx L87-112 static declarations; bodies at
//! L2047-2453).
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/Draft/
//!         Draft_Modification_1.cxx
//!
//! GAP carriers (plan §0.6; consumed by SmartParameter / Perform):
//! - GeomInt_IntSS (TKBool/GeomInt) — the surface-surface intersection with
//!   p-curve outputs; Perform is deferred (IsDone=false), the same reduction
//!   as feat/loc_ope_split_drafts_b.rs GeomIntIntSS (which lacks the
//!   HasLineOnS1/S2 + LineOnS1/S2 + TolReached3d surface consumed here).
//! - ProjLib_HCompProjectedCurve + Approx_CurveOnSurface +
//!   Adaptor3d_CurveOnSurface (TKGeomAlgo/TKTopAlgo) — the p-curve
//!   approximation tail of SmartParameter; the carriers keep the OCCT
//!   construction surface and the not-performed failure output.
//! - Extrema_ExtCS (TKGeomBase) — the curve-surface extrema fallback of
//!   Perform (cxx L1561); IsDone=false keeps the OCCT !IsDone() branch.
//!
//! Pure-math re-hosts: ElCLib::Parameter overloads (ElCLib.cxx L1192-1273),
//! Standard_Real Epsilon, gp::NormalizeAngle, gp_Dir::AngleWithRef, the
//! gp_Trsf rotation and the Geom/GP transformed variants (the OCCT
//! transformations applied through the rcad analytic types).

use glam::{DVec2, DVec3};
use super::draft_modification::ShapeIndexedMap;
use rcad_kernel::base::extrema::ExtPC2d;
use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_ext_pc::ExtremaExtPC;
use rcad_kernel::base::geom_proj_lib as geom_proj_lib;
use rcad_kernel::base::proj_lib::geom_adaptor_curve::GeomCurveAdaptor;
use rcad_kernel::base::int_ana::{
    intersect_line_plane, intersect_plane_cone_intana, intersect_plane_cylinder_intana,
    intersect_plane_plane_intana, PlnConResult, PlnCylResult, PlnPlnResult,
};
use rcad_kernel::base::geom2d_convert::{approx_curve_to_bspline, compose_curves_to_bspline};
use rcad_kernel::base::geom_api::project_on_surf::ProjectPointOnSurf;
use rcad_kernel::geom::{
    BezierCurve2, Circle3, ConicalSurface, CylindricalSurface, Curve2d, Curve2dEval, Curve3,
    CurveEval, Ellipse3, Hyperbola2d, Hyperbola3, Line2d, Line3, Parabola2d, Parabola3, Plane,
    Surface3, SurfaceEval, TrimmedCurve2, TrimmedCurve3,
};
use rcad_kernel::math::gp::{Ax1, Trsf, TrsfForm, GP_RESOLUTION};
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{GeomAbsShape, Orientation, ShapeType};

use crate::bop::int_tools::bean_face_intersector::{BRepAdaptorCurve, BRepAdaptorSurface};

use super::draft_edge_info::DraftEdgeInfo;
use super::draft_face_info::DraftFaceInfo;
use super::draft_modification::{brep_tool_continuity, brep_tool_pnt};
use super::draft_vertex_info::DraftVertexInfo;

use crate::brep_algo::tool::{brep_tool_curve, brep_tool_tolerance, explorer};
use crate::geomalgo::geom_api_project_point_on_curve::GeomAPIProjectPointOnCurve;

// ===========================================================================
// Reduced-kernel carriers (GAP annotations live in the module header).
// ===========================================================================

/// OCCT IntAna_ResultTypeInt (IntAna_ResultTypeInt.hxx) — the members consumed
/// by the Draft package.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum IntAnaTypeInter {
    Empty,
    Line,
    Circle,
    Ellipse,
    Other,
}

/// OCCT IntAna_QuadQuadGeo (IntAna_QuadQuadGeo.hxx) — the plane/plane,
/// plane/cylinder and plane/cone intersections consumed by Perform,
/// Propagate and NewSurface, carried over the rcad_kernel base::int_ana
/// vehicles.
#[derive(Debug, Clone)]
pub(crate) struct IntAnaQuadQuadGeo {
    done: bool,
    typ: IntAnaTypeInter,
    line: Option<Line3>,
    circle: Option<Circle3>,
    ellipse: Option<Ellipse3>,
}

impl Default for IntAnaQuadQuadGeo {
    fn default() -> Self {
        // OCCT IntAna_QuadQuadGeo default state: not done.
        IntAnaQuadQuadGeo {
            done: false,
            typ: IntAnaTypeInter::Other,
            line: None,
            circle: None,
            ellipse: None,
        }
    }
}

impl IntAnaQuadQuadGeo {
    /// OCCT IntAna_QuadQuadGeo(P1, P2, TolAng, TolTol) — the two-plane
    /// constructor; parallel/coincident planes give the Empty type with
    /// IsDone() true.
    pub(crate) fn new_plane_plane(p1: &Plane, p2: &Plane) -> Self {
        let mut r = IntAnaQuadQuadGeo {
            done: true,
            typ: IntAnaTypeInter::Empty,
            line: None,
            circle: None,
            ellipse: None,
        };
        if let PlnPlnResult::Line(l) = intersect_plane_plane_intana(p1, p2) {
            r.typ = IntAnaTypeInter::Line;
            r.line = Some(l);
        }
        r
    }

    /// OCCT IntAna_QuadQuadGeo::Perform(P, Cyl, TolAng, TolTol).
    pub(crate) fn perform_plane_cylinder(
        &mut self,
        the_p: &Plane,
        the_cyl: &CylindricalSurface,
        _the_tolang: f64,
        _the_toltol: f64,
    ) {
        *self = IntAnaQuadQuadGeo {
            done: true,
            typ: IntAnaTypeInter::Empty,
            line: None,
            circle: None,
            ellipse: None,
        };
        match intersect_plane_cylinder_intana(the_p, the_cyl) {
            PlnCylResult::Circle(c) => {
                self.typ = IntAnaTypeInter::Circle;
                self.circle = Some(c);
            }
            PlnCylResult::Ellipse(e) => {
                self.typ = IntAnaTypeInter::Ellipse;
                self.ellipse = Some(e);
            }
            _ => {}
        }
    }

    /// OCCT IntAna_QuadQuadGeo::Perform(P, Cone, TolAng, TolTol).
    pub(crate) fn perform_plane_cone(
        &mut self,
        the_p: &Plane,
        the_cone: &ConicalSurface,
        _the_tolang: f64,
        _the_toltol: f64,
    ) {
        *self = IntAnaQuadQuadGeo {
            done: true,
            typ: IntAnaTypeInter::Empty,
            line: None,
            circle: None,
            ellipse: None,
        };
        match intersect_plane_cone_intana(the_p, the_cone) {
            PlnConResult::Circle(c) => {
                self.typ = IntAnaTypeInter::Circle;
                self.circle = Some(c);
            }
            PlnConResult::Ellipse(e) => {
                self.typ = IntAnaTypeInter::Ellipse;
                self.ellipse = Some(e);
            }
            _ => {}
        }
    }

    /// OCCT IntAna_QuadQuadGeo::IsDone().
    pub(crate) fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT IntAna_QuadQuadGeo::TypeInter().
    pub(crate) fn type_inter(&self) -> IntAnaTypeInter {
        self.typ
    }

    /// OCCT IntAna_QuadQuadGeo::Line(N).
    pub(crate) fn line(&self, the_n: usize) -> Line3 {
        let _ = the_n;
        self.line.expect("IntAna_QuadQuadGeo::Line")
    }

    /// OCCT IntAna_QuadQuadGeo::Circle(N).
    pub(crate) fn circle(&self, the_n: usize) -> Circle3 {
        let _ = the_n;
        self.circle.expect("IntAna_QuadQuadGeo::Circle")
    }

    /// OCCT IntAna_QuadQuadGeo::Ellipse(N).
    pub(crate) fn ellipse(&self, the_n: usize) -> Ellipse3 {
        let _ = the_n;
        self.ellipse.expect("IntAna_QuadQuadGeo::Ellipse")
    }
}

/// OCCT IntAna_IntConicQuad (IntAna_IntConicQuad.hxx) — the line/plane
/// intersection consumed by InternalAdd and NewCurve, carried over the
/// rcad_kernel base::int_ana line-plane vehicle (a parallel line gives
/// IsDone() with NbPoints()==0, the OCCT empty contract).
#[derive(Debug, Clone)]
pub(crate) struct IntAnaIntConicQuad {
    done: bool,
    points: Vec<DVec3>,
}

impl IntAnaIntConicQuad {
    /// OCCT IntAna_IntConicQuad(L, Q, Tolang).
    pub(crate) fn new(the_l: &Line3, the_q: &Plane, _the_tolang: f64) -> Self {
        match intersect_line_plane(the_l, the_q) {
            Some(ip) => IntAnaIntConicQuad {
                done: true,
                points: vec![ip.point],
            },
            None => IntAnaIntConicQuad {
                done: true,
                points: Vec::new(),
            },
        }
    }

    /// OCCT IntAna_IntConicQuad::IsDone().
    pub(crate) fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT IntAna_IntConicQuad::NbPoints().
    pub(crate) fn nb_points(&self) -> usize {
        self.points.len()
    }

    /// OCCT IntAna_IntConicQuad::Point(N).
    pub(crate) fn point(&self, the_n: usize) -> DVec3 {
        self.points[the_n - 1]
    }
}

/// OCCT GeomInt_IntSS (TKBool/GeomInt) — the surface-surface intersection
/// consumed by Draft_Modification::Perform; GAP: the Perform body is deferred
/// (IsDone=false), the same reduction as the feat::loc_ope_split_drafts_b
/// GeomIntIntSS carrier, extended with the HasLineOnS1/S2 + LineOnS1/S2 +
/// TolReached3d surface Draft consumes.
#[derive(Debug, Clone)]
pub(crate) struct GeomIntIntSS {
    my_done: bool,                        // OCCT: myDone
    my_s1: Option<Surface3>,              // OCCT: myS1 (carried argument)
    my_s2: Option<Surface3>,              // OCCT: myS2 (carried argument)
    my_lines: Vec<Curve3>,                // OCCT: myLines
    my_lines_on_s1: Vec<Option<Curve2d>>, // OCCT: myLinesOnS1
    my_lines_on_s2: Vec<Option<Curve2d>>, // OCCT: myLinesOnS2
    my_tol_reached3d: f64,                // OCCT: myTolReached3d
}

impl Default for GeomIntIntSS {
    fn default() -> Self {
        Self::new()
    }
}

impl GeomIntIntSS {
    /// OCCT GeomInt_IntSS::GeomInt_IntSS().
    pub(crate) fn new() -> Self {
        GeomIntIntSS {
            my_done: false,
            my_s1: None,
            my_s2: None,
            my_lines: Vec::new(),
            my_lines_on_s1: Vec::new(),
            my_lines_on_s2: Vec::new(),
            my_tol_reached3d: 0.0,
        }
    }

    /// OCCT GeomInt_IntSS::Perform(S1, S2, Tol, Deflection, AppS1, AppS2) —
    /// GAP: deferred until the GeomInt_IntSS body lands (IsDone stays false
    /// and Perform takes the OCCT `!i2s.IsDone() || NbLines() <= 0` branch,
    /// cxx L972-977).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn perform(
        &mut self,
        the_s1: &Surface3,
        the_s2: &Surface3,
        _the_tol: f64,
        _the_deflection: bool,
        _the_app_s1: bool,
        _the_app_s2: bool,
    ) {
        self.my_s1 = Some(the_s1.clone());
        self.my_s2 = Some(the_s2.clone());
        // deferred: myDone stays false until the GeomInt_IntSS body lands.
    }

    /// OCCT GeomInt_IntSS::IsDone().
    pub(crate) fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT GeomInt_IntSS::NbLines().
    pub(crate) fn nb_lines(&self) -> usize {
        self.my_lines.len()
    }

    /// OCCT GeomInt_IntSS::Line(I).
    pub(crate) fn line(&self, the_i: usize) -> Curve3 {
        self.my_lines[the_i - 1].clone()
    }

    /// OCCT GeomInt_IntSS::HasLineOnS1(I).
    pub(crate) fn has_line_on_s1(&self, the_i: usize) -> bool {
        self.my_lines_on_s1
            .get(the_i - 1)
            .map(|c| c.is_some())
            .unwrap_or(false)
    }

    /// OCCT GeomInt_IntSS::HasLineOnS2(I).
    pub(crate) fn has_line_on_s2(&self, the_i: usize) -> bool {
        self.my_lines_on_s2
            .get(the_i - 1)
            .map(|c| c.is_some())
            .unwrap_or(false)
    }

    /// OCCT GeomInt_IntSS::LineOnS1(I).
    pub(crate) fn line_on_s1(&self, the_i: usize) -> Curve2d {
        self.my_lines_on_s1[the_i - 1]
            .clone()
            .expect("GeomInt_IntSS::LineOnS1")
    }

    /// OCCT GeomInt_IntSS::LineOnS2(I).
    pub(crate) fn line_on_s2(&self, the_i: usize) -> Curve2d {
        self.my_lines_on_s2[the_i - 1]
            .clone()
            .expect("GeomInt_IntSS::LineOnS2")
    }

    /// OCCT GeomInt_IntSS::TolReached3d().
    pub(crate) fn tol_reached_3d(&self) -> f64 {
        self.my_tol_reached3d
    }
}

/// OCCT ProjLib_HCompProjectedCurve (TKTopAlgo/ProjLib) — GAP carrier for the
/// SmartParameter tail (the compiled projection of the curve-on-surface on
/// S2): Bounds/NbIntervals carry the GAP panic until the translation lands.
pub(crate) struct ProjLibHCompProjectedCurve;

impl ProjLibHCompProjectedCurve {
    /// OCCT ProjLib_HCompProjectedCurve(S2, Cons, Tol3, Tol2).
    pub(crate) fn new(
        _the_hsur2: &BRepAdaptorSurface,
        _the_hcons: &BRepAdaptorCurve,
        _the_tol3: f64,
        _the_tol2: f64,
    ) -> Self {
        ProjLibHCompProjectedCurve
    }

    /// OCCT ProjLib_HCompProjectedCurve::Bounds(WhichFirst, Udeb, Ufin) — GAP.
    pub(crate) fn bounds(&self, _the_which: i32) -> (f64, f64) {
        panic!("GAP: ProjLib_HCompProjectedCurve::Bounds (TKTopAlgo/ProjLib not translated)");
    }

    /// OCCT ProjLib_HCompProjectedCurve::NbIntervals(S) — GAP.
    pub(crate) fn nb_intervals(&self, _the_s: GeomAbsShape) -> usize {
        panic!("GAP: ProjLib_HCompProjectedCurve::NbIntervals (TKTopAlgo/ProjLib not translated)");
    }
}

/// OCCT Approx_CurveOnSurface (TKGeomAlgo/Approx) — GAP carrier for the
/// SmartParameter tail: Perform stays not-done and Curve2d/Curve3d return the
/// OCCT null output.
pub(crate) struct ApproxCurveOnSurface;

impl ApproxCurveOnSurface {
    /// OCCT Approx_CurveOnSurface(HProjector, HSurf, Udeb, Ufin, Tol3d).
    pub(crate) fn new(
        _the_hprojector: &ProjLibHCompProjectedCurve,
        _the_hsurf: &BRepAdaptorSurface,
        _the_udeb: f64,
        _the_ufin: f64,
        _the_tol3d: f64,
    ) -> Self {
        ApproxCurveOnSurface
    }

    /// OCCT Approx_CurveOnSurface::Perform(MaxSegments, MaxDegree, Continue,
    /// Only3d, Only2d).
    pub(crate) fn perform(
        &mut self,
        _the_max_segments: usize,
        _the_max_degree: usize,
        _the_continuity: GeomAbsShape,
        _the_only3d: bool,
        _the_only2d: bool,
    ) {
        // GAP: deferred until the Approx_CurveOnSurface body lands.
    }

    /// OCCT Approx_CurveOnSurface::Curve2d() — the null output of a
    /// not-performed approximation.
    pub(crate) fn curve2d(&self) -> Option<Curve2d> {
        None
    }

    /// OCCT Approx_CurveOnSurface::Curve3d().
    pub(crate) fn curve3d(&self) -> Option<Curve3> {
        None
    }
}

/// OCCT Extrema_ExtCS (TKGeomBase/Extrema) — GAP carrier for the Perform
/// vertex fallback (cxx L1561): IsDone=false keeps the OCCT `!IsDone()`
/// VertexRecomputation branch.
pub(crate) struct ExtremaExtCS;

impl ExtremaExtCS {
    /// OCCT Extrema_ExtCS(C, S, TolC, TolS).
    pub(crate) fn new(
        _the_c: &BRepAdaptorCurve,
        _the_s: &BRepAdaptorSurface,
        _the_tolc: f64,
        _the_tols: f64,
    ) -> Self {
        ExtremaExtCS
    }

    /// OCCT Extrema_ExtCS::IsDone().
    pub(crate) fn is_done(&self) -> bool {
        false
    }

    /// OCCT Extrema_ExtCS::NbExt().
    pub(crate) fn nb_ext(&self) -> usize {
        0
    }
}

/// OCCT GeomConvert_CompCurveToBSplineCurve (TKTopAlgo/GeomConvert) — GAP
/// carrier for the Perform BSpline glueing (cxx L1281-1290): the
/// constructor stores the first BSpline; Add is deferred until the
/// GeomConvert batch lands (the branch is unreachable while the
/// GeomInt_IntSS body is deferred, since the glueing runs only on i2s
/// BSpline results).
pub(crate) struct GeomConvertCompCurveToBSplineCurve {
    my_bspline: rcad_kernel::geom::BSplineCurve3, // the stored first curve
}

impl GeomConvertCompCurveToBSplineCurve {
    /// OCCT GeomConvert_CompCurveToBSplineCurve(BSpline).
    pub(crate) fn new(the_bspline: &rcad_kernel::geom::BSplineCurve3) -> Self {
        GeomConvertCompCurveToBSplineCurve {
            my_bspline: the_bspline.clone(),
        }
    }

    /// OCCT GeomConvert_CompCurveToBSplineCurve::Add(BSpline, Tol, After) —
    /// GAP.
    pub(crate) fn add(
        &mut self,
        _the_bspline: &rcad_kernel::geom::BSplineCurve3,
        _the_tol: f64,
        _the_after: bool,
    ) {
        panic!(
            "GAP: GeomConvert_CompCurveToBSplineCurve::Add (TKTopAlgo/GeomConvert not translated)"
        );
    }

    /// OCCT GeomConvert_CompCurveToBSplineCurve::BSplineCurve().
    pub(crate) fn bspline_curve(&self) -> rcad_kernel::geom::BSplineCurve3 {
        self.my_bspline.clone()
    }
}

// ===========================================================================
// Pure-math re-hosts.
// ===========================================================================

/// OCCT Standard_Real Epsilon(theValue) (Standard_Real.hxx L238-247) — the
/// distance to the nearest representable double.
pub(crate) fn standard_real_epsilon(the_value: f64) -> f64 {
    if the_value == 0.0 {
        f64::EPSILON
    } else {
        f64::EPSILON * the_value.abs()
    }
}

/// OCCT gp::NormalizeAngle(A) (gp.hxx) — the single-step normalization to
/// (-PI, PI].
pub(crate) fn normalize_angle(the_a: f64) -> f64 {
    if the_a > std::f64::consts::PI {
        the_a - 2.0 * std::f64::consts::PI
    } else if the_a < -std::f64::consts::PI {
        the_a + 2.0 * std::f64::consts::PI
    } else {
        the_a
    }
}

/// OCCT gp_Dir::AngleWithRef(V, Ref) (gp_Dir.cxx) — the signed angle between
/// `the_xdir` and `the_v` in the reference direction `the_ref`
/// (pure-math re-host).
pub(crate) fn dir_angle_with_ref(the_xdir: DVec3, the_v: DVec3, the_ref: DVec3) -> f64 {
    let c = the_xdir.dot(the_v);
    let s = the_ref.dot(the_xdir.cross(the_v));
    normalize_angle(c.atan2(s))
}

/// OCCT ElCLib::LineParameter(L, P) (ElCLib.cxx L1192-1196).
pub(crate) fn elclib_line_parameter(the_l: &Line3, the_p: DVec3) -> f64 {
    (the_p - the_l.origin).dot(the_l.direction)
}

/// OCCT ElCLib::CircleParameter(Pos, P) (ElCLib.cxx L1199-1224).
pub(crate) fn elclib_circle_parameter(the_c: &Circle3, the_p: DVec3) -> f64 {
    let a_vec = the_p - the_c.center;
    if a_vec.length_squared() < GP_RESOLUTION {
        // coinciding points -> infinite number of parameters
        return 0.0;
    }
    // OCCT: aVProj = dir.CrossCrossed(aVec, dir) = dir ^ (aVec ^ dir).
    let a_v_proj = the_c.normal.cross(the_c.normal.cross(a_vec));
    if a_v_proj.length_squared() < GP_RESOLUTION {
        return 0.0;
    }
    // OCCT: Teta = XDirection.AngleWithRef(aVProj, dir); normalizeAngle.
    dir_angle_with_ref(the_c.x_dir, a_v_proj, the_c.normal)
}

/// OCCT ElCLib::EllipseParameter(Pos, MajorRadius, MinorRadius, P)
/// (ElCLib.cxx L1226-1251).
pub(crate) fn elclib_ellipse_parameter(the_e: &Ellipse3, the_p: DVec3) -> f64 {
    let op = the_p - the_e.center;
    let xaxis = the_e.major_dir;
    // OCCT: yaxis = YDirection — the right-handed N ^ X.
    let yaxis = the_e.normal.cross(xaxis);
    let ny = op.dot(yaxis);
    let nx = op.dot(xaxis);
    if nx.abs() <= GP_RESOLUTION && ny.abs() <= GP_RESOLUTION {
        // The point P is on the axis of the ellipse.
        return 0.0;
    }
    let om = xaxis * nx + yaxis * (ny * (the_e.major_radius / the_e.minor_radius));
    dir_angle_with_ref(xaxis, om, the_e.normal)
}

/// OCCT ElCLib::HyperbolaParameter(Pos, MajorRadius, MinorRadius, P)
/// (ElCLib.cxx L1253-1267).
pub(crate) fn elclib_hyperbola_parameter(the_h: &Hyperbola3, the_p: DVec3) -> f64 {
    // OCCT: sht = (P - Loc).Dot(YDirection) / MinorRadius.
    let ydir = the_h.normal.cross(the_h.major_dir);
    let sht = (the_p - the_h.center).dot(ydir) / the_h.semi_minor;
    sht.asinh()
}

/// OCCT ElCLib::ParabolaParameter(Pos, P) (ElCLib.cxx L1269-1273).
pub(crate) fn elclib_parabola_parameter(the_pb: &Parabola3, the_p: DVec3) -> f64 {
    // OCCT: (P - Loc).Dot(YDirection) — the rcad parabola frame carries the
    // focal direction (the OCCT XDirection); Y = N ^ X.
    let ydir = the_pb.normal.cross(the_pb.axis_dir);
    (the_p - the_pb.vertex).dot(ydir)
}

/// OCCT gp_Trsf::SetRotation(A1, Ang) (gp_Trsf.cxx L120-181 via gp_Mat) —
/// the rotation matrix about `the_axe` by `the_angle` (pure-math re-host).
pub(crate) fn gp_trsf_rotation(the_axe: &Ax1, the_angle: f64) -> Trsf {
    let (x, y, z) = (the_axe.direction.x, the_axe.direction.y, the_axe.direction.z);
    let c = the_angle.cos();
    let s = the_angle.sin();
    let t = 1.0 - c;
    Trsf {
        matrix: [
            [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
            [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
            [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
        ],
        loc: DVec3::ZERO,
        scale: 1.0,
        form: TrsfForm::Rotation,
    }
}

/// OCCT Geom_Plane / gp_Pln::Rotated(A1, Ang) (Geom_Plane.cxx via gp_Trsf) —
/// every member of the rcad Plane frame transforms like the OCCT gp_Ax3.
pub(crate) fn geom_plane_rotated(the_pl: &Plane, the_axe: &Ax1, the_theta: f64) -> Plane {
    let t = gp_trsf_rotation(the_axe, the_theta);
    Plane {
        origin: t.apply(the_pl.origin),
        normal: t.transform_dir(the_pl.normal),
        u_dir: t.transform_dir(the_pl.u_dir),
        v_dir: t.transform_dir(the_pl.v_dir),
    }
}

/// OCCT Geom_Curve::Rotated(A1, Ang) (Geom_Geometry.cxx) — the curve
/// transform through the rotation gp_Trsf, applied per analytic member
/// exactly as the OCCT geometry classes do.
pub(crate) fn geom_curve_rotated(the_c: &Curve3, the_axe: &Ax1, the_theta: f64) -> Curve3 {
    let t = gp_trsf_rotation(the_axe, the_theta);
    let pt = |p: DVec3| t.apply(p);
    let dir = |v: DVec3| t.transform_dir(v);
    match the_c {
        Curve3::Line(l) => Curve3::Line(Line3 {
            origin: pt(l.origin),
            direction: dir(l.direction),
        }),
        Curve3::Circle(c) => Curve3::Circle(Circle3 {
            center: pt(c.center),
            normal: dir(c.normal),
            x_dir: dir(c.x_dir),
            y_dir: dir(c.y_dir),
            radius: c.radius,
        }),
        Curve3::Ellipse(e) => Curve3::Ellipse(Ellipse3 {
            center: pt(e.center),
            normal: dir(e.normal),
            major_dir: dir(e.major_dir),
            major_radius: e.major_radius,
            minor_radius: e.minor_radius,
        }),
        Curve3::Hyperbola(h) => Curve3::Hyperbola(Hyperbola3 {
            center: pt(h.center),
            normal: dir(h.normal),
            major_dir: dir(h.major_dir),
            semi_major: h.semi_major,
            semi_minor: h.semi_minor,
        }),
        Curve3::Parabola(p) => Curve3::Parabola(Parabola3 {
            vertex: pt(p.vertex),
            normal: dir(p.normal),
            axis_dir: dir(p.axis_dir),
            focal_param: p.focal_param,
        }),
        Curve3::BSpline(b) => {
            let mut b = b.clone();
            b.control_points = b.control_points.iter().map(|p| pt(*p)).collect();
            Curve3::BSpline(b)
        }
        Curve3::Bezier(b) => {
            let mut b = b.clone();
            b.control_points = b.control_points.iter().map(|p| pt(*p)).collect();
            Curve3::Bezier(b)
        }
        Curve3::Trimmed(tr) => Curve3::Trimmed(TrimmedCurve3 {
            curve: Box::new(geom_curve_rotated(&tr.curve, the_axe, the_theta)),
            first: tr.first,
            last: tr.last,
        }),
        Curve3::Offset(o) => {
            let mut o = o.clone();
            o.basis = Box::new(geom_curve_rotated(&o.basis, the_axe, the_theta));
            o.offset_dir = dir(o.offset_dir);
            Curve3::Offset(o)
        }
        Curve3::CircularHelix(h) => {
            let mut h = h.clone();
            h.origin = pt(h.origin);
            h.axis = dir(h.axis);
            h.ref_dir = dir(h.ref_dir);
            Curve3::CircularHelix(h)
        }
        Curve3::SineWave(w) => {
            let mut w = w.clone();
            w.origin = pt(w.origin);
            w.baseline_dir = dir(w.baseline_dir);
            w.amplitude_dir = dir(w.amplitude_dir);
            Curve3::SineWave(w)
        }
    }
}

/// OCCT gp_Circ::Translate(V) — the center translation (gp_Circ.cxx).
pub(crate) fn gp_circ_translate(the_c: &Circle3, the_v: DVec3) -> Circle3 {
    let mut c = *the_c;
    c.center = the_c.center + the_v;
    c
}

/// OCCT Geom_Curve::Reverse() — the reversal via SetTrim(RevP(Last),
/// RevP(First)) (Geom_TrimmedCurve.cxx): the basis is reversed when the new
/// range comes in decreasing order.
pub(crate) fn geom_curve_reverse(the_c: &Curve3) -> Curve3 {
    match the_c {
        Curve3::Line(l) => Curve3::Line(Line3 {
            origin: l.origin,
            direction: -l.direction,
        }),
        Curve3::Circle(c) => Curve3::Circle(Circle3 {
            center: c.center,
            normal: -c.normal,
            x_dir: c.x_dir,
            y_dir: c.y_dir,
            radius: c.radius,
        }),
        Curve3::Ellipse(e) => Curve3::Ellipse(Ellipse3 {
            center: e.center,
            normal: -e.normal,
            major_dir: e.major_dir,
            major_radius: e.major_radius,
            minor_radius: e.minor_radius,
        }),
        Curve3::Hyperbola(h) => Curve3::Hyperbola(Hyperbola3 {
            center: h.center,
            normal: -h.normal,
            major_dir: h.major_dir,
            semi_major: h.semi_major,
            semi_minor: h.semi_minor,
        }),
        Curve3::Parabola(p) => Curve3::Parabola(Parabola3 {
            vertex: p.vertex,
            normal: -p.normal,
            axis_dir: p.axis_dir,
            focal_param: p.focal_param,
        }),
        Curve3::BSpline(b) => Curve3::BSpline(b.reversed()),
        Curve3::Bezier(b) => {
            let mut b = b.clone();
            b.control_points.reverse();
            Curve3::Bezier(b)
        }
        Curve3::Trimmed(t) => Curve3::Trimmed(TrimmedCurve3 {
            curve: Box::new(geom_curve_reverse(&t.curve)),
            first: CurveEval::reversed_parameter(t.curve.as_ref(), t.last),
            last: CurveEval::reversed_parameter(t.curve.as_ref(), t.first),
        }),
        Curve3::Offset(_) | Curve3::CircularHelix(_) | Curve3::SineWave(_) => {
            // The swept/parametric members never appear in the Draft flows
            // (plane/cylinder/cone faces); the OCCT classes have the same
            // reversal, the rcad analytic carriers are not needed yet.
            panic!("geom_curve_reverse: unsupported curve member in the Draft flows");
        }
    }
}

/// OCCT Geom2d_Curve::Reverse() — the 2D reversal (Geom2d_Curve.cxx; the
/// conic reversals negate the Y direction per the gp_Ax22d frame).
pub(crate) fn geom_curve2d_reverse(the_c: &Curve2d) -> Curve2d {
    match the_c {
        Curve2d::Line(l) => Curve2d::Line(Line2d {
            origin: l.origin,
            direction: -l.direction,
        }),
        Curve2d::Circle(c) => Curve2d::Circle(rcad_kernel::geom::Circle2d {
            center: c.center,
            x_dir: c.x_dir,
            y_dir: -c.y_dir,
            radius: c.radius,
        }),
        Curve2d::Ellipse(e) => Curve2d::Ellipse(rcad_kernel::geom::Ellipse2d {
            center: e.center,
            major_dir: e.major_dir,
            minor_dir: -e.minor_dir,
            major_radius: e.major_radius,
            minor_radius: e.minor_radius,
        }),
        Curve2d::Parabola(p) => Curve2d::Parabola(Parabola2d {
            origin: p.origin,
            axis_dir: -p.axis_dir,
            focal_param: p.focal_param,
        }),
        Curve2d::Hyperbola(h) => Curve2d::Hyperbola(Hyperbola2d {
            center: h.center,
            major_dir: -h.major_dir,
            semi_major: h.semi_major,
            semi_minor: h.semi_minor,
        }),
        Curve2d::BSpline(b) => Curve2d::BSpline(bspline2_reversed(b)),
        Curve2d::Bezier(b) => {
            let mut b = b.clone();
            b.control_points.reverse();
            Curve2d::Bezier(b)
        }
        Curve2d::Trimmed(t) => Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(geom_curve2d_reverse(&t.curve)),
            t_min: Curve2dEval::reversed_parameter(t.curve.as_ref(), t.t_max),
            t_max: Curve2dEval::reversed_parameter(t.curve.as_ref(), t.t_min),
        }),
        _ => {
            // See geom_curve_reverse: never produced in the Draft flows.
            panic!("geom_curve2d_reverse: unsupported curve member in the Draft flows");
        }
    }
}

/// OCCT Geom2d_BSplineCurve::Reverse — the knot reflection + pole/weight
/// reversion (the bspline_ops::reduced vehicle for the 2D curve).
fn bspline2_reversed(b: &rcad_kernel::geom::BSplineCurve2) -> rcad_kernel::geom::BSplineCurve2 {
    let kfirst = b.knots[0];
    let klast = b.knots[b.knots.len() - 1];
    let knots: Vec<f64> = b.knots.iter().rev().map(|k| kfirst + klast - k).collect();
    let control_points: Vec<DVec2> = b.control_points.iter().rev().copied().collect();
    let mut weights = b.weights.clone();
    weights.reverse();
    rcad_kernel::geom::BSplineCurve2 {
        degree: b.degree,
        knots,
        control_points,
        weights,
    }
}

/// OCCT Geom2dAPI_ProjectPointOnCurve (Geom2dAPI_ProjectPointOnCurve.hxx) —
/// the nearest projection of a point on a 2D curve, carried over
/// Extrema_ExtPC2d.
pub(crate) struct Geom2dAPIProjectPointOnCurve {
    my_nb_points: usize,
    my_lower_distance: f64,
    my_nearest: DVec2,
}

impl Geom2dAPIProjectPointOnCurve {
    /// OCCT Geom2dAPI_ProjectPointOnCurve(P, C).
    pub(crate) fn new(the_p: DVec2, the_c: &Curve2d) -> Self {
        let dom = the_c.default_domain();
        let ext = ExtPC2d::new(the_p, the_c, CONFUSION, dom[0], dom[1]);
        if ext.is_done() && ext.nb_ext() >= 1 {
            let (mut best, mut best_d) = (1usize, ext.square_distance(1));
            for i in 2..=ext.nb_ext() {
                let d = ext.square_distance(i);
                if d < best_d {
                    best_d = d;
                    best = i;
                }
            }
            Geom2dAPIProjectPointOnCurve {
                my_nb_points: ext.nb_ext(),
                my_lower_distance: best_d.sqrt(),
                my_nearest: ext.point(best).point,
            }
        } else {
            Geom2dAPIProjectPointOnCurve {
                my_nb_points: 0,
                my_lower_distance: f64::MAX,
                my_nearest: DVec2::ZERO,
            }
        }
    }

    /// OCCT Geom2dAPI_ProjectPointOnCurve::NbPoints().
    pub(crate) fn nb_points(&self) -> usize {
        self.my_nb_points
    }

    /// OCCT Geom2dAPI_ProjectPointOnCurve::LowerDistance().
    pub(crate) fn lower_distance(&self) -> f64 {
        self.my_lower_distance
    }

    /// OCCT Geom2dAPI_ProjectPointOnCurve::NearestPoint().
    #[allow(dead_code)]
    pub(crate) fn nearest_point(&self) -> DVec2 {
        self.my_nearest
    }
}

// ===========================================================================
// Draft_Modification_1.cxx file-local statics.
// ===========================================================================

/// OCCT static Orientation(S, F) (Draft_Modification_1.cxx L2383-2399) — the
/// cumulated orientation of F inside S.
pub(crate) fn orientation(the_s: &Shape, the_f: &Shape) -> Orientation {
    // OCCT L2388-2389: TopExp_Explorer expl; expl.Init(S, TopAbs_FACE);
    for a_cur in explorer(the_s, ShapeType::Face, ShapeType::Shape) {
        if a_cur.is_same(the_f) {
            // OCCT L2392-2394: return expl.Current().Orientation();
            return a_cur.orientation;
        }
    }
    // OCCT L2398: return TopAbs_FORWARD;
    Orientation::Forward
}

/// OCCT static FindRotation(Pl, Oris, Direction, Angle, NeutralPlane, Axe,
/// theta) (Draft_Modification_1.cxx L2403-2453).
pub(crate) fn find_rotation(
    the_pl: &Plane,
    the_oris: Orientation,
    the_direction: DVec3,
    the_angle: f64,
    the_neutral_plane: &Plane,
    the_axe: &mut Ax1,
    the_theta: &mut f64,
) -> bool {
    // OCCT L2411: IntAna_QuadQuadGeo i2pl(Pl, NeutralPlane, Angular, Confusion);
    let i2pl = IntAnaQuadQuadGeo::new_plane_plane(the_pl, the_neutral_plane);

    // OCCT L2413: if (i2pl.IsDone() && i2pl.TypeInter() == IntAna_Line).
    if i2pl.is_done() && i2pl.type_inter() == IntAnaTypeInter::Line {
        // OCCT L2415: gp_Lin li = i2pl.Line(1);
        let li = i2pl.line(1);
        // OCCT L2417-2418: nx = li.Direction(); ny = Pl.Axis().Direction().Crossed(nx);
        let nx = li.direction;
        let ny = the_pl.normal.cross(nx);
        // OCCT L2419: double a = Direction.Dot(nx);
        let a = the_direction.dot(nx);
        // OCCT L2420: if (std::abs(a) <= 1 - Precision::Angular()).
        if a.abs() <= 1.0 - rcad_kernel::precision::ANGULAR {
            // OCCT L2422-2423.
            let mut b = the_direction.dot(ny);
            let mut c = the_direction.dot(the_pl.normal);
            // OCCT L2424: direct = Pl.Position().Direct() — the sense of the
            // rcad plane frame is (U ^ V) . N > 0.
            let direct = the_pl.u_dir.cross(the_pl.v_dir).dot(the_pl.normal) > 0.0;
            // OCCT L2425-2429.
            if (direct && the_oris == Orientation::Reversed)
                || (!direct && the_oris == Orientation::Forward)
            {
                b = -b;
                c = -c;
            }
            // OCCT L2430-2432.
            let denom = (1.0 - a * a).sqrt();
            let sina = the_angle.sin();
            if denom > sina.abs() {
                // OCCT L2434-2436.
                let phi = (b / denom).atan2(c / denom);
                let theta0 = (sina / denom).acos();
                let mut theta = theta0 - phi;
                // OCCT L2437-2440.
                if theta.cos() < 0.0 {
                    theta = -theta0 - phi;
                }
                // OCCT L2442-2446: while (std::abs(theta) > M_PI).
                while theta.abs() > std::f64::consts::PI {
                    theta += std::f64::consts::PI * if theta < 0.0 { 1.0 } else { -1.0 };
                }
                // OCCT L2447-2448: Axe = li.Position(); return true;
                *the_axe = Ax1::new(li.origin, li.direction);
                *the_theta = theta;
                return true;
            }
        }
    }
    // OCCT L2452: return false;
    false
}

/// OCCT static Parameter(C, P, done) (Draft_Modification_1.cxx L2180-2270).
pub(crate) fn parameter(the_c: &Curve3, the_p: DVec3, done: &mut i32) -> f64 {
    *done = 0;
    // OCCT L2183-2189: cbase = C; if (ctyp == Geom_TrimmedCurve)
    // { cbase = BasisCurve(); ctyp = cbase->DynamicType(); }
    let cbase: &Curve3 = match the_c {
        Curve3::Trimmed(t) => t.curve.as_ref(),
        _ => the_c,
    };
    let mut param;
    match cbase {
        // OCCT L2191-2194.
        Curve3::Line(l) => {
            param = elclib_line_parameter(l, the_p);
        }
        // OCCT L2195-2200.
        Curve3::Circle(c) => {
            let mut p = elclib_circle_parameter(c, the_p);
            let two_pi = 2.0 * std::f64::consts::PI;
            // OCCT: if (std::abs(2.*M_PI - param) <= Epsilon(2.*M_PI)) param = 0.;
            if (two_pi - p).abs() <= standard_real_epsilon(two_pi) {
                p = 0.0;
            }
            param = p;
        }
        // OCCT L2203-2210.
        Curve3::Ellipse(e) => {
            let mut p = elclib_ellipse_parameter(e, the_p);
            let two_pi = 2.0 * std::f64::consts::PI;
            if (two_pi - p).abs() <= standard_real_epsilon(two_pi) {
                p = 0.0;
            }
            param = p;
        }
        // OCCT L2211-2214.
        Curve3::Parabola(pb) => {
            param = elclib_parabola_parameter(pb, the_p);
        }
        // OCCT L2215-2218.
        Curve3::Hyperbola(h) => {
            param = elclib_hyperbola_parameter(h, the_p);
        }
        // OCCT L2219-2268: the generic branch.
        _ => {
            // OCCT L2221-2222: GeomAdaptor_Curve TheCurve(C);
            // Extrema_ExtPC myExtPC(P, TheCurve) — the two-arg ctor over the
            // full domain, the default theTolF is 1.0e-10.
            let dom = the_c.default_domain();
            let a_adaptor = GeomCurveAdaptor::new(the_c.clone());
            let a_tool = CurveToolHandle::for_curve3(the_c, &a_adaptor, &a_adaptor);
            let my_ext_pc = ExtremaExtPC::new_point_curve(the_p, &a_tool, 1.0e-10);
            // OCCT L2223-2226.
            if !my_ext_pc.is_done() {
                panic!("Standard_Failure: Draft_Modification_1::Parameter: ExtremaPC not done.");
            }
            if my_ext_pc.nb_ext() >= 1 {
                // OCCT L2227-2240.
                let mut dist2_min = my_ext_pc.square_distance(1);
                let mut jmin = 1usize;
                for j in 2..=my_ext_pc.nb_ext() {
                    let dist2 = my_ext_pc.square_distance(j);
                    if dist2 < dist2_min {
                        dist2_min = dist2;
                        jmin = j;
                    }
                }
                param = my_ext_pc.point(jmin).param;
            } else {
                // OCCT L2242-2257: myExtPC.TrimmedSquareDistances(dist1_2,
                // dist2_2, p1b, p2b).
                let (dist1_2, dist2_2, _p1b, _p2b) = my_ext_pc.trimmed_square_distances();
                if dist1_2 < dist2_2 {
                    *done = -1;
                    param = dom[0];
                } else {
                    *done = 1;
                    param = dom[1];
                }
            }

            // OCCT L2259-2267.
            if cbase.is_periodic() {
                let per = cbase.default_domain()[1] - cbase.default_domain()[0];
                // OCCT: Precision::Parametric(P) = Parametric(P, Confusion()).
                let tolp = rcad_kernel::precision::parametric_default(CONFUSION);
                if (per - param).abs() <= tolp {
                    param = 0.0;
                }
            }
        }
    }
    // OCCT L2269: return param;
    param
}

/// OCCT static SmartParameter(Einf, EdgeTol, Pnt, sign, S1, S2)
/// (Draft_Modification_1.cxx L2274-2379).
pub(crate) fn smart_parameter(
    the_einf: &mut DraftEdgeInfo,
    the_edge_tol: f64,
    the_pnt: DVec3,
    the_sign: i32,
    the_s1: &Surface3,
    the_s2: &Surface3,
) -> f64 {
    let tol = CONFUSION; // OCCT L2282: constexpr double Tol = Precision::Confusion();
    let _ = the_edge_tol; // OCCT L2283: double Etol = EdgeTol (kept for form).

    // OCCT L2285-2297.
    if the_einf.first_pc().is_none() {
        let the_curve = the_einf
            .geometry()
            .cloned()
            .expect("SmartParameter: null Geometry");
        let cdom = the_curve.default_domain();
        let sdom = the_s1.default_domain();
        let pcu1 = geom_proj_lib::curve2d(
            &the_curve, cdom[0], cdom[1], the_s1, sdom[0], sdom[1], sdom[2], sdom[3],
        );
        *the_einf.change_first_pc() = pcu1;
    }
    // OCCT L2298-2307.
    if the_einf.second_pc().is_none() {
        let the_curve = the_einf
            .geometry()
            .cloned()
            .expect("SmartParameter: null Geometry");
        let cdom = the_curve.default_domain();
        let sdom = the_s2.default_domain();
        let pcu2 = geom_proj_lib::curve2d(
            &the_curve, cdom[0], cdom[1], the_s2, sdom[0], sdom[1], sdom[2], sdom[3],
        );
        *the_einf.change_second_pc() = pcu2;
    }

    // OCCT L2309-2311: GeomAPI_ProjectPointOnSurf Projector(Pnt, S1);
    // LowerDistanceParameters(U, V).
    let projector = ProjectPointOnSurf::new_point(the_pnt, the_s1);
    let (u, v) = projector.lower_distance_parameters();

    // OCCT L2313-2317.
    let mut new_c2d = the_einf.first_pc().cloned().expect("null FirstPC");
    if let Curve2d::Trimmed(t) = &new_c2d {
        new_c2d = t.curve.as_ref().clone();
    }

    // OCCT L2319-2321.
    let p2d = DVec2::new(u, v);
    let projector2d = Geom2dAPIProjectPointOnCurve::new(p2d, &new_c2d);
    if projector2d.nb_points() == 0 || projector2d.lower_distance() > tol {
        // OCCT L2323-2352: the 2-pole Bezier patch on the BSpline.
        // OCCT L2324-2331: the Geom2dConvert vehicle (exact passthrough for
        // BSplines, approximate conversion otherwise).
        let mut b_curve = match &new_c2d {
            Curve2d::BSpline(b) => b.clone(),
            other => approx_curve_to_bspline(other, tol)
                .expect("Geom2dConvert::CurveToBSplineCurve"),
        };
        if the_sign == -1 {
            // OCCT L2332-2341: PntArray(1) = P2d; PntArray(2) = Pole(1);
            // Concat.Add(Patch, Tol, false).
            let patch = BezierCurve2 {
                control_points: vec![p2d, b_curve.control_points[0]],
                weights: vec![1.0, 1.0],
            };
            let composed = compose_curves_to_bspline(&[
                Curve2d::Bezier(patch),
                Curve2d::BSpline(b_curve.clone()),
            ])
            .expect("Geom2dConvert_CompCurveToBSplineCurve");
            b_curve = composed;
        } else {
            // OCCT L2342-2351: Pole(NbPoles), P2d; Concat.Add(Patch, Tol, true).
            let last = b_curve.control_points.len() - 1;
            let patch = BezierCurve2 {
                control_points: vec![b_curve.control_points[last], p2d],
                weights: vec![1.0, 1.0],
            };
            let composed = compose_curves_to_bspline(&[
                Curve2d::BSpline(b_curve.clone()),
                Curve2d::Bezier(patch),
            ])
            .expect("Geom2dConvert_CompCurveToBSplineCurve");
            b_curve = composed;
        }
        new_c2d = Curve2d::BSpline(b_curve);
    }
    // OCCT L2354: Einf.ChangeFirstPC() = NewC2d;
    *the_einf.change_first_pc() = Some(new_c2d.clone());

    // OCCT L2355-2368: the approximation tail (GAP carriers — see the module
    // header): ProjLib_HCompProjectedCurve over the curve-on-surface, then
    // Approx_CurveOnSurface.
    let h_cons = BRepAdaptorCurve::new(
        the_einf
            .geometry()
            .cloned()
            .expect("SmartParameter: null Geometry"),
    );
    let h_sur2 = BRepAdaptorSurface::new(the_s2.clone());
    let h_projector = ProjLibHCompProjectedCurve::new(&h_sur2, &h_cons, tol, tol);
    // OCCT L2362-2364: Bounds(1, Udeb, Ufin); MaxSeg = 20 + NbIntervals(C3).
    let (udeb, ufin) = h_projector.bounds(1);
    let max_seg = 20 + h_projector.nb_intervals(GeomAbsShape::C3);
    let mut appr = ApproxCurveOnSurface::new(&h_projector, &h_sur2, udeb, ufin, tol);
    appr.perform(max_seg, 10, GeomAbsShape::C1, false, false);
    // OCCT L2367-2369.
    *the_einf.change_second_pc() = appr.curve2d();
    *the_einf.change_geometry() = appr.curve3d();
    the_einf.set_new_geometry(true);

    // OCCT L2371-2378: the null Curve3d dereference of the OCCT failure
    // output maps to the GAP panic.
    if the_sign == -1 {
        the_einf
            .geometry()
            .expect("SmartParameter: null Curve3d (GAP: Approx_CurveOnSurface)")
            .default_domain()[0]
    } else {
        the_einf
            .geometry()
            .expect("SmartParameter: null Curve3d (GAP: Approx_CurveOnSurface)")
            .default_domain()[1]
    }
}

/// OCCT static Choose(theFMap, theEMap, Vtx, Vinf, AC, AS)
/// (Draft_Modification_1.cxx L2047-2176).  The OCCT in-out adaptor carriers
/// (AC/AS loaded inside) are carried as the tuple return.
pub(crate) fn choose(
    the_fmap: &ShapeIndexedMap<DraftFaceInfo>,
    the_emap: &mut ShapeIndexedMap<DraftEdgeInfo>,
    the_vtx: &Shape,
    the_vinf: &mut DraftVertexInfo,
) -> Option<(BRepAdaptorCurve, BRepAdaptorSurface)> {
    // OCCT L2055: gp_Vec tgref; (uninitialized in the OCCT text).
    let tgref;
    the_vinf.init_edge_iterator();

    // OCCT L2058-2076: find a regular edge with null SecondFace.
    while the_vinf.more_edge() {
        let e1 = the_vinf.edge();
        let einf1 = the_emap.find_from_key(&e1);
        if einf1.second_face().is_null() {
            break;
        } else {
            let te = brep_tool_continuity(&e1, einf1.first_face(), einf1.second_face());
            if te >= GeomAbsShape::G1 {
                break;
            }
        }
        the_vinf.next_edge();
    }
    // OCCT L2077-2080: if (!Vinf.MoreEdge()) Vinf.InitEdgeIterator(); — take
    // the first edge.
    if !the_vinf.more_edge() {
        the_vinf.init_edge_iterator();
    }

    // OCCT L2082-2084.
    let e_ref = the_vinf.edge();
    let einf_geom = the_emap.find_from_key(&e_ref).geometry().cloned();
    // OCCT L2086: AC.Load(Einf.Geometry()); — None carries the OCCT null
    // raise.
    let ac = BRepAdaptorCurve::new(einf_geom.expect("Choose: null Geometry"));

    // OCCT L2088-2091.
    let (c_raw, _f, _l) = brep_tool_curve(&e_ref).expect("BRep_Tool::Curve(Eref)");
    // OCCT: C = C->Transformed(Loc.Transformation()) — identity.
    let c = c_raw;
    // OCCT L2093-2095: int done; double param = Parameter(C, BRep_Tool::Pnt(Vtx), done);
    let mut done = 0i32;
    let param = parameter(&c, brep_tool_pnt(the_vtx), &mut done);
    // OCCT L2096-2105.
    let prm;
    if done != 0 {
        let einf = the_emap.find_from_key(&e_ref);
        let s1 = the_fmap
            .find_from_key(einf.first_face())
            .geometry()
            .cloned()
            .expect("Choose: null S1");
        let s2 = the_fmap
            .find_from_key(einf.second_face())
            .geometry()
            .cloned()
            .expect("Choose: null S2");
        let tol = brep_tool_tolerance(&e_ref);
        let pnt = brep_tool_pnt(the_vtx);
        let einf = the_emap.change_from_key(&e_ref);
        prm = smart_parameter(einf, tol, pnt, done, &s1, &s2);
    } else {
        prm = param;
    }
    // OCCT L2106: C->D1(prm, ptbid, tgref);
    tgref = c.derivative_at(prm);

    // OCCT L2108-2145: find a non-tangent edge.
    the_vinf.init_edge_iterator();
    while the_vinf.more_edge() {
        let edg = the_vinf.edge();
        if !edg.is_same(&e_ref) {
            let (ef, es) = {
                let einfo = the_emap.find_from_key(&edg);
                (einfo.first_face().clone(), einfo.second_face().clone())
            };
            if !es.is_null() && brep_tool_continuity(&edg, &ef, &es) <= GeomAbsShape::C0 {
                let (c2_raw, _f2, _l2) = brep_tool_curve(&edg).expect("BRep_Tool::Curve(Edg)");
                let c2 = c2_raw; // identity location re-application
                // OCCT L2122-2135.
                let mut anewdone = 0i32;
                let anewparam = parameter(&c2, brep_tool_pnt(the_vtx), &mut anewdone);
                let prm2;
                if anewdone != 0 {
                    let s1 = the_fmap
                        .find_from_key(&ef)
                        .geometry()
                        .cloned()
                        .expect("Choose: null S1");
                    let s2 = the_fmap
                        .find_from_key(&es)
                        .geometry()
                        .cloned()
                        .expect("Choose: null S2");
                    let tol = brep_tool_tolerance(&edg);
                    let pnt = brep_tool_pnt(the_vtx);
                    let einfo = the_emap.change_from_key(&edg);
                    prm2 = smart_parameter(einfo, tol, pnt, anewdone, &s1, &s2);
                } else {
                    prm2 = anewparam;
                }
                // OCCT L2136-2137: gp_Vec tg; C->D1(prm, ptbid, tg);
                let tg = c2.derivative_at(prm2);
                // OCCT L2138: if (tg.CrossMagnitude(tgref) > Confusion()) break;
                if tg.cross(tgref).length() > CONFUSION {
                    break;
                }
            }
        }
        the_vinf.next_edge();
    }
    // OCCT L2146-2149.
    if !the_vinf.more_edge() {
        return None;
    }

    // OCCT L2151-2174: load AS with the face of the found edge that is not
    // shared with Einf.
    let e2 = the_vinf.edge();
    let einf2 = the_emap.find_from_key(&e2).clone();
    let einf = the_emap.find_from_key(&e_ref).clone();
    let as_surface = if !einf.second_face().is_null() {
        if einf2.first_face().is_same(einf.first_face())
            || einf2.first_face().is_same(einf.second_face())
        {
            the_fmap.find_from_key(einf2.second_face()).geometry().cloned()
        } else {
            the_fmap.find_from_key(einf2.first_face()).geometry().cloned()
        }
    } else if einf2.first_face().is_same(einf.first_face()) {
        the_fmap.find_from_key(einf2.second_face()).geometry().cloned()
    } else {
        the_fmap.find_from_key(einf2.first_face()).geometry().cloned()
    };
    // OCCT: AS.Load(...) — None carries the OCCT null-handle raise.
    let has = BRepAdaptorSurface::new(as_surface.expect("Choose: null AS surface"));
    // OCCT L2175: return true;
    Some((ac, has))
}
