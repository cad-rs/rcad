//! OCCT ChFi2d_AnaFilletAlgo — 1:1 translation of the analytical 2D fillet
//! algorithm (two conditional edges + a plane -> circular arc transition).
//!
//! Sources:
//!   - ChFi2d_AnaFilletAlgo.hxx (class shell, L25-141)
//!   - ChFi2d_AnaFilletAlgo.cxx (L1-1101)
//!
//! OCCT class `ChFi2d_AnaFilletAlgo` has no base class; rcad keeps it a plain
//! struct with the same member names.
//!
//! Architecture differences (rcad TShape pool vs OCCT handle graphs):
//!   - OCCT reads geometry straight from the shape handles; rcad stores the
//!     shapes in a BRep pool. The algo owns an internal pool `my_brep` as
//!     the carrier for the constructed edges (see the struct note); input
//!     shapes are read through their own TShape handles.
//!   - OCCT exceptions (Standard_TypeMismatch / Standard_Failure) map to
//!     panics carrying the same messages.
//!   - gp_Pnt/gp_Vec map to glam DVec3/DVec2; gp_Lin2d -> [`Line2d`],
//!     gp_Circ2d -> [`Circle2d`], gp_Pln -> [`Plane`] (missing primitive
//!     methods carry their OCCT anchors below).

use std::f64::consts::FRAC_PI_2;
use std::f64::consts::PI;

use glam::DVec2;
use glam::DVec3;

use rcad_kernel::base::extrema::closest_point_on_curve_with_range;
use rcad_kernel::base::extrema::extrema_curve_curve;
use rcad_kernel::base::int_ana2d::AnaIntersection2d;
use rcad_kernel::core::precision::ANGULAR;
use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::core::precision::SQUARE_CONFUSION;
use rcad_kernel::core::precision::is_negative_infinite_value;
use rcad_kernel::core::precision::is_positive_infinite_value;
use rcad_kernel::math::el::elslib_plane_value;
use rcad_kernel::geom::CurveEval as _;
use rcad_kernel::geom::{Circle2d, Circle3, Curve3, Line2d, Plane, Surface3, TrimmedCurve3};
use rcad_kernel::topo::topods::BRepTool as _;
use rcad_kernel::topo::topods::{BRep, Orientation, Shape};

// =========================================================================
// gp primitives missing from rcad (OCCT anchors inline).
// =========================================================================

/// OCCT gp::Resolution() (gp.hxx L60: `return RealSmall();`) — the smallest
/// positive double (DBL_MIN).
#[inline]
fn gp_resolution() -> f64 {
    f64::MIN_POSITIVE
}

/// OCCT gp_Vec2d::Angle (gp_Vec2d.cxx L47-89) — signed angle in (-PI, PI]
/// with the acos/asin precision branching.
fn gp_vec2d_angle(a: DVec2, b: DVec2) -> f64 {
    let a_norm = a.length();
    let b_norm = b.length();
    if a_norm <= gp_resolution() || b_norm <= gp_resolution() {
        // OCCT throws gp_VectorWithNullMagnitude; the callers guard against
        // null vectors before calling Angle, so reaching here is a defect.
        panic!("gp_Vec2d::Angle: null magnitude vector");
    }

    let d = a_norm * b_norm;
    let cosinus = a.dot(b) / d;
    let sinus = (a.x * b.y - a.y * b.x) / d;

    // M_SQRT1_2 — above 45 degrees arccos gives the best precision,
    // otherwise arcsin.
    const COS_45_DEG: f64 = std::f64::consts::FRAC_1_SQRT_2;

    if cosinus > -COS_45_DEG && cosinus < COS_45_DEG {
        if sinus > 0.0 {
            cosinus.acos()
        } else {
            -cosinus.acos()
        }
    } else if cosinus > 0.0 {
        sinus.asin()
    } else if sinus > 0.0 {
        PI - sinus.asin()
    } else {
        -PI - sinus.asin()
    }
}

/// OCCT gp_Vec2d::Rotate(angle) (gp_Vec2d.cxx): X' = X*cos - Y*sin,
/// Y' = X*sin + Y*cos.
fn gp_vec2d_rotate(v: DVec2, angle: f64) -> DVec2 {
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    DVec2::new(
        v.x * cos_a - v.y * sin_a,
        v.x * sin_a + v.y * cos_a,
    )
}

/// OCCT gp_Vec::AngleWithRef(Other, VRef) -> gp_Dir::AngleWithRef
/// (gp_Dir.cxx): unsigned angle in [0, PI], signed by the cross product
/// against `vref` (negative when (a ^ b) . vref < 0).
fn gp_dir_angle_with_ref(a: DVec3, b: DVec3, vref: DVec3) -> f64 {
    let xyz = a.cross(b);
    let cosinus = a.dot(b);
    let sinus = xyz.length();
    let ang = if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        cosinus.acos()
    } else if cosinus < 0.0 {
        PI - sinus.asin()
    } else {
        sinus.asin()
    };
    if xyz.dot(vref) >= 0.0 {
        ang
    } else {
        -ang
    }
}

/// OCCT gp_Circ2d::Distance(P) (gp_Circ2d.hxx L274-284):
/// `|radius - |P - Loc||`.
fn gp_circ2d_distance(p: DVec2, circ: &Circle2d) -> f64 {
    let mut d = circ.radius - (p - circ.center).length();
    if d < 0.0 {
        d = -d;
    }
    d
}

/// OCCT static IsEqual(gp_Pnt, gp_Pnt) (ChFi2d_AnaFilletAlgo.cxx L72-75) —
/// equality through the square distance vs Precision::SquareConfusion().
fn point_is_equal_pnt(p1: DVec3, p2: DVec3) -> bool {
    p1.distance_squared(p2) < SQUARE_CONFUSION
}

/// OCCT static IsEqual(gp_Pnt2d, gp_Pnt2d) (ChFi2d_AnaFilletAlgo.cxx L77-80).
fn point_is_equal_pnt2d(p1: DVec2, p2: DVec2) -> bool {
    p1.distance_squared(p2) < SQUARE_CONFUSION
}

/// OCCT ProjLib::Project(const gp_Pln&, const gp_Pnt&) (ProjLib.cxx L50-56)
/// = ElSLib::Parameters(Pl, P) = the plane-frame UV coordinates.
/// rcad has no plane-point ProjLib yet (gap), ElSLib::PlaneParameters exists
/// in rcad-kernel::math::el.
fn projlib_project_plane_point(plane: &Plane, p: DVec3) -> DVec2 {
    let d = p - plane.origin;
    DVec2::new(d.dot(plane.u_dir), d.dot(plane.v_dir))
}

/// OCCT ElSLib::Value(u, v, Pl) — 3D point of plane parameters.
fn elslib_plane_point(u: f64, v: f64, plane: &Plane) -> DVec3 {
    elslib_plane_value(u, v, plane.origin, plane.u_dir, plane.v_dir)
}

/// OCCT TopoDS_Shape::Reversed() — the orientation composed with
/// TopAbs_REVERSED (FORWARD <-> REVERSED swap, INTERNAL/EXTERNAL kept).
fn topods_reversed(s: &Shape) -> Shape {
    let mut r = s.clone();
    r.orientation = match s.orientation {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        other => other,
    };
    r
}

/// OCCT gp_Pln::Position().Direct() — directness of the plane's Ax3.
/// Architecture difference: rcad Plane stores the frame directions
/// explicitly, so the directness test is the triple product sign.
fn plane_position_direct(plane: &Plane) -> bool {
    plane.normal.cross(plane.u_dir).dot(plane.v_dir) > 0.0
}

/// OCCT TopExp::Vertices(E, Vfirst, Vlast, CumOri = Standard_True): the
/// returned vertices take the edge orientation into account. rcad edges
/// store the canonical (FORWARD-order) ends on the TShape; CumOri swaps them
/// for REVERSED edges (architecture difference).
fn topexp_vertices_cum_ori(e: &Shape) -> (Shape, Shape) {
    let ed = e.as_edge().expect("TopExp::Vertices: not an edge");
    if e.orientation == Orientation::Reversed {
        (ed.last.clone(), ed.first.clone())
    } else {
        (ed.first.clone(), ed.last.clone())
    }
}

// =========================================================================
// OCCT GeomAbs_CurveType subset used by BRepAdaptor_Curve::GetType().
// =========================================================================
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GeomAbsCurveType {
    Line,
    Circle,
    Other,
}

/// OCCT BRepAdaptor_Curve stand-in: the edge's 3D curve plus its parameter
/// range (rcad edges carry both on the TShape; OCCT reads them through
/// BRep_Tool::Curve with the Location applied — rcad `edge_curve_world`).
struct BRepAdaptorCurve {
    curve: Curve3,
    range: [f64; 2],
}

impl BRepAdaptorCurve {
    /// OCCT BRepAdaptor_Curve(E).
    fn new(brep: &BRep, e: &Shape) -> Self {
        let (curve, range) = brep
            .edge_curve_world(e)
            .expect("BRepAdaptor_Curve: edge without 3d curve");
        BRepAdaptorCurve { curve, range }
    }

    /// OCCT GetType().
    fn get_type(&self) -> GeomAbsCurveType {
        match self.curve {
            Curve3::Line(_) => GeomAbsCurveType::Line,
            Curve3::Circle(_) => GeomAbsCurveType::Circle,
            _ => GeomAbsCurveType::Other,
        }
    }

    /// OCCT FirstParameter().
    fn first_parameter(&self) -> f64 {
        self.range[0]
    }

    /// OCCT LastParameter().
    fn last_parameter(&self) -> f64 {
        self.range[1]
    }

    /// OCCT Value(t).
    fn value(&self, t: f64) -> DVec3 {
        self.curve.point_at(t)
    }

    /// OCCT Circle().
    fn circle(&self) -> Circle3 {
        match self.curve {
            Curve3::Circle(c) => c,
            _ => panic!("BRepAdaptor_Curve::Circle: not a circle"),
        }
    }

    fn curve(&self) -> &Curve3 {
        &self.curve
    }
}

// =========================================================================
// BRepBuilderAPI / ShapeAnalysis stand-ins (architecture differences).
// =========================================================================

/// OCCT BRepBuilderAPI_MakeWire stand-in: rcad builds the wire directly on
/// the BRep pool (BRep::add_twire, the BRepPrim_Builder path). OCCT
/// MakeWire additionally checks connectivity and reports IsDone; rcad
/// returns None only when a member is not an edge TShape.
fn brepbuilderapi_make_wire(brep: &mut BRep, edges: Vec<Shape>) -> Option<Shape> {
    for e in &edges {
        if e.as_edge().is_none() {
            return None;
        }
    }
    Some(brep.add_twire(edges))
}

/// OCCT BRepBuilderAPI_MakeFace(const gp_Pln&) stand-in: an infinite planar
/// face (natural restriction, no wires). IsDone is always true in OCCT.
fn brepbuilderapi_make_face_plane(brep: &mut BRep, plane: &Plane) -> Option<Shape> {
    Some(brep.add_tface(
        Some(Surface3::Plane(plane.clone())),
        Shape::null(),
        Vec::new(),
        None,
        None,
        Vec::new(),
        true,
    ))
}

/// OCCT ShapeAnalysis_Wire(W, F, Preci).CheckSelfIntersection().
///
/// OCCT GAP: rcad has no ShapeAnalysis_Wire equivalent yet
/// (the shhealing::shape_analysis module only exposes index-based whole-BRep
/// reports, not a (wire, face) pair query). Reported to the main agent as an
/// infrastructure gap; until filled, the check reports "no self
/// intersection" (the OCCT `false` branch).
fn shape_analysis_wire_check_self_intersection(
    _brep: &BRep,
    _wire: &Shape,
    _face: &Shape,
    _preci: f64,
) -> bool {
    false
}

/// OCCT BRepLib_MakeEdge::Project(C, V, p) (BRepLib_MakeEdge.cxx L55-101):
/// vertex parameter by projection onto the curve. First the distances at the
/// curve ends are verified (within the vertex tolerance); otherwise the
/// extrema are computed and the minimal one accepted when within Eps2.
/// rcad: closest_point_on_curve_with_range (Extrema_ExtPC equivalent);
/// its single returned projection is the minimal one.
fn breplib_make_edge_project(brep: &BRep, curve: &Curve3, v: &Shape, p: &mut f64) -> bool {
    let eps2 = {
        let t = brep.vertex_tolerance(v);
        t * t
    };
    let point = brep.vertex_position(v);

    // Afin de faire les extremas, on verifie les distances en bout.
    let dom = curve.default_domain();
    let p1 = curve.point_at(dom[0]);
    let p2 = curve.point_at(dom[1]);
    let d1 = p1.distance_squared(point);
    let d2 = p2.distance_squared(point);
    if d1 < d2 && d1 <= eps2 {
        *p = dom[0];
        return true;
    } else if d2 < d1 && d2 <= eps2 {
        *p = dom[1];
        return true;
    }

    // Sinon, on calcule les extremas.
    let ext = closest_point_on_curve_with_range(curve, point, 64, dom[0], dom[1]);
    let dist2 = ext.distance * ext.distance;
    if dist2 <= eps2 {
        *p = ext.param;
        return true;
    }
    false
}

/// OCCT ElCLib::AdjustPeriodic(UFirst, ULast, Preci, U1, U2)
/// (ElCLib.cxx): adjusts U1 into [UFirst, ULast] and U2 = U1 + adjusted
/// (U2 - U1) so that U2 > U1 within the period.
fn elclib_adjust_periodic(ufirst: f64, ulast: f64, preci: f64, u1: &mut f64, u2: &mut f64) {
    if !ufirst.is_finite() || !ulast.is_finite() {
        *u1 = ufirst;
        *u2 = ulast;
        return;
    }

    let period = ulast - ufirst;
    if period < ulast.abs() * f64::EPSILON {
        *u1 = ufirst;
        *u2 = ulast;
        return;
    }

    *u1 -= ((*u1 - ufirst) / period).floor() * period;
    if ulast - *u1 < preci {
        *u1 -= period;
    }
    *u2 -= ((*u2 - *u1) / period).floor() * period;
    if *u2 - *u1 < preci {
        *u2 += period;
    }
}

/// OCCT BRepLib_MakeEdge::Init(C, VV1, VV2, pp1, pp2) — "this one really
/// makes the job" (BRepLib_MakeEdge.cxx L602-789): kill trimmed curves,
/// periodic adjust / reordonate, compute points on the curve, closed-curve
/// handling, build the edge (vertices FORWARD / REVERSED, B.Range).
/// rcad: BRep::add_tvertex + add_tedge; the vertex orientation slots and the
/// closed/degenerated flags are carried by the edge TShape (architecture
/// difference).
fn breplib_make_edge_build(
    brep: &mut BRep,
    cc: &Curve3,
    vv1: Option<&Shape>,
    vv2: Option<&Shape>,
    pp1: f64,
    pp2: f64,
) -> Option<Shape> {
    // Kill trimmed curves.
    let mut c = cc.clone();
    while let Curve3::Trimmed(tc) = c.clone() {
        c = tc.basis_curve().clone();
    }

    // Check parameters.
    let mut p1 = pp1;
    let mut p2 = pp2;
    let cf = c.default_domain()[0];
    let cl = c.default_domain()[1];
    let epsilon = CONFUSION * 0.01; // Precision::PConfusion()
    let periodic = match c {
        Curve3::Circle(_) | Curve3::Ellipse(_) => true,
        _ => false,
    };

    let (mut v1, mut v2) = if periodic {
        // Adjust in period.
        elclib_adjust_periodic(cf, cl, epsilon, &mut p1, &mut p2);
        (vv1.cloned(), vv2.cloned())
    } else {
        // Reordonate.
        let (r1, r2) = if p1 < p2 {
            (vv1.cloned(), vv2.cloned())
        } else {
            let x = p1;
            p1 = p2;
            p2 = x;
            (vv2.cloned(), vv1.cloned())
        };

        // Check range.
        if cf - p1 > epsilon || p2 - cl > epsilon {
            return None; // BRepLib_ParameterOutOfRange
        }

        // Check punctuality.
        if p2 - p1 <= gp_resolution() {
            return None; // BRepLib_LineThroughIdenticPoints
        }
        (r1, r2)
    };

    // Compute points on the curve.
    let p1inf = is_negative_infinite_value(p1);
    let p2inf = is_positive_infinite_value(p2);
    let point1 = if !p1inf { c.point_at(p1) } else { DVec3::ZERO };
    let point2 = if !p2inf { c.point_at(p2) } else { DVec3::ZERO };

    let preci = CONFUSION; // BRepLib::Precision()

    // Check for closed curve.
    let mut closed = false;
    let mut degenerated = false;
    if !p1inf && !p2inf {
        closed = point1.distance(point2) <= preci;
    }

    // Check if the vertices are on the curve.
    if closed {
        if v1.is_none() && v2.is_none() {
            let v = brep.add_tvertex(point1);
            v1 = Some(v.clone());
            v2 = Some(v);
        } else if v1.is_none() {
            v1 = v2.clone();
        } else if v2.is_none() {
            v2 = v1.clone();
        } else {
            let uv1 = v1.clone().unwrap_or_else(Shape::null);
            let uv2 = v2.clone().unwrap_or_else(Shape::null);
            if !uv1.is_same(&uv2) {
                return None; // BRepLib_DifferentPointsOnClosedCurve
            } else if point1.distance(brep.vertex_position(&uv1))
                > preci.max(brep.vertex_tolerance(&uv1))
            {
                return None; // BRepLib_DifferentPointsOnClosedCurve
            } else {
                let pm = c.point_at(0.5 * (p1 + p2));
                if point1.distance(pm) < preci {
                    degenerated = true;
                }
            }
        }
    } else {
        // not closed
        if p1inf {
            if v1.is_some() {
                return None; // BRepLib_PointWithInfiniteParameter
            }
        } else {
            match v1.clone() {
                None => v1 = Some(brep.add_tvertex(point1)),
                Some(uv1) => {
                    if point1.distance(brep.vertex_position(&uv1))
                        > preci.max(brep.vertex_tolerance(&uv1))
                    {
                        return None; // BRepLib_DifferentsPointAndParameter
                    }
                }
            }
        }

        if p2inf {
            if v2.is_some() {
                return None; // BRepLib_PointWithInfiniteParameter
            }
        } else {
            match v2.clone() {
                None => v2 = Some(brep.add_tvertex(point2)),
                Some(uv2) => {
                    if point2.distance(brep.vertex_position(&uv2))
                        > preci.max(brep.vertex_tolerance(&uv2))
                    {
                        return None; // BRepLib_DifferentsPointAndParameter
                    }
                }
            }
        }
    }

    // V1.Orientation(TopAbs_FORWARD); V2.Orientation(TopAbs_REVERSED);
    // B.MakeEdge(E, C, preci); B.Add(E, V1); B.Add(E, V2);
    // B.Range(E, p1, p2); B.Degenerated(E, degenerated); E.Closed(closed).
    // rcad add_tedge stores the canonical ends and derives the degenerated
    // flag; the OCCT closed flag is implicit in coinciding range points.
    let e = brep.add_tedge(
        Some(c),
        v1.unwrap_or_else(Shape::null),
        v2.unwrap_or_else(Shape::null),
        [p1, p2],
    );
    let _ = (closed, degenerated);
    Some(e)
}

/// OCCT BRepBuilderAPI_MakeEdge(C, P1, P2) (BRepLib_MakeEdge.cxx L522-540):
/// vertices created at the points (BRepLib::Precision() tolerance), then
/// Init(C, V1, V2) projects them onto the curve for the parameters.
fn brepbuilderapi_make_edge_curve_point_point(
    brep: &mut BRep,
    c: &Curve3,
    p1: DVec3,
    p2: DVec3,
) -> Option<Shape> {
    let tol = CONFUSION; // BRepLib::Precision()
    let v1 = brep.add_tvertex(p1);
    let v2 = if p1.distance(p2) < tol {
        v1.clone()
    } else {
        brep.add_tvertex(p2)
    };

    // Init(C, V1, V2): try projecting the vertices on the curve.
    let mut pp1 = 0.0;
    let mut pp2 = 0.0;
    if !breplib_make_edge_project(brep, c, &v1, &mut pp1) {
        return None; // BRepLib_PointProjectionFailed
    }
    if !breplib_make_edge_project(brep, c, &v2, &mut pp2) {
        return None; // BRepLib_PointProjectionFailed
    }
    breplib_make_edge_build(brep, c, Some(&v1), Some(&v2), pp1, pp2)
}

/// OCCT BRepBuilderAPI_MakeEdge(C, p1, p2) (BRepLib_MakeEdge.cxx L512-517):
/// parameter-range constructor with null vertices.
fn brepbuilderapi_make_edge_curve_params(
    brep: &mut BRep,
    c: &Curve3,
    p1: f64,
    p2: f64,
) -> Option<Shape> {
    breplib_make_edge_build(brep, c, None, None, p1, p2)
}

// =========================================================================
// OCCT Standard_TypeMismatch / Standard_Failure architecture note: OCCT
// throws these from Init and the constructors; rcad panics with the same
// messages (no exception channel in the rcad API surface).
// =========================================================================

// =========================================================================
// OCCT ChFi2d_AnaFilletAlgo class (ChFi2d_AnaFilletAlgo.hxx L25-141).
// =========================================================================
//
// Architecture difference: OCCT BRepBuilderAPI_* builders create new TShapes
// in the ambient shape graph; rcad shapes live in a BRep pool. The algo
// therefore owns an internal pool `my_brep` where Perform/Cut build the
// fillet and shrinked edges (the MakeEdge / MakeWire / MakeFace stand-in
// carrier). The Init step only reads the input shapes.
#[derive(Debug, Clone)]
pub struct ChFi2dAnaFilletAlgo {
    /// Internal TShape pool for the constructed edges (rcad-only carrier,
    /// see the architecture note above).
    my_brep: BRep,
    /// Plane (hxx L111).
    plane: Plane,

    /// Left neighbour (hxx L114-123).
    e1: Shape,
    segment1: bool,
    x11: f64,
    y11: f64,
    x12: f64,
    y12: f64,
    xc1: f64,
    yc1: f64,
    radius1: f64,
    cw1: bool,

    /// Right neighbour (hxx L126-135).
    e2: Shape,
    segment2: bool,
    x21: f64,
    y21: f64,
    x22: f64,
    y22: f64,
    xc2: f64,
    yc2: f64,
    radius2: f64,
    cw2: bool,

    /// Fillet (result) (hxx L138-140).
    fillet: Shape,
    shrinke1: Shape,
    shrinke2: Shape,
}

impl ChFi2dAnaFilletAlgo {
    /// OCCT static isCW (ChFi2d_AnaFilletAlgo.cxx L40-69) — computes the
    /// CW || CCW flag of an adaptor curve on a circle.
    fn is_cw(ac: &BRepAdaptorCurve) -> bool {
        let f = ac.first_parameter();
        let l = ac.last_parameter();
        let circle = ac.circle(); // occ::handle<Geom_Circle>
        let start = ac.value(f);
        let end = ac.value(l);
        let center = circle.center; // AC.Circle().Location()
        let plane_dir = circle.normal; // AC.Circle().Position().Direction()

        // Get point on circle at half angle.
        let m = circle.point_at(0.5 * (f + l)); // circle->D0(0.5 * (f + l), m)

        // Compare angles between vectors to middle point and to the end
        // point.
        let startv = start - center;
        let endv = end - center;
        let middlev = m - center;
        let mut middlea = gp_dir_angle_with_ref(startv, middlev, plane_dir);
        while middlea < 0.0 {
            middlea += 2.0 * PI;
        }
        let mut enda = gp_dir_angle_with_ref(startv, endv, plane_dir);
        while enda < 0.0 {
            enda += 2.0 * PI;
        }

        let is_cw = middlea > enda;
        is_cw
    }

    /// OCCT ChFi2d_AnaFilletAlgo() — an empty constructor (cxx L84-104).
    /// Use the method init() to initialize the class.
    pub fn new() -> Self {
        ChFi2dAnaFilletAlgo {
            my_brep: BRep::new(),
            plane: Plane {
                origin: DVec3::ZERO,
                normal: DVec3::Z,
                u_dir: DVec3::X,
                v_dir: DVec3::Y,
            },
            e1: Shape::null(),
            segment1: false,
            x11: 0.0,
            y11: 0.0,
            x12: 0.0,
            y12: 0.0,
            xc1: 0.0,
            yc1: 0.0,
            radius1: 0.0,
            cw1: false,
            e2: Shape::null(),
            segment2: false,
            x21: 0.0,
            y21: 0.0,
            x22: 0.0,
            y22: 0.0,
            xc2: 0.0,
            yc2: 0.0,
            radius2: 0.0,
            cw2: false,
            fillet: Shape::null(),
            shrinke1: Shape::null(),
            shrinke2: Shape::null(),
        }
    }

    /// OCCT ChFi2d_AnaFilletAlgo(theWire, thePlane) (cxx L110-132) — a
    /// constructor. It expects a wire consisting of two edges of type (any
    /// combination of): segment, arc of circle.
    pub fn new_wire(the_wire: &Shape, the_plane: &Plane) -> Self {
        let mut algo = Self::new();
        algo.plane = the_plane.clone();
        algo.init_wire(the_wire, the_plane);
        algo
    }

    /// OCCT ChFi2d_AnaFilletAlgo(theEdge1, theEdge2, thePlane)
    /// (cxx L138-163) — a constructor. It expects two edges having a common
    /// point of type: segment, arc of circle.
    pub fn new_edges(the_edge1: &Shape, the_edge2: &Shape, the_plane: &Plane) -> Self {
        let mut algo = Self::new();
        algo.plane = the_plane.clone();
        // Make a wire consisting of two edges.
        algo.init_edges(the_edge1, the_edge2, the_plane);
        algo
    }

    /// OCCT Init(theWire, thePlane) (cxx L166-253) — initializes the class
    /// by a wire consisting of two edges.
    ///
    /// Architecture difference: OCCT throws Standard_TypeMismatch /
    /// Standard_Failure; rcad panics with the same message (no exception
    /// channel in the rcad API surface).
    pub fn init_wire(&mut self, the_wire: &Shape, the_plane: &Plane) {
        self.plane = the_plane.clone();

        // OCCT L169: TopoDS_Iterator itr(theWire) — iterate the wire's
        // sub-shapes (rcad TWireData::edges).
        let wire_edges: &[Shape] = match the_wire.as_wire() {
            Some(wd) => &wd.edges,
            None => &[],
        };
        for itr_value in wire_edges {
            if self.e1.is_null() {
                // TopoDS::Edge(itr.Value()).
                self.e1 = itr_value.clone();
            } else if self.e2.is_null() {
                self.e2 = itr_value.clone();
            }
        }
        if self.e1.is_null() || self.e2.is_null() {
            // throw Standard_TypeMismatch(...)
            panic!("The algorithm expects a wire consisting of two linear or circular edges.");
        }

        // Left neighbour.
        let ac1 = BRepAdaptorCurve::new(&self.my_brep, &self.e1);
        if ac1.get_type() != GeomAbsCurveType::Line && ac1.get_type() != GeomAbsCurveType::Circle {
            // throw Standard_TypeMismatch(...)
            panic!("A segment or an arc of circle is expected.");
        }

        let (v1, v2) = topexp_vertices_cum_ori(&self.e1);
        if v1.is_null() || v2.is_null() {
            // throw Standard_Failure(...)
            panic!("An infinite edge.");
        }

        let big_p1 = self.my_brep.vertex_position(&v1);
        let big_p2 = self.my_brep.vertex_position(&v2);
        let p1 = projlib_project_plane_point(the_plane, big_p1);
        let p2 = projlib_project_plane_point(the_plane, big_p2);
        self.x11 = p1.x;
        self.y11 = p1.y;
        self.x12 = p2.x;
        self.y12 = p2.y;

        self.segment1 = true;
        if ac1.get_type() == GeomAbsCurveType::Circle {
            self.segment1 = false;
            let c = ac1.circle();

            let loc = projlib_project_plane_point(the_plane, c.center);
            self.xc1 = loc.x;
            self.yc1 = loc.y;

            self.radius1 = c.radius;
            self.cw1 = Self::is_cw(&ac1);
        }

        // Right neighbour.
        let ac2 = BRepAdaptorCurve::new(&self.my_brep, &self.e2);
        if ac2.get_type() != GeomAbsCurveType::Line && ac2.get_type() != GeomAbsCurveType::Circle {
            // throw Standard_TypeMismatch(...)
            panic!("A segment or an arc of circle is expected.");
        }

        let (v1, v2) = topexp_vertices_cum_ori(&self.e2);
        if v1.is_null() || v2.is_null() {
            // throw Standard_Failure(...)
            panic!("An infinite edge.");
        }

        let big_p1 = self.my_brep.vertex_position(&v1);
        let big_p2 = self.my_brep.vertex_position(&v2);
        let p1 = projlib_project_plane_point(the_plane, big_p1);
        let p2 = projlib_project_plane_point(the_plane, big_p2);
        self.x21 = p1.x;
        self.y21 = p1.y;
        self.x22 = p2.x;
        self.y22 = p2.y;

        self.segment2 = true;
        if ac2.get_type() == GeomAbsCurveType::Circle {
            self.segment2 = false;
            let c = ac2.circle();

            let loc = projlib_project_plane_point(the_plane, c.center);
            self.xc2 = loc.x;
            self.yc2 = loc.y;

            self.radius2 = c.radius;
            self.cw2 = Self::is_cw(&ac2);
        }
    }

    /// OCCT Init(theEdge1, theEdge2, thePlane) (cxx L256-326) — initializes
    /// the class by two edges.
    pub fn init_edges(&mut self, the_edge1: &Shape, the_edge2: &Shape, the_plane: &Plane) {
        // Make a wire consisting of two edges.

        // Get common point.
        let (v11, v12) = topexp_vertices_cum_ori(the_edge1);
        let (v21, v22) = topexp_vertices_cum_ori(the_edge2);
        if v11.is_null() || v12.is_null() || v21.is_null() || v22.is_null() {
            // throw Standard_Failure(...)
            panic!("An infinite edge.");
        }

        let p11 = self.my_brep.vertex_position(&v11);
        let p12 = self.my_brep.vertex_position(&v12);
        let p21 = self.my_brep.vertex_position(&v21);
        let p22 = self.my_brep.vertex_position(&v22);

        let pcommon;
        if point_is_equal_pnt(p11, p21) || point_is_equal_pnt(p11, p22) {
            pcommon = p11;
        } else if point_is_equal_pnt(p12, p21) || point_is_equal_pnt(p12, p22) {
            pcommon = p12;
        } else {
            // throw Standard_Failure(...)
            panic!("The edges have no common point.");
        }

        // Reverse the edges in case of need (to construct a wire).
        let (mut is1st_reversed, mut is2nd_reversed) = (false, false);
        if point_is_equal_pnt(pcommon, p11) {
            is1st_reversed = true;
        } else if point_is_equal_pnt(pcommon, p22) {
            is2nd_reversed = true;
        }

        // Make a wire (BRepBuilderAPI_MakeWire mkWire).
        let mut mk_wire_edges: Vec<Shape> = Vec::new();
        if is1st_reversed {
            mk_wire_edges.push(topods_reversed(the_edge1));
        } else {
            mk_wire_edges.push(the_edge1.clone());
        }
        if is2nd_reversed {
            mk_wire_edges.push(topods_reversed(the_edge2));
        } else {
            mk_wire_edges.push(the_edge2.clone());
        }
        let mk_wire = brepbuilderapi_make_wire(&mut self.my_brep, mk_wire_edges);
        let w = match mk_wire {
            Some(w) => w,
            // throw Standard_Failure("Can't make a wire.")
            None => panic!("Can't make a wire."),
        };

        self.init_wire(&w, the_plane);
    }

    /// OCCT Perform(radius) (cxx L329-571) — calculates a fillet.
    pub fn perform(&mut self, radius: f64) -> bool {
        let mut b_ret = false;
        if self.e1.is_null() || self.e2.is_null() || radius < CONFUSION {
            return b_ret;
        }

        // Fillet definition.
        let mut xc = 0.0;
        let mut yc = 0.0;
        let mut start = 0.0;
        let mut end = 0.0; // parameters on neighbours
        let mut xstart = f64::MAX;
        let mut ystart = f64::MAX; // point on left neighbour
        let mut xend = f64::MAX;
        let mut yend = f64::MAX; // point on right neighbour
        let mut cw = false;

        // Analytical algorithm works for non-intersecting arcs only.
        // Check arcs on self-intersection.
        let mut is_cut = false;
        if !self.segment1 || !self.segment2 {
            let mk_wire =
                brepbuilderapi_make_wire(&mut self.my_brep, vec![self.e1.clone(), self.e2.clone()]);
            if let Some(w) = mk_wire {
                let mk_face = brepbuilderapi_make_face_plane(&mut self.my_brep, &self.plane);
                if let Some(f) = mk_face {
                    if shape_analysis_wire_check_self_intersection(&self.my_brep, &w, &f, CONFUSION)
                    {
                        // Cut the edges at the point of intersection.
                        is_cut = true;
                        // OCCT Cut(plane, e1, e2) passes the member edges by
                        // reference; rcad works on clones and writes them
                        // back (borrow rules).
                        let mut e1 = self.e1.clone();
                        let mut e2 = self.e2.clone();
                        let the_plane = self.plane.clone();
                        let cut_ok = self.cut(&the_plane, &mut e1, &mut e2);
                        self.e1 = e1;
                        self.e2 = e2;
                        if !cut_ok {
                            return false;
                        }
                    }
                }
            }
        } // a case of segment - segment

        // Choose the case.
        let ac1 = BRepAdaptorCurve::new(&self.my_brep, &self.e1);
        let ac2 = BRepAdaptorCurve::new(&self.my_brep, &self.e2);
        if self.segment1 && self.segment2 {
            b_ret =
                self.segment_fillet_segment(radius, &mut xc, &mut yc, &mut cw, &mut start, &mut end);
        } else if self.segment1 && !self.segment2 {
            b_ret = self.segment_fillet_arc(
                radius, &mut xc, &mut yc, &mut cw, &mut start, &mut end, &mut xend, &mut yend,
            );
        } else if !self.segment1 && self.segment2 {
            b_ret = self.arc_fillet_segment(
                radius, &mut xc, &mut yc, &mut cw, &mut start, &mut end, &mut xstart, &mut ystart,
            );
        } else if !self.segment1 && !self.segment2 {
            b_ret = self.arc_fillet_arc(radius, &mut xc, &mut yc, &mut cw, &mut start, &mut end);
        }

        if !b_ret {
            return false;
        }

        // Invert the fillet for left-handed plane.
        if !plane_position_direct(&self.plane) {
            cw = !cw;
        }

        // Construct a fillet.
        // Make circle.
        let mut center = elslib_plane_point(xc, yc, &self.plane); // ElSLib::Value(xc, yc, plane)
        let normal = self.plane.normal; // plane.Position().Direction()
        let mut circ = Circle3::new(
            center,
            if cw { -normal } else { normal },
            radius,
        ); // gp_Circ circ(gp_Ax2(center, cw ? -normal : normal), radius)

        // Fillet may only shrink a neighbour edge, it can't prolongate it.
        let delta1 = ac1.last_parameter() - ac1.first_parameter();
        let delta2 = ac2.last_parameter() - ac2.first_parameter();
        if !is_cut && (start > delta1 || end > delta2) {
            // Check a case when a neighbour edge almost disappears:
            // try to reduce the fillet radius for a little (1.e-5 mm).
            let little = 100.0 * CONFUSION;
            let d1 = (start - delta1).abs();
            let d2 = (end - delta2).abs();
            if d1 < little || d2 < little {
                if self.segment1 && self.segment2 {
                    b_ret = self.segment_fillet_segment(
                        radius - little,
                        &mut xc,
                        &mut yc,
                        &mut cw,
                        &mut start,
                        &mut end,
                    );
                } else if self.segment1 && !self.segment2 {
                    b_ret = self.segment_fillet_arc(
                        radius - little,
                        &mut xc,
                        &mut yc,
                        &mut cw,
                        &mut start,
                        &mut end,
                        &mut xend,
                        &mut yend,
                    );
                } else if !self.segment1 && self.segment2 {
                    b_ret = self.arc_fillet_segment(
                        radius - little,
                        &mut xc,
                        &mut yc,
                        &mut cw,
                        &mut start,
                        &mut end,
                        &mut xstart,
                        &mut ystart,
                    );
                } else if !self.segment1 && !self.segment2 {
                    b_ret = self.arc_fillet_arc(
                        radius - little,
                        &mut xc,
                        &mut yc,
                        &mut cw,
                        &mut start,
                        &mut end,
                    );
                }
                if b_ret {
                    // Invert the fillet for left-handed planes.
                    if !plane_position_direct(&self.plane) {
                        cw = !cw;
                    }

                    // Make the circle again.
                    center = elslib_plane_point(xc, yc, &self.plane);
                    circ.center = center; // circ.SetLocation(center)
                    circ.radius = radius - little; // circ.SetRadius(radius - little)
                } else {
                    return false;
                }
            } else {
                return false;
            }
        }
        if b_ret {
            // start: (xstart, ystart) - pstart.
            let pstart;
            if xstart != f64::MAX {
                pstart = elslib_plane_point(xstart, ystart, &self.plane);
            } else if self.e1.orientation == Orientation::Forward {
                pstart = ac1.value(ac1.last_parameter() - start);
            } else {
                pstart = ac1.value(ac1.first_parameter() + start);
            }
            // end: (xend, yend) -> pend.
            let pend;
            if xend != f64::MAX {
                pend = elslib_plane_point(xend, yend, &self.plane);
            } else if self.e2.orientation == Orientation::Forward {
                pend = ac2.value(ac2.first_parameter() + end);
            } else {
                pend = ac2.value(ac2.last_parameter() - end);
            }

            // Make arc.
            let mk_edge = brepbuilderapi_make_edge_curve_point_point(
                &mut self.my_brep,
                &Curve3::Circle(circ.clone()),
                pstart,
                pend,
            );
            b_ret = mk_edge.is_some();
            if b_ret {
                self.fillet = mk_edge.unwrap();

                // Limit the neighbours.
                // Left neighbour.
                let mut p1;
                let mut p2;
                self.shrinke1 = Shape::null(); // shrinke1.Nullify()
                if self.e1.orientation == Orientation::Forward {
                    p1 = ac1.value(ac1.first_parameter());
                    p2 = pstart;
                } else {
                    p1 = pstart;
                    p2 = ac1.value(ac1.last_parameter());
                }
                if self.segment1 {
                    let mk_segment1 = brepbuilderapi_make_edge_curve_point_point(
                        &mut self.my_brep,
                        ac1.curve(),
                        p1,
                        p2,
                    );
                    if let Some(e) = mk_segment1 {
                        self.shrinke1 = e;
                    }
                } else {
                    let mk_circ1 = brepbuilderapi_make_edge_curve_point_point(
                        &mut self.my_brep,
                        ac1.curve(),
                        p1,
                        p2,
                    );
                    if let Some(e) = mk_circ1 {
                        self.shrinke1 = e;
                    }
                }

                // Right neighbour.
                self.shrinke2 = Shape::null(); // shrinke2.Nullify()
                if self.e1.orientation == Orientation::Forward {
                    // OCCT keeps the e1 test here (cxx L537) — kept as is.
                    p1 = pend;
                    p2 = ac2.value(ac2.last_parameter());
                } else {
                    p1 = ac2.value(ac2.first_parameter());
                    p2 = pend;
                }
                if self.segment2 {
                    let mk_segment2 = brepbuilderapi_make_edge_curve_point_point(
                        &mut self.my_brep,
                        ac2.curve(),
                        p1,
                        p2,
                    );
                    if let Some(e) = mk_segment2 {
                        self.shrinke2 = e;
                    }
                } else {
                    let mk_circ2 = brepbuilderapi_make_edge_curve_point_point(
                        &mut self.my_brep,
                        ac2.curve(),
                        p1,
                        p2,
                    );
                    if let Some(e) = mk_circ2 {
                        self.shrinke2 = e;
                    }
                }

                b_ret = !self.shrinke1.is_null() && !self.shrinke2.is_null();
            } // fillet edge is done
        } // shrinking is good

        b_ret
    }

    /// OCCT Result(theE1, theE2) (cxx L574-579) — retrieves a result (fillet
    /// and shrinked neighbours). OCCT returns a const reference to `fillet`;
    /// rcad returns a clone (architecture difference).
    pub fn result(&self, the_e1: &mut Shape, the_e2: &mut Shape) -> Shape {
        *the_e1 = self.shrinke1.clone();
        *the_e2 = self.shrinke2.clone();
        self.fillet.clone()
    }

    /// OCCT SegmentFilletSegment (cxx L587-640) — WW5 method to compute
    /// fillet. It returns a constructed fillet definition: center point
    /// (xc, yc), point on the 1st segment (start), point on the 2nd segment
    /// (end), is the arc of fillet clockwise (cw = true) or
    /// counterclockwise (cw = false).
    fn segment_fillet_segment(
        &self,
        radius: f64,
        xc: &mut f64,
        yc: &mut f64,
        cw: &mut bool,
        start: &mut f64,
        end: &mut f64,
    ) -> bool {
        // Make normalized vectors at p12.
        let p11 = DVec2::new(self.x11, self.y11);
        let p12 = DVec2::new(self.x12, self.y12);
        let p22 = DVec2::new(self.x22, self.y22);

        // Check length of segments.
        if point_is_equal_pnt2d(p12, p11) || point_is_equal_pnt2d(p12, p22) {
            return false;
        }

        // Make vectors.
        let mut v1 = p11 - p12; // gp_Vec2d v1(p12, p11)
        let mut v2 = p22 - p12; // gp_Vec2d v2(p12, p22)
        v1 = v1.normalize();
        v2 = v2.normalize();

        // Make bisectrissa.
        let mut bisec = 0.5 * (v1 + v2);

        // Check bisectrissa.
        if bisec.length_squared() < SQUARE_CONFUSION {
            return false;
        }

        // Normalize the bisectrissa.
        bisec = bisec.normalize();

        // Angle at bisectrissa.
        let beta = gp_vec2d_angle(v1, bisec);

        // Length along the bisectrissa till the center of fillet.
        let l_len = radius / beta.abs().sin();

        // Center point of fillet.
        let pc = p12 + l_len * bisec; // p12.Translated(L * bisec)
        *xc = pc.x;
        *yc = pc.y;

        // Shrinking length along segments.
        *start = (l_len * l_len - radius * radius).sqrt();
        *end = *start;

        // Orientation of fillet.
        *cw = beta > 0.0;
        true
    }

    /// OCCT SegmentFilletArc (cxx L643-788) — a function constructs a fillet
    /// between a segment and an arc.
    #[allow(clippy::too_many_arguments)]
    fn segment_fillet_arc(
        &self,
        radius: f64,
        xc: &mut f64,
        yc: &mut f64,
        cw: &mut bool,
        start: &mut f64,
        end: &mut f64,
        xend: &mut f64,
        yend: &mut f64,
    ) -> bool {
        // Make a line parallel to the segment at the side of center point of
        // fillet. This side may be defined through making a bisectrissa for
        // vectors at p12 (or p21).

        // Make 2D points.
        let p12 = DVec2::new(self.x12, self.y12);
        let p11 = DVec2::new(self.x11, self.y11);
        let pc2 = DVec2::new(self.xc2, self.yc2);

        // Check length of segment.
        if p11.distance_squared(p12) < gp_resolution() {
            return false;
        }

        // Make 2D vectors.
        let mut v1 = p11 - p12; // gp_Vec2d v1(p12, p11)
        let mut v2 = pc2 - p12; // gp_Vec2d v2(p12, pc2)

        // Rotate the arc vector to become tangential at p21.
        if self.cw2 {
            v2 = gp_vec2d_rotate(v2, FRAC_PI_2);
        } else {
            v2 = gp_vec2d_rotate(v2, -FRAC_PI_2);
        }

        // If vectors coincide (segment and arc are tangent),
        // the algorithm doesn't work...
        let mut angle = gp_vec2d_angle(v1, v2);
        if angle.abs() < ANGULAR {
            return false;
        }

        // Make a bissectrisa of vectors at p12.
        v2 = v2.normalize();
        v1 = v1.normalize();
        let mut bisec = 0.5 * (v1 + v2);

        // If segment and arc look in opposite direction,
        // no fillet is possible.
        if bisec.length_squared() < gp_resolution() {
            return false;
        }

        // Define an appropriate point to choose center of fillet.
        bisec = bisec.normalize();
        let nearp = p12 + radius * bisec; // p12.Translated(radius * bisec)
        let nearl = Line2d::new(p12, bisec);

        // Make a line parallel to segment and
        // passing near the "near" point.
        let mut d1 = v1; // gp_Vec2d d1(v1)
        let mut line = Line2d::new(p11, -d1);
        d1 = gp_vec2d_rotate(d1, FRAC_PI_2);
        line = line.translate(radius * d1);
        if line.distance(nearp) > radius {
            line = line.translate(-2.0 * radius * d1);
        }

        // Make a circle of radius of the arc +/- fillet radius.
        // gp_Ax2d axes(pc2, gp::DX2d()); gp_Circ2d circ(axes, radius2 + radius)
        let mut circ = Circle2d::new(pc2, self.radius2 + radius);
        if self.radius2 > radius && gp_circ2d_distance(nearp, &circ) > radius {
            circ.radius = self.radius2 - radius; // circ.SetRadius(...)
        }

        // Calculate intersection of the line and the circle.
        let mut intersector = AnaIntersection2d::new();
        intersector.perform_lin_circ(&line, &circ); // IntAna2d_AnaIntersection(line, circ)
        if !intersector.is_done() || intersector.nb_points() == 0 {
            return false;
        }

        // Find center point of fillet.
        let mut min_dist = f64::MAX; // DBL_MAX
        for i in 1..=intersector.nb_points() {
            let intp = intersector.point(i);
            let p = intp.value();

            let d = nearl.distance(p);
            if d < min_dist {
                min_dist = d;
                *xc = p.x;
                *yc = p.y;
            }
        }

        // Shrink of segment.
        let pc = DVec2::new(*xc, *yc);
        let l2 = pc.distance_squared(p12);
        let rf2 = radius * radius;
        *start = (l2 - rf2).sqrt();

        // Shrink of arc.
        let pcc = pc - pc2; // gp_Vec2d pcc(pc2, pc)
        *end = gp_vec2d_angle(p12 - pc2, pcc).abs();

        // Duplicate the information on shrink the arc:
        // calculate a point on the arc coinciding with the end of fillet.
        line.origin = pc2; // line.SetLocation(pc2)
        line.direction = pcc.normalize(); // line.SetDirection(pcc)
        circ.center = pc2; // circ.SetLocation(pc2)
        circ.radius = self.radius2; // circ.SetRadius(radius2)
        intersector.perform_lin_circ(&line, &circ);
        if !intersector.is_done() || intersector.nb_points() == 0 {
            return false;
        }

        *xend = f64::MAX; // DBL_MAX
        *yend = f64::MAX;
        for i in 1..=intersector.nb_points() {
            let intp = intersector.point(i);
            let p = intp.value();

            let d2 = p.distance_squared(pc);
            if (d2 - rf2).abs() < CONFUSION {
                *xend = p.x;
                *yend = p.y;
                break;
            }
        }

        // Orientation of the fillet.
        angle = gp_vec2d_angle(v1, v2);
        *cw = angle > 0.0;
        true
    }

    /// OCCT ArcFilletSegment (cxx L791-936) — a function constructs a fillet
    /// between an arc and a segment.
    #[allow(clippy::too_many_arguments)]
    fn arc_fillet_segment(
        &self,
        radius: f64,
        xc: &mut f64,
        yc: &mut f64,
        cw: &mut bool,
        start: &mut f64,
        end: &mut f64,
        xstart: &mut f64,
        ystart: &mut f64,
    ) -> bool {
        // Make a line parallel to the segment at the side of center point of
        // fillet. This side may be defined through making a bisectrissa for
        // vectors at p12 (or p21).

        // Make 2D points.
        let p12 = DVec2::new(self.x12, self.y12);
        let p22 = DVec2::new(self.x22, self.y22);
        let pc1 = DVec2::new(self.xc1, self.yc1);

        // Check length of segment.
        if p12.distance_squared(p22) < gp_resolution() {
            return false;
        }

        // Make 2D vectors.
        let mut v1 = pc1 - p12; // gp_Vec2d v1(p12, pc1)
        let mut v2 = p22 - p12; // gp_Vec2d v2(p12, p22)

        // Rotate the arc vector to become tangential at p21.
        if self.cw1 {
            v1 = gp_vec2d_rotate(v1, -FRAC_PI_2);
        } else {
            v1 = gp_vec2d_rotate(v1, FRAC_PI_2);
        }

        // If vectors coincide (segment and arc are tangent),
        // the algorithm doesn't work...
        let mut angle = gp_vec2d_angle(v1, v2);
        if angle.abs() < ANGULAR {
            return false;
        }

        // Make a bisectrissa of vectors at p12.
        v1 = v1.normalize();
        v2 = v2.normalize();
        let mut bisec = 0.5 * (v1 + v2);

        // If segment and arc look in opposite direction,
        // no fillet is possible.
        if bisec.length_squared() < gp_resolution() {
            return false;
        }

        // Define an appropriate point to choose center of fillet.
        bisec = bisec.normalize();
        let near_point = p12 + radius * bisec; // p12.Translated(radius * bisec)
        let near_line = Line2d::new(p12, bisec);

        // Make a line parallel to segment and
        // passing near the "near" point.
        let mut a_d2_vec = v2; // gp_Vec2d aD2Vec(v2)
        let mut line = Line2d::new(p22, -a_d2_vec);
        a_d2_vec = gp_vec2d_rotate(a_d2_vec, FRAC_PI_2);
        line = line.translate(radius * a_d2_vec);
        if line.distance(near_point) > radius {
            line = line.translate(-2.0 * radius * a_d2_vec);
        }

        // Make a circle of radius of the arc +/- fillet radius.
        // gp_Ax2d axes(pc1, gp::DX2d()); gp_Circ2d circ(axes, radius1 + radius)
        let mut circ = Circle2d::new(pc1, self.radius1 + radius);
        if self.radius1 > radius && gp_circ2d_distance(near_point, &circ) > radius {
            circ.radius = self.radius1 - radius; // circ.SetRadius(...)
        }

        // Calculate intersection of the line and the big circle.
        let mut intersector = AnaIntersection2d::new();
        intersector.perform_lin_circ(&line, &circ);
        if !intersector.is_done() || intersector.nb_points() == 0 {
            return false;
        }

        // Find center point of fillet.
        let mut min_dist = f64::MAX; // DBL_MAX
        for i in 1..=intersector.nb_points() {
            let intp = intersector.point(i);
            let p = intp.value();

            let d = near_line.distance(p);
            if d < min_dist {
                min_dist = d;
                *xc = p.x;
                *yc = p.y;
            }
        }

        // Shrink of segment.
        let pc = DVec2::new(*xc, *yc);
        let l2 = pc.distance_squared(p12);
        let rf2 = radius * radius;
        *end = (l2 - rf2).sqrt();

        // Shrink of arc.
        let pcc = pc - pc1; // gp_Vec2d pcc(pc1, pc)
        *start = gp_vec2d_angle(p12 - pc1, pcc).abs();

        // Duplicate the information on shrink the arc:
        // calculate a point on the arc coinciding with the start of fillet.
        line.origin = pc1; // line.SetLocation(pc1)
        line.direction = pcc.normalize(); // line.SetDirection(pcc)
        circ.center = pc1; // circ.SetLocation(pc1)
        circ.radius = self.radius1; // circ.SetRadius(radius1)
        intersector.perform_lin_circ(&line, &circ);
        if !intersector.is_done() || intersector.nb_points() == 0 {
            return false;
        }

        *xstart = f64::MAX; // DBL_MAX
        *ystart = f64::MAX;
        for i in 1..=intersector.nb_points() {
            let intp = intersector.point(i);
            let p = intp.value();

            let d2 = p.distance_squared(pc);
            if (d2 - rf2).abs() < SQUARE_CONFUSION {
                *xstart = p.x;
                *ystart = p.y;
                break;
            }
        }

        // Orientation of the fillet.
        angle = gp_vec2d_angle(v2, v1);
        *cw = angle < 0.0;
        true
    }

    /// OCCT ArcFilletArc (cxx L944-1048) — WW5 method to compute fillet:
    /// arc - arc. It returns a constructed fillet definition: center point
    /// (xc, yc), shrinking parameter of the 1st circle (start), shrinking
    /// parameter of the 2nd circle (end), if the arc of fillet clockwise
    /// (cw = true) or counterclockwise (cw = false).
    fn arc_fillet_arc(
        &self,
        radius: f64,
        xc: &mut f64,
        yc: &mut f64,
        cw: &mut bool,
        start: &mut f64,
        end: &mut f64,
    ) -> bool {
        // Make points.
        let pc1 = DVec2::new(self.xc1, self.yc1);
        let pc2 = DVec2::new(self.xc2, self.yc2);
        let p12 = DVec2::new(self.x12, self.y12);

        // Make vectors at p12.
        let mut v1 = p12 - pc1; // gp_Vec2d v1(pc1, p12)
        let mut v2 = p12 - pc2; // gp_Vec2d v2(pc2, p12)

        // Rotate the vectors so that they are tangent to circles at p12.
        if self.cw1 {
            v1 = gp_vec2d_rotate(v1, FRAC_PI_2);
        } else {
            v1 = gp_vec2d_rotate(v1, -FRAC_PI_2);
        }
        if self.cw2 {
            v2 = gp_vec2d_rotate(v2, -FRAC_PI_2);
        } else {
            v2 = gp_vec2d_rotate(v2, FRAC_PI_2);
        }

        // Make a "check" point for choosing an offset circle.
        v1 = v1.normalize();
        v2 = v2.normalize();
        let bisec = 0.5 * (v1 + v2);
        if bisec.length_squared() < gp_resolution() {
            return false;
        }

        let checkp = p12 + radius * bisec; // p12.Translated(radius * bisec)
        let _checkl = Line2d::new(p12, bisec); // gp_Lin2d checkl(p12, bisec)

        // Make two circles of radius r1 +/- r and r2 +/- r
        // with center point equal to pc1 and pc2.
        // Arc 1.
        // gp_Ax2d axes(pc1, gp::DX2d()); gp_Circ2d c1(axes, radius1 + radius)
        let mut c1 = Circle2d::new(pc1, self.radius1 + radius);
        if self.radius1 > radius && gp_circ2d_distance(checkp, &c1) > radius {
            c1.radius = self.radius1 - radius; // c1.SetRadius(...)
        }
        // Arc 2.
        // axes.SetLocation(pc2); gp_Circ2d c2(axes, radius2 + radius)
        let mut c2 = Circle2d::new(pc2, self.radius2 + radius);
        if self.radius2 > radius && gp_circ2d_distance(checkp, &c2) > radius {
            c2.radius = self.radius2 - radius; // c2.SetRadius(...)
        }

        // Calculate an intersection point of these two circles
        // and choose the one closer to the "check" point.
        let mut intersector = AnaIntersection2d::new();
        intersector.perform_circ_circ(&c1, &c2); // IntAna2d_AnaIntersection(c1, c2)
        if !intersector.is_done() || intersector.nb_points() == 0 {
            return false;
        }

        // Find center point of fillet.
        let mut pc = DVec2::ZERO; // gp_Pnt2d pc;
        let mut min_dist = f64::MAX; // DBL_MAX
        for i in 1..=intersector.nb_points() {
            let intp = intersector.point(i);
            let p = intp.value();

            let d = checkp.distance_squared(p);
            if d < min_dist {
                min_dist = d;
                pc = p;
            }
        }
        *xc = pc.x;
        *yc = pc.y;

        // Orientation of fillet.
        let mut angle = gp_vec2d_angle(v1, v2);
        if angle.abs() < ANGULAR {
            angle = gp_vec2d_angle(pc1 - pc, pc2 - pc); // gp_Vec2d(pc, pc1).Angle(gp_Vec2d(pc, pc2))
            *cw = angle < 0.0;
        } else {
            *cw = angle > 0.0;
        }

        // Shrinking of circles.
        *start = gp_vec2d_angle(p12 - pc1, pc - pc1).abs();
        *end = gp_vec2d_angle(p12 - pc2, pc - pc2).abs();
        true
    }

    /// OCCT Cut (cxx L1051-1101) — cuts intersecting edges of a contour.
    fn cut(
        &mut self,
        the_plane: &Plane,
        the_e1: &mut Shape,
        the_e2: &mut Shape,
    ) -> bool {
        let mut p = DVec3::ZERO; // gp_Pnt p;
        let mut found = false;
        let mut param1 = 0.0;
        let mut param2 = 0.0;

        // OCCT L1057-1058: BRep_Tool::Curve(theE1, f1, l1) etc.
        let (c1, range1) = match self.my_brep.edge_curve_world(the_e1) {
            Some(x) => x,
            None => return false, // OCCT would carry a null curve handle
        };
        let f1 = range1[0];
        let l1 = range1[1];
        let (c2, range2) = match self.my_brep.edge_curve_world(the_e2) {
            Some(x) => x,
            None => return false,
        };
        let f2 = range2[0];
        let l2 = range2[1];

        // OCCT L1059: GeomAPI_ExtremaCurveCurve extrema(c1, c2, f1, l1, f2, l2)
        // Architecture stand-in: rcad extrema_curve_curve(c1, c2, samples)
        // works on the curves' natural domains; the OCCT parameter ranges are
        // applied by wrapping the curves into Geom_TrimmedCurve equivalents
        // (curve_domain reads default_domain of TrimmedCurve3).
        let tc1 = Curve3::Trimmed(TrimmedCurve3::new(c1.clone(), f1, l1));
        let tc2 = Curve3::Trimmed(TrimmedCurve3::new(c2.clone(), f2, l2));
        let extrema = extrema_curve_curve(&tc1, &tc2, 64);

        if !extrema.pairs.is_empty() {
            // extrema.NbExtrema()
            let nb = extrema.pairs.len();
            for i in 1..=nb {
                let d = extrema.pairs[i - 1].distance; // extrema.Distance(i)
                if d < CONFUSION {
                    param1 = extrema.pairs[i - 1].param1; // extrema.Parameters(i, param1, param2)
                    param2 = extrema.pairs[i - 1].param2;
                    if (l1 - param1).abs() > CONFUSION && (f2 - param2).abs() > CONFUSION {
                        found = true;
                        p = extrema.pairs[i - 1].point1; // extrema.Points(i, p, p)
                        break;
                    }
                }
            }
        }

        if found {
            // BRepBuilderAPI_MakeEdge mkEdge1(c1, f1, param1).
            let mk_edge1 = brepbuilderapi_make_edge_curve_params(&mut self.my_brep, &c1, f1, param1);
            if let Some(e1) = mk_edge1 {
                *the_e1 = e1;

                // BRepBuilderAPI_MakeEdge mkEdge2(c2, param2, l2).
                let mk_edge2 =
                    brepbuilderapi_make_edge_curve_params(&mut self.my_brep, &c2, param2, l2);
                if let Some(e2) = mk_edge2 {
                    *the_e2 = e2;

                    let p2d = projlib_project_plane_point(the_plane, p);
                    self.x12 = p2d.x;
                    self.y12 = p2d.y;
                    self.x21 = self.x12;
                    self.y21 = self.y12;
                    return true;
                }
            }
        }
        false
    }
}
