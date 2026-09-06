//! OCCT `BRepLib_MakeEdge2d` (TKTopAlgo/BRepLib) — 1:1 translation of
//! `BRepLib_MakeEdge2d.cxx` (L17-669), with:
//! - the `BRepLib_MakeEdge.hxx` L254-256 fields (`myError`, `myVertex1`,
//!   `myVertex2`) and the `BRepLib_MakeEdge()` default-ctor error init
//!   (`BRepLib_MakeEdge.cxx` L178-181),
//! - the `BRepLib_EdgeError.hxx` L20-31 error enum,
//! - the `ElCLib` helpers the file calls (`AdjustPeriodic` ElCLib.cxx L115-146,
//!   `LineParameter` L1276-1281, `CircleParameter` L1285-1292 with
//!   `normalizeAngle` L56-67 and `gp_Vec2d::Angle` gp_Vec2d.cxx L47-75),
//! - the package static plane `BRepLib::Plane()` (BRepLib.cxx L85, L131-145),
//!   living in the [`super`] module.
//!
//! The C++ constructor families over `gp_Lin2d` / `gp_Circ2d` / `gp_Elips2d` /
//! `gp_Hypr2d` / `gp_Parab2d` / `Geom2d_Curve` collapse onto the kernel
//! [`Curve2d`] enum: in Rust the caller wraps the gp-level curve into its
//! `Curve2d` variant and calls the matching `Init` overload method
//! (`init` / `init_params` / `init_points` / `init_vertices` /
//! `init_points_params` / `init` + explicit parameters).
//!
//! rcad data-model notes:
//! - OCCT builds the edge with a local `BRep_Builder` (standalone TShapes);
//!   rcad TShapes live in a [`BRep`], so the builder borrows the caller's
//!   `&mut BRep` arena (the same role the HLRBRep callers hand over).
//! - OCCT stores the 2D curve over `BRepLib::Plane()`
//!   (`B.UpdateEdge(E, C, Plane, TopLoc_Location(), preci)`); the rcad
//!   pcurve map is keyed by an owning face, so the edge carries the plane
//!   image of the 2D curve as its 3D curve — a parameterization-preserving
//!   map of every supported `Curve2d` family onto its `Curve3` counterpart
//!   (the same representation the HLRBRep consumers read back).

use crate::base::extrema::ExtPC2d;
use crate::core::precision::{
    is_negative_infinite_value, is_positive_infinite_value, CONFUSION, COMPUTATIONAL,
    INFINITE_VALUE,
};
use crate::geom::{
    Circle2d, Circle3, Curve2d, Curve2dEval, Curve3, Ellipse3, Hyperbola3, Line2d, Line3,
    Parabola3, Plane, SurfaceEval, BSplineCurve3, BezierCurve3,
};
use crate::topo::topods::{BRep, Orientation, Shape};
use glam::{DVec2, DVec3};

use super::plane;

/// OCCT `gp::Resolution()` = `RealSmall()` — the smallest positive double.
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

/// OCCT `BRepLib_EdgeError.hxx` L20-31 — errors that can occur at edge
/// construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeError {
    /// BRepLib_EdgeDone.
    EdgeDone,
    /// BRepLib_PointProjectionFailed.
    PointProjectionFailed,
    /// BRepLib_ParameterOutOfRange.
    ParameterOutOfRange,
    /// BRepLib_DifferentPointsOnClosedCurve.
    DifferentPointsOnClosedCurve,
    /// BRepLib_PointWithInfiniteParameter.
    PointWithInfiniteParameter,
    /// BRepLib_DifferentsPointAndParameter.
    DifferentsPointAndParameter,
    /// BRepLib_LineThroughIdenticPoints.
    LineThroughIdenticPoints,
}

/// OCCT `Point` static (BRepLib_MakeEdge2d.cxx L50-53) — make a 3d point on
/// the current plane: `BRepLib::Plane()->Value(P.X(), P.Y())`.
fn plane_point(p: DVec2) -> DVec3 {
    plane().point_at(p.x, p.y)
}

/// OCCT `Project(const TopoDS_Vertex& Ve)` static (cxx L60-66) — project a
/// vertex on the current plane: `ElSLib::Parameters(Plane()->Pln(), P)`.
fn project_vertex(v: &Shape) -> DVec2 {
    let p = match v.as_vertex() {
        Some(vd) => vd.point,
        None => DVec3::ZERO,
    };
    let pln: Plane = plane();
    let d = p - pln.origin;
    DVec2::new(d.dot(pln.u_dir), d.dot(pln.v_dir))
}

/// OCCT `ElCLib::LineParameter(const gp_Ax2d& L, const gp_Pnt2d& P)`
/// (ElCLib.cxx L1276-1281).
fn elclib_line_parameter(l: &Line2d, p: DVec2) -> f64 {
    (p - l.origin).dot(l.direction)
}

/// OCCT `normalizeAngle` (ElCLib.cxx L56-71) — normalize to [0, 2*PI],
/// preserving the closing seam value of exactly 2*PI.
fn normalize_angle(the_angle: &mut f64) {
    const PIPI: f64 = std::f64::consts::PI + std::f64::consts::PI;
    let negative_resolution = -COMPUTATIONAL;
    while *the_angle < negative_resolution {
        *the_angle += PIPI;
    }
    while *the_angle > PIPI * (1.0 + GP_RESOLUTION) {
        *the_angle -= PIPI;
    }
    if *the_angle < 0.0 {
        *the_angle = 0.0;
    }
}

/// OCCT `gp_Vec2d::Angle` (gp_Vec2d.cxx L47-75) — the angle from `a` to `b`
/// in (-PI, PI]. The 45-degree acos/asin precision split reduces to atan2 of
/// the same sine / cosine pair.
fn gp_vec2d_angle(a: DVec2, b: DVec2) -> f64 {
    let a_norm = a.length();
    let b_norm = b.length();
    if a_norm <= GP_RESOLUTION || b_norm <= GP_RESOLUTION {
        // OCCT throws gp_VectorWithNullMagnitude; the callers never pass a
        // null vector (the circle frame direction).
        return 0.0;
    }
    let sinus = (a.x * b.y - a.y * b.x) / (a_norm * b_norm);
    let cosinus = a.dot(b) / (a_norm * b_norm);
    sinus.atan2(cosinus)
}

/// OCCT `ElCLib::CircleParameter(const gp_Ax22d& Pos, const gp_Pnt2d& P)`
/// (ElCLib.cxx L1285-1292).
fn elclib_circle_parameter(pos: &Circle2d, p: DVec2) -> f64 {
    let v = p - pos.center;
    let mut teta = gp_vec2d_angle(pos.x_dir, v);
    // ((Pos.XDirection() ^ Pos.YDirection()) >= 0.0) ? Teta : -Teta.
    let cross = pos.x_dir.x * pos.y_dir.y - pos.x_dir.y * pos.y_dir.x;
    teta = if cross >= 0.0 { teta } else { -teta };
    normalize_angle(&mut teta);
    teta
}

/// OCCT `ElCLib::AdjustPeriodic` (ElCLib.cxx L115-146) — adjust `u1` / `u2`
/// into the period [`u_first`, `u_last`].
fn elclib_adjust_periodic(u_first: f64, u_last: f64, preci: f64, u1: &mut f64, u2: &mut f64) {
    if u_first.abs() >= 0.5 * INFINITE_VALUE || u_last.abs() >= 0.5 * INFINITE_VALUE {
        *u1 = u_first;
        *u2 = u_last;
        return;
    }

    let a_period = u_last - u_first;

    // OCCT L128-135: the Epsilon(ULast) overflow guard.
    if a_period < u_last.abs() * COMPUTATIONAL {
        *u1 = u_first;
        *u2 = u_last;
        return;
    }

    *u1 -= ((*u1 - u_first) / a_period).floor() * a_period;
    if u_last - *u1 < preci {
        *u1 -= a_period;
    }
    *u2 -= ((*u2 - *u1) / a_period).floor() * a_period;
    if *u2 - *u1 < preci {
        *u2 += a_period;
    }
}

/// OCCT `Project(const occ::handle<Geom2d_Curve>& C, const TopoDS_Vertex& V,
/// double& p)` static (cxx L70-108) — the parameter of the vertex projection
/// on the curve. Returns false when the extrema are not done
/// (`BRepLib_PointProjectionFailed`). The `p` out-parameter starts at the
/// caller-initialized value (OCCT leaves it uninitialized).
fn project_curve_vertex(c: &Curve2d, v: &Shape, p: &mut f64) -> bool {
    // gp_Pnt2d P = Project(V);
    let pp = project_vertex(v);
    match c {
        // OCCT L74-77: GeomAbs_Line.
        Curve2d::Line(l) => {
            *p = elclib_line_parameter(l, pp);
        }
        // OCCT L78-81: GeomAbs_Circle.
        Curve2d::Circle(circ) => {
            *p = elclib_circle_parameter(circ, pp);
        }
        // OCCT L83-105: Extrema_ExtPC2d over the adaptor's parameter range.
        _ => {
            let [uinf, usup] = c.default_domain();
            let mut extrema = ExtPC2d::new(pp, c, CONFUSION, uinf, usup);
            if extrema.is_done() {
                let n = extrema.nb_ext();
                let mut d2 = f64::MAX; // RealLast
                for i in 1..=n {
                    let dd2 = extrema.square_distance(i);
                    if dd2 < d2 {
                        d2 = dd2;
                        *p = extrema.point(i).param;
                    }
                }
            } else {
                return false;
            }
        }
    }
    true
}

/// The plane image of the 2D support curve — the `B.UpdateEdge(E, C,
/// BRepLib::Plane(), TopLoc_Location(), preci)` representation (see the
/// module notes). The map preserves the parameterization of every supported
/// family; the remaining `Curve2d` families have no `Curve3` counterpart.
fn plane_image(c: &Curve2d) -> Option<Curve3> {
    let pln = plane();
    let p3 = |p: DVec2| pln.point_at(p.x, p.y);
    let d3 = |d: DVec2| pln.u_dir * d.x + pln.v_dir * d.y;
    match c {
        Curve2d::Line(l) => Some(Curve3::Line(Line3::new(p3(l.origin), d3(l.direction)))),
        Curve2d::Circle(c2) => Some(Curve3::Circle(Circle3 {
            center: p3(c2.center),
            normal: pln.normal,
            x_dir: d3(c2.x_dir),
            y_dir: d3(c2.y_dir),
            radius: c2.radius,
        })),
        Curve2d::Ellipse(e2) => Some(Curve3::Ellipse(Ellipse3 {
            center: p3(e2.center),
            normal: pln.normal,
            major_dir: d3(e2.major_dir),
            major_radius: e2.major_radius,
            minor_radius: e2.minor_radius,
        })),
        Curve2d::Hyperbola(h2) => Some(Curve3::Hyperbola(Hyperbola3 {
            center: p3(h2.center),
            normal: pln.normal,
            major_dir: d3(h2.major_dir),
            semi_major: h2.semi_major,
            semi_minor: h2.semi_minor,
        })),
        Curve2d::Parabola(p2) => Some(Curve3::Parabola(Parabola3 {
            vertex: p3(p2.origin),
            normal: pln.normal,
            axis_dir: d3(p2.axis_dir),
            focal_param: p2.focal_param,
        })),
        Curve2d::Bezier(b) => Some(Curve3::Bezier(BezierCurve3 {
            control_points: b.control_points.iter().map(|p| p3(*p)).collect(),
            weights: b.weights.clone(),
        })),
        Curve2d::BSpline(b) => Some(Curve3::BSpline(BSplineCurve3 {
            degree: b.degree,
            knots: b.knots.clone(),
            control_points: b.control_points.iter().map(|p| p3(*p)).collect(),
            weights: b.weights.clone(),
            is_periodic: false,
        })),
        _ => None,
    }
}

/// OCCT `BRepLib_MakeEdge2d` — the 2D edge builder on the current plane.
pub struct MakeEdge2d<'a> {
    /// The rcad stand-in for the OCCT local `BRep_Builder` (the TShape store).
    brep: &'a mut BRep,
    /// OCCT `BRepLib_MakeEdge.hxx` L254 — `BRepLib_EdgeError myError;` with
    /// the `BRepLib_MakeEdge()` default-ctor init (MakeEdge.cxx L178-181):
    /// `BRepLib_PointProjectionFailed`. Note the OCCT flow only ever sets
    /// `myError` on failure — `Done()` does not touch it — so `error()`
    /// returns this value after a successful build too (the OCCT behavior).
    my_error: EdgeError,
    /// OCCT `myVertex1` (hxx L255).
    my_vertex1: Shape,
    /// OCCT `myVertex2` (hxx L256).
    my_vertex2: Shape,
    /// OCCT `BRepLib_MakeShape.hxx` L74 — `TopoDS_Shape myShape;`.
    my_shape: Shape,
    /// OCCT `BRepLib_Command` done flag (set by `Done()`).
    is_done: bool,
}

impl<'a> MakeEdge2d<'a> {
    /// OCCT `BRepLib_MakeEdge()` default-ctor state (MakeEdge.cxx L178-181):
    /// `myError(BRepLib_PointProjectionFailed)`.
    pub fn new(brep: &'a mut BRep) -> Self {
        MakeEdge2d {
            brep,
            my_error: EdgeError::PointProjectionFailed,
            my_vertex1: Shape::null(),
            my_vertex2: Shape::null(),
            my_shape: Shape::null(),
            is_done: false,
        }
    }

    /// OCCT ctor (cxx L112-125): `BRepLib_MakeEdge2d(V1, V2)` — the line
    /// through the two projected vertices.
    pub fn from_vertices(brep: &'a mut BRep, vv1: &Shape, vv2: &Shape) -> Self {
        let mut m = MakeEdge2d::new(brep);
        let p1 = project_vertex(vv1);
        let p2 = project_vertex(vv2);
        let l = p1.distance(p2);
        if l <= GP_RESOLUTION {
            m.my_error = EdgeError::LineThroughIdenticPoints;
            return m;
        }
        let gl = Curve2d::Line(Line2d {
            origin: p1,
            direction: (p2 - p1).normalize_or_zero(),
        });
        // Init(GL, V1, V2, 0, l).
        m.init_core(&gl, vv1, vv2, 0.0, l);
        m
    }

    /// OCCT ctor (cxx L129-140): `BRepLib_MakeEdge2d(P1, P2)` — the line
    /// through the two plane points.
    pub fn from_points(brep: &'a mut BRep, p1: DVec2, p2: DVec2) -> Self {
        let mut m = MakeEdge2d::new(brep);
        let l = p1.distance(p2);
        if l <= GP_RESOLUTION {
            m.my_error = EdgeError::LineThroughIdenticPoints;
            return m;
        }
        let gl = Curve2d::Line(Line2d {
            origin: p1,
            direction: (p2 - p1).normalize_or_zero(),
        });
        // Init(GL, P1, P2, 0, l).
        m.init_points_params(&gl, p1, p2, 0.0, l);
        m
    }

    /// OCCT `Init(C)` (cxx L370-373).
    pub fn init(&mut self, c: &Curve2d) {
        let [p1, p2] = c.default_domain();
        self.init_params(c, p1, p2);
    }

    /// OCCT `Init(C, p1, p2)` (cxx L377-383).
    pub fn init_params(&mut self, c: &Curve2d, p1: f64, p2: f64) {
        // TopoDS_Vertex V1, V2 — null.
        let v1 = Shape::null();
        let v2 = Shape::null();
        self.init_core(c, &v1, &v2, p1, p2);
    }

    /// OCCT `Init(C, P1, P2)` (cxx L387-403).
    pub fn init_points(&mut self, c: &Curve2d, p1: DVec2, p2: DVec2) {
        // BRep_Builder B;
        let v1 = self.make_vertex(plane_point(p1));
        let v2 = if p1.distance(p2) < CONFUSION {
            v1.clone()
        } else {
            self.make_vertex(plane_point(p2))
        };
        // OCCT L402: Init(C, V1, V2) — the projecting Init.
        self.init_vertices(c, &v1, &v2);
    }

    /// OCCT `Init(C, V1, V2)` (cxx L407-435) — try projecting the vertices on
    /// the curve.
    pub fn init_vertices(&mut self, c: &Curve2d, vv1: &Shape, vv2: &Shape) {
        let mut p1 = 0.0; // the OCCT uninitialized double, 0-initialized
        let mut p2 = 0.0;

        if vv1.is_null() {
            p1 = c.default_domain()[0];
        } else if !project_curve_vertex(c, vv1, &mut p1) {
            self.my_error = EdgeError::PointProjectionFailed;
            return;
        }
        if vv2.is_null() {
            p2 = c.default_domain()[1];
        } else if !project_curve_vertex(c, vv2, &mut p2) {
            self.my_error = EdgeError::PointProjectionFailed;
            return;
        }

        self.init_core(c, vv1, vv2, p1, p2);
    }

    /// OCCT `Init(C, P1, P2, p1, p2)` (cxx L439-459).
    pub fn init_points_params(&mut self, c: &Curve2d, p1: DVec2, p2: DVec2, pp1: f64, pp2: f64) {
        let v1 = self.make_vertex(plane_point(p1));
        let v2 = if p1.distance(p2) < CONFUSION {
            v1.clone()
        } else {
            self.make_vertex(plane_point(p2))
        };
        self.init_core(c, &v1, &v2, pp1, pp2);
    }

    /// OCCT `BRepLib_MakeVertex` (BRepLib_MakeVertex.cxx) — `B.MakeVertex(V,
    /// Point(P), Precision::Confusion())` through the arena.  MakeVertex is
    /// always a new TShape (the position-quantized [`BRep::add_tvertex`]
    /// registry has no OCCT counterpart and must not serve the
    /// OCCT-translated builders).
    fn make_vertex(&mut self, p: DVec3) -> Shape {
        self.brep.add_tvertex_unique(p)
    }

    /// OCCT `Init(C, VV1, VV2, pp1, pp2)` (cxx L466-632) — this one really
    /// makes the job.
    pub fn init_core(&mut self, cc: &Curve2d, vv1: &Shape, vv2: &Shape, pp1: f64, pp2: f64) {
        // kill trimmed curves (cxx L472-479).
        let mut c: Curve2d = cc.clone();
        while let Curve2d::Trimmed(ct) = c.clone() {
            c = (*ct.curve).clone();
        }

        // check parameters (cxx L482-487).
        let mut p1 = pp1;
        let mut p2 = pp2;
        let [cf, cl] = c.default_domain();
        let epsilon = CONFUSION;
        let periodic = c.is_periodic();

        let mut v1: Shape;
        let mut v2: Shape;
        if periodic {
            // adjust in period (cxx L493-495).
            elclib_adjust_periodic(cf, cl, epsilon, &mut p1, &mut p2);
            v1 = vv1.clone();
            v2 = vv2.clone();
        } else {
            // reordonate (cxx L499-512).
            if p1 < p2 {
                v1 = vv1.clone();
                v2 = vv2.clone();
            } else {
                v2 = vv1.clone();
                v1 = vv2.clone();
                let x = p1;
                p1 = p2;
                p2 = x;
            }

            // check range (cxx L515-519).
            if (cf - p1 > epsilon) || (p2 - cl > epsilon) {
                self.my_error = EdgeError::ParameterOutOfRange;
                return;
            }
        }

        // compute points on the curve (cxx L522-533).
        let p1inf = is_negative_infinite_value(p1);
        let p2inf = is_positive_infinite_value(p2);
        let mut pt1 = DVec2::ZERO;
        let mut pt2 = DVec2::ZERO;
        if !p1inf {
            pt1 = c.point_at(p1);
        }
        if !p2inf {
            pt2 = c.point_at(p2);
        }

        let preci = CONFUSION;

        // check for closed curve (cxx L539-543).
        let mut closed = false;
        if !p1inf && !p2inf {
            closed = pt1.distance(pt2) <= preci;
        }

        // check if the vertices are on the curve (cxx L546-612).
        if closed {
            if v1.is_null() && v2.is_null() {
                v1 = self.make_vertex(plane_point(pt1));
                v2 = v1.clone();
            } else if v1.is_null() {
                v1 = v2.clone();
            } else if v2.is_null() {
                v2 = v1.clone();
            } else if !v1.is_same(&v2) {
                self.my_error = EdgeError::DifferentPointsOnClosedCurve;
                return;
            } else if plane_point(pt1).distance(vertex_point(&v1)) > preci {
                self.my_error = EdgeError::DifferentPointsOnClosedCurve;
                return;
            }
        } else {
            // not closed (cxx L577-612).
            if p1inf {
                if !v1.is_null() {
                    self.my_error = EdgeError::PointWithInfiniteParameter;
                    return;
                }
            } else {
                let p = plane_point(pt1);
                if v1.is_null() {
                    v1 = self.make_vertex(p);
                }
            }

            if p2inf {
                if !v2.is_null() {
                    self.my_error = EdgeError::PointWithInfiniteParameter;
                    return;
                }
            } else {
                let p = plane_point(pt2);
                if v2.is_null() {
                    v2 = self.make_vertex(p);
                }
            }
        }

        // OCCT L614-617.
        v1.orientation = Orientation::Forward;
        v2.orientation = Orientation::Reversed;
        self.my_vertex1 = v1.clone();
        self.my_vertex2 = v2.clone();

        // OCCT L619-630: B.MakeEdge(E); B.UpdateEdge(E, C, Plane, Loc, preci);
        // B.Add(E, V1); B.Add(E, V2); B.Range(E, p1, p2) — through the arena
        // (the plane image stands in for the pcurve representation; the
        // vertex orientations are carried on the stored shapes).
        let curve3d = plane_image(&c);
        self.my_shape = self.brep.add_tedge(curve3d, v1, v2, [p1, p2]);
        // OCCT L631: Done().
        self.is_done = true;
    }

    /// OCCT `IsDone()` (BRepLib_Command).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// OCCT `Error()` (cxx L636-639).
    pub fn error(&self) -> EdgeError {
        self.my_error
    }

    /// OCCT `Edge()` (cxx L643-646) — `TopoDS::Edge(Shape())` with the
    /// NotDone check of `BRepLib_MakeShape::Shape()`.
    pub fn edge(&self) -> Shape {
        assert!(self.is_done, "StdFail_NotDone: BRepLib_MakeEdge2d::Edge");
        self.my_shape.clone()
    }

    /// OCCT `Vertex1()` (cxx L650-654) — `Check()` then the field.
    pub fn vertex1(&self) -> Shape {
        assert!(
            self.is_done,
            "StdFail_NotDone: BRepLib_MakeEdge2d::Vertex1"
        );
        self.my_vertex1.clone()
    }

    /// OCCT `Vertex2()` (cxx L658-662).
    pub fn vertex2(&self) -> Shape {
        assert!(
            self.is_done,
            "StdFail_NotDone: BRepLib_MakeEdge2d::Vertex2"
        );
        self.my_vertex2.clone()
    }
}

/// `BRep_Tool::Pnt(V)` through the vertex data.
fn vertex_point(v: &Shape) -> DVec3 {
    match v.as_vertex() {
        Some(vd) => vd.point,
        None => DVec3::ZERO,
    }
}

/// The edge data view of a built edge (the TEdgeData behind the TShape).
macro_rules! edge_data {
    ($e:expr) => {
        match &*$e.data {
            TShape::Edge(ed) => ed,
            _ => panic!("not an edge"),
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{Ellipse2d, Hyperbola2d, Parabola2d, BSplineCurve2, BezierCurve2};
    use crate::math::bspl::de_boor_2d;
    use crate::topo::topods::TShape;
    use std::f64::consts::TAU;

    /// The plane-image 3D point of a 2D plane point (z = 0 for the default
    /// XOY plane).
    fn img(p: DVec2) -> DVec3 {
        DVec3::new(p.x, p.y, 0.0)
    }

    /// OCCT anchor (cxx L129-140 + L348-355 + L439-459 + L466-632): the line
    /// through two plane points — the vertices stand at the trimmed ends, the
    /// range is [0, l], the tolerance is Precision::Confusion, and the vertex
    /// parameters match the range ends.
    #[test]
    fn make_edge2d_line_anchor() {
        let mut brep = BRep::new();
        let p1 = DVec2::new(1.0, 2.0);
        let p2 = DVec2::new(4.0, 6.0);
        let m = MakeEdge2d::from_points(&mut brep, p1, p2);
        assert!(m.is_done());
        let e = m.edge();

        // BRep_Tool::Tolerance(E) == Precision::Confusion(); BRep_Tool::Range.
        let ed = edge_data!(e);
        assert_eq!(ed.tolerance, CONFUSION);
        let l = p1.distance(p2);
        assert_eq!(ed.range[0], 0.0);
        assert!((ed.range[1] - l).abs() < 1e-12);

        // the vertices at the ends; the first is FORWARD, the last REVERSED.
        let v1 = m.vertex1();
        let v2 = m.vertex2();
        assert_eq!(v1.orientation, Orientation::Forward);
        assert_eq!(v2.orientation, Orientation::Reversed);
        assert!(vertex_point(&v1).distance(img(p1)) <= CONFUSION);
        assert!(vertex_point(&v2).distance(img(p2)) <= CONFUSION);

        // BRep_Tool::Parameter(V, E) — the range ends.
        assert_eq!(ed.vertex_params[&v1.ptr_id()], 0.0);
        assert!((ed.vertex_params[&v2.ptr_id()] - l).abs() < 1e-12);
    }

    /// OCCT anchor (cxx L186-190 + L377-383): the circle arc between two
    /// parameters — vertices on the circle, range [p1, p2].
    #[test]
    fn make_edge2d_circle_anchor() {
        let mut brep = BRep::new();
        let c2 = Circle2d {
            center: DVec2::new(1.0, -1.0),
            x_dir: DVec2::X,
            y_dir: DVec2::Y,
            radius: 2.0,
        };
        let c = Curve2d::Circle(c2);
        let (p1, p2) = (0.5, 1.5);
        let mut m = MakeEdge2d::new(&mut brep);
        m.init_params(&c, p1, p2);
        assert!(m.is_done());
        let e = m.edge();

        let ed = edge_data!(e);
        assert_eq!(ed.range, [p1, p2]);
        assert_eq!(ed.tolerance, CONFUSION);
        let v1 = m.vertex1();
        let v2 = m.vertex2();
        let q1 = DVec2::new(
            c2.center.x + c2.radius * p1.cos(),
            c2.center.y + c2.radius * p1.sin(),
        );
        let q2 = DVec2::new(
            c2.center.x + c2.radius * p2.cos(),
            c2.center.y + c2.radius * p2.sin(),
        );
        assert!(vertex_point(&v1).distance(img(q1)) <= CONFUSION);
        assert!(vertex_point(&v2).distance(img(q2)) <= CONFUSION);
        assert_eq!(ed.vertex_params[&v1.ptr_id()], p1);
        assert_eq!(ed.vertex_params[&v2.ptr_id()], p2);
    }

    /// OCCT anchor (cxx L220-224): the ellipse arc — the parameterization
    /// center + a*cos*major_dir + b*sin*minor_dir, vertices at the ends.
    #[test]
    fn make_edge2d_ellipse_anchor() {
        let mut brep = BRep::new();
        let e2 = Ellipse2d {
            center: DVec2::ZERO,
            major_dir: DVec2::X,
            minor_dir: DVec2::Y,
            major_radius: 3.0,
            minor_radius: 1.5,
        };
        let c = Curve2d::Ellipse(e2);
        let (p1, p2) = (0.25, 2.0);
        let mut m = MakeEdge2d::new(&mut brep);
        m.init_params(&c, p1, p2);
        assert!(m.is_done());
        let e = m.edge();

        let minor = DVec2::new(-e2.major_dir.y, e2.major_dir.x);
        let q = |t: f64| {
            e2.center
                + e2.major_dir * (e2.major_radius * t.cos())
                + minor * (e2.minor_radius * t.sin())
        };
        let ed = edge_data!(e);
        assert_eq!(ed.range, [p1, p2]);
        let v1 = m.vertex1();
        let v2 = m.vertex2();
        assert!(vertex_point(&v1).distance(img(q(p1))) <= CONFUSION);
        assert!(vertex_point(&v2).distance(img(q(p2))) <= CONFUSION);
        assert_eq!(ed.vertex_params[&v1.ptr_id()], p1);
        assert_eq!(ed.vertex_params[&v2.ptr_id()], p2);
    }

    /// OCCT anchor (cxx L254-258): the hyperbola arc — X(t) = a*cosh(t),
    /// Y(t) = b*sinh(t); the range check passes against the infinite domain.
    #[test]
    fn make_edge2d_hyperbola_anchor() {
        let mut brep = BRep::new();
        let h2 = Hyperbola2d {
            center: DVec2::new(1.0, 2.0),
            major_dir: DVec2::X,
            semi_major: 2.0,
            semi_minor: 1.0,
        };
        let c = Curve2d::Hyperbola(h2);
        let (p1, p2) = (-0.5, 0.8);
        let mut m = MakeEdge2d::new(&mut brep);
        m.init_params(&c, p1, p2);
        assert!(m.is_done());
        let e = m.edge();

        let minor = DVec2::new(-h2.major_dir.y, h2.major_dir.x);
        let q = |t: f64| {
            h2.center
                + h2.major_dir * (h2.semi_major * t.cosh())
                + minor * (h2.semi_minor * t.sinh())
        };
        let ed = edge_data!(e);
        assert_eq!(ed.range, [p1, p2]);
        assert!(vertex_point(&m.vertex1()).distance(img(q(p1))) <= CONFUSION);
        assert!(vertex_point(&m.vertex2()).distance(img(q(p2))) <= CONFUSION);
    }

    /// OCCT anchor (cxx L288-292): the parabola arc — X(t) = t^2/(2p),
    /// Y(t) = t relative to the apex frame.
    #[test]
    fn make_edge2d_parabola_anchor() {
        let mut brep = BRep::new();
        let p2d = Parabola2d {
            origin: DVec2::new(-1.0, 0.5),
            axis_dir: DVec2::X,
            focal_param: 2.0,
        };
        let c = Curve2d::Parabola(p2d);
        let (a1, a2) = (-1.0, 2.0);
        let mut m = MakeEdge2d::new(&mut brep);
        m.init_params(&c, a1, a2);
        assert!(m.is_done());
        let e = m.edge();

        let perp = DVec2::new(-p2d.axis_dir.y, p2d.axis_dir.x);
        let q = |t: f64| p2d.origin + (t * t / (2.0 * p2d.focal_param)) * p2d.axis_dir + t * perp;
        let ed = edge_data!(e);
        assert_eq!(ed.range, [a1, a2]);
        assert!(vertex_point(&m.vertex1()).distance(img(q(a1))) <= CONFUSION);
        assert!(vertex_point(&m.vertex2()).distance(img(q(a2))) <= CONFUSION);
    }

    /// OCCT anchor (cxx L321-326 + L370-373): the Bezier edge over the whole
    /// curve ([0, 1]); the poles image through the plane.
    #[test]
    fn make_edge2d_bezier_anchor() {
        let mut brep = BRep::new();
        let poles = [
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 2.0),
            DVec2::new(2.0, -1.0),
            DVec2::new(3.0, 1.0),
        ];
        let b = BezierCurve2 {
            control_points: poles.to_vec(),
            weights: vec![1.0; 4],
        };
        let c = Curve2d::Bezier(b);
        let mut m = MakeEdge2d::new(&mut brep);
        m.init(&c);
        assert!(m.is_done());
        let e = m.edge();

        // the De Casteljau value of the cubic at t.
        let value = |t: f64| -> DVec2 {
            let mut pts: Vec<DVec2> = poles.to_vec();
            while pts.len() > 1 {
                let mut next = Vec::with_capacity(pts.len() - 1);
                for k in 0..pts.len() - 1 {
                    next.push(pts[k] * (1.0 - t) + pts[k + 1] * t);
                }
                pts = next;
            }
            pts[0]
        };
        let ed = edge_data!(e);
        assert_eq!(ed.range, [0.0, 1.0]);
        assert!(vertex_point(&m.vertex1()).distance(img(value(0.0))) <= CONFUSION);
        assert!(vertex_point(&m.vertex2()).distance(img(value(1.0))) <= CONFUSION);
        // the plane image of the poles is the stored 3D curve.
        match &ed.curve {
            Some(Curve3::Bezier(b3)) => {
                for (i, p) in poles.iter().enumerate() {
                    assert!(b3.control_points[i].distance(img(*p)) <= CONFUSION);
                }
            }
            _ => panic!("no bezier image"),
        }
    }

    /// OCCT anchor (cxx L348-355): the BSpline edge over a sub-range of the
    /// clamped cubic domain.
    #[test]
    fn make_edge2d_bspline_anchor() {
        let mut brep = BRep::new();
        // a clamped cubic over [0, 2]: knots {0,0,0,0,2,2,2,2}, 4 poles.
        let poles = [
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 3.0),
            DVec2::new(2.0, -2.0),
            DVec2::new(4.0, 0.0),
        ];
        let knots = vec![0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2.0];
        let b = BSplineCurve2 {
            degree: 3,
            knots: knots.clone(),
            control_points: poles.to_vec(),
            weights: vec![1.0; 4],
        };
        let c = Curve2d::BSpline(b);
        let (p1, p2) = (0.25, 1.75);
        let mut m = MakeEdge2d::new(&mut brep);
        m.init_params(&c, p1, p2);
        assert!(m.is_done());
        let e = m.edge();

        let ed = edge_data!(e);
        assert_eq!(ed.range, [p1, p2]);
        assert!(vertex_point(&m.vertex1())
            .distance(img(de_boor_2d(3, &knots, &poles, &[1.0; 4], p1)))
            <= CONFUSION);
        assert!(vertex_point(&m.vertex2())
            .distance(img(de_boor_2d(3, &knots, &poles, &[1.0; 4], p2)))
            <= CONFUSION);
    }

    /// OCCT anchor (cxx L117-121): two identical plane points raise
    /// BRepLib_LineThroughIdenticPoints.
    #[test]
    fn make_edge2d_identic_points_error() {
        let mut brep = BRep::new();
        let p = DVec2::new(1.0, 2.0);
        let m = MakeEdge2d::from_points(&mut brep, p, p);
        assert!(!m.is_done());
        assert_eq!(m.error(), EdgeError::LineThroughIdenticPoints);
    }

    /// OCCT anchor (cxx L515-519): parameters outside the (finite) curve
    /// range raise BRepLib_ParameterOutOfRange.
    #[test]
    fn make_edge2d_parameter_out_of_range_error() {
        let mut brep = BRep::new();
        let b = BezierCurve2 {
            control_points: vec![DVec2::ZERO, DVec2::X, DVec2::Y, DVec2::ONE],
            weights: vec![1.0; 4],
        };
        let c = Curve2d::Bezier(b);
        let mut m = MakeEdge2d::new(&mut brep);
        m.init_params(&c, -0.5, 0.5);
        assert!(!m.is_done());
        assert_eq!(m.error(), EdgeError::ParameterOutOfRange);
    }

    /// OCCT anchor (cxx L539-573): the full circle is a closed curve — the
    /// two null vertices collapse to one shared vertex (V2 = V1) at the seam.
    #[test]
    fn make_edge2d_closed_circle_shared_vertex() {
        let mut brep = BRep::new();
        let c2 = Circle2d {
            center: DVec2::ZERO,
            x_dir: DVec2::X,
            y_dir: DVec2::Y,
            radius: 1.0,
        };
        let c = Curve2d::Circle(c2);
        let mut m = MakeEdge2d::new(&mut brep);
        m.init_params(&c, 0.0, TAU);
        assert!(m.is_done());
        let v1 = m.vertex1();
        let v2 = m.vertex2();
        assert!(v1.is_same(&v2));
        // the shared vertex sits at the seam point (radius, 0).
        assert!(vertex_point(&v1).distance(img(DVec2::new(1.0, 0.0))) <= CONFUSION);
        // the edge keeps the shared vertex at both ends.
        let e_tmp = m.edge();
        let ed = edge_data!(e_tmp);
        assert!(ed.first.is_same(&v1));
        assert!(ed.last.is_same(&v2));
    }

    /// OCCT anchor (cxx L563-572): two distinct vertices at the seam of a
    /// closed curve raise BRepLib_DifferentPointsOnClosedCurve.
    #[test]
    fn make_edge2d_different_points_on_closed_curve_error() {
        let mut brep = BRep::new();
        let c2 = Circle2d {
            center: DVec2::ZERO,
            x_dir: DVec2::X,
            y_dir: DVec2::Y,
            radius: 1.0,
        };
        // two distinct vertices within Confusion of the seam point.
        let va = brep.add_tvertex(DVec3::new(1.0, 0.0, 0.0));
        let vb = brep.add_tvertex(DVec3::new(1.0 + 5e-8, 0.0, 0.0));
        assert!(!va.is_same(&vb));
        let c = Curve2d::Circle(c2);
        let mut m = MakeEdge2d::new(&mut brep);
        m.init_core(&c, &va, &vb, 0.0, 1e-9);
        assert!(!m.is_done());
        assert_eq!(m.error(), EdgeError::DifferentPointsOnClosedCurve);
    }

    /// OCCT anchor (cxx L579-585): a vertex with an infinite parameter raises
    /// BRepLib_PointWithInfiniteParameter.
    #[test]
    fn make_edge2d_point_with_infinite_parameter_error() {
        let mut brep = BRep::new();
        let va = brep.add_tvertex(DVec3::new(1.0, 0.0, 0.0));
        let c = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
        let mut m = MakeEdge2d::new(&mut brep);
        m.init_core(&c, &va, &Shape::null(), -INFINITE_VALUE, INFINITE_VALUE);
        assert!(!m.is_done());
        assert_eq!(m.error(), EdgeError::PointWithInfiniteParameter);
    }

    /// OCCT anchor (cxx L499-512): p2 < p1 reordonates the edge — the range
    /// is [p2, p1] and the first vertex sits at the smaller parameter.
    #[test]
    fn make_edge2d_reordonate_swapped_parameters() {
        let mut brep = BRep::new();
        let b = BezierCurve2 {
            control_points: vec![DVec2::ZERO, DVec2::X, DVec2::Y, DVec2::ONE],
            weights: vec![1.0; 4],
        };
        let c = Curve2d::Bezier(b);
        let mut m = MakeEdge2d::new(&mut brep);
        m.init_params(&c, 0.8, 0.2);
        assert!(m.is_done());
        let e_tmp = m.edge();
        let ed = edge_data!(e_tmp);
        assert_eq!(ed.range, [0.2, 0.8]);
        assert_eq!(ed.vertex_params[&m.vertex1().ptr_id()], 0.2);
        assert_eq!(ed.vertex_params[&m.vertex2().ptr_id()], 0.8);
    }

    /// OCCT anchor (cxx L407-435 + L74-77): the line projection of a vertex
    /// (Init(C, V1, V2)) picks the projected parameter; the null second
    /// vertex keeps the infinite open end.
    #[test]
    fn make_edge2d_init_vertices_line_project() {
        let mut brep = BRep::new();
        let va = brep.add_tvertex(DVec3::new(2.0, 0.0, 0.0));
        let c = Curve2d::Line(Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::X,
        });
        let mut m = MakeEdge2d::new(&mut brep);
        m.init_vertices(&c, &va, &Shape::null());
        assert!(m.is_done());
        // p1 = 2 (the projection), p2 = +inf: the edge carries one vertex and
        // the [2, +inf] range.
        let e_tmp = m.edge();
        let ed = edge_data!(e_tmp);
        assert_eq!(ed.range[0], 2.0);
        assert_eq!(ed.range[1], f64::INFINITY);
    }
}

