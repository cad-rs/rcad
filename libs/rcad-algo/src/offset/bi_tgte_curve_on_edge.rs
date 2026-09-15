// OCCT BiTgte_CurveOnEdge.cxx L1-301 (+ BiTgte_CurveOnEdge.hxx L42-119) —
// 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BiTgte/
//         BiTgte_CurveOnEdge.cxx / .hxx
//
// OCCT inheritance chain (hxx L42): BiTgte_CurveOnEdge -> Adaptor3d_Curve.
// Architecture difference: rcad has no Adaptor3d_Curve class hierarchy —
// the adaptor is carried as a plain struct with the OCCT method names; the
// methods that throw Standard_NotImplemented in OCCT keep the same panic.
// GeomAbs_CurveType is carried by the rcad proj_lib CurveType stand-in.

use glam::DVec3;
use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::geom::{Circle3, Curve3, Line3, TrimmedCurve3};
use rcad_kernel::topo_shape::Shape;

use super::brep_offset_make_simple_offset::edge_curve_of;
use super::brep_offset_offset::GeomAbsShapeKind;

use rcad_kernel::core::precision::{ANGULAR as PRECISION_ANGULAR, CONFUSION as PRECISION_CONFUSION};

// ---------------------------------------------------------------------------
// OCCT GeomAPI_ProjectPointOnCurve (TKTopAlgo/GeomAPI) — the real 1:1 body
// lives in crate::geomalgo::geom_api_project_point_on_curve (arch. diff. #12
// of brep_offset_offset.rs retired); the local GAP carrier is deleted.
// ---------------------------------------------------------------------------
pub use crate::geomalgo::geom_api_project_point_on_curve::GeomAPIProjectPointOnCurve;

/// OCCT Geom_Curve::ResD1 (Geom_Curve.hxx L62-66) — the D0/D1 result carrier.
pub struct ResD1 {
    /// OCCT: Point.
    pub point: DVec3,
    /// OCCT: D1.
    pub d1: DVec3,
}

/// OCCT Geom_Curve::ResD2 (Geom_Curve.hxx L68-72) — the D2 result carrier
/// (Point, D1, D2).
pub struct ResD2 {
    /// OCCT: Point.
    pub point: DVec3,
    /// OCCT: D1.
    pub d1: DVec3,
    /// OCCT: D2.
    pub d2: DVec3,
}

/// OCCT Geom_Curve::ResD3 (Geom_Curve.hxx L74-78) — the D3 result carrier
/// (Point, D1, D2, D3).
pub struct ResD3 {
    /// OCCT: Point.
    pub point: DVec3,
    /// OCCT: D1.
    pub d1: DVec3,
    /// OCCT: D2.
    pub d2: DVec3,
    /// OCCT: D3.
    pub d3: DVec3,
}

/// OCCT gp_Dir::Angle (gp_Dir.cxx L27-50) — the angle in [0, PI]: arccos of
/// the dot above 45 degrees, arcsin of the cross magnitude below.
fn dir_angle(coord: DVec3, other_coord: DVec3) -> f64 {
    // OCCT L36: Cosinus = coord.Dot(Other.coord).
    let cosinus = coord.dot(other_coord);
    // OCCT L37-39: above 45 degrees arccos gives the best precision.
    if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        return cosinus.acos();
    }
    // OCCT L40-49: below 45 degrees, arcsin on the crossed magnitude.
    let sinus = coord.cross(other_coord).length();
    if cosinus < 0.0 {
        std::f64::consts::PI - sinus.asin()
    } else {
        sinus.asin()
    }
}

/// OCCT gp_Dir::IsEqual (gp_Dir.hxx L164-167) — Angle(Other) <= tolerance.
fn dir_is_equal(coord: DVec3, other_coord: DVec3, angular_tolerance: f64) -> bool {
    dir_angle(coord, other_coord) <= angular_tolerance
}

/// OCCT gp_Ax1::IsCoaxial (gp_Ax1.cxx L30-44) — the coaxiality test of the
/// (location, direction) axis pair.
fn ax1_is_coaxial(
    the_loc1: DVec3,
    the_dir1: DVec3,
    the_loc2: DVec3,
    the_dir2: DVec3,
    angular_tolerance: f64,
    linear_tolerance: f64,
) -> bool {
    // OCCT L34-37: XYZ1 = (loc - Other.loc) ^ Other.vdir; D1 = |XYZ1|.
    let xyz1 = (the_loc1 - the_loc2).cross(the_dir2);
    let d1 = xyz1.length();
    // OCCT L38-41: XYZ2 = (Other.loc - loc) ^ vdir; D2 = |XYZ2|.
    let xyz2 = (the_loc2 - the_loc1).cross(the_dir1);
    let d2 = xyz2.length();
    // OCCT L42-43: vdir.IsEqual(Other.vdir, AngularTolerance)
    //              && D1 <= LinearTolerance && D2 <= LinearTolerance.
    dir_is_equal(the_dir1, the_dir2, angular_tolerance)
        && d1 <= linear_tolerance
        && d2 <= linear_tolerance
}

// ---------------------------------------------------------------------------
// OCCT class (BiTgte_CurveOnEdge.hxx L42-119, cxx L40-301).
// ---------------------------------------------------------------------------

/// OCCT BiTgte_CurveOnEdge (BiTgte_CurveOnEdge.hxx L42-119).
#[derive(Clone)]
pub struct BiTgteCurveOnEdge {
    my_edge: Shape,               // OCCT: myEdge (hxx L109)
    my_eon_f: Shape,              // OCCT: myEonF (hxx L110)
    my_curv: Option<Curve3>,      // OCCT: myCurv (hxx L111)
    my_conf: Option<Curve3>,      // OCCT: myConF (hxx L112)
    my_type: CurveType,           // OCCT: myType (hxx L113)
    my_circ: Option<Circle3>,     // OCCT: myCirc (hxx L114)
}

impl Default for BiTgteCurveOnEdge {
    fn default() -> Self {
        Self::new()
    }
}

impl BiTgteCurveOnEdge {
    /// OCCT BiTgte_CurveOnEdge::BiTgte_CurveOnEdge() (cxx L42-45).
    pub fn new() -> Self {
        BiTgteCurveOnEdge {
            my_edge: Shape::null(),
            my_eon_f: Shape::null(),
            my_curv: None,
            my_conf: None,
            my_type: CurveType::Other,
            my_circ: None,
        }
    }

    /// OCCT BiTgte_CurveOnEdge::BiTgte_CurveOnEdge(EonF, Edge) (cxx L49-55).
    pub fn with_edges(the_eon_f: &Shape, the_edge: &Shape) -> Self {
        let mut res = BiTgteCurveOnEdge {
            my_edge: the_edge.clone(),
            my_eon_f: the_eon_f.clone(),
            my_curv: None,
            my_conf: None,
            my_type: CurveType::Other,
            my_circ: None,
        };
        res.init(the_eon_f, the_edge);
        res
    }

    /// OCCT BiTgte_CurveOnEdge::ShallowCopy() (cxx L59-71) — the member-wise
    /// copy (architecture difference: value semantics, no handle indirection).
    pub fn shallow_copy(&self) -> Self {
        // OCCT L61: aCopy = new BiTgte_CurveOnEdge().
        let mut a_copy = BiTgteCurveOnEdge::new();
        // OCCT L63-68: myEdge, myEonF, myCurv, myConF, myType, myCirc.
        a_copy.my_edge = self.my_edge.clone();
        a_copy.my_eon_f = self.my_eon_f.clone();
        a_copy.my_curv = self.my_curv.clone();
        a_copy.my_conf = self.my_conf.clone();
        a_copy.my_type = self.my_type;
        a_copy.my_circ = self.my_circ.clone();
        a_copy
    }

    /// OCCT BiTgte_CurveOnEdge::Init(EonF, Edge) (cxx L75-102).
    pub fn init(&mut self, eon_f: &Shape, edge: &Shape) {
        // OCCT L80-81: myCurv = BRep_Tool::Curve(myEdge, f, l);
        // myCurv = new Geom_TrimmedCurve(myCurv, f, l).
        if let Some((c, f, l)) = edge_curve_of(edge) {
            self.my_curv = Some(Curve3::Trimmed(TrimmedCurve3 {
                curve: Box::new(c),
                first: f,
                last: l,
            }));
        }
        self.my_edge = edge.clone();

        // OCCT L83-85: myConF = BRep_Tool::Curve(myEonF, f, l);
        // myConF = new Geom_TrimmedCurve(myConF, f, l).
        if let Some((c, f, l)) = edge_curve_of(eon_f) {
            self.my_conf = Some(Curve3::Trimmed(TrimmedCurve3 {
                curve: Box::new(c),
                first: f,
                last: l,
            }));
        }
        self.my_eon_f = eon_f.clone();

        // peut on generer un cercle de rayon nul
        // OCCT L88-89: GeomAdaptor_Curve Curv(myCurv); ConF(myConF) — the
        // adaptor flattens Geom_TrimmedCurve to its basis curve.
        let curv = self.my_curv.as_ref().map(|cc| match cc {
            Curve3::Trimmed(tc) => tc.curve.as_ref(),
            _ => cc,
        });
        let conf = self.my_conf.as_ref().map(|cc| match cc {
            Curve3::Trimmed(tc) => tc.curve.as_ref(),
            _ => cc,
        });

        // OCCT L91-101.
        self.my_type = CurveType::Other;
        if let (Some(Curve3::Line(line)), Some(Curve3::Circle(circle))) = (curv, conf) {
            // OCCT L94-95: a1 = Curv.Line().Position(); a2 =
            // ConF.Circle().Axis().
            let line3: &Line3 = line;
            let a1 = (line3.origin, line3.direction);
            let a2 = (circle.center, circle.normal);
            if ax1_is_coaxial(
                a1.0,
                a1.1,
                a2.0,
                a2.1,
                PRECISION_ANGULAR,
                PRECISION_CONFUSION,
            ) {
                // OCCT L98-99: myType = GeomAbs_Circle;
                // myCirc = gp_Circ(ConF.Circle().Position(), 0.).
                self.my_type = CurveType::Circle;
                let mut circ = *circle;
                circ.radius = 0.0;
                self.my_circ = Some(circ);
            }
        }
    }

    /// OCCT BiTgte_CurveOnEdge::FirstParameter() (cxx L106-109).
    pub fn first_parameter(&self) -> f64 {
        // myConF->FirstParameter().
        match &self.my_conf {
            Some(Curve3::Trimmed(tc)) => tc.first,
            _ => panic!("BiTgte_CurveOnEdge: null myConF"),
        }
    }

    /// OCCT BiTgte_CurveOnEdge::LastParameter() (cxx L113-116).
    pub fn last_parameter(&self) -> f64 {
        // myConF->LastParameter().
        match &self.my_conf {
            Some(Curve3::Trimmed(tc)) => tc.last,
            _ => panic!("BiTgte_CurveOnEdge: null myConF"),
        }
    }

    /// OCCT BiTgte_CurveOnEdge::Continuity() (cxx L120-123).
    pub fn continuity(&self) -> GeomAbsShapeKind {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::NbIntervals(S) (cxx L127-130).
    pub fn nb_intervals(&self, _s: GeomAbsShapeKind) -> i32 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::Intervals(T, S) (cxx L134-137).
    pub fn intervals(&self, _t: &mut [f64], _s: GeomAbsShapeKind) {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::Trim(F, L, Tol) (cxx L141-146).
    pub fn trim(&self, _first: f64, _last: f64, _tol: f64) -> Curve3 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::IsClosed() (cxx L150-153).
    pub fn is_closed(&self) -> bool {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::IsPeriodic() (cxx L157-160).
    pub fn is_periodic(&self) -> bool {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::Period() (cxx L164-167).
    pub fn period(&self) -> f64 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::EvalD0(theU) (cxx L171-177).
    pub fn eval_d0(&self, the_u: f64) -> DVec3 {
        // OCCT L173-176: GeomAPI_ProjectPointOnCurve aProjector;
        // aP = myConF->Value(theU); aProjector.Init(aP, myCurv);
        // return aProjector.NearestPoint().
        let my_conf = self.my_conf.as_ref().expect("BiTgte_CurveOnEdge: null myConF");
        let my_curv = self.my_curv.as_ref().expect("BiTgte_CurveOnEdge: null myCurv");
        let a_p = use_point_at(my_conf, the_u);
        let mut a_projector = GeomAPIProjectPointOnCurve::new();
        a_projector.init_point_curve(a_p, my_curv);
        a_projector.nearest_point()
    }

    /// OCCT BiTgte_CurveOnEdge::EvalD1(theU) (cxx L181-184).
    pub fn eval_d1(&self, _the_u: f64) -> ResD1 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::EvalD2(theU) (cxx L188-191).
    pub fn eval_d2(&self, _the_u: f64) -> ResD2 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::EvalD3(theU) (cxx L195-198).
    pub fn eval_d3(&self, _the_u: f64) -> ResD3 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::EvalDN(theU, theN) (cxx L202-205).
    pub fn eval_dn(&self, _the_u: f64, _the_n: i32) -> DVec3 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::Resolution(theR3d) (cxx L209-212).
    pub fn resolution(&self, _the_r3d: f64) -> f64 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::GetType() (cxx L216-219).
    pub fn get_type(&self) -> CurveType {
        self.my_type
    }

    /// OCCT BiTgte_CurveOnEdge::Line() (cxx L223-226).
    pub fn line(&self) -> Line3 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::Circle() (cxx L230-238).
    pub fn circle(&self) -> Circle3 {
        if self.my_type != CurveType::Circle {
            // OCCT L232-235: throw Standard_NoSuchObject.
            panic!("NoSuchObject: BiTgte_CurveOnEdge::Circle");
        }
        self.my_circ.expect("BiTgte_CurveOnEdge: null myCirc")
    }

    /// OCCT BiTgte_CurveOnEdge::Ellipse() (cxx L242-245).
    pub fn ellipse(&self) -> rcad_kernel::geom::Ellipse3 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::Hyperbola() (cxx L249-252).
    pub fn hyperbola(&self) -> rcad_kernel::geom::Hyperbola3 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::Parabola() (cxx L256-259).
    pub fn parabola(&self) -> rcad_kernel::geom::Parabola3 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::Degree() (cxx L263-266).
    pub fn degree(&self) -> i32 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::IsRational() (cxx L270-273).
    pub fn is_rational(&self) -> bool {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::NbPoles() (cxx L277-280).
    pub fn nb_poles(&self) -> i32 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::NbKnots() (cxx L284-287).
    pub fn nb_knots(&self) -> i32 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::Bezier() (cxx L291-294).
    pub fn bezier(&self) -> rcad_kernel::geom::BezierCurve3 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }

    /// OCCT BiTgte_CurveOnEdge::BSpline() (cxx L298-301).
    pub fn bspline(&self) -> rcad_kernel::geom::BSplineCurve3 {
        panic!("NotImplemented: BiTgte_CurveOnEdge");
    }
}

/// OCCT myConF->Value(theU) — the 3d curve evaluation via the rcad
/// CurveEval trait (the trait-object dispatch of Geom_Curve::Value).
fn use_point_at(c: &Curve3, the_u: f64) -> DVec3 {
    use rcad_kernel::geom::CurveEval;
    c.point_at(the_u)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use rcad_kernel::topods::{tshape_flags, Orientation, TEdgeData, TShape};

    /// An edge TShape carrying (curve, range) — the BRep_Tool::Curve /
    /// BRep_Tool::Range stand-in input of edge_curve_of.
    fn edge_shape(curve: Option<Curve3>, range: [f64; 2]) -> Shape {
        Shape {
            data: Arc::new(TShape::Edge(TEdgeData {
                my_shapes: Vec::new(),
                flags: tshape_flags::DEFAULT,
                curve,
                first: Shape::null(),
                last: Shape::null(),
                range,
                degenerated: false,
                pcurves: indexmap::IndexMap::new(),
                representations: Vec::new(),
                vertex_params: std::collections::HashMap::new(),
                tolerance: 0.0,
                same_parameter: false,
                same_range: false,
            })),
            index: usize::MAX,
            location: 0,
            orientation: Orientation::Forward,
        }
    }

    /// OCCT gp_Lin(gp_Pnt(5, 0, 0), gp_Dir(0, 0, 1)) — the X = 5 Z-axis line
    /// edge on [0, 10].
    fn line_edge() -> Shape {
        edge_shape(
            Some(Curve3::Line(Line3::new(
                DVec3::new(5.0, 0.0, 0.0),
                DVec3::new(0.0, 0.0, 1.0),
            ))),
            [0.0, 10.0],
        )
    }

    /// OCCT gp_Circ(gp_Ax2(gp_Pnt(5, 0, 10), N), 2.5) — the blend circle
    /// edge on [0, 2*PI].
    fn circle_edge(normal: DVec3, center: DVec3) -> Shape {
        edge_shape(
            Some(Curve3::Circle(Circle3::new(center, normal, 2.5))),
            [0.0, std::f64::consts::TAU],
        )
    }

    #[test]
    fn coaxial_line_and_circle_yields_zero_radius_circle() {
        // OCCT Init L91-101: GeomAbs_Line on a Z axis + GeomAbs_Circle around
        // the same axis -> myType = GeomAbs_Circle,
        // myCirc = gp_Circ(ConF.Circle().Position(), 0.).
        let c =
            BiTgteCurveOnEdge::with_edges(&circle_edge(DVec3::Z, DVec3::new(5.0, 0.0, 10.0)), &line_edge());
        assert_eq!(c.get_type(), CurveType::Circle);
        let circ = c.circle();
        assert_eq!(circ.radius, 0.0);
        assert_eq!(circ.center, DVec3::new(5.0, 0.0, 10.0));
        assert_eq!(circ.normal, DVec3::Z);
        // OCCT L106-116: FirstParameter/LastParameter delegate to myConF.
        assert_eq!(c.first_parameter(), 0.0);
        assert_eq!(c.last_parameter(), std::f64::consts::TAU);
        // OCCT EvalD0 L171-177: project myConF->Value(u) onto myCurv — the
        // circle point at u = 0 is (7.5, 0, 10); its projection on the axis
        // line is the (5, 0, 10) point of the circle plane.
        let p = c.eval_d0(0.0);
        assert!((p - DVec3::new(5.0, 0.0, 10.0)).length() < 1e-12);
    }

    #[test]
    fn shallow_copy_preserves_the_members() {
        // OCCT ShallowCopy (cxx L59-71): the member-wise copy.
        let c =
            BiTgteCurveOnEdge::with_edges(&circle_edge(DVec3::Z, DVec3::new(5.0, 0.0, 10.0)), &line_edge());
        let copy = c.shallow_copy();
        assert_eq!(copy.get_type(), c.get_type());
        assert_eq!(copy.first_parameter(), c.first_parameter());
        assert_eq!(copy.last_parameter(), c.last_parameter());
        let a = c.circle();
        let b = copy.circle();
        assert_eq!(a.radius, b.radius);
        assert_eq!(a.center, b.center);
        assert_eq!(a.normal, b.normal);
    }

    #[test]
    fn anti_parallel_axis_is_not_coaxial() {
        // OCCT gp_Dir::IsEqual (gp_Dir.hxx L164) rejects the PI angle; the
        // former |v1 x v2| stand-in wrongly accepted the anti-parallel axis
        // (zero cross product).
        let c = BiTgteCurveOnEdge::with_edges(
            &circle_edge(-DVec3::Z, DVec3::new(5.0, 0.0, 10.0)),
            &line_edge(),
        );
        assert_eq!(c.get_type(), CurveType::Other);
    }

    #[test]
    fn offset_axis_is_not_coaxial() {
        // OCCT gp_Ax1.cxx L42-43: the crossed momentum D1 = 5 exceeds
        // Precision::Confusion.
        let c = BiTgteCurveOnEdge::with_edges(
            &circle_edge(DVec3::Z, DVec3::new(0.0, 0.0, 10.0)),
            &line_edge(),
        );
        assert_eq!(c.get_type(), CurveType::Other);
    }

    #[test]
    fn dir_angle_matches_the_three_branch_form() {
        // OCCT gp_Dir.cxx L27-50: 60 degrees -> the arccos branch.
        let a = DVec3::new(1.0, 0.0, 0.0);
        let b = DVec3::new(0.5, 0.5_f64.sqrt(), 0.0);
        assert!((dir_angle(a, b) - std::f64::consts::FRAC_PI_3).abs() < 1e-15);
        // Anti-parallel: PI (the asin branch with cosinus < 0).
        assert!((dir_angle(a, -a) - std::f64::consts::PI).abs() < 1e-15);
        assert!(dir_is_equal(a, a, PRECISION_ANGULAR));
        assert!(!dir_is_equal(a, -a, PRECISION_ANGULAR));
    }
}
