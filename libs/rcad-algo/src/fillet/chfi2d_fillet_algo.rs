//! ChFi2d_FilletAlgo — OCCT TKFillet 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKFillet/ChFi2d/
//!   - ChFi2d_FilletAlgo.hxx (L55-228: ChFi2d_FilletAlgo + FilletPoint)
//!   - ChFi2d_FilletAlgo.cxx (L36-977)
//!
//! Algorithm (hxx L29-54): creates a fillet arc tangent to two edges; the
//! initial edges must lie on a plane and be connected by an end/start
//! point.  The search is an iterative Newton method over the parameter of
//! the first curve; the fillet function value is
//! |projection of the arc center onto curve2|² - radius².
//!
//! Architecture differences (Rust vs C++ data model):
//!   - TopoDS_Edge / TopoDS_Wire -> `rcad_kernel::topo::topo_shape::Shape`
//!     handles.  OCCT BRep_Tool::Curve reads the 3D curve + range from the
//!     edge TShape (with the shape Location applied); rcad reads the
//!     `Arc<TShape>` payload directly (curve + range), identical for the
//!     identity-location shapes this sketcher-level API operates on.
//!   - `occ::handle<Geom2d_Curve>` -> `Curve2d` value (Option for the null
//!     pre-Init state); `occ::handle<Geom_Plane>` -> `Plane` value.
//!   - BRepBuilderAPI_MakeEdge has no standalone rcad translation; the
//!     fillet / trimmed result edges are created as new TShapes appended
//!     to `my_brep` (an owned empty BRep table).
//!   - Geom2dAPI_ProjectPointOnCurve / GeomAPI_ProjectPointOnCurve are
//!     bridged below over the OCCT-aligned rcad extrema kernels
//!     (`Extrema_ExtPC2d` equivalent `ExtPC2d` and
//!     `closest_point_on_curve_range`).
//!   - OCCT `gp_Vec2d::Angle / IsOpposite` and `gp_Vec::Angle` are bridged
//!     by free functions translated from gp_Vec2d.cxx / gp_Vec.hxx.

use glam::{DVec2, DVec3};
use rcad_kernel::base::extrema::ExtPC2d;
use rcad_kernel::base::geom_api::project::closest_point_on_curve_range;
use rcad_kernel::base::geom_proj_lib::curve2d_simple;
use rcad_kernel::core::precision::{ANGULAR, CONFUSION, SQUARE_CONFUSION};
use rcad_kernel::geom::{
    Circle3, Curve2d, Curve2dEval as _, Curve3, CurveEval as _, Line2d, Plane, Surface3,
};
use rcad_kernel::math::el::{elclib_circle_value, elslib_plane_parameters, elslib_plane_value};
use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topods;

use crate::geomalgo::inter_cc::InterCurveCurve;

/// OCCT gp::Resolution() = RealSmall() = DBL_MIN (gp.hxx L60).
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

/// OCCT gp_Vec2d::Angle (gp_Vec2d.cxx L47-87) — the angular value between
/// two 2D vectors in [-PI, PI], positive in the trigonometric sense.
/// Above 45 degrees arccos gives the best precision; otherwise arcsin.
fn gp_vec2d_angle(the_v: DVec2, the_other: DVec2) -> f64 {
    let a_norm = the_v.length();
    let an_other_norm = the_other.length();
    if a_norm <= GP_RESOLUTION || an_other_norm <= GP_RESOLUTION {
        // OCCT: throw gp_VectorWithNullMagnitude();
        panic!("gp_Vec2d::Angle - vector with null magnitude");
    }

    let a_d = a_norm * an_other_norm;
    let a_cosinus = the_v.dot(the_other) / a_d;
    let a_sinus = (the_v.x * the_other.y - the_v.y * the_other.x) / a_d;

    // Use 1/sqrt(2) (approximately 0.7071067811865476) for better
    // readability and precision (OCCT: M_SQRT1_2).
    let a_cos_45_deg = std::f64::consts::FRAC_1_SQRT_2;

    if a_cosinus > -a_cos_45_deg && a_cosinus < a_cos_45_deg {
        // For angles near +/-90 degrees, use acos for better precision
        if a_sinus > 0.0 {
            a_cosinus.acos()
        } else {
            -a_cosinus.acos()
        }
    } else {
        // For angles near 0 degrees or +/-180 degrees, use asin for better precision
        if a_cosinus > 0.0 {
            a_sinus.asin()
        } else if a_sinus > 0.0 {
            std::f64::consts::PI - a_sinus.asin()
        } else {
            -std::f64::consts::PI - a_sinus.asin()
        }
    }
}

/// OCCT gp_Vec2d::IsOpposite (gp_Vec2d.hxx L364-368) —
/// `PI - abs(Angle(theOther)) <= theAngularTolerance`.
fn gp_vec2d_is_opposite(the_v: DVec2, the_other: DVec2, the_angular_tolerance: f64) -> bool {
    let an_ang = gp_vec2d_angle(the_v, the_other).abs();
    std::f64::consts::PI - an_ang <= the_angular_tolerance
}

/// OCCT gp_Vec::Angle (gp_Vec.hxx L488-494 -> gp_Dir::Angle) — the angular
/// value between two 3D vectors in [0, PI].
fn gp_vec3_angle(the_v: DVec3, the_other: DVec3) -> f64 {
    // OCCT: gp_VectorWithNullMagnitude_Raise_if(Modulus() <= gp::Resolution())
    if the_v.length() <= GP_RESOLUTION || the_other.length() <= GP_RESOLUTION {
        panic!("gp_Vec::Angle - vector with null magnitude");
    }
    the_v.angle_between(the_other)
}

/// OCCT Geom2d_Curve::Period() — architecture bridge: rcad `Curve2d` has no
/// explicit period; a periodic curve's period equals its natural domain
/// length (Geom2d_Line / Geom2d_Circle default domains).
fn geom2d_curve_period(the_curve: &Curve2d) -> f64 {
    let a_domain = the_curve.default_domain();
    a_domain[1] - a_domain[0]
}

/// OCCT BRepAdaptor_Curve::GetType() — the kind of the underlying 3D curve,
/// mapped to the GeomAbs_CurveType enumeration values (Line=0, Circle=1,
/// Ellipse=2, Hyperbola=3, Parabola=4, Bezier=5, BSpline=6, Offset=7,
/// Other=8) so that the `<` comparison in Init() keeps OCCT ordering
/// semantics.
fn brep_adaptor_curve_type(the_curve: Option<&Curve3>) -> i32 {
    match the_curve {
        Some(Curve3::Line(_)) => 0,
        Some(Curve3::Circle(_)) => 1,
        Some(Curve3::Ellipse(_)) => 2,
        Some(Curve3::Hyperbola(_)) => 3,
        Some(Curve3::Parabola(_)) => 4,
        Some(Curve3::Bezier(_)) => 5,
        Some(Curve3::BSpline(_)) => 6,
        Some(Curve3::Offset(_)) => 7,
        _ => 8, // OCCT: GeomAbs_OtherCurve
    }
}

// =========================================================================
// OCCT Geom2dAPI_ProjectPointOnCurve (Geom2dAPI_ProjectPointOnCurve.hxx /
// .cxx).  Architecture bridge: OCCT wraps Extrema_ExtPC2d; rcad has the
// OCCT-aligned `Extrema_ExtPC2d` equivalent `rcad_kernel::base::extrema::
// ExtPC2d` (extrema sorted by ascending distance), so this bridge keeps the
// Geom2dAPI call surface used by the algorithm above (NbPoints / Point /
// Parameter / Distance / LowerDistanceParameter / NearestPoint / Perform /
// Init) on top of the ExtPC2d engine.
// =========================================================================
#[derive(Debug, Clone)]
pub struct Geom2dAPIProjectPointOnCurve {
    /// OCCT myCurve — the projected curve (owned clone of the handle).
    my_curve: Option<Curve2d>,
    /// OCCT myU1, myU2 — the parameter range used by Perform.
    my_u1: f64,
    my_u2: f64,
    /// Extrema results per point: (parameter, point, distance), ordered by
    /// ascending distance (Geom2dAPI / Extrema_ExtPC2d ordering).
    my_points: Vec<(f64, DVec2, f64)>,
}

impl Geom2dAPIProjectPointOnCurve {
    /// OCCT ctor (P, C) — projection on the curve's own parametric bounds.
    pub fn new(the_p: DVec2, the_curve: &Curve2d) -> Self {
        let a_domain = the_curve.default_domain();
        Self::new_with_range(the_p, the_curve, a_domain[0], a_domain[1])
    }

    /// OCCT ctor (P, C, U1, U2) — projection restricted to [U1, U2].
    pub fn new_with_range(the_p: DVec2, the_curve: &Curve2d, the_u1: f64, the_u2: f64) -> Self {
        let mut a_proj = Geom2dAPIProjectPointOnCurve {
            my_curve: Some(the_curve.clone()),
            my_u1: the_u1,
            my_u2: the_u2,
            my_points: Vec::new(),
        };
        a_proj.perform(the_p);
        a_proj
    }

    /// OCCT Init(P, C, U1, U2) — re-initialization and immediate Perform.
    pub fn init(&mut self, the_p: DVec2, the_curve: &Curve2d, the_u1: f64, the_u2: f64) {
        self.my_curve = Some(the_curve.clone());
        self.my_u1 = the_u1;
        self.my_u2 = the_u2;
        self.perform(the_p);
    }

    /// OCCT Perform(P) — computes the distances from P to the curve on
    /// [myU1, myU2] (Extrema_ExtPC2d in OCCT).
    pub fn perform(&mut self, the_p: DVec2) {
        self.my_points.clear();
        let the_curve = self
            .my_curve
            .as_ref()
            .expect("Geom2dAPI_ProjectPointOnCurve - curve is null");
        let an_ext_pc = ExtPC2d::new(the_p, the_curve, CONFUSION, self.my_u1, self.my_u2);
        if an_ext_pc.is_done() {
            for n in 1..=an_ext_pc.nb_ext() {
                let a_pon = an_ext_pc.point(n);
                self.my_points
                    .push((a_pon.param, a_pon.point, an_ext_pc.square_distance(n).sqrt()));
            }
        }
    }

    /// OCCT NbPoints().
    pub fn nb_points(&self) -> usize {
        self.my_points.len()
    }

    /// OCCT Point(N) — 1-indexed.
    pub fn point(&self, the_n: usize) -> DVec2 {
        self.my_points
            .get(the_n - 1)
            .expect("Geom2dAPI_ProjectPointOnCurve::Point")
            .1
    }

    /// OCCT Parameter(N) — 1-indexed.
    pub fn parameter(&self, the_n: usize) -> f64 {
        self.my_points
            .get(the_n - 1)
            .expect("Geom2dAPI_ProjectPointOnCurve::Parameter")
            .0
    }

    /// OCCT Distance(N) — 1-indexed.
    pub fn distance(&self, the_n: usize) -> f64 {
        self.my_points
            .get(the_n - 1)
            .expect("Geom2dAPI_ProjectPointOnCurve::Distance")
            .2
    }

    /// OCCT LowerDistance() — the smallest distance.
    pub fn lower_distance(&self) -> f64 {
        self.my_points
            .first()
            .expect("Geom2dAPI_ProjectPointOnCurve::LowerDistance")
            .2
    }

    /// OCCT LowerDistanceParameter() — the parameter of the nearest point.
    pub fn lower_distance_parameter(&self) -> f64 {
        self.my_points
            .first()
            .expect("Geom2dAPI_ProjectPointOnCurve::LowerDistanceParameter")
            .0
    }

    /// OCCT NearestPoint().
    pub fn nearest_point(&self) -> DVec2 {
        self.my_points
            .first()
            .expect("Geom2dAPI_ProjectPointOnCurve::NearestPoint")
            .1
    }
}

// =========================================================================
// OCCT GeomAPI_ProjectPointOnCurve (3D).  Architecture bridge over the
// kernel range-restricted projection `closest_point_on_curve_range`
// (Extrema_ExtPElC analytic for Line / Circle, sampling + Newton
// otherwise).  Keeps the GeomAPI call surface used by Result():
// LowerDistanceParameter / NearestPoint / Perform.
// =========================================================================
#[derive(Debug, Clone)]
pub struct GeomAPIProjectPointOnCurve {
    /// OCCT myCurve — the projected curve (owned clone of the handle).
    my_curve: Option<Curve3>,
    /// OCCT results of the last Perform.
    my_param: f64,
    my_point: DVec3,
    my_distance: f64,
}

impl GeomAPIProjectPointOnCurve {
    /// OCCT ctor (P, C) — projection on the curve's own parametric bounds.
    pub fn new(the_p: DVec3, the_curve: &Curve3) -> Self {
        let mut a_proj = GeomAPIProjectPointOnCurve {
            my_curve: Some(the_curve.clone()),
            my_param: 0.0,
            my_point: DVec3::ZERO,
            my_distance: 0.0,
        };
        a_proj.perform(the_p);
        a_proj
    }

    /// OCCT Perform(P) — recomputes the projection.
    pub fn perform(&mut self, the_p: DVec3) {
        let the_curve = self
            .my_curve
            .as_ref()
            .expect("GeomAPI_ProjectPointOnCurve - curve is null");
        let a_domain = the_curve.default_domain();
        let a_proj = closest_point_on_curve_range(the_curve, the_p, a_domain[0], a_domain[1], 64);
        self.my_param = a_proj.param;
        self.my_point = a_proj.point;
        self.my_distance = a_proj.distance;
    }

    /// OCCT LowerDistanceParameter().
    pub fn lower_distance_parameter(&self) -> f64 {
        self.my_param
    }

    /// OCCT NearestPoint().
    pub fn nearest_point(&self) -> DVec3 {
        self.my_point
    }

    /// OCCT LowerDistance().
    pub fn lower_distance(&self) -> f64 {
        self.my_distance
    }
}

// =========================================================================
// OCCT FilletPoint (ChFi2d_FilletAlgo.hxx L138-228, .cxx L761-977) —
// Private class. Corresponds to the point on the first curve, computed
// fillet function and derivative on it.
//
// OCCT NCollection_Sequence is 1-based; the rcad Vec<T> carriers keep the
// OCCT 1-based indexing convention at the call sites (index-1 access).
// =========================================================================
#[derive(Debug, Clone)]
pub struct FilletPoint {
    /// Parameter on the first curve (start fillet point).
    my_param: f64,
    /// Parameter on the second curve (end fillet point).
    my_param2: f64,
    /// Values and derivative values of the fillet function.
    /// May be several if there are many projections on the second curve.
    my_v: Vec<f64>,
    my_d: Vec<f64>,
    /// Center of the fillet arc (OCCT gp_Pnt2d default = (0, 0)).
    my_center: DVec2,
    /// Flags for storage the validity of solutions. Indexes corresponds to
    /// indexes in sequences myV, myD.
    my_valid: Vec<bool>,
    my_near: Vec<usize>,
}

impl FilletPoint {
    /// OCCT ChFi2d_FilletAlgo.cxx L761-765 — creates a point on a first
    /// curve by parameter on this curve.
    pub fn new(the_param: f64) -> Self {
        FilletPoint {
            my_param: the_param,
            my_param2: 0.0,
            my_v: Vec::new(),
            my_d: Vec::new(),
            my_center: DVec2::ZERO,
            my_valid: Vec::new(),
            my_near: Vec::new(),
        }
    }

    /// OCCT hxx L147 — changes the point position by changing point
    /// parameter on the first curve.
    pub fn set_param(&mut self, the_param: f64) {
        self.my_param = the_param;
    }

    /// OCCT hxx L150 — returns the point parameter on the first curve.
    pub fn get_param(&self) -> f64 {
        self.my_param
    }

    /// OCCT hxx L153 — returns number of found values of function in this
    /// point.
    pub fn get_nb_values(&self) -> usize {
        self.my_v.len()
    }

    /// OCCT hxx L156 — returns value of function in this point.
    pub fn get_value(&self, the_index: usize) -> f64 {
        self.my_v[the_index - 1]
    }

    /// OCCT hxx L159 — returns derivatives of function in this point.
    pub fn get_diff(&self, the_index: usize) -> f64 {
        self.my_d[the_index - 1]
    }

    /// OCCT hxx L162 — returns true if function is valid (radiuses vectors
    /// of fillet do not intersect any curve).
    pub fn is_valid(&self, the_index: usize) -> bool {
        self.my_valid[the_index - 1]
    }

    /// OCCT hxx L165 — returns the index of the nearest value.
    pub fn get_near(&self, the_index: usize) -> usize {
        self.my_near[the_index - 1]
    }

    /// OCCT hxx L168 — defines the parameter of the projected point on the
    /// second curve.
    pub fn set_param2(&mut self, the_param2: f64) {
        self.my_param2 = the_param2;
    }

    /// OCCT hxx L171 — returns the parameter of the projected point on the
    /// second curve.
    pub fn get_param2(&self) -> f64 {
        self.my_param2
    }

    /// OCCT hxx L174 — center of the fillet.
    pub fn set_center(&mut self, the_point: DVec2) {
        self.my_center = the_point;
    }

    /// OCCT hxx L177 — center of the fillet.
    pub fn get_center(&self) -> DVec2 {
        self.my_center
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L767-781 — appends value of the function.
    pub fn append_value(&mut self, the_value: f64, the_valid: bool) {
        for a in 1..=self.my_v.len() {
            if the_value < self.my_v[a - 1] {
                // OCCT: myV.InsertBefore(a, theValue); myValid.InsertBefore(a, theValid);
                self.my_v.insert(a - 1, the_value);
                self.my_valid.insert(a - 1, the_valid);
                return;
            }
        }
        self.my_v.push(the_value);
        self.my_valid.push(the_valid);
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L783-829 — computes difference between
    /// this point and the given. Stores difference in myD.
    pub fn calculate_diff(&mut self, the_point: &FilletPoint) -> bool {
        let a_diffs_set = !self.my_d.is_empty();
        let a_dx = the_point.get_param() - self.my_param;
        let mut a_dy: f64 = 0.0;
        if the_point.my_v.len() == self.my_v.len() {
            // absolutely the same points
            for a in 1..=self.my_v.len() {
                a_dy = the_point.my_v[a - 1] - self.my_v[a - 1];
                if a_diffs_set {
                    // OCCT: myD.SetValue(a, aDY / aDX);
                    self.my_d[a - 1] = a_dy / a_dx;
                } else {
                    self.my_d.push(a_dy / a_dx);
                }
            }
            return true;
        }
        // between the diffeerent points searching for nearest analogs
        for a in 1..=self.my_v.len() {
            for b in 1..=the_point.my_v.len() {
                if b == 1 || (the_point.my_v[b - 1] - self.my_v[a - 1]).abs() < a_dy.abs() {
                    a_dy = the_point.my_v[b - 1] - self.my_v[a - 1];
                }
            }
            if a_diffs_set {
                if (a_dy / a_dx).abs() < self.my_d[a - 1].abs() {
                    self.my_d[a - 1] = a_dy / a_dx;
                }
            } else {
                self.my_d.push(a_dy / a_dx);
            }
        } // for

        false
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L831-942 — filters out the values and
    /// leaves the most optimal one.
    pub fn filter_points(&mut self, the_point: &mut FilletPoint) {
        // OCCT: NCollection_Sequence<double> aDiffs;
        let mut a_diffs: Vec<f64> = Vec::new();
        let a_dx = the_point.get_param() - self.my_param;
        let mut a: usize = 1;
        while a <= self.my_v.len() {
            // searching for near point from thePoint
            let mut a_near: usize = 0;
            let mut a_diff: f64 = a_dx * 10000.0;
            let a_y = self.my_v[a - 1] + self.my_d[a - 1] * a_dx;
            for b in 1..=the_point.my_v.len() {
                // calculate hypothesis value of the Y2 with the constant
                // first and second derivative
                let a_y2 = a_y + a_dx * (the_point.my_d[b - 1] - self.my_d[a - 1]) / 2.0;
                if a_near == 0 || (a_y2 - the_point.my_v[b - 1]).abs() < a_diff.abs() {
                    a_near = b;
                    a_diff = a_y2 - the_point.my_v[b - 1];
                }
            } // for b...

            if a_near != 0 {
                // the same sign at the same sides of the interval
                if self.my_v[a - 1] * the_point.my_v[a_near - 1] > 0.0 {
                    if self.my_v[a - 1] * self.my_d[a - 1] > 0.0 {
                        if self.my_d[a - 1].abs() > CONFUSION {
                            a_near = 0;
                        }
                    } else {
                        if self.my_v[a - 1].abs() > the_point.my_v[a_near - 1].abs() {
                            if the_point.my_v[a_near - 1] * the_point.my_d[a_near - 1] < 0.0
                                && the_point.my_d[a_near - 1].abs() > CONFUSION
                            {
                                a_near = 0;
                            }
                        }
                    }
                }
            } // if aNear

            if a_near != 0 {
                if self.my_v[a - 1] * the_point.my_v[a_near - 1] > 0.0 {
                    if (self.my_v[a - 1] + self.my_d[a - 1] * a_dx) * self.my_v[a - 1] > CONFUSION
                        && (the_point.my_v[a_near - 1] + the_point.my_d[a_near - 1] * a_dx)
                            * the_point.my_v[a_near - 1]
                            > CONFUSION
                    {
                        a_near = 0;
                    }
                }
            } // if aNear

            if a_near != 0 {
                if (a_diff / a_dx).abs() > 1.0e+7 {
                    a_near = 0;
                }
            }

            if a_near == 0 {
                // there is no near: remove it from the list
                self.my_v.remove(a - 1);
                self.my_d.remove(a - 1);
                self.my_valid.remove(a - 1);
                a -= 1;
            } else {
                let mut a_found = false;
                for b in 1..=self.my_near.len() {
                    if self.my_near[b - 1] == a_near {
                        if a_diffs[b - 1].abs() < a_diff.abs() {
                            // return this 'near'
                            a_found = true;
                            self.my_v.remove(a - 1);
                            self.my_d.remove(a - 1);
                            self.my_valid.remove(a - 1);
                            a -= 1;
                            break;
                        } else {
                            // remove the old 'near'
                            self.my_v.remove(b - 1);
                            self.my_d.remove(b - 1);
                            self.my_valid.remove(b - 1);
                            self.my_near.remove(b - 1);
                            a_diffs.remove(b - 1);
                            a -= 1;
                            break;
                        }
                    }
                } // for b...
                if !a_found {
                    self.my_near.push(a_near);
                    a_diffs.push(a_diff);
                }
            } // else
            a += 1;
        } // for a...
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L944-955 — returns a pointer to created
    /// copy of the point.
    /// warning: this is not the full copy! Copies only: myParam, myV, myD, myValid
    pub fn copy(&self) -> FilletPoint {
        let mut a_copy = FilletPoint::new(self.my_param);
        for a in 1..=self.my_v.len() {
            a_copy.my_v.push(self.my_v[a - 1]);
            a_copy.my_d.push(self.my_d[a - 1]);
            a_copy.my_valid.push(self.my_valid[a - 1]);
        }
        a_copy
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L957-969 — returns the index of the
    /// solution or zero if there is no solution.
    pub fn has_solution(&self, the_radius: f64) -> usize {
        for a in 1..=self.my_v.len() {
            // OCCT: abs(sqrt(abs(abs(myV.Value(a)) + theRadius * theRadius)) - theRadius)
            if ((self.my_v[a - 1].abs() + the_radius * the_radius).abs().sqrt() - the_radius)
                .abs()
                < CONFUSION
            {
                return a;
            }
        }
        0
    }

    /// OCCT hxx L195-209 — for debug only.
    pub fn lower_value(&self) -> f64 {
        let mut a_result_index: usize = 0;
        // OCCT reads an uninitialized aValue when myV is empty; the rcad
        // carrier initializes it to 0.0.
        let mut a_value: f64 = 0.0;
        for a in (1..=self.my_v.len()).rev() {
            if a_result_index == 0 || a_value.abs() > self.my_v[a - 1].abs() {
                a_result_index = a;
                a_value = self.my_v[a - 1];
            }
        }
        a_value
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L971-977 — removes the found value by the
    /// given index.
    pub fn remove(&mut self, the_index: usize) {
        self.my_v.remove(the_index - 1);
        self.my_d.remove(the_index - 1);
        self.my_valid.remove(the_index - 1);
        self.my_near.remove(the_index - 1);
    }
}

// =========================================================================
// OCCT ChFi2d_FilletAlgo (ChFi2d_FilletAlgo.hxx L55-136, .cxx L36-759).
// =========================================================================

/// OCCT ChFi2d_FilletAlgo — algorithm that creates fillet edge: arc tangent
/// to two edges in the start and in the end vertices. Initial edges must be
/// located on the plane and must be connected by the end or start points
/// (shared vertices are not obligatory). Created fillet arc is created with
/// the given radius, that is useful in sketcher applications.
///
/// The algorithm is iterative that allows to create fillet on any curves
/// of initial edges, that supports projection of point and C2 continuous.
/// Principles of algorithm can de reduced to the Newton method
/// (see ChFi2d_FilletAlgo.hxx L36-54).
#[derive(Debug, Clone)]
pub struct ChFi2dFilletAlgo {
    /// OCCT: TopoDS_Edge myEdge1, myEdge2 — initial edges where the fillet
    /// must be computed.
    pub my_edge1: Shape,
    pub my_edge2: Shape,
    /// OCCT: occ::handle<Geom_Plane> myPlane — plane where fillet arc must
    /// be created (None = OCCT null handle before Init).
    pub my_plane: Option<Plane>,
    /// OCCT: occ::handle<Geom2d_Curve> myCurve1, myCurve2 — underlying
    /// curves of the initial edges.
    pub my_curve1: Option<Curve2d>,
    pub my_curve2: Option<Curve2d>,
    /// OCCT: double myStart1, myEnd1, myStart2, myEnd2, myRadius — start
    /// and end parameters of curves of initial edges.
    pub my_start1: f64,
    pub my_end1: f64,
    pub my_start2: f64,
    pub my_end2: f64,
    pub my_radius: f64,
    /// OCCT: NCollection_List<double> myResultParams — list of params where
    /// roots were found.
    pub my_result_params: Vec<f64>,
    /// OCCT: NCollection_Sequence<int> myResultOrientation — sequence of 0
    /// or 1: position of the fillet relatively to the first curve.
    pub my_result_orientation: Vec<i32>,
    /// OCCT: bool myStartSide — position of the fillet relatively to the
    /// first curve.
    pub my_start_side: bool,
    /// OCCT: bool myEdgesExchnged — are initial edges where exchanged in
    /// the beginning: to make first edge more simple and minimize number of
    /// iterations.
    pub my_edges_exchnged: bool,
    /// OCCT: int myDegreeOfRecursion — number to avoid infinity recursion:
    /// indicates how deep the recursion is performed.
    pub my_degree_of_recursion: i32,
    /// rcad architecture: Result() creates the fillet arc and the trimmed
    /// edges via the BRepBuilderAPI_MakeEdge bridge; new TShapes live in a
    /// BRep table, so the algorithm owns one (OCCT: created in the global
    /// handle graph).
    pub my_brep: topods::BRep,
}

impl ChFi2dFilletAlgo {
    /// OCCT ChFi2d_FilletAlgo.cxx L36-46 — an empty constructor of the
    /// fillet algorithm. Call a method Init() to initialize the algorithm
    /// before calling of a Perform() method.
    pub fn new() -> Self {
        ChFi2dFilletAlgo {
            my_edge1: Shape::null(),
            my_edge2: Shape::null(),
            my_plane: None,
            my_curve1: None,
            my_curve2: None,
            my_start1: 0.0,
            my_end1: 0.0,
            my_start2: 0.0,
            my_end2: 0.0,
            my_radius: 0.0,
            my_result_params: Vec::new(),
            my_result_orientation: Vec::new(),
            my_start_side: false,
            my_edges_exchnged: false,
            my_degree_of_recursion: 0,
            my_brep: topods::BRep::new(),
        }
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L48-59 — a constructor of a fillet
    /// algorithm: accepts a wire consisting of two edges in a plane.
    pub fn new_wire(the_wire: &Shape, the_plane: &Plane) -> Self {
        let mut an_algo = ChFi2dFilletAlgo::new();
        an_algo.init_wire(the_wire, the_plane);
        an_algo
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L61-76 — a constructor of a fillet
    /// algorithm: accepts two edges in a plane (myEdge1 / myEdge2 set in
    /// the initializer list).
    pub fn new_edges(the_edge1: &Shape, the_edge2: &Shape, the_plane: &Plane) -> Self {
        let mut an_algo = ChFi2dFilletAlgo::new();
        an_algo.my_edge1 = the_edge1.clone();
        an_algo.my_edge2 = the_edge2.clone();
        an_algo.init_edges(the_edge1, the_edge2, the_plane);
        an_algo
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L78-103 — initializes a fillet algorithm:
    /// accepts a wire consisting of two edges in a plane.
    pub fn init_wire(&mut self, the_wire: &Shape, the_plane: &Plane) {
        let mut the_edge1 = Shape::null();
        let mut the_edge2 = Shape::null();
        // OCCT TopoDS_Iterator itr(theWire); for (; itr.More(); itr.Next())
        // Architecture bridge: the wire children are the TWireData::edges
        // list read through the Arc<TShape> payload.
        let a_children = the_wire.as_wire().expect("TopoDS::Wire").edges.clone();
        for a_shape in &a_children {
            if the_edge1.is_null() {
                // OCCT: theEdge1 = TopoDS::Edge(itr.Value());
                the_edge1 = a_shape.clone();
            } else if the_edge2.is_null() {
                // OCCT: theEdge2 = TopoDS::Edge(itr.Value());
                the_edge2 = a_shape.clone();
            } else {
                break;
            }
        }
        if the_edge1.is_null() || the_edge2.is_null() {
            // OCCT: throw Standard_ConstructionError(
            //   "The fillet algorithms expects a wire consisting of two edges.");
            panic!("The fillet algorithms expects a wire consisting of two edges.");
        }
        self.init_edges(&the_edge1, &the_edge2, the_plane);
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L105-161 — initializes a fillet algorithm:
    /// accepts two edges in a plane.
    pub fn init_edges(&mut self, the_edge1: &Shape, the_edge2: &Shape, the_plane: &Plane) {
        // OCCT L109: myPlane = new Geom_Plane(thePlane);
        self.my_plane = Some(*the_plane);

        // OCCT L111
        self.my_edges_exchnged = false;

        // OCCT L113-115: BRepAdaptor_Curve aBAC1(theEdge1); aBAC2(theEdge2);
        // if (aBAC1.GetType() < aBAC2.GetType())
        let a_bac1_type =
            brep_adaptor_curve_type(self.my_edge1.as_edge().and_then(|a_ed| a_ed.curve.as_ref()));
        let a_bac2_type =
            brep_adaptor_curve_type(self.my_edge2.as_edge().and_then(|a_ed| a_ed.curve.as_ref()));
        if a_bac1_type < a_bac2_type {
            // first curve must be more complicated
            self.my_edge1 = the_edge2.clone();
            self.my_edge2 = the_edge1.clone();
            self.my_edges_exchnged = true;
        } else {
            self.my_edge1 = the_edge1.clone();
            self.my_edge2 = the_edge2.clone();
        }

        // OCCT L127-128: BRep_Tool::Curve(myEdge1, myStart1, myEnd1);
        let a_ed1 = self
            .my_edge1
            .as_edge()
            .expect("ChFi2d_FilletAlgo::Init - myEdge1 is not an edge");
        let a_curve1 = a_ed1
            .curve
            .clone()
            .expect("ChFi2d_FilletAlgo::Init - myEdge1 has no 3D curve");
        self.my_start1 = a_ed1.range[0];
        self.my_end1 = a_ed1.range[1];
        let a_ed2 = self
            .my_edge2
            .as_edge()
            .expect("ChFi2d_FilletAlgo::Init - myEdge2 is not an edge");
        let a_curve2 = a_ed2
            .curve
            .clone()
            .expect("ChFi2d_FilletAlgo::Init - myEdge2 has no 3D curve");
        self.my_start2 = a_ed2.range[0];
        self.my_end2 = a_ed2.range[1];

        // OCCT L130-131: myCurve1 = GeomProjLib::Curve2d(aCurve1, myStart1,
        // myEnd1, myPlane);
        let a_plane = self.my_plane.expect("myPlane");
        self.my_curve1 = curve2d_simple(
            &a_curve1,
            self.my_start1,
            self.my_end1,
            &Surface3::Plane(a_plane),
        );
        self.my_curve2 = curve2d_simple(
            &a_curve2,
            self.my_start2,
            self.my_end2,
            &Surface3::Plane(a_plane),
        );

        // OCCT L133-136
        while self
            .my_curve1
            .as_ref()
            .expect("myCurve1")
            .is_periodic()
            && self.my_start1 >= self.my_end1
        {
            self.my_end1 += geom2d_curve_period(self.my_curve1.as_ref().expect("myCurve1"));
        }
        // OCCT L137-140
        while self
            .my_curve2
            .as_ref()
            .expect("myCurve2")
            .is_periodic()
            && self.my_start2 >= self.my_end2
        {
            self.my_end2 += geom2d_curve_period(self.my_curve2.as_ref().expect("myCurve2"));
        }

        // OCCT L142-160
        if a_bac1_type == a_bac2_type {
            if self.my_end2 - self.my_start2 < self.my_end1 - self.my_start1 {
                // first curve must be parametrically shorter
                let an_edge = self.my_edge1.clone();
                self.my_edge1 = self.my_edge2.clone();
                self.my_edge2 = an_edge;
                let a_curve = self.my_curve1.clone();
                self.my_curve1 = self.my_curve2.clone();
                self.my_curve2 = a_curve;
                let mut a = self.my_start1;
                self.my_start1 = self.my_start2;
                self.my_start2 = a;
                a = self.my_end1;
                self.my_end1 = self.my_end2;
                self.my_end2 = a;
                self.my_edges_exchnged = true;
            }
        }
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L350-413 — constructs a fillet edge.
    /// Returns true, if at least one result was found.
    pub fn perform(&mut self, the_radius: f64) -> bool {
        self.my_degree_of_recursion = 0;
        self.my_result_params.clear();
        self.my_result_orientation.clear();

        // OCCT L358-375: Geom2dAdaptor_Curve aGAC(myCurve1);
        // switch (aGAC.GetType()) — number of Newton steps from the curve
        // complexity.
        let a_nb_steps: f64 = match self.my_curve1.as_ref().expect("myCurve1") {
            Curve2d::Line(_) => 2.0, // OCCT: GeomAbs_Line
            Curve2d::Circle(_) => 4.0, // OCCT: GeomAbs_Circle
            Curve2d::Ellipse(_) => 5.0, // OCCT: GeomAbs_Ellipse
            Curve2d::BSpline(a_bs) => {
                // OCCT: 2 + aGAC.Degree() * aGAC.NbPoles();
                2.0 + (a_bs.degree * a_bs.control_points.len()) as f64
            }
            _ => 100.0, // unknown: maximum
        };

        self.my_radius = the_radius;
        // OCCT L379-381
        let a_step = (self.my_end1 - self.my_start1) / a_nb_steps;
        let a_d_step = 1.0e-4 * a_step;

        // OCCT L383-410:
        // for (aCycle = 2, myStartSide = false; aCycle; myStartSide = !myStartSide, aCycle--)
        let mut a_cycle: i32 = 2;
        self.my_start_side = false;
        while a_cycle != 0 {
            let mut a_left: Option<FilletPoint> = None;

            // OCCT: for (aParam = myStart1 + aStep;
            //      aParam < myEnd1 || std::abs(myEnd1 - aParam) < Precision::Confusion();
            //      aParam += aStep)
            let mut a_param = self.my_start1 + a_step;
            while a_param < self.my_end1 || (self.my_end1 - a_param).abs() < CONFUSION {
                if a_left.is_none() {
                    let mut a_left_point = FilletPoint::new(a_param - a_step);
                    self.fill_point(&mut a_left_point, a_param);
                    self.fill_diff(&mut a_left_point, a_d_step, true);
                    a_left = Some(a_left_point);
                }

                let mut a_right = FilletPoint::new(a_param);
                self.fill_point(&mut a_right, a_param - a_step);
                self.fill_diff(&mut a_right, a_d_step, false);

                // OCCT: aLeft->FilterPoints(aRight);
                let a_left_ref = a_left.as_mut().expect("aLeft");
                a_left_ref.filter_points(&mut a_right);
                self.perform_newton(a_left_ref, &mut a_right);

                // OCCT: delete aLeft; aLeft = aRight;
                a_left = Some(a_right);

                a_param += a_step;
            } // for
              // OCCT L409: delete aLeft;
            drop(a_left);

            // OCCT for-increment: myStartSide = !myStartSide, aCycle--
            self.my_start_side = !self.my_start_side;
            a_cycle -= 1;
        } // for

        !self.my_result_params.is_empty()
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L566-590 — returns number of possible
    /// solutions. <thePoint> chooses a particular fillet in case of several
    /// fillets may be constructed (for example, a circle intersecting a
    /// segment in 2 points). Put the intersecting (or common) point of the
    /// edges.
    pub fn nb_results(&mut self, the_point: DVec3) -> i32 {
        let a_plane = self.my_plane.expect("myPlane");
        // OCCT L569-572: ElSLib::PlaneParameters(myPlane->Pln().Position(),
        // thePoint, aX, aY); aTargetPoint2d.SetCoord(aX, aY);
        let (a_x, a_y) = elslib_plane_parameters(the_point, a_plane.origin, a_plane.u_dir, a_plane.v_dir);
        let a_target_point2d = DVec2::new(a_x, a_y);
        // Kept: the OCCT source computes the target point but does not use
        // it in NbResults either.
        let _ = a_target_point2d;

        // iterate through all possible solutions.
        // OCCT: NCollection_List<double>::Iterator anIter(myResultParams);
        let a_params = self.my_result_params.clone();
        let mut i: usize = 1;
        let mut nb: i32 = 0;
        for &an_iter_value in a_params.iter() {
            self.my_start_side = self.my_result_orientation[i - 1] != 0;
            let mut a_point = FilletPoint::new(an_iter_value);
            self.fill_point(&mut a_point, an_iter_value + 1.0);
            if a_point.has_solution(self.my_radius) != 0 {
                nb += 1;
            }
            // OCCT: delete aPoint;
            drop(a_point);
            i += 1;
        } // for

        nb
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L594-759 — returns result (fillet edge,
    /// modified edge1, modified edge2), nearest to the given point
    /// <thePoint>.  OCCT default argument: iSolution = -1.
    // OCCT L672-673: aPointOnCircle is reassigned (p0) and never read
    // afterward; the OCCT form is kept, the warning suppressed.
    #[allow(unused_assignments)]
    pub fn result(
        &mut self,
        the_point: DVec3,
        the_edge1: &mut Shape,
        the_edge2: &mut Shape,
        i_solution: i32,
    ) -> Shape {
        let mut a_result = Shape::null();
        let a_plane = self.my_plane.expect("myPlane");
        // OCCT L600-603: ElSLib::PlaneParameters + aTargetPoint2d.SetCoord
        let (a_x, a_y) = elslib_plane_parameters(the_point, a_plane.origin, a_plane.u_dir, a_plane.v_dir);
        let a_target_point2d = DVec2::new(a_x, a_y);

        // choose the nearest circle
        // OCCT L606-643
        let mut a_distance: f64 = 0.0;
        let mut a_nearest: Option<FilletPoint> = None;
        let mut i_sol: i32 = 1;
        let a_params = self.my_result_params.clone();
        let mut a: usize = 1;
        for &an_iter_value in a_params.iter() {
            self.my_start_side = self.my_result_orientation[a - 1] != 0;
            let mut a_point = FilletPoint::new(an_iter_value);
            self.fill_point(&mut a_point, an_iter_value + 1.0);
            if a_point.has_solution(self.my_radius) == 0 {
                // OCCT: delete aPoint; continue;
                continue;
            }
            let mut a_p = f64::MAX; // OCCT: aP = DBL_MAX;
            if i_solution == -1 {
                a_p = (a_point.get_center().distance(a_target_point2d) - self.my_radius).abs();
            } else if i_solution == i_sol {
                a_p = 0.0;
            }
            if a_nearest.is_none() || a_p < a_distance {
                a_nearest = Some(a_point);
                a_distance = a_p;
            }
            // else: OCCT delete aPoint (dropped)
            if i_solution == i_sol {
                break;
            }
            i_sol += 1;
            a += 1;
        } // for

        // OCCT L645-648
        let a_nearest = match a_nearest {
            Some(a_point) => a_point,
            None => return a_result,
        };

        // create circle edge
        // OCCT L650-655: gp_Pnt aCenter = ElSLib::PlaneValue(...);
        // occ::handle<Geom_Circle> aCircle = new Geom_Circle(gp_Ax2(aCenter,
        // myPlane->Pln().Axis().Direction()), myRadius);
        let a_nearest_center = a_nearest.get_center();
        let a_center = elslib_plane_value(
            a_nearest_center.x,
            a_nearest_center.y,
            a_plane.origin,
            a_plane.u_dir,
            a_plane.v_dir,
        );
        let a_circle = Circle3::new(a_center, a_plane.normal, self.my_radius);
        let a_circle_curve = Curve3::Circle(a_circle);

        // OCCT L656-660
        let a_curve1 = self.my_curve1.clone().expect("myCurve1");
        let a_curve2 = self.my_curve2.clone().expect("myCurve2");
        let a_point2d1 = a_curve1.point_at(a_nearest.get_param());
        let a_point2d2 = a_curve2.point_at(a_nearest.get_param2());
        let a_point1 = elslib_plane_value(
            a_point2d1.x,
            a_point2d1.y,
            a_plane.origin,
            a_plane.u_dir,
            a_plane.v_dir,
        );
        let a_point2 = elslib_plane_value(
            a_point2d2.x,
            a_point2d2.y,
            a_plane.origin,
            a_plane.u_dir,
            a_plane.v_dir,
        );

        // OCCT L662-664: GeomAPI_ProjectPointOnCurve aProj(thePoint, aCircle);
        let mut a_proj = GeomAPIProjectPointOnCurve::new(the_point, &a_circle_curve);
        let mut a_target_param = a_proj.lower_distance_parameter();
        let mut a_point_on_circle = a_proj.nearest_point();

        // There is a bug in Open CASCADE in calculation of nearest point to
        // a circle near the parameter 0.0 Therefore I check this extrema
        // point manually: (OCCT L666-673)
        // OCCT: gp_Pnt p0 = ElCLib::Value(0.0, aCircle->Circ());
        let p0 = elclib_circle_value(
            0.0,
            a_circle.center,
            a_circle.x_dir,
            a_circle.y_dir,
            a_circle.radius,
        );
        if p0.distance(the_point) < a_point_on_circle.distance(the_point) {
            a_target_param = 0.0;
            a_point_on_circle = p0;
        }

        // OCCT L675-684
        a_proj.perform(a_point1);
        let a_param1 = a_proj.lower_distance_parameter();
        a_proj.perform(a_point2);
        let a_param2 = a_proj.lower_distance_parameter();
        let mut a_is_out = (a_param1 < a_target_param && a_param2 < a_target_param)
            || (a_param1 > a_target_param && a_param2 > a_target_param);
        if a_param1 > a_param2 {
            a_is_out = !a_is_out;
        }

        // OCCT L685-688: BRepBuilderAPI_MakeEdge aBuilder(aCircle->Circ(),
        // aIsOut ? aParam2 : aParam1, aIsOut ? aParam1 : aParam2);
        // aResult = aBuilder.Edge();
        let a_builder_first = if a_is_out { a_param2 } else { a_param1 };
        let a_builder_last = if a_is_out { a_param1 } else { a_param2 };
        a_result = self.make_circle_edge(&a_circle, a_builder_first, a_builder_last);

        // divide edges
        // OCCT L690-694: aCurve = BRep_Tool::Curve(myEdge1, aStart, anEnd);
        // gp_Vec aDir; aCurve->D1(aNearest->getParam(), aPoint1, aDir);
        let a_ed1 = self
            .my_edge1
            .as_edge()
            .expect("ChFi2d_FilletAlgo::Result - myEdge1 is not an edge");
        let a_curve = a_ed1
            .curve
            .clone()
            .expect("ChFi2d_FilletAlgo::Result - myEdge1 has no 3D curve");
        let mut a_start = a_ed1.range[0];
        let mut an_end = a_ed1.range[1];
        // OCCT: aCurve->D1(aNearest->getParam(), aPoint1, aDir);
        let a_dir = a_curve.derivative_at(a_nearest.get_param());

        // OCCT L696-697: gp_Vec aCircleDir; aCircle->D1(aParam1, aPoint1, aCircleDir);
        let a_circle_dir = a_circle.derivative_at(a_param1);

        // OCCT L699-706
        let a_cond = if gp_vec3_angle(a_circle_dir, a_dir) > std::f64::consts::PI / 2.0 {
            !a_is_out
        } else {
            a_is_out
        };
        if a_cond {
            a_start = a_nearest.get_param();
        } else {
            an_end = a_nearest.get_param();
        }

        // Check the case when start and end are identical. This happens
        // when the edge decreases to size 0. Old ww5 allows such
        // cases. So we are again bug compatible (OCCT L708-714)
        if (a_start - an_end).abs() < CONFUSION {
            an_end = a_start + CONFUSION;
        }
        // Divide edge (OCCT L715-724)
        // OCCT: BRepBuilderAPI_MakeEdge aDivider1(aCurve, aStart, anEnd);
        let a_divider1 = self.make_curve_edge(&a_curve, a_start, an_end);
        if self.my_edges_exchnged {
            *the_edge2 = a_divider1;
        } else {
            *the_edge1 = a_divider1;
        }

        // OCCT L726-738
        let a_ed2 = self
            .my_edge2
            .as_edge()
            .expect("ChFi2d_FilletAlgo::Result - myEdge2 is not an edge");
        let a_curve = a_ed2
            .curve
            .clone()
            .expect("ChFi2d_FilletAlgo::Result - myEdge2 has no 3D curve");
        let mut a_start = a_ed2.range[0];
        let mut an_end = a_ed2.range[1];
        // OCCT: aCurve->D1(aNearest->getParam2(), aPoint2, aDir);
        let a_dir = a_curve.derivative_at(a_nearest.get_param2());

        // OCCT L729: aCircle->D1(aParam2, aPoint2, aCircleDir);
        let a_circle_dir = a_circle.derivative_at(a_param2);

        // OCCT L731-738
        let a_cond = if gp_vec3_angle(a_circle_dir, a_dir) > std::f64::consts::PI / 2.0 {
            a_is_out
        } else {
            !a_is_out
        };
        if a_cond {
            a_start = a_nearest.get_param2();
        } else {
            an_end = a_nearest.get_param2();
        }

        // Check the case when start and end are identical. This happens
        // when the edge decreases to size 0. Old ww5 allows such
        // cases. So we are again bug compatible (OCCT L740-746)
        if (a_start - an_end).abs() < CONFUSION {
            an_end = a_start + CONFUSION;
        }
        // OCCT L747-755
        // OCCT: BRepBuilderAPI_MakeEdge aDivider2(aCurve, aStart, anEnd);
        let a_divider2 = self.make_curve_edge(&a_curve, a_start, an_end);
        if self.my_edges_exchnged {
            *the_edge1 = a_divider2;
        } else {
            *the_edge2 = a_divider2;
        }

        // OCCT L757: delete aNearest;
        drop(a_nearest);
        a_result
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L258-334 — FillPoint: computes the value
    /// the function in the current point. <theLimit> is end parameter of
    /// the segment.
    fn fill_point(&mut self, the_point: &mut FilletPoint, the_limit: f64) {
        let a_curve1 = self.my_curve1.clone().expect("myCurve1");
        let a_curve2 = self.my_curve2.clone().expect("myCurve2");

        // on the intersection point
        let mut a_valid;
        let mut a_step = CONFUSION;
        // OCCT: gp_Pnt2d aCenter, aPoint; — default gp_Pnt2d is (0, 0)
        let mut a_center = DVec2::ZERO; // center of fillet and point on curve1
        let mut a_point = DVec2::ZERO;
        let mut a_param = the_point.get_param();
        if the_limit < a_param {
            a_step = -a_step;
        }
        // OCCT: for (aValid = false; !aValid; aParam += aStep)
        a_valid = false;
        while !a_valid {
            if (a_param - a_step - the_limit) * (a_param - the_limit) <= 0.0 {
                break; // limit was exceeded
            }
            a_step *= 2.0;
            // OCCT: myCurve1->D1(aParam, aPoint, aVec);
            a_point = a_curve1.point_at(a_param);
            let a_vec = a_curve1.derivative_at(a_param);
            if a_vec.length_squared() < CONFUSION {
                // OCCT: continue; -> the for-increment runs
                a_param += a_step;
                continue;
            }

            // OCCT L284-287: perpendicular vector to the tangent, scaled by
            // the fillet radius; the side is chosen by myStartSide.
            let a_sign = if self.my_start_side { -1.0 } else { 1.0 };
            let mut a_perp = DVec2::new(a_sign * a_vec.y, -a_sign * a_vec.x);
            a_perp = a_perp.normalize(); // OCCT: aPerp.Normalize();
            a_perp *= self.my_radius; // OCCT: aPerp.Multiply(myRadius);
            a_center = a_point + a_perp; // OCCT: aCenter = aPoint.Translated(aPerp);

            // OCCT L289-295
            let a_proj_int = Geom2dAPIProjectPointOnCurve::new_with_range(
                a_point,
                &a_curve2,
                self.my_start2,
                self.my_end2,
            );
            if a_proj_int.nb_points() == 0
                || a_point.distance(a_proj_int.nearest_point()) > CONFUSION
            {
                a_valid = true;
                break;
            }
            // for-increment
            a_param += a_step;
        }
        if a_valid {
            the_point.set_param(a_param);
            the_point.set_center(a_center);
            a_valid = !is_radius_intersected(
                &a_curve2,
                self.my_start2,
                self.my_end2,
                a_point,
                a_center,
                true,
            );
        }

        // OCCT L304: Geom2dAPI_ProjectPointOnCurve aProj(aCenter, myCurve2);
        let a_proj = Geom2dAPIProjectPointOnCurve::new(a_center, &a_curve2);
        let a_nb = a_proj.nb_points();
        for a in (1..=a_nb).rev() {
            if a_point.distance_squared(a_proj.point(a)) < CONFUSION {
                continue;
            }

            let mut a_valid2 = a_valid;
            if a_valid2 {
                a_valid2 = !is_radius_intersected(
                    &a_curve1,
                    self.my_start1,
                    self.my_end1,
                    a_center,
                    a_proj.point(a),
                    false,
                );
            }

            // checking the right parameter
            let mut a_param_proj = a_proj.parameter(a);
            while a_curve2.is_periodic() && a_param_proj < self.my_start2 {
                a_param_proj += geom2d_curve_period(&a_curve2);
            }

            let d = a_proj.distance(a);
            the_point.append_value(
                d * d - self.my_radius * self.my_radius,
                a_param_proj >= self.my_start2 && a_param_proj <= self.my_end2 && a_valid2,
            );
            if (d - self.my_radius).abs() < CONFUSION {
                the_point.set_param2(a_param_proj);
            }
        }
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L336-348 — FillDiff: computes the
    /// derivative value of the function in the current point.
    /// <theDiffStep> is small step for approximate derivative computation;
    /// <theFront> is direction of the step: from or reversed.
    fn fill_diff(&mut self, the_point: &mut FilletPoint, the_diff_step: f64, the_front: bool) {
        let a_delta = if the_front { the_diff_step } else { -the_diff_step };
        let mut a_diff = FilletPoint::new(the_point.get_param() + a_delta);
        self.fill_point(&mut a_diff, a_delta * 999.0);
        if !the_point.calculate_diff(&a_diff) {
            a_diff.set_param(the_point.get_param() - a_delta);
            self.fill_point(&mut a_diff, -a_delta * 999.0);
            the_point.calculate_diff(&a_diff);
        }
        // OCCT L347: delete aDiff;
        drop(a_diff);
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L415-463 — ProcessPoint: splits segment
    /// by the parameter and calls Newton method for both segments. It
    /// supplies recursive iterations of the Newton methods calls
    /// (PerformNewton calls this function and this calls Newton two times).
    fn process_point(
        &mut self,
        the_left: &FilletPoint,
        the_right: &mut FilletPoint,
        mut the_parameter: f64,
    ) -> bool {
        if the_parameter >= the_left.get_param() && the_parameter < the_right.get_param() {
            let a_dx = the_right.get_param() - the_left.get_param();
            if the_parameter - the_left.get_param() < a_dx / 100.0 {
                the_parameter = the_left.get_param() + a_dx / 100.0;
            }
            if the_right.get_param() - the_parameter < a_dx / 100.0 {
                the_parameter = the_right.get_param() - a_dx / 100.0;
            }

            // Protection on infinite loops.
            self.my_degree_of_recursion += 1;
            let mut diffx = 0.001 * a_dx;
            if self.my_degree_of_recursion > 100 {
                diffx *= 10.0;
                if self.my_degree_of_recursion > 1000 {
                    diffx *= 10.0;
                    if self.my_degree_of_recursion > 3000 {
                        return true;
                    }
                }
            }

            // OCCT L447-450: FilletPoint* aPoint1 = theLeft->Copy();
            // FilletPoint* aPoint2 = new FilletPoint(theParameter);
            let mut a_point1 = the_left.copy();
            let mut a_point2 = FilletPoint::new(the_parameter);
            self.fill_point(&mut a_point2, a_point1.get_param());
            self.fill_diff(&mut a_point2, diffx, true);

            // OCCT L452-455
            a_point1.filter_points(&mut a_point2);
            self.perform_newton(&mut a_point1, &mut a_point2);
            a_point2.filter_points(the_right);
            self.perform_newton(&mut a_point2, the_right);

            // OCCT L457-458: delete aPoint1; delete aPoint2;
            drop(a_point1);
            drop(a_point2);
            return true;
        }

        false
    }

    /// OCCT ChFi2d_FilletAlgo.cxx L465-564 — PerformNewton: using Newton
    /// methods computes optimal point, that can be root of the function
    /// taking into account two input points, functions value and
    /// derivatives. Performs iteration until root is found or failed to
    /// find root. Stores roots in myResultParams.
    fn perform_newton(&mut self, the_left: &mut FilletPoint, the_right: &mut FilletPoint) {
        // check the left: if this is solution store it and remove it from
        // the list of researching points of theLeft (OCCT L468-479)
        let a = the_left.has_solution(self.my_radius);
        if a != 0 {
            if the_left.is_valid(a) {
                self.my_result_params.push(the_left.get_param());
                // OCCT: myResultOrientation.Append(myStartSide); (bool -> int)
                self.my_result_orientation
                    .push(if self.my_start_side { 1 } else { 0 });
            }
            return;
        }

        let a_dx = the_right.get_param() - the_left.get_param();
        if a_dx < 1.0e-6 * CONFUSION {
            let a = the_right.has_solution(self.my_radius);
            if a != 0 && the_right.is_valid(a) {
                self.my_result_params.push(the_right.get_param());
                self.my_result_orientation
                    .push(if self.my_start_side { 1 } else { 0 });
            }
            return;
        }
        // OCCT: for (a = 1; a <= theLeft->getNBValues(); a++)
        let mut a: usize = 1;
        while a <= the_left.get_nb_values() {
            let a_near = the_left.get_near(a);

            let a_a = (the_right.get_diff(a_near) - the_left.get_diff(a)) / a_dx;
            let a_b = the_left.get_diff(a) - a_a * the_left.get_param();
            let a_c = the_left.get_value(a) - the_left.get_diff(a) * the_left.get_param()
                + a_a * the_left.get_param() * the_left.get_param() / 2.0;
            let mut a_det = a_b * a_b - 2.0 * a_a * a_c;

            if a_a.abs() < CONFUSION {
                // linear case
                if a_b.abs() > 10e-20 {
                    let a_x0 = -a_c / a_b; // use extremum
                    if a_x0 > the_left.get_param() && a_x0 < the_right.get_param() {
                        self.process_point(the_left, the_right, a_x0);
                    }
                } else {
                    // linear division otherwise
                    self.process_point(the_left, the_right, the_left.get_param() + a_dx / 2.0);
                }
            } else {
                if a_b.abs() > (a_det * 1000000.0).abs() {
                    // possible floating point operations accuracy errors
                    // linear division otherwise
                    self.process_point(the_left, the_right, the_left.get_param() + a_dx / 2.0);
                } else {
                    if a_det > 0.0 {
                        // two solutions
                        a_det = a_det.sqrt();
                        let mut a_res = self.process_point(the_left, the_right, (-a_b + a_det) / a_a);
                        if !a_res {
                            a_res = self.process_point(the_left, the_right, (-a_b - a_det) / a_a);
                        }
                        if !a_res {
                            // linear division otherwise
                            self.process_point(
                                the_left,
                                the_right,
                                the_left.get_param() + a_dx / 2.0,
                            );
                        }
                    } else {
                        let a_x0 = -a_b / a_a; // use extremum
                        if a_x0 > the_left.get_param() && a_x0 < the_right.get_param() {
                            self.process_point(the_left, the_right, a_x0);
                        } else {
                            // linear division otherwise
                            self.process_point(
                                the_left,
                                the_right,
                                the_left.get_param() + a_dx / 2.0,
                            );
                        }
                    }
                }
            }
            a += 1;
        } // for
    }

    /// OCCT BRepBuilderAPI_MakeEdge(gp_Circ, p1, p2) — architecture bridge:
    /// the fillet arc edge is a new TShape appended to `my_brep`; vertices
    /// are created at the circle points of the two parameters.
    fn make_circle_edge(&mut self, the_circ: &Circle3, the_p1: f64, the_p2: f64) -> Shape {
        let a_v1 = self.my_brep.add_tvertex(elclib_circle_value(
            the_p1,
            the_circ.center,
            the_circ.x_dir,
            the_circ.y_dir,
            the_circ.radius,
        ));
        let a_v2 = self.my_brep.add_tvertex(elclib_circle_value(
            the_p2,
            the_circ.center,
            the_circ.x_dir,
            the_circ.y_dir,
            the_circ.radius,
        ));
        self.my_brep
            .add_tedge(Some(Curve3::Circle(*the_circ)), a_v1, a_v2, [the_p1, the_p2])
    }

    /// OCCT BRepBuilderAPI_MakeEdge(Curve, p1, p2) — architecture bridge:
    /// the trimmed edge is a new TShape appended to `my_brep`; vertices are
    /// created at the curve points of the two parameters.
    fn make_curve_edge(&mut self, the_curve: &Curve3, the_p1: f64, the_p2: f64) -> Shape {
        let a_v1 = self.my_brep.add_tvertex(the_curve.point_at(the_p1));
        let a_v2 = self.my_brep.add_tvertex(the_curve.point_at(the_p2));
        self.my_brep
            .add_tedge(Some(the_curve.clone()), a_v1, a_v2, [the_p1, the_p2])
    }
}

/// OCCT ChFi2d_FilletAlgo.cxx L163-256 — static IsRadiusIntersected: this
/// function returns true if linear segment from start point of the fillet
/// arc to the end point is intersected by the first or second curve: in
/// this case fillet is invalid.
fn is_radius_intersected(
    the_curve: &Curve2d,
    the_curve_min: f64,
    the_curve_max: f64,
    the_start: DVec2,
    the_end: DVec2,
    the_start_connected: bool,
) -> bool {
    // Check the given start and end if they are identical. If yes
    // return false (OCCT L175-178)
    if the_start.distance_squared(the_end) < SQUARE_CONFUSION {
        return false;
    }
    // OCCT L179: line from theStart towards theEnd
    // occ::handle<Geom2d_Line> line = new Geom2d_Line(theStart, gp_Dir2d(gp_Vec2d(theStart, theEnd)));
    let a_line = Curve2d::Line(Line2d::new(the_start, the_end - the_start));
    // OCCT L180: Geom2dAPI_InterCurveCurve anInter(theCurve, line, Precision::Confusion());
    let an_inter = InterCurveCurve::new_curves(Some(the_curve), Some(&a_line), CONFUSION);

    // OCCT L183-207: cross-intersection points
    for a in (1..=an_inter.nb_points()).rev() {
        let a_point = an_inter.point(a);
        // OCCT L186: Geom2dAPI_ProjectPointOnCurve aProjInt(aPoint, theCurve,
        // theCurveMin, theCurveMax);
        let a_proj_int = Geom2dAPIProjectPointOnCurve::new_with_range(
            a_point,
            the_curve,
            the_curve_min,
            the_curve_max,
        );
        if a_proj_int.nb_points() < 1 || a_proj_int.lower_distance_parameter() > CONFUSION {
            continue; // point is not on edge
        }

        if a_point.distance(the_start) < CONFUSION {
            if !the_start_connected {
                return true;
            }
        }
        if a_point.distance(the_end) < CONFUSION {
            return true;
        }
        // OCCT L203: gp_Vec2d(aPoint, theStart).IsOpposite(gp_Vec2d(aPoint, theEnd),
        // Precision::Angular())
        if gp_vec2d_is_opposite(the_start - a_point, the_end - a_point, ANGULAR) {
            return true;
        }
    }

    // OCCT L208-254: tangential segments. anInter.Segment() is not
    // implemented in OCCT ("bug in OCC", L211) — the curve endpoint values
    // are probed instead.
    for _a in (1..=an_inter.nb_segments()).rev() {
        // OCCT L212: aPoint = aCurve->Value(aCurve->FirstParameter());
        let mut a_point = the_curve.point_at(the_curve.default_domain()[0]);

        let mut a_proj_int = Geom2dAPIProjectPointOnCurve::new_with_range(
            a_point,
            the_curve,
            the_curve_min,
            the_curve_max,
        );
        if a_proj_int.nb_points() > 0 && a_proj_int.lower_distance_parameter() < CONFUSION {
            // point is on edge
            if a_point.distance(the_start) < CONFUSION {
                if !the_start_connected {
                    return true;
                }
            }
            if a_point.distance(the_end) < CONFUSION {
                return true;
            }
            if gp_vec2d_is_opposite(the_start - a_point, the_end - a_point, ANGULAR) {
                return true;
            }
        }
        // OCCT L233: aPoint = aCurve->Value(aCurve->LastParameter());
        a_point = the_curve.point_at(the_curve.default_domain()[1]);

        // OCCT L235: aProjInt.Init(aPoint, theCurve, theCurveMin, theCurveMax);
        a_proj_int.init(a_point, the_curve, the_curve_min, the_curve_max);
        if a_proj_int.nb_points() > 0 && a_proj_int.lower_distance_parameter() < CONFUSION {
            // point is on edge
            if a_point.distance(the_start) < CONFUSION {
                if !the_start_connected {
                    return true;
                }
            }
            if a_point.distance(the_end) < CONFUSION {
                return true;
            }
            if gp_vec2d_is_opposite(the_start - a_point, the_end - a_point, ANGULAR) {
                return true;
            }
        }
    }
    false
}
