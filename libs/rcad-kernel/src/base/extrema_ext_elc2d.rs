//! OCCT Extrema_ExtElC2d (TKGeomBase/Extrema/Extrema_ExtElC2d.hxx L35-99 and
//! Extrema_ExtElC2d.cxx L32-543) — all the extremum distances between two
//! elementary 2D curves (line-line, line-circle, line-ellipse, line-hyperbola,
//! line-parabola, circle-circle, circle-ellipse, circle-hyperbola,
//! circle-parabola).
//!
//! The circle-ellipse / circle-hyperbola / circle-parabola constructors run
//! `Extrema_ExtPElC2d`, translated in [`super::extrema_ext_p_elc2d`].
//!
//! Also carries the ElCLib 2D evaluation helpers the two 2D Extrema bodies
//! consume (`ElCLib::Value` over the gp_Lin2d / gp_Circ2d / gp_Elips2d /
//! gp_Hypr2d / gp_Parab2d position payloads, ElCLib.cxx L521-591) plus the
//! `gp_Dir2d::Angle` / `gp_Dir2d::IsParallel` / `gp_Lin2d::SquareDistance`
//! primitives (gp_Dir2d.cxx L26-56, gp_Dir2d.hxx L422-436, gp_Lin2d.hxx
//! L240-247).
//!
//! gp payload mapping (the rcad `geom` 2D types): gp_Lin2d -> [`Line2d`]
//! (origin + unit direction); gp_Circ2d -> [`Circle2d`] (center + x_dir +
//! y_dir + radius); gp_Elips2d -> [`Ellipse2d`] (major_dir = XDirection,
//! minor_dir = YDirection); gp_Hypr2d -> [`Hyperbola2d`] (semi_major =
//! MajorRadius, semi_minor = MinorRadius, y_dir = the kernel 90-degree
//! turn of major_dir, the same convention the geom evaluator uses);
//! gp_Parab2d -> [`Parabola2d`] (axis_dir = the X direction of the OCCT
//! position, y_dir = the kernel 90-degree turn, focal_param = the OCCT
//! parameter `2 * Focal`, so OCCT `Focal()` reads `focal_param / 2`).

use glam::DVec2;

use crate::base::extrema::POnCurve2d;
use crate::base::extrema_ext_elc::{elclib_circle2d_parameter, REAL_EPSILON, REAL_FIRST, REAL_LAST};
use crate::base::extrema_ext_p_elc2d::ExtremaExtPElC2d;
use crate::core::precision::{ANGULAR, CONFUSION, SQUARE_CONFUSION};
use crate::geom::{Circle2d, Ellipse2d, Hyperbola2d, Line2d, Parabola2d};

// =============================================================================
// ElCLib 2D evaluation (ElCLib.cxx L521-591) and gp 2D primitives
// =============================================================================

/// OCCT ElCLib::LineValue(U, Pos) / ElCLib::Value(U, gp_Lin2d)
/// (ElCLib.cxx L521-528).
pub(crate) fn elclib2d_line_value(u: f64, lin: &Line2d) -> DVec2 {
    lin.origin + u * lin.direction
}

/// OCCT ElCLib::CircleValue(U, Pos, Radius) / ElCLib::Value(U, gp_Circ2d)
/// (ElCLib.cxx L530-541).
pub(crate) fn elclib2d_circle_value(u: f64, circ: &Circle2d) -> DVec2 {
    circ.center + circ.radius * (u.cos() * circ.x_dir + u.sin() * circ.y_dir)
}

/// OCCT ElCLib::EllipseValue(U, Pos, MajorRadius, MinorRadius) /
/// ElCLib::Value(U, gp_Elips2d) (ElCLib.cxx L543-557).
pub(crate) fn elclib2d_ellipse_value(u: f64, elips: &Ellipse2d) -> DVec2 {
    elips.center
        + elips.major_radius * u.cos() * elips.major_dir
        + elips.minor_radius * u.sin() * elips.minor_dir
}

/// The Y direction of the gp_Hypr2d position — the kernel 2D curves store no
/// explicit second direction; the geom evaluator turns major_dir by 90
/// degrees (geom/eval.rs `Conic2dEval for Hyperbola2d`).
pub(crate) fn hyperbola2d_y_dir(hypr: &Hyperbola2d) -> DVec2 {
    crate::geom::turn_2d(hypr.major_dir)
}

/// OCCT ElCLib::HyperbolaValue(U, Pos, MajorRadius, MinorRadius) /
/// ElCLib::Value(U, gp_Hypr2d) (ElCLib.cxx L559-573).
pub(crate) fn elclib2d_hyperbola_value(u: f64, hypr: &Hyperbola2d) -> DVec2 {
    hypr.center
        + hypr.semi_major * u.cosh() * hypr.major_dir
        + hypr.semi_minor * u.sinh() * hyperbola2d_y_dir(hypr)
}

/// The Y direction of the gp_Parab2d position — the kernel 2D curves store no
/// explicit second direction; the geom evaluator turns axis_dir by 90 degrees
/// (geom/eval.rs `Conic2dEval for Parabola2d`).
pub(crate) fn parabola2d_y_dir(parab: &Parabola2d) -> DVec2 {
    crate::geom::turn_2d(parab.axis_dir)
}

/// OCCT ElCLib::ParabolaValue(U, Pos, Focal) / ElCLib::Value(U, gp_Parab2d)
/// (ElCLib.cxx L575-591): `A1 = U * U / (4.0 * Focal)`.
///
/// The rcad [`Parabola2d`] carries the OCCT parameter `p = 2 * Focal`
/// (geom `focal_param`), so the OCCT `Focal` of the position reads
/// `focal_param / 2` and `U * U / (4 * Focal)` becomes
/// `U * U / (2 * focal_param)` — the identical evaluation the geom
/// `point_at` runs.
pub(crate) fn elclib2d_parabola_value(u: f64, parab: &Parabola2d) -> DVec2 {
    let a_focal = parab.focal_param / 2.0;
    if a_focal.abs() <= crate::base::extrema_ext_elc::GP_RESOLUTION {
        // ElCLib.cxx L577-583: the degenerate-focal arm answers along XDir.
        return parab.origin + u * parab.axis_dir;
    }
    let a1 = u * u / (4.0 * a_focal);
    parab.origin + a1 * parab.axis_dir + u * parabola2d_y_dir(parab)
}

/// OCCT gp_Dir2d::Angle (gp_Dir2d.cxx L26-56) — the angle in [-PI, PI], the
/// arccos above 45 degrees and the arcsin below, with the sign of the sine.
fn dir2d_angle(a: DVec2, other: DVec2) -> f64 {
    let a_cosinus = a.dot(other);
    let a_sinus = a.x * other.y - a.y * other.x;
    if a_cosinus > -0.70710678118655 && a_cosinus < 0.70710678118655 {
        if a_sinus > 0.0 {
            a_cosinus.acos()
        } else {
            -a_cosinus.acos()
        }
    } else if a_cosinus > 0.0 {
        a_sinus.abs().asin()
    } else {
        let a_val = std::f64::consts::PI - a_sinus.abs().asin();
        if a_sinus > 0.0 {
            a_val
        } else {
            -a_val
        }
    }
}

/// OCCT gp_Dir2d::IsParallel (gp_Dir2d.hxx L422-436): the absolute angle is
/// within `AngularTolerance` of 0 or of PI.
pub(crate) fn dir2d_is_parallel(a: DVec2, other: DVec2, angular_tolerance: f64) -> bool {
    let mut an_ang = dir2d_angle(a, other);
    if an_ang < 0.0 {
        an_ang = -an_ang;
    }
    an_ang <= angular_tolerance || std::f64::consts::PI - an_ang <= angular_tolerance
}

/// OCCT gp_Lin2d::SquareDistance(gp_Pnt2d) (gp_Lin2d.hxx L240-247): the
/// squared cross of (P - Location) with the direction.
pub(crate) fn line2d_square_distance(lin: &Line2d, p: DVec2) -> f64 {
    let a_coord = p - lin.origin;
    let a_d = a_coord.x * lin.direction.y - a_coord.y * lin.direction.x;
    a_d * a_d
}

// =============================================================================
// Extrema_ExtElC2d (hxx L35-99, cxx L32-543)
// =============================================================================

/// OCCT Extrema_ExtElC2d (hxx L35-99) — all the distances between two
/// elementary 2D curves.
#[derive(Debug, Clone)]
pub struct ExtremaExtElC2d {
    /// hxx L94: bool myDone.
    my_done: bool,
    /// hxx L95: bool myIsPar.
    my_is_par: bool,
    /// hxx L96: int myNbExt.
    my_nb_ext: usize,
    /// hxx L97: double mySqDist[8].
    my_sq_dist: [f64; 8],
    /// hxx L98: Extrema_POnCurv2d myPoint[8][2].
    my_point: [[POnCurve2d; 2]; 8],
}

/// The OCCT default-constructed Extrema_POnCurv2d (parameter 0, point (0,0)).
pub(crate) fn default_p_on_curve2d() -> POnCurve2d {
    POnCurve2d {
        param: 0.0,
        point: DVec2::ZERO,
    }
}

impl ExtremaExtElC2d {
    /// OCCT Extrema_ExtElC2d(const gp_Lin2d& C1, const gp_Lin2d& C2, AngTol)
    /// (cxx L45-101).
    pub fn line_line(the_c1: &Line2d, the_c2: &Line2d, _ang_tol: f64) -> Self {
        let mut this = ExtremaExtElC2d {
            my_done: false,
            my_is_par: false,
            my_nb_ext: 0,
            my_sq_dist: [REAL_LAST; 8],
            my_point: std::array::from_fn(|_| [default_p_on_curve2d(), default_p_on_curve2d()]),
        };

        // cxx L68-69: gp_Vec2d D1(C1.Direction()), D2(C2.Direction()).
        let a_d1 = the_c1.direction;
        let a_d2 = the_c2.direction;
        // cxx L70: D1.IsParallel(D2, Precision::Angular()).
        if dir2d_is_parallel(a_d1, a_d2, ANGULAR) {
            // cxx L72-75: the parallel arm — distance from C1.Location() to C2.
            this.my_is_par = true;
            this.my_sq_dist[0] = line2d_square_distance(the_c2, the_c1.origin);
            this.my_nb_ext = 1;
        } else {
            // cxx L79: gp_Vec2d aP1P2(C1.Location(), C2.Location()) = P2 - P1.
            let a_p1p2 = the_c2.origin - the_c1.origin;
            // cxx L86: aDelim = 1 / (D1 ^ D2).
            let a_delim = 1.0 / (a_d1.x * a_d2.y - a_d1.y * a_d2.x);
            // cxx L88-89.
            let a_param1 = (a_p1p2.x * a_d2.y - a_p1p2.y * a_d2.x) * a_delim;
            let a_param2 = -(a_d1.x * a_p1p2.y - a_d1.y * a_p1p2.x) * a_delim;

            // cxx L91-92.
            let a_p1 = elclib2d_line_value(a_param1, the_c1);
            let a_p2 = elclib2d_line_value(a_param2, the_c2);

            // cxx L94-97.
            this.my_sq_dist[this.my_nb_ext] = 0.0;
            this.my_point[this.my_nb_ext][0] = POnCurve2d {
                param: a_param1,
                point: a_p1,
            };
            this.my_point[this.my_nb_ext][1] = POnCurve2d {
                param: a_param2,
                point: a_p2,
            };
            this.my_nb_ext = 1;
        }

        // cxx L100.
        this.my_done = true;
        this
    }

    /// OCCT Extrema_ExtElC2d(const gp_Lin2d& C1, const gp_Circ2d& C2, Tol)
    /// (cxx L104-167).
    pub fn line_circle(the_c1: &Line2d, the_c2: &Circle2d, _tol: f64) -> Self {
        let mut this = ExtremaExtElC2d {
            my_done: false,
            my_is_par: false,
            my_nb_ext: 0,
            my_sq_dist: [REAL_LAST; 8],
            my_point: std::array::from_fn(|_| [default_p_on_curve2d(), default_p_on_curve2d()]),
        };

        // cxx L126-132: T1 in the reference of the circle.
        let a_d = the_c1.direction;
        let a_x2 = the_c2.x_dir;
        let a_y2 = the_c2.y_dir;
        let a_dx = a_d.dot(a_x2);
        let a_dy = a_d.dot(a_y2);
        let a_o1 = the_c1.origin;

        // cxx L137-149.
        let mut a_teta = [0.0f64; 2];
        if a_dy.abs() <= REAL_EPSILON {
            a_teta[0] = std::f64::consts::FRAC_PI_2;
        } else {
            a_teta[0] = (-a_dx / a_dy).atan();
        }
        a_teta[1] = a_teta[0] + std::f64::consts::PI;
        if a_teta[0] < 0.0 {
            a_teta[0] += std::f64::consts::PI + std::f64::consts::PI;
        }

        // cxx L151-157.
        let a_p2 = elclib2d_circle_value(a_teta[0], the_c2);
        let a_u1 = (a_p2 - a_o1).dot(a_d);
        let a_p1 = elclib2d_line_value(a_u1, the_c1);
        this.my_sq_dist[this.my_nb_ext] = a_p1.distance_squared(a_p2);
        this.my_point[this.my_nb_ext][0] = POnCurve2d {
            param: a_u1,
            point: a_p1,
        };
        this.my_point[this.my_nb_ext][1] = POnCurve2d {
            param: a_teta[0],
            point: a_p2,
        };
        this.my_nb_ext += 1;

        // cxx L159-165.
        let a_p2 = elclib2d_circle_value(a_teta[1], the_c2);
        let a_u1 = (a_p2 - a_o1).dot(a_d);
        let a_p1 = elclib2d_line_value(a_u1, the_c1);
        this.my_sq_dist[this.my_nb_ext] = a_p1.distance_squared(a_p2);
        this.my_point[this.my_nb_ext][0] = POnCurve2d {
            param: a_u1,
            point: a_p1,
        };
        this.my_point[this.my_nb_ext][1] = POnCurve2d {
            param: a_teta[1],
            point: a_p2,
        };
        this.my_nb_ext += 1;

        // cxx L166.
        this.my_done = true;
        this
    }

    /// OCCT Extrema_ExtElC2d(const gp_Lin2d& C1, const gp_Elips2d& C2)
    /// (cxx L170-222).
    pub fn line_ellipse(the_c1: &Line2d, the_c2: &Ellipse2d) -> Self {
        let mut this = ExtremaExtElC2d {
            // cxx L172-175: myDone = true; myIsPar = false; myDone = false;
            // myNbExt = 0.
            my_done: true,
            my_is_par: false,
            my_nb_ext: 0,
            my_sq_dist: [REAL_LAST; 8],
            my_point: std::array::from_fn(|_| [default_p_on_curve2d(), default_p_on_curve2d()]),
        };
        this.my_done = false;

        // cxx L182-190: T1 in the reference of the ellipse.
        let a_d = the_c1.direction;
        let a_x2 = the_c2.major_dir;
        let a_y2 = the_c2.minor_dir;
        let a_dx = a_d.dot(a_x2);
        let a_dy = a_d.dot(a_y2);
        let a_r1 = the_c2.major_radius;
        let a_r2 = the_c2.minor_radius;
        let a_o1 = the_c1.origin;

        // cxx L192-205.
        let mut a_teta = [0.0f64; 2];
        if a_dy.abs() <= REAL_EPSILON {
            a_teta[0] = std::f64::consts::FRAC_PI_2;
        } else {
            a_teta[0] = (-a_dx * a_r2 / (a_dy * a_r1)).atan();
        }
        a_teta[1] = a_teta[0] + std::f64::consts::PI;
        if a_teta[0] < 0.0 {
            a_teta[0] += std::f64::consts::PI + std::f64::consts::PI;
        }

        // cxx L206-212.
        let a_p2 = elclib2d_ellipse_value(a_teta[0], the_c2);
        let a_u1 = (a_p2 - a_o1).dot(a_d);
        let a_p1 = elclib2d_line_value(a_u1, the_c1);
        this.my_sq_dist[this.my_nb_ext] = a_p1.distance_squared(a_p2);
        this.my_point[this.my_nb_ext][0] = POnCurve2d {
            param: a_u1,
            point: a_p1,
        };
        this.my_point[this.my_nb_ext][1] = POnCurve2d {
            param: a_teta[0],
            point: a_p2,
        };
        this.my_nb_ext += 1;

        // cxx L214-220.
        let a_p2 = elclib2d_ellipse_value(a_teta[1], the_c2);
        let a_u1 = (a_p2 - a_o1).dot(a_d);
        let a_p1 = elclib2d_line_value(a_u1, the_c1);
        this.my_sq_dist[this.my_nb_ext] = a_p1.distance_squared(a_p2);
        this.my_point[this.my_nb_ext][0] = POnCurve2d {
            param: a_u1,
            point: a_p1,
        };
        this.my_point[this.my_nb_ext][1] = POnCurve2d {
            param: a_teta[1],
            point: a_p2,
        };
        this.my_nb_ext += 1;

        // cxx L221.
        this.my_done = true;
        this
    }

    /// OCCT Extrema_ExtElC2d(const gp_Lin2d& C1, const gp_Hypr2d& C2)
    /// (cxx L226-269).
    pub fn line_hyperbola(the_c1: &Line2d, the_c2: &Hyperbola2d) -> Self {
        let mut this = ExtremaExtElC2d {
            my_done: false,
            my_is_par: false,
            my_nb_ext: 0,
            my_sq_dist: [REAL_LAST; 8],
            my_point: std::array::from_fn(|_| [default_p_on_curve2d(), default_p_on_curve2d()]),
        };

        // cxx L237-242: T1 in the reference of the hyperbola.
        let a_d = the_c1.direction;
        let a_x2 = the_c2.major_dir;
        let a_y2 = hyperbola2d_y_dir(the_c2);
        let a_dx = a_d.dot(a_x2);
        let a_dy = a_d.dot(a_y2);

        // cxx L244-249: the degenerate arms return with myDone == false.
        let a_r = the_c2.semi_major;
        let a_r_minor = the_c2.semi_minor;
        if a_dy.abs() < REAL_EPSILON {
            return this;
        }
        if (a_r - a_r_minor * a_dx / a_dy).abs() < REAL_EPSILON {
            return this;
        }

        // cxx L255-259: v2 = ...; U2 starts at 0.0 and only the v2 > 0 arm
        // takes the logarithm.
        let a_v2 = (a_r + a_r_minor * a_dx / a_dy) / (a_r - a_r_minor * a_dx / a_dy);
        let mut a_u2 = 0.0f64;
        if a_v2 > 0.0 {
            a_u2 = a_v2.sqrt().ln();
        }
        let a_p2 = elclib2d_hyperbola_value(a_u2, the_c2);

        // cxx L262-268.
        let a_u1 = (a_p2 - the_c1.origin).dot(a_d);
        let a_p1 = elclib2d_line_value(a_u1, the_c1);
        this.my_sq_dist[this.my_nb_ext] = a_p1.distance_squared(a_p2);
        this.my_point[this.my_nb_ext][0] = POnCurve2d {
            param: a_u1,
            point: a_p1,
        };
        this.my_point[this.my_nb_ext][1] = POnCurve2d {
            param: a_u2,
            point: a_p2,
        };
        this.my_nb_ext += 1;
        this.my_done = true;
        this
    }

    /// OCCT Extrema_ExtElC2d(const gp_Lin2d& C1, const gp_Parab2d& C2)
    /// (cxx L273-307).
    pub fn line_parabola(the_c1: &Line2d, the_c2: &Parabola2d) -> Self {
        let mut this = ExtremaExtElC2d {
            my_done: false,
            my_is_par: false,
            my_nb_ext: 0,
            my_sq_dist: [REAL_LAST; 8],
            my_point: std::array::from_fn(|_| [default_p_on_curve2d(), default_p_on_curve2d()]),
        };

        // cxx L284-289: T1 in the reference of the parabola; MirrorAxis
        // direction is the X direction of the OCCT position, Axis().YAxis()
        // direction the Y direction.
        let a_d = the_c1.direction;
        let a_x2 = the_c2.axis_dir;
        let a_y2 = parabola2d_y_dir(the_c2);
        let a_dx = a_d.dot(a_x2);
        let a_dy = a_d.dot(a_y2);

        // cxx L291: double U1, U2, P = C2.Parameter() — the OCCT parameter
        // reads 2 * Focal, which is the rcad focal_param.
        let a_p = the_c2.focal_param;

        // cxx L293-296: the degenerate arm returns with myDone == false.
        if a_dy.abs() < REAL_EPSILON {
            return this;
        }

        // cxx L297-298.
        let a_u2 = a_dx * a_p / a_dy;
        let a_p2 = elclib2d_parabola_value(a_u2, the_c2);

        // cxx L300-306.
        let a_u1 = (a_p2 - the_c1.origin).dot(a_d);
        let a_p1 = elclib2d_line_value(a_u1, the_c1);
        this.my_sq_dist[this.my_nb_ext] = a_p1.distance_squared(a_p2);
        this.my_point[this.my_nb_ext][0] = POnCurve2d {
            param: a_u1,
            point: a_p1,
        };
        this.my_point[this.my_nb_ext][1] = POnCurve2d {
            param: a_u2,
            point: a_p2,
        };
        this.my_nb_ext += 1;
        this.my_done = true;
        this
    }

    /// OCCT Extrema_ExtElC2d(const gp_Circ2d& C1, const gp_Circ2d& C2)
    /// (cxx L311-366).
    pub fn circle_circle(the_c1: &Circle2d, the_c2: &Circle2d) -> Self {
        let mut this = ExtremaExtElC2d {
            // cxx L313-320: myIsPar = false; myDone = false; myNbExt = 0;
            // myDone = true; the RealLast fill.
            my_done: false,
            my_is_par: false,
            my_nb_ext: 0,
            my_sq_dist: [REAL_LAST; 8],
            my_point: std::array::from_fn(|_| [default_p_on_curve2d(), default_p_on_curve2d()]),
        };
        this.my_done = true;

        // cxx L322-323.
        let a_o1 = the_c1.center;
        let a_o2 = the_c2.center;

        // cxx L325-335: the concentric arm.
        let a_do1o2 = a_o2 - a_o1;
        let a_sq_d_centers = a_do1o2.length_squared();
        if a_sq_d_centers < SQUARE_CONFUSION {
            this.my_is_par = true;
            this.my_nb_ext = 1;
            this.my_done = true;
            let a_dr = the_c1.radius - the_c2.radius;
            this.my_sq_dist[0] = a_dr * a_dr;
            return this;
        }

        // cxx L337-352.
        let a_r1 = the_c1.radius;
        let a_r2 = the_c2.radius;
        let a_o1o2 = a_do1o2 / a_sq_d_centers.sqrt();

        let a_p1 = [a_o1 + a_r1 * a_o1o2, a_o1 - a_r1 * a_o1o2];
        let a_usol1 = [
            elclib_circle2d_parameter(the_c1, a_p1[0]),
            elclib_circle2d_parameter(the_c1, a_p1[1]),
        ];
        let a_p2 = [a_o2 + a_r2 * a_o1o2, a_o2 - a_r2 * a_o1o2];
        let a_usol2 = [
            elclib_circle2d_parameter(the_c2, a_p2[0]),
            elclib_circle2d_parameter(the_c2, a_p2[1]),
        ];

        // cxx L354-365.
        for a_no_sol in 0..=1usize {
            let a_u1 = a_usol1[a_no_sol];
            for a_kk in 0..=1usize {
                let a_u2 = a_usol2[a_kk];
                this.my_sq_dist[this.my_nb_ext] = a_p2[a_kk].distance_squared(a_p1[a_no_sol]);
                this.my_point[this.my_nb_ext][0] = POnCurve2d {
                    param: a_u1,
                    point: a_p1[a_no_sol],
                };
                this.my_point[this.my_nb_ext][1] = POnCurve2d {
                    param: a_u2,
                    point: a_p2[a_kk],
                };
                this.my_nb_ext += 1;
            }
        }
        this
    }

    /// OCCT Extrema_ExtElC2d(const gp_Circ2d& C1, const gp_Elips2d& C2)
    /// (cxx L370-406).
    pub fn circle_ellipse(the_c1: &Circle2d, the_c2: &Ellipse2d) -> Self {
        let mut this = ExtremaExtElC2d {
            my_done: false,
            my_is_par: false,
            my_nb_ext: 0,
            my_sq_dist: [REAL_LAST; 8],
            my_point: std::array::from_fn(|_| [default_p_on_curve2d(), default_p_on_curve2d()]),
        };

        // cxx L382: Extrema_ExtPElC2d ExtElip(C1.Location(), C2,
        // Precision::Confusion(), 0.0, 2.0 * M_PI).
        let an_ext_elip = ExtremaExtPElC2d::point_ellipse(
            the_c1.center,
            the_c2,
            CONFUSION,
            0.0,
            std::f64::consts::PI + std::f64::consts::PI,
        );

        // cxx L384-405.
        if an_ext_elip.is_done() {
            for an_i in 1..=an_ext_elip.nb_ext() {
                // cxx L388-392: Extrema_ExtPElC2d ExtCirc(ExtElip.Point(i).Value(),
                // C1, Precision::Confusion(), 0.0, 2.0 * M_PI).
                let an_ext_circ = ExtremaExtPElC2d::point_circle(
                    an_ext_elip.point(an_i).point,
                    the_c1,
                    CONFUSION,
                    0.0,
                    std::f64::consts::PI + std::f64::consts::PI,
                );
                if an_ext_circ.is_done() {
                    for a_j in 1..=an_ext_circ.nb_ext() {
                        this.my_sq_dist[this.my_nb_ext] = an_ext_circ.square_distance(a_j);
                        this.my_point[this.my_nb_ext][0] = an_ext_circ.point(a_j);
                        this.my_point[this.my_nb_ext][1] = an_ext_elip.point(an_i);
                        this.my_nb_ext += 1;
                    }
                }
                this.my_done = true;
            }
        }
        this
    }

    /// OCCT Extrema_ExtElC2d(const gp_Circ2d& C1, const gp_Hypr2d& C2)
    /// (cxx L410-446).
    pub fn circle_hyperbola(the_c1: &Circle2d, the_c2: &Hyperbola2d) -> Self {
        let mut this = ExtremaExtElC2d {
            my_done: false,
            my_is_par: false,
            my_nb_ext: 0,
            my_sq_dist: [REAL_LAST; 8],
            my_point: std::array::from_fn(|_| [default_p_on_curve2d(), default_p_on_curve2d()]),
        };

        // cxx L422: Extrema_ExtPElC2d ExtHyp(C1.Location(), C2,
        // Precision::Confusion(), RealFirst(), RealLast()).
        let an_ext_hyp = ExtremaExtPElC2d::point_hyperbola(
            the_c1.center,
            the_c2,
            CONFUSION,
            REAL_FIRST,
            REAL_LAST,
        );

        // cxx L424-445.
        if an_ext_hyp.is_done() {
            for an_i in 1..=an_ext_hyp.nb_ext() {
                // cxx L428-432.
                let an_ext_circ = ExtremaExtPElC2d::point_circle(
                    an_ext_hyp.point(an_i).point,
                    the_c1,
                    CONFUSION,
                    0.0,
                    std::f64::consts::PI + std::f64::consts::PI,
                );
                if an_ext_circ.is_done() {
                    for a_j in 1..=an_ext_circ.nb_ext() {
                        this.my_sq_dist[this.my_nb_ext] = an_ext_circ.square_distance(a_j);
                        this.my_point[this.my_nb_ext][0] = an_ext_circ.point(a_j);
                        this.my_point[this.my_nb_ext][1] = an_ext_hyp.point(an_i);
                        this.my_nb_ext += 1;
                    }
                }
                this.my_done = true;
            }
        }
        this
    }

    /// OCCT Extrema_ExtElC2d(const gp_Circ2d& C1, const gp_Parab2d& C2)
    /// (cxx L450-486).
    pub fn circle_parabola(the_c1: &Circle2d, the_c2: &Parabola2d) -> Self {
        let mut this = ExtremaExtElC2d {
            my_done: false,
            my_is_par: false,
            my_nb_ext: 0,
            my_sq_dist: [REAL_LAST; 8],
            my_point: std::array::from_fn(|_| [default_p_on_curve2d(), default_p_on_curve2d()]),
        };

        // cxx L462: Extrema_ExtPElC2d ExtParab(C1.Location(), C2,
        // Precision::Confusion(), RealFirst(), RealLast()).
        let an_ext_parab = ExtremaExtPElC2d::point_parabola(
            the_c1.center,
            the_c2,
            CONFUSION,
            REAL_FIRST,
            REAL_LAST,
        );

        // cxx L464-485.
        if an_ext_parab.is_done() {
            for an_i in 1..=an_ext_parab.nb_ext() {
                // cxx L468-472.
                let an_ext_circ = ExtremaExtPElC2d::point_circle(
                    an_ext_parab.point(an_i).point,
                    the_c1,
                    CONFUSION,
                    0.0,
                    std::f64::consts::PI + std::f64::consts::PI,
                );
                if an_ext_circ.is_done() {
                    for a_j in 1..=an_ext_circ.nb_ext() {
                        this.my_sq_dist[this.my_nb_ext] = an_ext_circ.square_distance(a_j);
                        this.my_point[this.my_nb_ext][0] = an_ext_circ.point(a_j);
                        this.my_point[this.my_nb_ext][1] = an_ext_parab.point(an_i);
                        this.my_nb_ext += 1;
                    }
                }
                this.my_done = true;
            }
        }
        this
    }

    /// OCCT IsDone() (cxx L490-493).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT IsParallel() (cxx L497-504).
    pub fn is_parallel(&self) -> bool {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_is_par
    }

    /// OCCT NbExt() (cxx L508-516).
    pub fn nb_ext(&self) -> usize {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_nb_ext
    }

    /// OCCT SquareDistance(N) (cxx L520-528) — 1-based.
    pub fn square_distance(&self, n: usize) -> f64 {
        if n < 1 || n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.my_sq_dist[n - 1]
    }

    /// OCCT Points(N, P1, P2) (cxx L532-540) — 1-based.
    pub fn points(&self, n: usize, p1: &mut POnCurve2d, p2: &mut POnCurve2d) {
        if n < 1 || n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        *p1 = self.my_point[n - 1][0].clone();
        *p2 = self.my_point[n - 1][1].clone();
    }
}
