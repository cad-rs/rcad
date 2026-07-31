//! 2D analytic intersection algorithms (IntAna2d package).
//!
//! Provides analytical intersection between 2D geometric elements:
//! lines, circles, and general conics (ellipse, parabola, hyperbola).
//!
//! OCCT TKGeomBase IntAna2d package: IntAna2d_AnaIntersection, IntAna2d_Conic,
//! IntAna2d_IntPoint, IntAna2d_Outils.

use glam::DVec2;

use crate::geom::{Circle2d, Ellipse2d, Hyperbola2d, Line2d, Parabola2d, Point2, Vec2};

// ============================================================================
// Tolerance constants (IntAna2d — machine precision level)
// ============================================================================

/// Angular tolerance (1e-12 rad).
const TOL_ANG: f64 = 1e-12;
/// Confusion tolerance — point coincidence (1e-10).
const TOL_CONF: f64 = 1e-10;
/// Polynomial root tolerance (1e-12).
const TOL_ROOT: f64 = 1e-12;
/// Zero discriminant threshold for polynomial solvers.
const TOL_DISC: f64 = 1e-14;

// ============================================================================
// IntAna2d_Conic
// ============================================================================

/// A conic defined by its implicit quadratic equation:
///
/// ```text
/// A·x² + B·y² + 2·C·x·y + 2·D·x + 2·E·y + F = 0
/// ```
///
/// OCCT: `IntAna2d_Conic`.
#[derive(Debug, Clone, Copy)]
pub struct Conic2d {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
    pub e: f64,
    pub f: f64,
}

impl Conic2d {
    /// Construct a conic from a 2D line.
    ///
    /// A line in implicit form: the squared distance form is degenerate,
    /// but we can represent a line as `(n·P - d)² = 0` or expand as a pair
    /// of coincident lines. OCCT represents a line as its linear equation
    /// embedded in the quadratic form: A=B=C=0, giving 2*(D*x + E*y) + F = 0.
    /// To detect line-conic intersections, we store the line coefficients
    /// directly: D and E are the normal components, F is the constant term.
    pub fn from_line(line: &Line2d) -> Self {
        // The line: n·(P - origin) = 0, where n = (-dir.y, dir.x) is the normal.
        // n.x * (x - ox) + n.y * (y - oy) = 0
        // n.x * x + n.y * y - (n.x*ox + n.y*oy) = 0
        let n = Vec2::new(-line.direction.y, line.direction.x);
        let d_const = -(n.x * line.origin.x + n.y * line.origin.y);
        // Embed in quadratic: A=B=C=0, D=n.x/2, E=n.y/2, F=d_const
        Conic2d {
            a: 0.0,
            b: 0.0,
            c: 0.0,
            d: n.x * 0.5,
            e: n.y * 0.5,
            f: d_const,
        }
    }

    /// Construct a conic from a 2D circle.
    ///
    /// Circle: (x - cx)² + (y - cy)² = R²
    ///       = x² + y² - 2*cx*x - 2*cy*y + (cx² + cy² - R²) = 0
    /// So: A=1, B=1, C=0, D=-cx, E=-cy, F=cx²+cy²-R²
    pub fn from_circle(circle: &Circle2d) -> Self {
        let cx = circle.center.x;
        let cy = circle.center.y;
        let r = circle.radius;
        Conic2d {
            a: 1.0,
            b: 1.0,
            c: 0.0,
            d: -cx,
            e: -cy,
            f: cx * cx + cy * cy - r * r,
        }
    }

    /// Construct a conic from an ellipse.
    ///
    /// Parametric: P(θ) = center + major_dir * a*cos(θ) + minor_dir * b*sin(θ)
    /// In local frame (X = major_dir, Y = minor_dir):
    ///   (x - cx)²/a² + (y - cy)²/b² = 1
    /// Expanded to implicit form.
    pub fn from_ellipse(ellipse: &Ellipse2d) -> Self {
        // Local frame: X = major_dir, Y = minor_dir = (-major_dir.y, major_dir.x)
        let a = ellipse.major_radius;
        let b = ellipse.minor_radius;
        let cx = ellipse.center.x;
        let cy = ellipse.center.y;
        let cos_a = ellipse.major_dir.x;
        let sin_a = ellipse.major_dir.y;

        if a.abs() < TOL_CONF || b.abs() < TOL_CONF {
            return Conic2d {
                a: 0.0,
                b: 0.0,
                c: 0.0,
                d: 0.0,
                e: 0.0,
                f: 0.0,
            };
        }

        // Transform implicit form back to global coordinates
        // The ellipse in local frame (u,v): u²/a² + v²/b² = 1
        // where u = (x-cx)*cos_a + (y-cy)*sin_a
        //       v = -(x-cx)*sin_a + (y-cy)*cos_a
        let a2 = a * a;
        let b2 = b * b;
        let cos2 = cos_a * cos_a;
        let sin2 = sin_a * sin_a;
        let cs = cos_a * sin_a;

        let aa = cos2 / a2 + sin2 / b2;
        let bb = sin2 / a2 + cos2 / b2;
        let cc = cs / a2 - cs / b2;
        let dd = -(aa * cx + cc * cy);
        let ee = -(cc * cx + bb * cy);
        let ff = aa * cx * cx + 2.0 * cc * cx * cy + bb * cy * cy - 1.0;

        Conic2d {
            a: aa,
            b: bb,
            c: cc,
            d: dd,
            e: ee,
            f: ff,
        }
    }

    /// Construct a conic from a parabola.
    ///
    /// Parabola in local frame: y² = 2*p*x (where p = focal_param)
    /// In standard position: y² - 2*p*x = 0
    /// General position: rotated and translated.
    pub fn from_parabola(parabola: &Parabola2d) -> Self {
        let p = parabola.focal_param;
        let ox = parabola.origin.x;
        let oy = parabola.origin.y;
        let cos_a = parabola.axis_dir.x;
        let sin_a = parabola.axis_dir.y;

        if p.abs() < TOL_CONF {
            return Conic2d {
                a: 0.0,
                b: 0.0,
                c: 0.0,
                d: 0.0,
                e: 0.0,
                f: 0.0,
            };
        }

        // In local frame (u,v): v² = 2*p*u
        // u = (x-ox)*cos_a + (y-oy)*sin_a
        // v = -(x-ox)*sin_a + (y-oy)*cos_a
        //
        // v² - 2*p*u = 0
        let cos2 = cos_a * cos_a;
        let sin2 = sin_a * sin_a;
        let cs = cos_a * sin_a;

        // v² = ( -(x-ox)*sin + (y-oy)*cos )²
        //    = (x-ox)²*sin² - 2*(x-ox)*(y-oy)*sin*cos + (y-oy)²*cos²
        // 2*p*u = 2*p*((x-ox)*cos + (y-oy)*sin)
        Conic2d {
            a: sin2,
            b: cos2,
            c: -cs,
            d: -(sin2 * ox - cs * oy + p * cos_a),
            e: -(cos2 * oy - cs * ox + p * sin_a),
            f: sin2 * ox * ox - 2.0 * cs * ox * oy + cos2 * oy * oy - 2.0 * p * (ox * cos_a + oy * sin_a),
        }
    }

    /// Construct a conic from a hyperbola.
    ///
    /// Hyperbola in local frame: u²/a² - v²/b² = 1
    /// where u = (x-cx)*cos + (y-cy)*sin, v = -(x-cx)*sin + (y-cy)*cos
    pub fn from_hyperbola(hyperbola: &Hyperbola2d) -> Self {
        let a = hyperbola.semi_major;
        let b = hyperbola.semi_minor;
        let cx = hyperbola.center.x;
        let cy = hyperbola.center.y;
        let cos_a = hyperbola.major_dir.x;
        let sin_a = hyperbola.major_dir.y;

        if a.abs() < TOL_CONF || b.abs() < TOL_CONF {
            return Conic2d {
                a: 0.0,
                b: 0.0,
                c: 0.0,
                d: 0.0,
                e: 0.0,
                f: 0.0,
            };
        }

        let a2 = a * a;
        let b2 = b * b;
        let cos2 = cos_a * cos_a;
        let sin2 = sin_a * sin_a;
        let cs = cos_a * sin_a;

        // u²/a² - v²/b² = 1
        let aa = cos2 / a2 - sin2 / b2;
        let bb = sin2 / a2 - cos2 / b2;
        let cc = cs / a2 + cs / b2;
        let dd = -(aa * cx + cc * cy);
        let ee = -(cc * cx + bb * cy);
        let ff = aa * cx * cx + 2.0 * cc * cx * cy + bb * cy * cy - 1.0;

        Conic2d {
            a: aa,
            b: bb,
            c: cc,
            d: dd,
            e: ee,
            f: ff,
        }
    }

    /// Evaluate the conic function at (x, y).
    ///
    /// OCCT: `Value(X, Y)` — returns F(x, y).
    pub fn value(&self, x: f64, y: f64) -> f64 {
        self.a * x * x + self.b * y * y + 2.0 * self.c * x * y + 2.0 * self.d * x + 2.0 * self.e * y + self.f
    }

    /// Evaluate the gradient at (x, y).
    ///
    /// OCCT: `Grad(X, Y)` — returns (∂F/∂x, ∂F/∂y).
    pub fn grad(&self, x: f64, y: f64) -> DVec2 {
        DVec2::new(
            2.0 * (self.a * x + self.c * y + self.d),
            2.0 * (self.b * y + self.c * x + self.e),
        )
    }

    /// Evaluate the function and gradient at (x, y).
    ///
    /// OCCT: `ValAndGrad(X, Y, Val, Grd)`.
    pub fn val_and_grad(&self, x: f64, y: f64) -> (f64, DVec2) {
        let val = self.value(x, y);
        let grd = self.grad(x, y);
        (val, grd)
    }

    /// Returns the six coefficients (A, B, C, D, E, F).
    ///
    /// OCCT: `Coefficients(A, B, C, D, E, F)`.
    pub fn coefficients(&self) -> (f64, f64, f64, f64, f64, f64) {
        (self.a, self.b, self.c, self.d, self.e, self.f)
    }

    /// Compute coefficients in a new coordinate system.
    ///
    /// OCCT: `NewCoefficients(A, B, C, D, E, F, Axe)`.
    pub fn new_coefficients(&self, origin: DVec2, u_dir: DVec2) -> (f64, f64, f64, f64, f64, f64) {
        let cos_a = u_dir.x;
        let sin_a = u_dir.y;
        let cos2 = cos_a * cos_a;
        let sin2 = sin_a * sin_a;
        let cs = cos_a * sin_a;

        // Transform coordinates: x = ox + u*cos_a + v*(-sin_a)
        //                       y = oy + u*sin_a + v*cos_a
        // Substitute into A*x² + B*y² + 2*C*x*y + 2*D*x + 2*E*y + F
        let ox = origin.x;
        let oy = origin.y;

        // Coefficients after translation to origin, then rotation
        let aa = self.a * cos2 + self.b * sin2 + 2.0 * self.c * cs;
        let bb = self.a * sin2 + self.b * cos2 - 2.0 * self.c * cs;
        let cc = (self.b - self.a) * cs + self.c * (cos2 - sin2);
        let dd = self.a * ox * cos_a + self.b * oy * sin_a + self.c * (ox * sin_a + oy * cos_a) + self.d * cos_a + self.e * sin_a;
        let ee = -self.a * ox * sin_a + self.b * oy * cos_a + self.c * (-ox * cos_a + oy * sin_a) - self.d * sin_a + self.e * cos_a;
        let ff = self.value(ox, oy);

        (aa, bb, cc, dd, ee, ff)
    }
}

// ============================================================================
// IntAna2d_IntPoint
// ============================================================================

/// A 2D intersection point between two curves, storing parametric values.
///
/// OCCT: `IntAna2d_IntPoint`.
#[derive(Debug, Clone, Copy)]
pub struct IntPoint2d {
    /// The 2D coordinates of the intersection point.
    point: Point2,
    /// Parameter on the first curve.
    param1: f64,
    /// Parameter on the second curve (undefined if second is implicit).
    param2: f64,
    /// True when the second curve is defined by an implicit equation.
    second_is_implicit: bool,
}

impl IntPoint2d {
    /// Create an intersection point between two parametric curves.
    ///
    /// OCCT: `IntAna2d_IntPoint(X, Y, U1, U2)`.
    pub fn new(x: f64, y: f64, u1: f64, u2: f64) -> Self {
        IntPoint2d {
            point: DVec2::new(x, y),
            param1: u1,
            param2: u2,
            second_is_implicit: false,
        }
    }

    /// Create an intersection point between a parametric curve and an implicit curve.
    ///
    /// OCCT: `IntAna2d_IntPoint(X, Y, U1)`.
    pub fn new_implicit(x: f64, y: f64, u1: f64) -> Self {
        IntPoint2d {
            point: DVec2::new(x, y),
            param1: u1,
            param2: 0.0,
            second_is_implicit: true,
        }
    }

    /// The geometric point.
    ///
    /// OCCT: `Value()`.
    pub fn value(&self) -> Point2 {
        self.point
    }

    /// Parameter on the first curve.
    ///
    /// OCCT: `ParamOnFirst()`.
    pub fn param_on_first(&self) -> f64 {
        self.param1
    }

    /// Parameter on the second curve.
    ///
    /// OCCT: `ParamOnSecond()`.
    /// Panics if the second curve is implicit.
    pub fn param_on_second(&self) -> f64 {
        if self.second_is_implicit {
            panic!("IntPoint2d: second curve is implicit, no param_on_second");
        }
        self.param2
    }

    /// Returns `true` when the second curve is implicit.
    ///
    /// OCCT: `SecondIsImplicit()`.
    pub fn second_is_implicit(&self) -> bool {
        self.second_is_implicit
    }

    /// Set values for a non-implicit point.
    ///
    /// OCCT: `SetValue(X, Y, U1, U2)`.
    pub fn set_value(&mut self, x: f64, y: f64, u1: f64, u2: f64) {
        self.point = DVec2::new(x, y);
        self.param1 = u1;
        self.param2 = u2;
        self.second_is_implicit = false;
    }

    /// Set values for an implicit point.
    ///
    /// OCCT: `SetValue(X, Y, U1)`.
    pub fn set_value_implicit(&mut self, x: f64, y: f64, u1: f64) {
        self.point = DVec2::new(x, y);
        self.param1 = u1;
        self.param2 = 0.0;
        self.second_is_implicit = true;
    }
}

// ============================================================================
// IntAna2d_AnaIntersection
// ============================================================================

/// Analytical intersection between 2D curves.
///
/// Supports line-line, line-circle, circle-circle, and conic-conic intersections.
///
/// OCCT: `IntAna2d_AnaIntersection`.
#[derive(Debug, Clone)]
pub struct AnaIntersection2d {
    done: bool,
    para: bool,
    iden: bool,
    empt: bool,
    nbp: usize,
    lpnt: [IntPoint2d; 4],
}

impl AnaIntersection2d {
    /// Empty constructor. `IsDone()` returns false.
    ///
    /// OCCT: default constructor.
    pub fn new() -> Self {
        AnaIntersection2d {
            done: false,
            para: false,
            iden: false,
            empt: false,
            nbp: 0,
            lpnt: [IntPoint2d::new(0.0, 0.0, 0.0, 0.0); 4],
        }
    }

    /// Intersect two lines.
    ///
    /// OCCT: `Perform(gp_Lin2d, gp_Lin2d)`.
    pub fn perform_lin_lin(&mut self, l1: &Line2d, l2: &Line2d) {
        self.done = true;

        // Check if lines are parallel
        let cross = l1.direction.x * l2.direction.y - l1.direction.y * l2.direction.x;
        if cross.abs() < TOL_ANG {
            // Parallel or coincident
            self.para = true;
            // Check if coincident: point of l1 on l2
            let d = l2.distance(l1.origin);
            if d < TOL_CONF {
                self.iden = true;
            } else {
                self.empt = true;
            }
            self.nbp = 0;
            return;
        }

        // Intersection of two lines: solve for parameters
        // l1: o1 + t1*d1, l2: o2 + t2*d2
        // o1 + t1*d1 = o2 + t2*d2
        // In 2D: t1 = cross(o2-o1, d2) / cross(d1, d2)
        //         t2 = cross(o2-o1, d1) / cross(d1, d2)
        let diff = l2.origin - l1.origin;
        let t1 = (diff.x * l2.direction.y - diff.y * l2.direction.x) / cross;
        let t2 = (diff.x * l1.direction.y - diff.y * l1.direction.x) / cross;
        let p = l1.origin + t1 * l1.direction;

        self.nbp = 1;
        self.lpnt[0] = IntPoint2d::new(p.x, p.y, t1, t2);
    }

    /// Intersect a line and a circle.
    ///
    /// OCCT: `Perform(gp_Lin2d, gp_Circ2d)`.
    pub fn perform_lin_circ(&mut self, line: &Line2d, circle: &Circle2d) {
        self.done = true;

        // Line: P = origin + t * direction
        // Circle: |P - center|² = R²
        // |origin + t*direction - center|² = R²
        // |d|²*t² + 2*d·(o-c)*t + |o-c|² - R² = 0
        let d = line.direction;
        let diff = line.origin - circle.center;
        let a = d.dot(d);
        let b = 2.0 * d.dot(diff);
        let c = diff.dot(diff) - circle.radius * circle.radius;

        solve_quadratic_intersection(a, b, c, &mut self.lpnt, &mut self.nbp, &mut self.para, &mut self.iden, &mut self.empt);

        // Compute params on circle for each intersection point
        if self.nbp > 0 {
            for i in 0..self.nbp {
                let p = self.lpnt[i].value();
                let d_c = p - circle.center;
                let theta = d_c.y.atan2(d_c.x)
                    - circle.x_dir.y.atan2(circle.x_dir.x);
                let two_pi = 2.0 * std::f64::consts::PI;
                let theta = theta.rem_euclid(two_pi);
                self.lpnt[i] = IntPoint2d::new(p.x, p.y, self.lpnt[i].param_on_first(), theta);
            }
        }
    }

    /// Intersect two circles.
    ///
    /// OCCT: `Perform(gp_Circ2d, gp_Circ2d)`.
    pub fn perform_circ_circ(&mut self, c1: &Circle2d, c2: &Circle2d) {
        self.done = true;

        let d = c2.center - c1.center;
        let dist = d.length();

        // Check concentric
        if dist < TOL_CONF {
            self.para = true;
            if (c1.radius - c2.radius).abs() < TOL_CONF {
                self.iden = true; // Same circle
            } else {
                self.empt = true; // Concentric, different radii
            }
            self.nbp = 0;
            return;
        }

        // Check if circles intersect
        let r1 = c1.radius;
        let r2 = c2.radius;
        if dist > r1 + r2 + TOL_CONF || dist < (r1 - r2).abs() - TOL_CONF {
            self.empt = true;
            self.nbp = 0;
            return;
        }

        // Check tangent
        if (dist - (r1 + r2)).abs() < TOL_CONF || (dist - (r1 - r2).abs()).abs() < TOL_CONF {
            // Tangent: single point along the line connecting centers
            let t = if r1 > r2 { r1 / dist } else { -r1 / dist };
            let p = c1.center + d * t;
            self.nbp = 1;
            self.lpnt[0] = IntPoint2d::new(p.x, p.y, 0.0, 0.0);
            // Compute params
            let theta1 = (p - c1.center).y.atan2((p - c1.center).x)
                - c1.x_dir.y.atan2(c1.x_dir.x);
            let theta2 = (p - c2.center).y.atan2((p - c2.center).x)
                - c2.x_dir.y.atan2(c2.x_dir.x);
            let two_pi = 2.0 * std::f64::consts::PI;
            self.lpnt[0] = IntPoint2d::new(p.x, p.y, theta1.rem_euclid(two_pi), theta2.rem_euclid(two_pi));
            return;
        }

        // Two intersection points: use radical line method
        // The radical line is: 2*(c2-c1)·P + |c1|² - |c2|² + r2² - r1² = 0
        // Intersect this line with either circle.
        let a = 2.0 * d.x;
        let b = 2.0 * d.y;
        let c_const = c1.center.dot(c1.center) - c2.center.dot(c2.center) + r2 * r2 - r1 * r1;

        // Find intersection of radical line with c1
        // The radical line is: a*x + b*y + c = 0 → y = -(a*x + c)/b when |b| > |a|
        // Or x = -(b*y + c)/a when |a| > |b|
        // Substitute into circle equation
        let (t_vals, n) = if b.abs() > a.abs() {
            // y = -(a*x + c)/b
            // (x-cx)² + (-(a*x+c)/b - cy)² = r²
            let slope = -a / b;
            let intercept = -c_const / b - c1.center.y;
            let dx = c1.center.x;
            let aa = 1.0 + slope * slope;
            let bb = 2.0 * (slope * intercept - dx);
            let cc = dx * dx + intercept * intercept - r1 * r1;
            solve_quadratic_raw(aa, bb, cc)
        } else {
            // x = -(b*y + c)/a
            let slope = -b / a;
            let intercept = -c_const / a - c1.center.x;
            let dy = c1.center.y;
            let aa = 1.0 + slope * slope;
            let bb = 2.0 * (slope * intercept - dy);
            let cc = dy * dy + intercept * intercept - r1 * r1;
            solve_quadratic_raw(aa, bb, cc)
        };

        self.nbp = n.min(4);
        let two_pi = 2.0 * std::f64::consts::PI;
        for i in 0..self.nbp {
            let t = t_vals[i];
            let p = if b.abs() > a.abs() {
                let y = -(a * t + c_const) / b;
                DVec2::new(t, y)
            } else {
                let x = -(b * t + c_const) / a;
                DVec2::new(x, t)
            };
            let theta1 = (p - c1.center).y.atan2((p - c1.center).x)
                - c1.x_dir.y.atan2(c1.x_dir.x);
            let theta2 = (p - c2.center).y.atan2((p - c2.center).x)
                - c2.x_dir.y.atan2(c2.x_dir.x);
            self.lpnt[i] = IntPoint2d::new(p.x, p.y, theta1.rem_euclid(two_pi), theta2.rem_euclid(two_pi));
        }

        if self.nbp > 1 {
            // Sort by param on first circle
            let points = &mut self.lpnt[..self.nbp];
            points.sort_by(|a, b| a.param_on_first().partial_cmp(&b.param_on_first()).unwrap_or(std::cmp::Ordering::Equal));
            // Deduplicate
            let mut j = 1;
            for i in 1..self.nbp {
                let d = (points[i].value() - points[j - 1].value()).length();
                if d > TOL_CONF {
                    if i != j {
                        points[j] = points[i];
                    }
                    j += 1;
                }
            }
            self.nbp = j;
        }
    }

    /// Intersect a line with a conic.
    ///
    /// OCCT: `Perform(gp_Lin2d, IntAna2d_Conic)`.
    pub fn perform_lin_conic(&mut self, line: &Line2d, conic: &Conic2d) {
        self.done = true;

        // Line: P = origin + t * direction
        // Conic: A*x² + B*y² + 2*C*x*y + 2*D*x + 2*E*y + F = 0
        // Substitute line into conic → quadratic in t
        let ox = line.origin.x;
        let oy = line.origin.y;
        let dx = line.direction.x;
        let dy = line.direction.y;

        let a = conic.a * dx * dx
            + conic.b * dy * dy
            + 2.0 * conic.c * dx * dy;
        let b = 2.0 * conic.a * ox * dx
            + 2.0 * conic.b * oy * dy
            + 2.0 * conic.c * (ox * dy + oy * dx)
            + 2.0 * conic.d * dx
            + 2.0 * conic.e * dy;
        let c = conic.a * ox * ox
            + conic.b * oy * oy
            + 2.0 * conic.c * ox * oy
            + 2.0 * conic.d * ox
            + 2.0 * conic.e * oy
            + conic.f;

        solve_quadratic_intersection(a, b, c, &mut self.lpnt, &mut self.nbp, &mut self.para, &mut self.iden, &mut self.empt);

        // Fill in the point coordinates
        for i in 0..self.nbp {
            let t = self.lpnt[i].param_on_first();
            let p = line.origin + t * line.direction;
            self.lpnt[i] = IntPoint2d::new_implicit(p.x, p.y, t);
        }
    }

    /// Intersect a circle with a conic.
    ///
    /// OCCT: `Perform(gp_Circ2d, IntAna2d_Conic)`.
    pub fn perform_circ_conic(&mut self, circle: &Circle2d, conic: &Conic2d) {
        self.done = true;
        self.intersect_circle_conic(circle, conic);
    }

    /// Intersect an ellipse with a conic.
    ///
    /// OCCT: `Perform(gp_Elips2d, IntAna2d_Conic)`.
    pub fn perform_ellipse_conic(&mut self, ellipse: &Ellipse2d, conic: &Conic2d) {
        self.done = true;
        // Convert ellipse to parametric form and substitute into conic
        self.intersect_ellipse_conic(ellipse, conic);
    }

    /// Intersect a parabola with a conic.
    ///
    /// OCCT: `Perform(gp_Parab2d, IntAna2d_Conic)`.
    pub fn perform_parabola_conic(&mut self, parabola: &Parabola2d, conic: &Conic2d) {
        self.done = true;
        // Parabola: P(t) = origin + axis_dir * (t²/(2p)) + perp_dir * t
        // Substitute into conic → quartic in t
        self.intersect_parabola_conic(parabola, conic);
    }

    /// Intersect a hyperbola with a conic.
    ///
    /// OCCT: `Perform(gp_Hypr2d, IntAna2d_Conic)`.
    pub fn perform_hyperbola_conic(&mut self, hyperbola: &Hyperbola2d, conic: &Conic2d) {
        self.done = true;
        // Hyperbola: P(t) = center + major_dir * a*cosh(t) + minor_dir * b*sinh(t)
        // Or use rational parameterization: x = a*(1+t²)/(1-t²), y = 2*b*t/(1-t²)
        self.intersect_hyperbola_conic(hyperbola, conic);
    }

    // ========================================================================
    // Internal: circle-conic intersection via tangent half-angle substitution
    // ========================================================================

    fn intersect_circle_conic(&mut self, circle: &Circle2d, conic: &Conic2d) {
        let cx = circle.center.x;
        let cy = circle.center.y;
        let r = circle.radius;
        let cos0 = circle.x_dir.x;
        let sin0 = circle.x_dir.y;

        // Circle in global coords:
        // x = cx + r*(cos0*cosθ - sin0*sinθ) = cx + r*cos(θ+φ) where φ = atan2(sin0, cos0)
        // y = cy + r*(sin0*cosθ + cos0*sinθ) = cy + r*sin(θ+φ)
        //
        // Using tangent half-angle: u = tan(θ/2)
        // cosθ = (1-u²)/(1+u²), sinθ = 2u/(1+u²)
        //
        // Substitute into conic implicit equation → quartic in u

        let a = conic.a;
        let b = conic.b;
        let cc = conic.c;
        let d = conic.d;
        let e = conic.e;
        let f = conic.f;

        // Precompute rotated circle frame
        let r_cos0 = r * cos0;
        let r_sin0 = r * sin0;

        // x(u) = cx + r_cos0*(1-u²)/(1+u²) - r_sin0*2u/(1+u²)
        //       = (cx*(1+u²) + r_cos0*(1-u²) - 2*r_sin0*u) / (1+u²)
        // y(u) = cy + r_sin0*(1-u²)/(1+u²) + r_cos0*2u/(1+u²)
        //       = (cy*(1+u²) + r_sin0*(1-u²) + 2*r_cos0*u) / (1+u²)
        //
        // Let denom = (1+u²). Then:
        // X(u) = X_num / denom  where X_num = (cx+r_cos0) + (-2*r_sin0)*u + (cx-r_cos0)*u²
        // Y(u) = Y_num / denom  where Y_num = (cy+r_sin0) + (2*r_cos0)*u + (cy-r_sin0)*u²
        //
        // F(X/denom, Y/denom) = 0
        // Multiply by denom²: F(X_num, Y_num) = 0
        // This is a quartic in u.

        let x0 = cx + r_cos0;
        let x1 = -2.0 * r_sin0;
        let x2 = cx - r_cos0;
        let y0 = cy + r_sin0;
        let y1 = 2.0 * r_cos0;
        let y2 = cy - r_sin0;

        // F(x, y) = A*x² + B*y² + 2C*x*y + 2D*x + 2E*y + F
        // x(u) = x0 + x1*u + x2*u²
        // y(u) = y0 + y1*u + y2*u²
        //
        // Expand to quartic in u.

        // Compute coefficients for x², y², xy, x, y
        let qa = |r: f64| -> [f64; 5] {
            // (r0 + r1*u + r2*u²)²
            // = r0² + 2*r0*r1*u + (2*r0*r2 + r1²)*u² + 2*r1*r2*u³ + r2²*u⁴
            [r * r, 0.0, 0.0, 0.0, 0.0]
        };

        let x2_coeff = quad_mul(x0, x1, x2, x0, x1, x2);
        let y2_coeff = quad_mul(y0, y1, y2, y0, y1, y2);
        let xy_coeff = quad_mul(x0, x1, x2, y0, y1, y2);
        // After multiplying F(X_num/denom, Y_num/denom) by denom², the linear
        // terms become X_num*denom and Y_num*denom (not just X_num/Y_num).
        let x_coeff = quad_mul(x0, x1, x2, 1.0, 0.0, 1.0);
        let y_coeff = quad_mul(y0, y1, y2, 1.0, 0.0, 1.0);

        // Constant term F → F*(1+u²)² = F*(1 + 2u² + u⁴)
        let f_const = [f, 0.0, 2.0 * f, 0.0, f];

        // Combine:
        let mut coeff = [0.0_f64; 5];
        for i in 0..5 {
            coeff[i] = a * x2_coeff[i]
                + b * y2_coeff[i]
                + 2.0 * cc * xy_coeff[i]
                + 2.0 * d * x_coeff[i]
                + 2.0 * e * y_coeff[i]
                + f_const[i];
        }

        // Reverse to descending degree order (highest first) to match the
        // solvers below (solve_quadratic_raw / solve_quartic_raw).
        coeff.reverse();

        // Normalize: drop leading zeros
        let mut start = 0;
        while start < 5 && coeff[start].abs() < TOL_DISC {
            start += 1;
        }

        if start >= 5 {
            // Degenerate case: no meaningful roots
            self.empt = true;
            self.nbp = 0;
            return;
        }

        // Check if the quartic is actually a quadratic or lower
        let degree = 4 - start;

        let (roots, n) = if degree <= 2 {
            if degree == 2 {
                solve_quadratic_raw(coeff[start], coeff[start + 1], coeff[start + 2])
            } else if degree == 1 {
                if coeff[start].abs() > TOL_DISC {
                    ([-coeff[start + 1] / coeff[start]; 4], 1)
                } else {
                    ([0.0; 4], 0)
                }
            } else {
                ([0.0; 4], 0)
            }
        } else {
            solve_quartic_raw(coeff[start], coeff[start + 1], coeff[start + 2], coeff[start + 3], coeff[start + 4])
        };

        // Convert root values (u = tan(θ/2)) back to (x, y, theta)
        let two_pi = 2.0 * std::f64::consts::PI;
        let mut pts = Vec::with_capacity(4);
        for i in 0..n {
            let u = roots[i];
            let denom = 1.0 + u * u;
            if denom.abs() < TOL_CONF {
                continue;
            }
            let x = (x0 + x1 * u + x2 * u * u) / denom;
            let y = (y0 + y1 * u + y2 * u * u) / denom;
            let theta = 2.0 * u.atan(); // θ = 2*arctan(u)
            let theta = theta.rem_euclid(two_pi);

            // Verify the point is actually on the conic
            if conic.value(x, y).abs() > TOL_CONF * 100.0 {
                continue;
            }

            pts.push((x, y, theta));
        }

        // Sort by theta
        pts.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

        // Deduplicate
        let mut j = 0;
        for i in 0..pts.len() {
            if j == 0 || (pts[i].0 - pts[j - 1].0).abs() > TOL_CONF || (pts[i].1 - pts[j - 1].1).abs() > TOL_CONF {
                if i != j {
                    pts[j] = pts[i];
                }
                j += 1;
            }
        }
        pts.truncate(j);

        self.nbp = pts.len().min(4);
        for i in 0..self.nbp {
            self.lpnt[i] = IntPoint2d::new_implicit(pts[i].0, pts[i].1, pts[i].2);
        }

        if self.nbp == 0 {
            self.empt = true;
        }
    }

    fn intersect_ellipse_conic(&mut self, ellipse: &Ellipse2d, conic: &Conic2d) {
        // Ellipse parametric: P(θ) = center + major_dir*a*cosθ + minor_dir*b*sinθ
        // Use tangent half-angle: u = tan(θ/2)
        // cosθ = (1-u²)/(1+u²), sinθ = 2u/(1+u²)
        //
        // x(u) = cx + (a*cos0)*(1-u²)/(1+u²) - (b*sin0)*2u/(1+u²)
        // y(u) = cy + (a*sin0)*(1-u²)/(1+u²) + (b*cos0)*2u/(1+u²)

        let cx = ellipse.center.x;
        let cy = ellipse.center.y;
        let a_r = ellipse.major_radius;
        let b_r = ellipse.minor_radius;
        let cos0 = ellipse.major_dir.x;
        let sin0 = ellipse.major_dir.y;

        let a_cos0 = a_r * cos0;
        let a_sin0 = a_r * sin0;
        let b_cos0 = b_r * cos0;
        let b_sin0 = b_r * sin0;

        let x0 = cx + a_cos0;
        let x1 = -2.0 * b_sin0;
        let x2 = cx - a_cos0;
        let y0 = cy + a_sin0;
        let y1 = 2.0 * b_cos0;
        let y2 = cy - a_sin0;

        self.solve_conic_param_quartic(conic, x0, x1, x2, y0, y1, y2);
    }

    fn intersect_parabola_conic(&mut self, parabola: &Parabola2d, conic: &Conic2d) {
        // Parabola in local frame: P(t) = origin + axis_dir*(t²/(2p)) + perp_dir*t
        // In global coords:
        // x(t) = ox + cos_a*(t²/(2p)) - sin_a*t
        // y(t) = oy + sin_a*(t²/(2p)) + cos_a*t
        //
        // This is a quadratic in t: x = x0 + x1*t + x2*t², y = y0 + y1*t + y2*t²
        // Substituting into conic → quartic in t

        let p = parabola.focal_param;
        let ox = parabola.origin.x;
        let oy = parabola.origin.y;
        let cos_a = parabola.axis_dir.x;
        let sin_a = parabola.axis_dir.y;

        let inv_2p = 1.0 / (2.0 * p);
        let x0 = ox;
        let x1 = -sin_a;
        let x2 = cos_a * inv_2p;
        let y0 = oy;
        let y1 = cos_a;
        let y2 = sin_a * inv_2p;

        self.solve_conic_param_quartic(conic, x0, x1, x2, y0, y1, y2);
    }

    fn intersect_hyperbola_conic(&mut self, hyperbola: &Hyperbola2d, conic: &Conic2d) {
        // Hyperbola rational parameterization:
        // u = tanh(t/2): x = cx + a*(1+u²)/(1-u²), y = cy + b*2u/(1-u²)
        // Using t-param: x = cx + a*cosh(t), y = cy + b*sinh(t)
        // For the rational form using u = tan(θ/2) substitution:
        // Use: x = a*(1+t²)/(1-t²), y = 2*b*t/(1-t²) centered at origin
        //
        // Actually, let's use the simple approach: parameterize the hyperbola in t.
        // x = cx + semi_major * cosh(t)
        // y = cy + semi_minor * sinh(t)
        // This is not rational, but we can use the transformation:
        // cosh(t) = (e^t + e^-t)/2, sinh(t) = (e^t - e^-t)/2
        //
        // For the general conic intersection, we use rational params:
        // x = cx + a*(1+t²)/(2t), y = cy + b*(1-t²)/(2t)
        // OR using the standard rational form:
        // x = a*(1+u²)/(1-u²), y = 2*b*u/(1-u²) (rectangle hyperbola)
        //
        // For real hyperbola: x²/a² - y²/b² = 1
        // Parametric: x = a*(1+t²)/(1-t²), y = 2*b*t/(1-t²)
        // where t = tanh(θ/2) ∈ (-1, 1)

        let cx = hyperbola.center.x;
        let cy = hyperbola.center.y;
        let a_r = hyperbola.semi_major;
        let b_r = hyperbola.semi_minor;
        let cos0 = hyperbola.major_dir.x;
        let sin0 = hyperbola.major_dir.y;

        // In local frame: u = a_r*(1+t²)/(1-t²), v = 2*b_r*t/(1-t²)
        // Then transform back: P = center + u*major_dir + v*minor_dir
        // minor_dir = (-sin0, cos0)

        // x(t) = cx + cos0*a_r*(1+t²)/(1-t²) - sin0*2*b_r*t/(1-t²)
        // y(t) = cy + sin0*a_r*(1+t²)/(1-t²) + cos0*2*b_r*t/(1-t²)
        //
        // Multiply numerator and denominator:
        // Let denom = 1-t²
        // X_num = cx*(1-t²) + a_r*cos0*(1+t²) - 2*b_r*sin0*t
        //       = (cx + a_r*cos0) + (-2*b_r*sin0)*t + (-cx + a_r*cos0)*t²
        // Y_num = cy*(1-t²) + a_r*sin0*(1+t²) + 2*b_r*cos0*t
        //       = (cy + a_r*sin0) + (2*b_r*cos0)*t + (-cy + a_r*sin0)*t²
        //
        // So x = X_num/(1-t²), y = Y_num/(1-t²)
        // F(X/denom, Y/denom) = 0 → F(X_num, Y_num) = 0 (multiply by denom²)
        // This gives a quartic in t.

        let x0 = cx + a_r * cos0;
        let x1 = -2.0 * b_r * sin0;
        let x2 = -cx + a_r * cos0;
        let y0 = cy + a_r * sin0;
        let y1 = 2.0 * b_r * cos0;
        let y2 = -cy + a_r * sin0;

        self.solve_conic_param_quartic(conic, x0, x1, x2, y0, y1, y2);
    }

    /// Generic solver for conic intersection where the first curve is
    /// parameterized as a quadratic rational: x(t) = (x0 + x1*t + x2*t²) / (1 + t²)
    /// or more generally as a polynomial: x(t) = x0 + x1*t + x2*t², y(t) = y0 + y1*t + y2*t²
    /// and substituted into the conic's implicit equation.
    fn solve_conic_param_quartic(
        &mut self,
        conic: &Conic2d,
        x0: f64, x1: f64, x2: f64,
        y0: f64, y1: f64, y2: f64,
    ) {
        let a = conic.a;
        let b = conic.b;
        let cc = conic.c;
        let d = conic.d;
        let e = conic.e;
        let f = conic.f;

        // x (t) = x0 + x1*t + x2*t²
        // y (t) = y0 + y1*t + y2*t²
        //
        // F(x(t), y(t)) = A*x² + B*y² + 2C*x*y + 2D*x + 2E*y + F

        let x2_c = quad_mul(x0, x1, x2, x0, x1, x2);
        let y2_c = quad_mul(y0, y1, y2, y0, y1, y2);
        let xy_c = quad_mul(x0, x1, x2, y0, y1, y2);
        let x_c = [x0, x1, x2, 0.0, 0.0];
        let y_c = [y0, y1, y2, 0.0, 0.0];

        let mut coeff = [0.0_f64; 5];
        for i in 0..5 {
            coeff[i] = a * x2_c[i] + b * y2_c[i] + 2.0 * cc * xy_c[i] + 2.0 * d * x_c[i] + 2.0 * e * y_c[i];
        }
        // Add constant F (multiply by 1 = (1+0*t+0*t²) for denom consistency)
        coeff[0] += f;

        // Solve the quartic
        let degree = drop_leading_zeros(&mut coeff);
        let (roots, n) = if degree <= 2 {
            solve_quadratic_raw(coeff[0], coeff[1], coeff[2])
        } else if degree == 3 {
            let rc = solve_cubic_raw(coeff[0], coeff[1], coeff[2], coeff[3]);
            (rc.0, rc.1)
        } else {
            solve_quartic_raw(coeff[0], coeff[1], coeff[2], coeff[3], coeff[4])
        };

        // Convert roots to (x, y)
        let mut pts = Vec::with_capacity(4);
        for i in 0..n {
            let t = roots[i];
            let x = x0 + x1 * t + x2 * t * t;
            let y = y0 + y1 * t + y2 * t * t;
            // Verify
            if conic.value(x, y).abs() > TOL_CONF * 100.0 {
                continue;
            }
            pts.push((x, y, t));
        }

        // Sort by t
        pts.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

        // Deduplicate
        let mut j = 0;
        for i in 0..pts.len() {
            if j == 0 || (pts[i].0 - pts[j - 1].0).abs() > TOL_CONF || (pts[i].1 - pts[j - 1].1).abs() > TOL_CONF {
                if i != j {
                    pts[j] = pts[i];
                }
                j += 1;
            }
        }
        pts.truncate(j);

        self.nbp = pts.len().min(4);
        for i in 0..self.nbp {
            self.lpnt[i] = IntPoint2d::new_implicit(pts[i].0, pts[i].1, pts[i].2);
        }

        if self.nbp == 0 {
            self.empt = true;
        }
    }

    // ========================================================================
    // Query methods
    // ========================================================================

    /// Returns `true` if the computation was performed.
    ///
    /// OCCT: `IsDone()`.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// Returns `true` when there is no intersection.
    ///
    /// OCCT: `IsEmpty()`.
    pub fn is_empty(&self) -> bool {
        if !self.done {
            panic!("AnaIntersection2d: not done");
        }
        self.nbp == 0 && !self.iden
    }

    /// Returns `true` if the elements are identical (coincident).
    ///
    /// OCCT: `IdenticalElements()`.
    pub fn identical_elements(&self) -> bool {
        if !self.done {
            panic!("AnaIntersection2d: not done");
        }
        self.iden
    }

    /// Returns `true` if the elements are parallel.
    ///
    /// OCCT: `ParallelElements()`.
    pub fn parallel_elements(&self) -> bool {
        if !self.done {
            panic!("AnaIntersection2d: not done");
        }
        self.para
    }

    /// Returns the number of intersection points.
    ///
    /// OCCT: `NbPoints()`.
    pub fn nb_points(&self) -> usize {
        if !self.done {
            panic!("AnaIntersection2d: not done");
        }
        self.nbp
    }

    /// Returns the N-th intersection point (1-indexed).
    ///
    /// OCCT: `Point(N)`.
    /// Panics if N is out of range (N < 1 or N > NbPoints).
    pub fn point(&self, n: usize) -> &IntPoint2d {
        if !self.done {
            panic!("AnaIntersection2d: not done");
        }
        if n == 0 || n > self.nbp {
            panic!("AnaIntersection2d: point index {} out of range (1..{})", n, self.nbp);
        }
        &self.lpnt[n - 1]
    }
}

impl Default for AnaIntersection2d {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Outils — polynomial root solvers and utility functions
// ============================================================================

/// Compute the product (p0 + p1*t + p2*t²) * (q0 + q1*t + q2*t²) as coefficients
/// of a quartic: result[0] + result[1]*t + result[2]*t² + result[3]*t³ + result[4]*t⁴
fn quad_mul(p0: f64, p1: f64, p2: f64, q0: f64, q1: f64, q2: f64) -> [f64; 5] {
    [
        p0 * q0,
        p0 * q1 + p1 * q0,
        p0 * q2 + p1 * q1 + p2 * q0,
        p1 * q2 + p2 * q1,
        p2 * q2,
    ]
}

/// Drop leading zeros from coefficients (degree 4 polynomial).
/// Returns the actual degree.
fn drop_leading_zeros(coeff: &mut [f64; 5]) -> usize {
    let mut start = 0;
    while start < 5 && coeff[start].abs() < TOL_DISC {
        start += 1;
    }
    if start > 0 && start < 5 {
        let d = 4 - start;
        for i in 0..=d {
            coeff[i] = coeff[start + i];
        }
        for i in d + 1..5 {
            coeff[i] = 0.0;
        }
    }
    4 - start
}

/// Solve quadratic ax² + bx + c = 0 and populate intersection result.
fn solve_quadratic_intersection(
    a: f64, b: f64, c: f64,
    lpnt: &mut [IntPoint2d; 4],
    nbp: &mut usize,
    para: &mut bool,
    iden: &mut bool,
    empt: &mut bool,
) {
    *nbp = 0;
    *para = false;
    *iden = false;
    *empt = false;

    if a.abs() < TOL_DISC {
        if b.abs() < TOL_DISC {
            if c.abs() < TOL_DISC {
                *iden = true; // All points satisfy (identity)
            } else {
                *empt = true; // No solution
            }
            return;
        }
        // Linear: b*t + c = 0
        let t = -c / b;
        if t.is_finite() {
            *nbp = 1;
            lpnt[0] = IntPoint2d::new(0.0, 0.0, t, 0.0);
        }
        return;
    }

    let disc = b * b - 4.0 * a * c;
    if disc < -TOL_DISC {
        *empt = true;
        return;
    }

    if disc.abs() < TOL_DISC {
        let t = -b / (2.0 * a);
        *nbp = 1;
        lpnt[0] = IntPoint2d::new(0.0, 0.0, t, 0.0);
        return;
    }

    let sqrt_disc = disc.sqrt();
    let t1 = (-b - sqrt_disc) / (2.0 * a);
    let t2 = (-b + sqrt_disc) / (2.0 * a);
    *nbp = 2;
    lpnt[0] = IntPoint2d::new(0.0, 0.0, t1.min(t2), 0.0);
    lpnt[1] = IntPoint2d::new(0.0, 0.0, t1.max(t2), 0.0);
}

/// Solve quadratic ax² + bx + c = 0, return raw roots.
fn solve_quadratic_raw(a: f64, b: f64, c: f64) -> ([f64; 4], usize) {
    if a.abs() < TOL_DISC {
        if b.abs() < TOL_DISC {
            return ([0.0; 4], 0);
        }
        return ([-c / b, 0.0, 0.0, 0.0], 1);
    }
    let disc = b * b - 4.0 * a * c;
    if disc < -TOL_DISC {
        return ([0.0; 4], 0);
    }
    if disc.abs() < TOL_DISC {
        return ([-b / (2.0 * a), 0.0, 0.0, 0.0], 1);
    }
    let sqrt_disc = disc.sqrt();
    let t1 = (-b - sqrt_disc) / (2.0 * a);
    let t2 = (-b + sqrt_disc) / (2.0 * a);
    if t1 < t2 {
        ([t1, t2, 0.0, 0.0], 2)
    } else {
        ([t2, t1, 0.0, 0.0], 2)
    }
}

/// Solve cubic a*x³ + b*x² + c*x + d = 0, return real roots.
fn solve_cubic_raw(a: f64, b: f64, c: f64, d: f64) -> ([f64; 4], usize) {
    if a.abs() < TOL_DISC {
        return solve_quadratic_raw(b, c, d);
    }

    // Normalize to monic: x³ + p*x² + q*x + r = 0
    let p = b / a;
    let q = c / a;
    let r = d / a;

    // Substitute x = t - p/3 to get depressed cubic: t³ + α*t + β = 0
    let alpha = q - p * p / 3.0;
    let beta = 2.0 * p * p * p / 27.0 - p * q / 3.0 + r;

    let disc = beta * beta / 4.0 + alpha * alpha * alpha / 27.0;
    let offset = p / 3.0;

    if disc > TOL_DISC {
        // One real root
        let sqrt_disc = disc.sqrt();
        let u = (-beta / 2.0 + sqrt_disc).cbrt();
        let v = (-beta / 2.0 - sqrt_disc).cbrt();
        ([u + v - offset, 0.0, 0.0, 0.0], 1)
    } else if disc < -TOL_DISC {
        // Three distinct real roots
        let radius = (-alpha * alpha * alpha / 27.0).sqrt();
        let theta = (-beta / (2.0 * radius)).acos() / 3.0;
        let cbrt_r = radius.cbrt();
        (
            [
                2.0 * cbrt_r * theta.cos() - offset,
                2.0 * cbrt_r * (theta + 2.0 * std::f64::consts::FRAC_PI_3).cos() - offset,
                2.0 * cbrt_r * (theta + 4.0 * std::f64::consts::FRAC_PI_3).cos() - offset,
                0.0,
            ],
            3,
        )
    } else {
        // One or two real roots (repeated)
        if alpha.abs() < TOL_DISC {
            ([0.0; 4], 0)
        } else {
            let u = (-beta / 2.0).cbrt();
            ([2.0 * u - offset, -u - offset, 0.0, 0.0], 2)
        }
    }
}

/// Solve quartic a*x⁴ + b*x³ + c*x² + d*x + e = 0 using Ferrari's method.
fn solve_quartic_raw(a: f64, b: f64, c: f64, d: f64, e: f64) -> ([f64; 4], usize) {
    if a.abs() < TOL_DISC {
        return solve_cubic_raw(b, c, d, e);
    }

    // Normalize to monic: x⁴ + p*x³ + q*x² + r*x + s = 0
    let p = b / a;
    let q = c / a;
    let r = d / a;
    let s = e / a;

    // Substitute x = y - p/4 to eliminate cubic term
    let p2 = p * p;
    let p3 = p2 * p;
    let p4 = p3 * p;

    let alpha = q - 3.0 * p2 / 8.0;
    let beta = r + p3 / 8.0 - p * q / 2.0;
    let gamma = s - 3.0 * p4 / 256.0 + p2 * q / 16.0 - p * r / 4.0;

    if beta.abs() < TOL_DISC {
        // Depressed quartic: y⁴ + α*y² + γ = 0
        // This is quadratic in y²
        let disc = alpha * alpha - 4.0 * gamma;
        if disc < -TOL_DISC {
            return ([0.0; 4], 0);
        }
        let mut roots = Vec::with_capacity(4);
        let y1_sq = if disc.abs() < TOL_DISC {
            -alpha / 2.0
        } else {
            let sqrt_disc = disc.sqrt();
            let s1 = (-alpha + sqrt_disc) / 2.0;
            let s2 = (-alpha - sqrt_disc) / 2.0;
            if s1 >= 0.0 {
                let y = s1.sqrt();
                roots.push(y - p / 4.0);
                roots.push(-y - p / 4.0);
            }
            s2
        };
        if y1_sq >= 0.0 {
            let y = y1_sq.sqrt();
            roots.push(y - p / 4.0);
            roots.push(-y - p / 4.0);
        }
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let n = roots.len().min(4);
        let mut result = [0.0_f64; 4];
        for i in 0..n {
            result[i] = roots[i];
        }
        return (result, n);
    }

    // General Ferrari method
    // Find a real root t of the resolvent cubic: t³ - α*t² - 4γ*t + (4αγ - β²) = 0
    let resolvent = solve_cubic_raw(1.0, -alpha, -4.0 * gamma, 4.0 * alpha * gamma - beta * beta);

    // Find a positive real root
    let t = resolvent.0.iter().copied()
        .find(|&r| r > TOL_ROOT)
        .unwrap_or_else(|| {
            resolvent.0.iter().copied()
                .find(|&r| r > -TOL_ROOT)
                .unwrap_or(0.0)
        });

    if t < -TOL_ROOT {
        return ([0.0; 4], 0);
    }

    let sqrt_t = t.max(0.0).sqrt();
    let alpha_plus_t2 = alpha + t;

    // The quartic factors into two quadratics:
    // y² + sqrt_t*y + (alpha + t)/2 + β/(2*sqrt_t) = 0  [for sqrt_t > 0]
    // y² - sqrt_t*y + (alpha + t)/2 - β/(2*sqrt_t) = 0
    // OR if sqrt_t is near zero:
    // The factor is just from the depressed form

    let mut roots = Vec::with_capacity(4);
    let eps = 1e-10;

    if sqrt_t > eps {
        let inv_2t = 1.0 / (2.0 * sqrt_t);
        let ap = (alpha_plus_t2) / 2.0;
        let bp = beta * inv_2t;
        let am = (alpha_plus_t2) / 2.0;
        let bm = -beta * inv_2t;

        // Solve first quadratic: y² + sqrt_t*y + ap + bp = 0
        let (q1, n1) = solve_quadratic_raw(1.0, sqrt_t, ap + bp);
        for i in 0..n1 {
            roots.push(q1[i] - p / 4.0);
        }

        // Solve second quadratic: y² - sqrt_t*y + am + bm = 0
        let (q2, n2) = solve_quadratic_raw(1.0, -sqrt_t, am + bm);
        for i in 0..n2 {
            roots.push(q2[i] - p / 4.0);
        }
    } else {
        // sqrt_t ≈ 0: use alternative approach
        let (q, n) = solve_quadratic_raw(1.0, 0.0, alpha_plus_t2);
        for i in 0..n {
            let y = q[i];
            // Solve y² + β/(2*y) = 0 or similar
            roots.push(y - p / 4.0);
            if y.abs() > eps {
                let y2 = -beta / (2.0 * y);
                roots.push(y2 - p / 4.0);
            }
        }
    }

    // Deduplicate and sort
    roots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    roots.dedup_by(|a, b| (*a - *b).abs() < TOL_ROOT);

    let n = roots.len().min(4);
    let mut result = [0.0_f64; 4];
    for i in 0..n {
        result[i] = roots[i];
    }
    (result, n)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_approx(a: f64, b: f64, eps: f64) {
        assert!((a - b).abs() < eps, "expected {} ≈ {}, eps={}", a, b, eps);
    }

    #[test]
    fn test_conic_from_line() {
        let line = Line2d::new(DVec2::ZERO, DVec2::X);
        let conic = Conic2d::from_line(&line);
        // Point (0, 1) is at distance 1 from line y=0
        assert_approx(conic.value(0.0, 1.0), 1.0, 1e-10); // should be on the line? No, (0,1) is off
        // Actually y=0 line has normal (0,1), so the value at (0,1) = 1
        // At (0,0): should be 0 (on the line)
        assert_approx(conic.value(0.0, 0.0), 0.0, 1e-10);
    }

    #[test]
    fn test_conic_from_circle() {
        let circle = Circle2d::new(DVec2::ZERO, 5.0);
        let conic = Conic2d::from_circle(&circle);
        // At (5, 0): should be 0
        assert_approx(conic.value(5.0, 0.0), 0.0, 1e-10);
        // At (0, 0): should be -R² = -25
        assert_approx(conic.value(0.0, 0.0), -25.0, 1e-10);
    }

    #[test]
    fn test_int_point_basic() {
        let p = IntPoint2d::new(1.0, 2.0, 0.5, 1.5);
        assert_approx(p.value().x, 1.0, 1e-15);
        assert_approx(p.param_on_first(), 0.5, 1e-15);
        assert_approx(p.param_on_second(), 1.5, 1e-15);
        assert!(!p.second_is_implicit());
    }

    #[test]
    fn test_intersection_line_line() {
        let l1 = Line2d::new(DVec2::ZERO, DVec2::X);
        let l2 = Line2d::new(DVec2::new(0.0, 1.0), DVec2::new(1.0, -1.0));
        let mut inter = AnaIntersection2d::new();
        inter.perform_lin_lin(&l1, &l2);
        assert!(inter.is_done());
        assert_eq!(inter.nb_points(), 1);
        let p = inter.point(1);
        assert_approx(p.value().x, 1.0, 1e-10);
        assert_approx(p.value().y, 0.0, 1e-10);
    }

    #[test]
    fn test_intersection_line_line_parallel() {
        let l1 = Line2d::new(DVec2::ZERO, DVec2::X);
        let l2 = Line2d::new(DVec2::new(0.0, 1.0), DVec2::X);
        let mut inter = AnaIntersection2d::new();
        inter.perform_lin_lin(&l1, &l2);
        assert!(inter.is_done());
        assert!(inter.parallel_elements());
        assert!(inter.is_empty());
    }

    #[test]
    fn test_intersection_line_circle() {
        let line = Line2d::new(DVec2::new(-2.0, 0.0), DVec2::X);
        let circle = Circle2d::new(DVec2::ZERO, 1.0);
        let mut inter = AnaIntersection2d::new();
        inter.perform_lin_circ(&line, &circle);
        assert!(inter.is_done());
        assert_eq!(inter.nb_points(), 2);
    }

    #[test]
    fn test_intersection_circle_circle() {
        let c1 = Circle2d::new(DVec2::ZERO, 5.0);
        // d = 3, |5 - 4| = 1 < 3 < 9 = 5 + 4 → two intersection points
        let c2 = Circle2d::new(DVec2::new(3.0, 0.0), 4.0);
        let mut inter = AnaIntersection2d::new();
        inter.perform_circ_circ(&c1, &c2);
        assert!(inter.is_done());
        assert_eq!(inter.nb_points(), 2);
    }

    #[test]
    fn test_intersection_line_conic() {
        let line = Line2d::new(DVec2::new(-3.0, 0.0), DVec2::X);
        let circle = Circle2d::new(DVec2::ZERO, 2.0);
        let conic = Conic2d::from_circle(&circle);
        let mut inter = AnaIntersection2d::new();
        inter.perform_lin_conic(&line, &conic);
        assert!(inter.is_done());
        assert_eq!(inter.nb_points(), 2);
        // Should hit circle at x = ±2
        assert_approx(inter.point(1).value().x.abs(), 2.0, 1e-7);
    }

    #[test]
    fn test_intersection_circle_conic() {
        let circle = Circle2d::new(DVec2::ZERO, 5.0);
        let conic = Conic2d::from_circle(&Circle2d::new(DVec2::new(3.0, 0.0), 4.0));
        let mut inter = AnaIntersection2d::new();
        inter.perform_circ_conic(&circle, &conic);
        assert!(inter.is_done());
        // Two circles should intersect at 2 points
        assert_eq!(inter.nb_points(), 2);
    }

    #[test]
    fn test_solve_quadratic() {
        // x² - 3x + 2 = 0 → x = 1, 2
        let (roots, n) = solve_quadratic_raw(1.0, -3.0, 2.0);
        assert_eq!(n, 2);
        assert_approx(roots[0], 1.0, 1e-12);
        assert_approx(roots[1], 2.0, 1e-12);
    }

    #[test]
    fn test_solve_cubic() {
        // (x-1)(x-2)(x-3) = x³ - 6x² + 11x - 6 = 0
        let (roots, n) = solve_cubic_raw(1.0, -6.0, 11.0, -6.0);
        assert_eq!(n, 3);
        let mut sorted = vec![roots[0], roots[1], roots[2]];
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_approx(sorted[0], 1.0, 1e-10);
        assert_approx(sorted[1], 2.0, 1e-10);
        assert_approx(sorted[2], 3.0, 1e-10);
    }
}
