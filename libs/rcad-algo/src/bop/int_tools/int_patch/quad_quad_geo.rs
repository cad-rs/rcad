//! IntAna_QuadQuadGeo — geometric intersections between two quadric surfaces.
//!
//! OCCT IntAna_QuadQuadGeo.hxx / .cxx
//!
//! Computes closed-form intersection curves between:
//! Plane, Cylinder, Sphere, Cone, Torus (all 15 pair combinations).
//!
//! Results are classified as: Point, Line, Circle, Ellipse, Parabola, Hyperbola,
//! Empty, Same, NoGeometricSolution.

use crate::topalgo::int_surf::quadric::Quadric;
use glam::DVec3;
use rcad_kernel::geom::{Circle3, Curve3, Ellipse3, Hyperbola3, Line3, Parabola3, Plane};

use super::int_quad_quad::any_perpendicular_axis;

/// OCCT Precision-based clamping minimum (IntAna_QuadQuadGeo InitTolerances).
const TOLERANCE_CLAMP_MIN: f64 = 1e-15;
/// Guard against dividing by a near-zero squared length.
const TOLERANCE_LEN_SQ_DIV_SAFE: f64 = 1e-30;

/// OCCT IntAna_ResultType.hxx
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnaResultType {
    Point,
    Line,
    Circle,
    PointAndCircle,
    Ellipse,
    Parabola,
    Hyperbola,
    Empty,
    Same,
    NoGeometricSolution,
}

/// IntAna_QuadQuadGeo
///
/// OCCT fields (L256-280):
///   done, nbint, typeres
///   pt1-4, dir1-4  — result points and directions (up to 4 solutions)
///   param1-4, param1bis, param2bis  — curve parameters (radius, angle, etc.)
///   myEPSILON_*  — internal tolerances
///   myCommonGen, myPChar  — common generator flag and characteristic point
pub struct QuadQuadGeo {
    // OCCT L256-258
    done: bool,
    nbint: i32,
    typeres: AnaResultType,
    // OCCT L259-266: result geometry (up to 4 points/directions)
    pt1: DVec3,
    pt2: DVec3,
    pt3: DVec3,
    pt4: DVec3,
    dir1: DVec3,
    dir2: DVec3,
    dir3: DVec3,
    dir4: DVec3,
    // OCCT L267-272: curve parameters
    param1: f64,
    param2: f64,
    param3: f64,
    param4: f64,
    param1bis: f64,
    param2bis: f64,
    // OCCT L273-278: internal tolerances
    my_epsilon_distance: f64,
    my_epsilon_angle_cone: f64,
    my_epsilon_mini_circle_radius: f64,
    my_epsilon_cylinder_delta_radius: f64,
    my_epsilon_cylinder_delta_distance: f64,
    my_epsilon_axes_para: f64,
    // OCCT L279-280
    my_common_gen: bool,
    my_p_char: DVec3,
}

impl QuadQuadGeo {
    pub fn new() -> Self {
        Self {
            done: false,
            nbint: 0,
            typeres: AnaResultType::Empty,
            pt1: DVec3::ZERO,
            pt2: DVec3::ZERO,
            pt3: DVec3::ZERO,
            pt4: DVec3::ZERO,
            dir1: DVec3::Z,
            dir2: DVec3::Z,
            dir3: DVec3::Z,
            dir4: DVec3::Z,
            param1: 0.0,
            param2: 0.0,
            param3: 0.0,
            param4: 0.0,
            param1bis: 0.0,
            param2bis: 0.0,
            my_epsilon_distance: 0.0,
            my_epsilon_angle_cone: 0.0,
            my_epsilon_mini_circle_radius: 0.0,
            my_epsilon_cylinder_delta_radius: 0.0,
            my_epsilon_cylinder_delta_distance: 0.0,
            my_epsilon_axes_para: 0.0,
            my_common_gen: false,
            my_p_char: DVec3::ZERO,
        }
    }

    /// OCCT L254: InitTolerances
    pub fn init_tolerances(&mut self) {
        self.my_epsilon_distance = 1e-10;
        self.my_epsilon_angle_cone = 1e-12;
        self.my_epsilon_mini_circle_radius = TOLERANCE_CLAMP_MIN;
        self.my_epsilon_cylinder_delta_radius = 1e-10;
        self.my_epsilon_cylinder_delta_distance = 1e-10;
        self.my_epsilon_axes_para = 1e-10;
    }

    // ---- Accessors (OCCT L214-250) ----

    pub fn is_done(&self) -> bool {
        self.done
    }
    pub fn type_inter(&self) -> AnaResultType {
        self.typeres
    }
    pub fn nb_solutions(&self) -> i32 {
        self.nbint
    }
    pub fn has_common_gen(&self) -> bool {
        self.my_common_gen
    }
    pub fn p_char(&self) -> DVec3 {
        self.my_p_char
    }

    pub fn point(&self, num: i32) -> DVec3 {
        match num {
            1 => self.pt1,
            2 => self.pt2,
            3 => self.pt3,
            _ => self.pt4,
        }
    }

    pub fn line(&self, num: i32) -> Line3 {
        match num {
            1 => Line3 {
                origin: self.pt1,
                direction: self.dir1,
            },
            2 => Line3 {
                origin: self.pt2,
                direction: self.dir2,
            },
            _ => Line3 {
                origin: self.pt3,
                direction: self.dir3,
            },
        }
    }

    /// OCCT IntAna_QuadQuadGeo::Circle(Index) — gp_Circ(DirToAx2(ptN, dirN),
    /// paramN).  The circle normal is `dirN`; the X/Y frame is re-derived
    /// perpendicular to it by OCCT's DirToAx2 (IntAna_QuadQuadGeo.cxx L237-257).
    /// (Previously the normal was mapped from `dir3`, which is wrong — dir3 is a
    /// plane-normal scratch slot, not the conic normal.)
    pub fn circle_n(&self, num: i32) -> Circle3 {
        let (pt, dir, param) = match num {
            2 => (self.pt2, self.dir2, self.param2),
            3 => (self.pt3, self.dir3, self.param3),
            4 => (self.pt4, self.dir4, self.param4),
            _ => (self.pt1, self.dir1, self.param1),
        };
        let normal = dir.normalize_or_zero();
        let x = normal.x;
        let y = normal.y;
        let z = normal.z;
        let ax = x.abs();
        let ay = y.abs();
        let az = z.abs();
        let v = if ax == 0.0 || (ax < ay && ax < az) {
            DVec3::new(0.0, -z, y)
        } else if ay == 0.0 || (ay < ax && ay < az) {
            DVec3::new(-z, 0.0, x)
        } else {
            DVec3::new(-y, x, 0.0)
        };
        let x_dir = v.normalize_or_zero();
        let y_dir = normal.cross(x_dir).normalize_or_zero();
        Circle3 {
            center: pt,
            normal,
            x_dir,
            y_dir,
            radius: param,
        }
    }

    /// OCCT IntAna_QuadQuadGeo::Circle(1) — single-solution shorthand.
    pub fn circle(&self) -> Circle3 {
        self.circle_n(1)
    }

    /// OCCT IntAna_QuadQuadGeo::Ellipse(1) — gp_Elips(gp_Ax2(pt1, dir1, dir2),
    /// R1, R2) with R1 = max(param1, param1bis), R2 = min(param1, param1bis).
    pub fn ellipse(&self) -> Ellipse3 {
        self.ellipse_n(1)
    }

    /// OCCT IntAna_QuadQuadGeo::Ellipse(Index) — the two ellipse solutions use
    /// pt1/param1/param1bis (Index 1, frame (dir1, dir2)) and
    /// pt2/param2/param2bis (Index 2, frame (dir2, dir1)) respectively.
    pub fn ellipse_n(&self, num: i32) -> Ellipse3 {
        let (center, dir_normal, dir_in_plane, p, p_bis) = if num == 2 {
            (self.pt2, self.dir2, self.dir1, self.param2, self.param2bis)
        } else {
            (self.pt1, self.dir1, self.dir2, self.param1, self.param1bis)
        };
        let normal = dir_normal.normalize_or_zero();
        let raw = dir_in_plane - normal * dir_in_plane.dot(normal);
        let major_dir = raw.normalize_or_zero();
        let (major, minor) = if p >= p_bis { (p, p_bis) } else { (p_bis, p) };
        Ellipse3 {
            center,
            normal,
            major_dir,
            major_radius: major,
            minor_radius: minor,
        }
    }

    /// OCCT IntAna_QuadQuadGeo::Parabola(1) — gp_Parab(gp_Ax2(pt1, dir1, dir2),
    /// param1).
    pub fn parabola(&self) -> Parabola3 {
        let normal = self.dir1.normalize_or_zero();
        let raw = self.dir2 - normal * self.dir2.dot(normal);
        let axis_dir = raw.normalize_or_zero();
        Parabola3 {
            vertex: self.pt1,
            normal,
            axis_dir,
            focal_param: self.param1,
        }
    }

    /// OCCT IntAna_QuadQuadGeo::Hyperbola(Index) — gp_Hypr(gp_Ax2(ptN, dir1,
    /// dir2), paramN, paramNbis).  The two branches use pt1/param1/param1bis and
    /// pt2/param2/param2bis respectively.
    pub fn hyperbola_n(&self, num: i32) -> Hyperbola3 {
        let (pt, param, param_bis) = match num {
            2 => (self.pt2, self.param2, self.param2bis),
            _ => (self.pt1, self.param1, self.param1bis),
        };
        let normal = self.dir1.normalize_or_zero();
        let raw = self.dir2 - normal * self.dir2.dot(normal);
        let major_dir = raw.normalize_or_zero();
        Hyperbola3 {
            center: pt,
            normal,
            major_dir,
            semi_major: param,
            semi_minor: param_bis,
        }
    }

    /// OCCT IntAna_QuadQuadGeo::Hyperbola(1) — single-solution shorthand.
    pub fn hyperbola(&self) -> Hyperbola3 {
        self.hyperbola_n(1)
    }

    // ---- Convert result to Curve3 for rcad integration ----
    pub fn to_curves(&self) -> Vec<Curve3> {
        if !self.done
            || self.typeres == AnaResultType::Empty
            || self.typeres == AnaResultType::NoGeometricSolution
        {
            return vec![];
        }
        let mut curves = Vec::new();
        match self.typeres {
            AnaResultType::Line => {
                for i in 1..=self.nbint {
                    let l = self.line(i);
                    curves.push(Curve3::Line(l));
                }
            }
            AnaResultType::Circle | AnaResultType::PointAndCircle => {
                // OCCT IntCoCo/IntCyCo/IntCySp/IntCoSp: one GLine per solution.
                for i in 1..=self.nbint {
                    curves.push(Curve3::Circle(self.circle_n(i)));
                }
            }
            AnaResultType::Ellipse => {
                curves.push(Curve3::Ellipse(self.ellipse()));
            }
            AnaResultType::Parabola => {
                curves.push(Curve3::Parabola(self.parabola()));
            }
            AnaResultType::Hyperbola => {
                // OCCT IntCoCo/IntCyCo/IntCoSp Hyperbola case: both branches.
                for i in 1..=self.nbint {
                    curves.push(Curve3::Hyperbola(self.hyperbola_n(i)));
                }
            }
            _ => {}
        }
        curves
    }

    // =====================================================================
    // OCCT Perform methods — one per quadric pair
    // =====================================================================

    /// OCCT IntAna_QuadQuadGeo.cxx L389-512: Perform(P1, P2, TolAng, Tol) — Plane/Plane
    pub fn perform_plane_plane(&mut self, p1: &Quadric, p2: &Quadric, tol_ang: f64, tol: f64) {
        self.init_tolerances();
        self.done = false;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;
        self.param2bis = 0.0;

        // OCCT L394: double A1,B1,C1,D1, A2,B2,C2,D2, dist1, dist2, aMVD
        let (a1, b1, c1, d1) = p1.plane_coeffs();
        let (a2, b2, c2, d2) = p2.plane_coeffs();

        // OCCT L402-404: gp_Vec aVN1(A1,B1,C1), aVN2(A2,B2,C2), vd = aVN1.Crossed(aVN2)
        let a_vn1 = DVec3::new(a1, b1, c1);
        let a_vn2 = DVec3::new(a2, b2, c2);
        let vd = a_vn1.cross(a_vn2);

        // OCCT L406-407: const gp_Pnt& aLocP1 = P1.Location(), aLocP2 = P2.Location()
        let a_loc_p1 = p1.axis_loc();
        let a_loc_p2 = p2.axis_loc();

        // OCCT L409-410: dist1 = A2*X + B2*Y + C2*Z + D2, dist2 = A1*X + B1*Y + C1*Z + D1
        let dist1 = a2 * a_loc_p1.x + b2 * a_loc_p1.y + c2 * a_loc_p1.z + d2;
        let dist2 = a1 * a_loc_p2.x + b1 * a_loc_p2.y + c1 * a_loc_p2.z + d1;

        // OCCT L412-417: if (aMVD <= TolAng) — normals collinear
        let a_mvd = vd.length();
        if a_mvd <= tol_ang {
            self.typeres = if dist1.abs() <= tol && dist2.abs() <= tol {
                AnaResultType::Same
            } else {
                AnaResultType::Empty
            };
            self.done = true;
            return;
        }

        // OCCT L420-447: compute intersection line
        let a_eps = 1e-16;
        let denom = a1 * a2 + b1 * b2 + c1 * c2;
        let denom2 = denom * denom;
        let mut ddenom = 1.0 - denom2;
        if ddenom.abs() <= a_eps {
            ddenom = a_eps;
        }

        let par1 = dist1 / ddenom;
        let par2 = -dist2 / ddenom;

        let inter1 = a_vn1.cross(vd);
        let inter2 = a_vn2.cross(vd);

        let x1 = a_loc_p1.x + par1 * inter1.x;
        let y1 = a_loc_p1.y + par1 * inter1.y;
        let z1 = a_loc_p1.z + par1 * inter1.z;
        let x2 = a_loc_p2.x + par2 * inter2.x;
        let y2 = a_loc_p2.y + par2 * inter2.y;
        let z2 = a_loc_p2.z + par2 * inter2.z;

        // OCCT L443-446: pt1 = midpoint, dir1 = vd normalized
        self.pt1 = DVec3::new((x1 + x2) * 0.5, (y1 + y2) * 0.5, (z1 + z2) * 0.5);
        self.dir1 = vd.normalize_or_zero();
        self.typeres = AnaResultType::Line;
        self.nbint = 1;

        // OCCT L458-509: refine origin when angle between planes is small
        let a_tresh_ang = 2e-6;
        let a_tresh_dist = 1e-12;
        if a_mvd < a_tresh_ang {
            let a_dist1 = a1 * self.pt1.x + b1 * self.pt1.y + c1 * self.pt1.z + d1;
            let a_dist2 = a2 * self.pt1.x + b2 * self.pt1.y + c2 * self.pt1.z + d2;
            if a_dist1.abs() > a_tresh_dist || a_dist2.abs() > a_tresh_dist {
                // OCCT L475-486: Perform(line along plane1 normal through pt1, plane1)
                let a_dn1 = a_vn1.normalize_or_zero();
                let a_pt1 = self.pt1; // copy before refinement
                // Line: P = a_pt1 + t * a_dn1, intersect with plane1: a_vn1 · P + d1 = 0
                let t_param1 = -(a1 * a_pt1.x + b1 * a_pt1.y + c1 * a_pt1.z + d1)
                    / (a1 * a_dn1.x + b1 * a_dn1.y + c1 * a_dn1.z).max(1e-16);
                let a_pnt1 = a_pt1 + t_param1 * a_dn1;

                // OCCT L489-507: line along dir1×norm1 through a_pnt1, intersect plane2
                let a_dl2 = self.dir1.cross(a_dn1);
                let a_l2_dir = a_dl2.normalize_or_zero();
                // Line: P = a_pnt1 + t * a_l2_dir, intersect with plane2
                let denom2_l = a2 * a_l2_dir.x + b2 * a_l2_dir.y + c2 * a_l2_dir.z;
                if denom2_l.abs() > 1e-16 {
                    let t_param2 = -(a2 * a_pnt1.x + b2 * a_pnt1.y + c2 * a_pnt1.z + d2) / denom2_l;
                    self.pt1 = a_pnt1 + t_param2 * a_l2_dir;
                }
            }
        }

        self.done = true;
    }

    /// OCCT L105-109: Perform(P, C, TolAng, Tol, H) — Plane/Cylinder
    ///
    /// Rewritten to match OCCT IntAna_QuadQuadGeo.cxx L543-722 algorithm.
    /// 1. Angle-dependent tolerance adjustment (H for short cylinders).
    /// 2. If axis parallel to plane: 0/1/2 lines (generatrices).
    /// 3. If not parallel: Circle or Ellipse at the piercing point.
    pub fn perform_plane_cylinder(
        &mut self,
        p: &Quadric,
        c: &Quadric,
        tol_ang: f64,
        tol: f64,
        h: f64,
    ) {
        self.init_tolerances();
        self.done = true;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;
        self.param2bis = 0.0;

        let (a, b, cc, d) = p.plane_coeffs();
        let normp = DVec3::new(a, b, cc);
        let radius = c.radius();
        let axis_loc = c.axis_loc();
        let axis_dir = c.axis_dir();
        let (ax, ay, az) = (axis_loc.x, axis_loc.y, axis_loc.z);

        // OCCT L554-565: signed distance from axis origin to plane
        let dist = a * ax + b * ay + cc * az + d;

        // OCCT L571-598: angle-dependent tolerance adjustment
        let mut tolang = tol_ang;
        let mut toltang = tol;
        let mut newparams = false;

        let n_len = normp.length();
        let d_len = axis_dir.length();
        let cos_da = normp.dot(axis_dir) / (n_len * d_len);
        let d_a = cos_da.acos().abs();

        if d_a > std::f64::consts::FRAC_PI_4 {
            let dang = d_a - std::f64::consts::FRAC_PI_2;
            let dangle = dang.abs();
            if dangle > tol_ang {
                let sinda = dangle.sin().abs();
                let dif = (sinda - tol).abs();
                if dif < tol || (h > 0.0 && sinda * h < 2.0 * tol) {
                    tolang = sinda * 2.0;
                    toltang = tol.max(sinda * h * 1.01);
                    newparams = true;
                }
            }
        }

        // OCCT L600-601: IntAna_IntConicQuad parallel check
        let dot_nd = normp.dot(axis_dir);
        let is_parallel = dot_nd.abs() / (n_len * d_len) <= tolang.sin();

        if is_parallel {
            // OCCT L603-683: parallel -> 0/1/2 lines
            let omega = DVec3::new(ax - dist * a, ay - dist * b, az - dist * cc);
            let abs_dist = dist.abs();

            if (abs_dist - radius).abs() < tol {
                self.nbint = 1;
                self.pt1 = omega;
                self.dir1 = if newparams {
                    let omega_xyz_trnsl = axis_loc + 100.0 * axis_dir;
                    let distt =
                        a * omega_xyz_trnsl.x + b * omega_xyz_trnsl.y + cc * omega_xyz_trnsl.z + d;
                    let omega1 = DVec3::new(
                        omega_xyz_trnsl.x - distt * a,
                        omega_xyz_trnsl.y - distt * b,
                        omega_xyz_trnsl.z - distt * cc,
                    );
                    (omega1 - omega).normalize_or_zero()
                } else {
                    axis_dir
                };
                self.typeres = AnaResultType::Line;
            } else if abs_dist < radius {
                self.nbint = 2;
                let axey = axis_dir.cross(normp).normalize_or_zero();
                let hh = (radius * radius - dist * dist).sqrt().max(0.0);
                self.pt1 = omega - hh * axey;
                self.pt2 = omega + hh * axey;

                if newparams {
                    let omega_xyz_trnsl = axis_loc + 100.0 * axis_dir;
                    let distt =
                        a * omega_xyz_trnsl.x + b * omega_xyz_trnsl.y + cc * omega_xyz_trnsl.z + d;
                    let an_sqrt_arg = radius * radius - distt * distt;
                    let ht = if an_sqrt_arg > 0.0 {
                        an_sqrt_arg.sqrt()
                    } else {
                        0.0
                    };
                    let omega1 = DVec3::new(
                        omega_xyz_trnsl.x - distt * a,
                        omega_xyz_trnsl.y - distt * b,
                        omega_xyz_trnsl.z - distt * cc,
                    );
                    self.dir1 = (omega1 - ht * axey - self.pt1).normalize_or_zero();
                    self.dir2 = (omega1 + ht * axey - self.pt2).normalize_or_zero();
                } else {
                    self.dir1 = axis_dir;
                    self.dir2 = axis_dir;
                }
                self.typeres = AnaResultType::Line;
            }
            // else: dist > radius -> keep Empty
        } else {
            // OCCT L685-719: not parallel -> circle or ellipse
            self.nbint = 1;

            // Piercing point of cylinder axis through plane
            let denom = normp.dot(axis_dir);
            let t_param = if denom.abs() > 1e-15 {
                -(normp.dot(axis_loc) + d) / denom
            } else {
                0.0
            };
            let pierce_pt = axis_loc + t_param * axis_dir;
            self.pt1 = pierce_pt;

            let axey = normp.cross(axis_dir);
            let sint = axey.length() / (n_len * d_len);

            if sint < tol / radius.max(1e-15) {
                // OCCT L695-703: Circle
                self.typeres = AnaResultType::Circle;
                let up = if axis_dir.x.abs() > 0.1 || axis_dir.y.abs() > 0.1 {
                    DVec3::Z
                } else {
                    DVec3::X
                };
                let x_dir = axis_dir.cross(up).cross(axis_dir).normalize_or_zero();
                self.dir1 = axis_dir.normalize_or_zero(); // circle X = cylinder axis
                self.dir2 = x_dir; // circle Y
                self.dir3 = normp.normalize_or_zero(); // circle normal = plane normal
                self.param1 = radius;
            } else {
                // OCCT L706-718: Ellipse
                self.typeres = AnaResultType::Ellipse;
                let cost = (axis_dir.dot(normp) / (d_len * n_len)).abs();
                let axex = axey.cross(normp);
                self.dir1 = normp.normalize_or_zero();
                self.dir2 = axex.normalize_or_zero();
                self.dir3 = axis_dir.normalize_or_zero();
                self.param1 = radius / cost.max(1e-15);
                self.param1bis = radius;
            }
        }

        self.done = true;
    }

    /// OCCT IntAna_QuadQuadGeo.cxx L980-1021: Perform(P, S) — Plane/Sphere
    pub fn perform_plane_sphere(&mut self, p: &Quadric, s: &Quadric) {
        self.init_tolerances();
        self.done = false;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;

        // OCCT L991-993: P.Coefficients(A,B,C,D); S.Location(X,Y,Z); radius = S.Radius()
        let (a, b, c_n, d) = p.plane_coeffs();
        let (x, y, z) = (s.axis_loc().x, s.axis_loc().y, s.axis_loc().z);
        let radius = s.radius();

        // OCCT L995: dist = A*X + B*Y + C*Z + D
        let dist = a * x + b * y + c_n * z + d;

        // OCCT L997-1004: if (Abs(Abs(dist)-radius) < Epsilon(radius)) — tangent
        if (dist.abs() - radius).abs() < radius.abs().max(1.0) * 1e-15 {
            self.nbint = 1;
            self.typeres = AnaResultType::Point;
            self.pt1 = DVec3::new(x - dist * a, y - dist * b, z - dist * c_n);
        // OCCT L1005-1018: else if (Abs(dist) < radius) — Circle
        } else if dist.abs() < radius {
            self.nbint = 1;
            self.typeres = AnaResultType::Circle;
            self.pt1 = DVec3::new(x - dist * a, y - dist * b, z - dist * c_n);
            self.dir1 = p.axis_dir();
            if !p.ax3_direct() {
                self.dir1 = -self.dir1;
            }
            self.dir2 = p.x_dir();
            self.param1 = (radius * radius - dist * dist).sqrt();
        }
        self.param2bis = 0.0;
        self.done = true;
    }

    // =====================================================================
    // Remaining 12 Perform methods
    // =====================================================================

    /// OCCT IntAna_QuadQuadGeo.cxx L752-922: Plane/Cone
    pub fn perform_plane_cone(&mut self, p: &Quadric, co: &Quadric, tol_ang: f64, tol: f64) {
        self.init_tolerances();
        self.done = false;
        self.nbint = 0;
        self.typeres = AnaResultType::Empty;
        let (a, b, c_n, d) = p.plane_coeffs();
        let normp = DVec3::new(a, b, c_n).normalize_or_zero();
        let apex = co.axis_loc();
        let axis_dir = co.axis_dir();
        let dist = a * apex.x + b * apex.y + c_n * apex.z + d;
        let semi_angle = co.semi_angle();
        let cosa = semi_angle.cos();
        let sina = semi_angle.sin().abs();
        let axey = normp.cross(axis_dir);
        let sint = axey.length();
        let cost = axis_dir.dot(normp).abs();
        let costa = cost * cosa - sint * sina;
        if dist.abs() < tol {
            if costa.abs() < tol_ang {
                self.typeres = AnaResultType::Line;
                self.nbint = 1;
                let p2 = apex + 10.0 * axis_dir;
                let d2 = a * p2.x + b * p2.y + c_n * p2.z + d;
                self.pt1 = apex;
                self.dir1 = p2 - d2 * normp - apex;
            } else if cost < sina {
                self.typeres = AnaResultType::Line;
                self.nbint = 2;
                let dh = (sina * sina - cost * cost).sqrt() / cosa;
                let xd = axey.cross(normp);
                self.pt1 = apex;
                self.pt2 = apex;
                self.dir1 = xd + dh * axey;
                self.dir2 = xd - dh * axey;
            } else {
                self.typeres = AnaResultType::Point;
                self.nbint = 1;
                self.pt1 = apex;
            }
        } else {
            let xv = axey.cross(normp);
            if cost < tol_ang {
                self.typeres = AnaResultType::Hyperbola;
                self.nbint = 2;
                self.pt1 = apex - dist * normp;
                self.pt2 = self.pt1;
                self.dir1 = normp;
                self.dir2 = xv.normalize_or_zero();
                self.param1 = (dist / semi_angle.tan()).abs();
                self.param2 = self.param1;
                self.param1bis = dist.abs();
                self.param2bis = dist.abs();
            } else {
                let denom = a * axis_dir.x + b * axis_dir.y + c_n * axis_dir.z;
                let centre = apex - dist * axis_dir / denom.max(1e-16);
                let distance = apex.distance(centre);
                if costa.abs() < tol_ang {
                    self.typeres = AnaResultType::Parabola;
                    self.nbint = 1;
                    let dc = distance / 2.0 / cosa;
                    let ax = xv.normalize_or_zero();
                    self.pt1 = centre - dc * ax;
                    self.dir1 = normp;
                    self.dir2 = ax;
                    self.param1 = dc * sina * sina;
                } else if sint < tol_ang {
                    self.typeres = AnaResultType::Circle;
                    self.nbint = 1;
                    self.pt1 = centre;
                    self.dir1 = axis_dir;
                    self.param1 = distance * semi_angle.tan().abs();
                } else if cost < sina {
                    self.typeres = AnaResultType::Hyperbola;
                    self.nbint = 2;
                    let ax = xv.normalize_or_zero();
                    let den = sina * sina - cost * cost;
                    let dc = sint * sina * sina * distance / den;
                    self.pt1 = centre - dc * ax;
                    self.pt2 = self.pt1;
                    self.dir1 = normp;
                    self.dir2 = ax;
                    self.param1 = cost * sina * cosa * distance / den;
                    self.param2 = self.param1;
                    self.param1bis = cost * sina * distance / den.sqrt();
                    self.param2bis = self.param1bis;
                } else {
                    self.typeres = AnaResultType::Ellipse;
                    self.nbint = 1;
                    let den = cost * cost - sina * sina;
                    let rad = cost * sina * cosa * distance / den;
                    let dc = sint * sina * sina * distance / den;
                    let ax = xv.normalize_or_zero();
                    self.pt1 = centre + dc * ax;
                    self.dir1 = normp;
                    self.dir2 = ax;
                    self.param1 = rad;
                    self.param1bis = cost * sina * distance / den.sqrt();
                }
            }
        }
        if self.typeres == AnaResultType::Ellipse && self.param1.abs() > 1.0E9 {
            self.done = false;
            return;
        }
        if self.typeres == AnaResultType::Hyperbola && self.param1.abs() > 2.0E6 {
            self.done = false;
            return;
        }
        self.done = true;
    }

    /// OCCT L1050: Cylinder/Cylinder
    pub fn perform_cylinder_cylinder(&mut self, c1: &Quadric, c2: &Quadric, tol: f64) {
        self.init_tolerances();
        self.done = true;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;
        let r1 = c1.radius();
        let r2 = c2.radius();
        let rm_r = if r1 > r2 { r1 - r2 } else { r2 - r1 };
        let rmax = if r1 > r2 { r1 } else { r2 };
        let rm_r_rel = rm_r / rmax;
        let a1 = c1.axis_dir();
        let p1 = c1.axis_loc();
        let a2 = c2.axis_dir();
        let p2t = c2.axis_loc();
        let cross = a1.cross(a2);
        let is_parallel = cross.length() <= self.my_epsilon_axes_para;
        let dist_a1a2 = if is_parallel {
            (p2t - p1).cross(a1.normalize_or_zero()).length()
        } else {
            (p2t - p1).dot(cross.normalize_or_zero()).abs()
        };
        if is_parallel {
            if dist_a1a2 <= tol {
                if rm_r <= tol {
                    self.typeres = AnaResultType::Same;
                } else {
                    self.typeres = AnaResultType::Empty;
                }
                return;
            }
            let dir_cyl = a1.normalize_or_zero();
            let p2 = p2t - (p2t - p1).dot(dir_cyl) * dir_cyl;
            let r1pr2 = r1 + r2;
            let dist = (p2 - p1).length();
            if dist > r1pr2 + tol {
                self.typeres = AnaResultType::Empty;
                return;
            }
            if (r1pr2 - dist) <= 1e-15 {
                self.typeres = AnaResultType::Line;
                self.nbint = 1;
                self.dir1 = dir_cyl;
                self.pt1 = p1 + (r1 / r1pr2) * (p2 - p1);
                return;
            }
            if dist > rm_r {
                self.typeres = AnaResultType::Line;
                self.nbint = 2;
                self.dir1 = dir_cyl;
                self.dir2 = dir_cyl;
                let a_cos = 0.5 * (r1 * r1 - r2 * r2 + dist * dist) / (r1 * dist);
                let a_sin2 = 1.0 - a_cos * a_cos;
                let is_tangent = (4.0 * r1 * r1 * a_sin2) < tol * tol;
                let dir_a1a2 = (p2 - p1) / dist;
                if is_tangent {
                    self.nbint = 1;
                    self.pt1 = p1 + dir_a1a2 * r1 * a_cos;
                } else {
                    let a_sin = a_sin2.sqrt();
                    let x_dir = c1.x_dir();
                    let y_dir = c1.y_dir();
                    let a_dx = dir_a1a2.dot(x_dir);
                    let a_dy = dir_a1a2.dot(y_dir);
                    self.pt1 = p1
                        + (a_dx * a_cos - a_dy * a_sin) * r1 * x_dir
                        + (a_dy * a_cos + a_dx * a_sin) * r1 * y_dir;
                    self.pt2 = p1
                        + (a_dx * a_cos + a_dy * a_sin) * r1 * x_dir
                        + (a_dy * a_cos - a_dx * a_sin) * r1 * y_dir;
                }
                return;
            }
            if dist > rm_r - tol {
                self.typeres = AnaResultType::Line;
                self.nbint = 1;
                self.dir1 = dir_cyl;
                let r1_rmr = if r1 < r2 { -(r1 / rm_r) } else { r1 / rm_r };
                self.pt1 = p1 + r1_rmr * (p2 - p1);
                return;
            }
            self.typeres = AnaResultType::Empty;
        } else {
            if rm_r_rel <= self.my_epsilon_cylinder_delta_radius {
                let dir1 = a1.normalize_or_zero();
                let dir2 = a2.normalize_or_zero();
                let cross_n = cross.normalize_or_zero();
                let t1 = (p2t - p1).cross(a2).dot(cross_n) / cross.length_squared().max(1e-16);
                let t2 = (p2t - p1).cross(a1).dot(cross_n) / cross.length_squared().max(1e-16);
                let pt_int = ((p1 + t1 * a1) + (p2t + t2 * a2)) * 0.5;
                let d_pt = (p1 + t1 * a1 - p2t - t2 * a2).length();
                if d_pt <= self.my_epsilon_distance {
                    self.typeres = AnaResultType::Ellipse;
                    self.nbint = 2;
                    self.pt1 = pt_int;
                    self.pt2 = pt_int;
                    let angle = dir1.dot(dir2).abs().acos();
                    let a_val = (0.5_f64 * (std::f64::consts::PI - angle)).sin().abs();
                    let b_val = (0.5_f64 * angle).sin().abs();
                    if a_val == 0.0 || b_val == 0.0 {
                        self.typeres = AnaResultType::Same;
                        return;
                    }
                    self.dir1 = (dir1 + dir2).normalize_or_zero();
                    self.dir2 = (dir1 - dir2).normalize_or_zero();
                    self.param2 = r1 / a_val;
                    self.param1 = r1 / b_val;
                    self.param2bis = r1;
                    self.param1bis = r1;
                    if self.param1 < self.param1bis {
                        std::mem::swap(&mut self.param1, &mut self.param1bis);
                    }
                    if self.param2 < self.param2bis {
                        std::mem::swap(&mut self.param2, &mut self.param2bis);
                    }
                    return;
                }
            }
            if (dist_a1a2 - r1 - r2).abs() < tol {
                self.typeres = AnaResultType::Point;
                self.nbint = 1;
                let n = cross.normalize_or_zero();
                let t1 = (p2t - p1).cross(a2).dot(n) / cross.length_squared().max(1e-16);
                self.pt1 = p1 + t1 * a1 + r1 * n;
            } else {
                self.typeres = AnaResultType::NoGeometricSolution;
            }
        }
    }

    /// OCCT IntAna_QuadQuadGeo.cxx L1324-1344: Cylinder/Cone (coaxial only)
    pub fn perform_cylinder_cone(&mut self, cyl: &Quadric, con: &Quadric, _tol: f64) {
        self.init_tolerances();
        self.done = true;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;
        let cyl_axis = cyl.axis_dir();
        let con_axis = con.axis_dir();
        let cross = cyl_axis.cross(con_axis);
        if cross.length() > self.my_epsilon_axes_para {
            self.typeres = AnaResultType::NoGeometricSolution;
            return;
        }
        let perp = (con.axis_loc() - cyl.axis_loc()).cross(cyl_axis).length();
        if perp > self.my_epsilon_distance {
            self.typeres = AnaResultType::NoGeometricSolution;
            return;
        }
        let dist = cyl.radius() / con.semi_angle().tan();
        let dir = cyl_axis.normalize_or_zero();
        let apex = con.axis_loc();
        self.pt1 = apex + dist * dir;
        self.pt2 = apex - dist * dir;
        self.dir1 = dir;
        self.dir2 = dir;
        self.param1 = cyl.radius();
        self.param2 = cyl.radius();
        self.nbint = 2;
        self.typeres = AnaResultType::Circle;
    }

    /// OCCT IntAna_QuadQuadGeo.cxx L1373-1405: Cylinder/Sphere
    pub fn perform_cylinder_sphere(&mut self, cyl: &Quadric, sph: &Quadric, _tol: f64) {
        self.init_tolerances();
        self.done = true;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;
        // OCCT L1377: AxeOperator + coaxial check
        let cyl_axis = cyl.axis_dir();
        let cyl_origin = cyl.axis_loc();
        let sph_center = sph.axis_loc();
        let sph_radius = sph.radius();
        let cyl_radius = cyl.radius();
        // Check if sphere center lies on cylinder axis
        let to_center = sph_center - cyl_origin;
        let proj = to_center.dot(cyl_axis);
        let perp = (to_center - proj * cyl_axis).length();
        if perp <= _tol.max(1e-15) {
            // OCCT L1378: coaxial — sphere center on cylinder axis
            if sph_radius < cyl_radius {
                // OCCT L1380-1383: sphere inside cylinder, no intersection
                self.typeres = AnaResultType::Empty;
            } else {
                // OCCT L1386-1398: 1 or 2 circles (parallel)
                let dist = (sph_radius * sph_radius - cyl_radius * cyl_radius).sqrt();
                let dir = cyl_axis.normalize_or_zero();
                self.dir1 = dir;
                self.dir2 = dir;
                self.typeres = AnaResultType::Circle;
                self.pt1 = sph_center + dist * dir;
                self.nbint = 1;
                self.param1 = cyl_radius;
                if dist > 1e-15 {
                    self.pt2 = sph_center - dist * dir;
                    self.param2 = cyl_radius;
                    self.nbint = 2;
                }
            }
        } else {
            // OCCT L1401-1404: not coaxial
            self.typeres = AnaResultType::NoGeometricSolution;
        }
    }

    /// OCCT IntAna_QuadQuadGeo.cxx L1433-1521+: Cone/Cone
    pub fn perform_cone_cone(&mut self, c1: &Quadric, c2: &Quadric, tol: f64) {
        self.init_tolerances();
        self.done = true;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;
        let tg1 = c1.semi_angle().tan();
        let tg2 = c2.semi_angle().tan();
        let a1 = c1.axis_dir();
        let a2 = c2.axis_dir();
        let apex1 = c1.axis_loc();
        let apex2 = c2.axis_loc();
        let a_da1a2 = apex1.distance_squared(apex2);
        let cross = a1.cross(a2);
        // OCCT AxeOperator::Same() = Parallel() && (Distance() < 1e-14) — the
        // axes must be BOTH parallel AND coincident.  Offset-parallel axes
        // (distance > 0) fall through to the parallel-plane or NoGeom branch.
        let parallel = cross.length() <= self.my_epsilon_axes_para;
        let dist_between = if parallel {
            (apex2 - apex1).cross(a1).length() / a1.length().max(1e-300)
        } else {
            f64::INFINITY
        };
        let is_same_axis = parallel && dist_between <= 1e-14;
        if is_same_axis {
            // OCCT L1478-1521: same axis
            let d = (apex2 - apex1).dot(a1);
            if (tg1 - tg2).abs() > self.my_epsilon_angle_cone {
                if d.abs() < 1e-10 {
                    self.typeres = AnaResultType::Point;
                    self.nbint = 1;
                    self.pt1 = apex1;
                    return;
                }
                let x1 = (d * tg2) / (tg1 + tg2);
                self.pt1 = apex1 + x1 * a1;
                self.param1 = (x1 * tg1).abs();
                let x2 = (d * tg2) / (tg2 - tg1);
                self.pt2 = apex1 + x2 * a1;
                self.param2 = (x2 * tg1).abs();
                self.dir1 = a1;
                self.dir2 = a1;
                self.nbint = 2;
                self.typeres = AnaResultType::Circle;
            } else {
                if d.abs() < 1e-10 {
                    self.typeres = AnaResultType::Same;
                    return;
                }
                let x = d * 0.5;
                self.pt1 = apex1 + x * a1;
                self.param1 = (x * tg1).abs();
                self.dir1 = a1;
                self.nbint = 1;
                self.typeres = AnaResultType::Circle;
            }
        } else if (tg1 - tg2).abs() < self.my_epsilon_angle_cone && parallel {
            // OCCT L1524-1605 (case 2): parallel axes with (nearly) equal
            // semi-angles.  The intersection of the two cones lies in a plane;
            // intersect that plane with cone 1.
            let dist_a1a2 = dist_between;
            let da1 = a1;
            let geom_apex1 = c1.axis_loc() - da1 * (c1.radius() / tg1.abs().max(1e-300));
            let geom_apex2 = c2.axis_loc() - da1 * (c2.radius() / tg2.abs().max(1e-300));
            let o1o2 = geom_apex2 - geom_apex1;
            let o1o2n = o1o2.normalize_or_zero();
            let o1o2_da1 = da1.dot(o1o2n);
            let o1_proj_a2 = o1o2n - o1o2_da1 * da1;
            let db1 = o1_proj_a2.normalize_or_zero();
            let y_o1o2 = o1o2.dot(da1);
            let abs_tg1 = tg1.abs();
            let x2 = (dist_a1a2 / abs_tg1.max(1e-300) - y_o1o2) * 0.5;
            let x1 = x2 + y_o1o2;
            let p1 = geom_apex1 + x1 * (da1 + abs_tg1 * db1);
            let m_o1o2 = (geom_apex1 + geom_apex2) * 0.5;
            let p1_m_o1o2 = m_o1o2 - p1;
            let da1_x_db1 = da1.cross(db1);
            let ortho_pln = da1_x_db1.cross(p1_m_o1o2.normalize_or_zero()).normalize_or_zero();
            let pln = rcad_kernel::geom::Plane {
                origin: p1,
                normal: ortho_pln,
                u_dir: any_perpendicular_axis(ortho_pln),
                v_dir: ortho_pln.cross(any_perpendicular_axis(ortho_pln)).normalize_or_zero(),
            };
            let mut inter_quad_pln = QuadQuadGeo::new();
            inter_quad_pln.perform_plane_cone(&Quadric::from_plane(&pln), c1, self.my_epsilon_angle_cone, tol);
            if inter_quad_pln.is_done() {
                match inter_quad_pln.type_inter() {
                    AnaResultType::Ellipse => {
                        self.typeres = AnaResultType::Ellipse;
                        self.pt1 = inter_quad_pln.pt1;
                        self.dir1 = inter_quad_pln.dir1;
                        self.dir2 = inter_quad_pln.dir2;
                        self.param1 = inter_quad_pln.param1;
                        self.param1bis = inter_quad_pln.param1bis;
                        self.nbint = 1;
                    }
                    AnaResultType::Circle => {
                        self.typeres = AnaResultType::Circle;
                        self.pt1 = inter_quad_pln.pt1;
                        self.dir1 = inter_quad_pln.dir1;
                        self.dir2 = inter_quad_pln.dir2;
                        self.param1 = inter_quad_pln.param1;
                        self.nbint = 1;
                    }
                    AnaResultType::Hyperbola => {
                        self.typeres = AnaResultType::Hyperbola;
                        self.pt1 = inter_quad_pln.pt1;
                        self.pt2 = inter_quad_pln.pt2;
                        self.dir1 = inter_quad_pln.dir1;
                        self.dir2 = inter_quad_pln.dir2;
                        self.param1 = inter_quad_pln.param1;
                        self.param2 = inter_quad_pln.param2;
                        self.param1bis = inter_quad_pln.param1bis;
                        self.param2bis = inter_quad_pln.param2bis;
                        self.nbint = 2;
                    }
                    AnaResultType::Line => {
                        self.typeres = AnaResultType::Line;
                        self.pt1 = inter_quad_pln.pt1;
                        self.pt2 = inter_quad_pln.pt2;
                        self.dir1 = inter_quad_pln.dir1;
                        self.dir2 = inter_quad_pln.dir2;
                        self.nbint = 2;
                    }
                    _ => {
                        self.typeres = AnaResultType::NoGeometricSolution;
                    }
                }
            } else {
                self.typeres = AnaResultType::NoGeometricSolution;
            }
        } else {
            self.typeres = AnaResultType::NoGeometricSolution;
        }
    }

    /// OCCT IntAna_QuadQuadGeo.cxx L1917-2033: Sphere/Cone (same axis only)
    pub fn perform_sphere_cone(&mut self, sph: &Quadric, con: &Quadric, tol: f64) {
        self.init_tolerances();
        self.done = true;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;
        let sph_center = sph.axis_loc();
        let sph_radius = sph.radius();
        let con_apex = con.axis_loc();
        let con_axis = con.axis_dir();
        let cross = con_axis.cross(sph_center - con_apex);
        if cross.length() > self.my_epsilon_distance {
            self.typeres = AnaResultType::NoGeometricSolution;
            return;
        }
        let d_apex_sph = sph_center.distance(con_apex);
        let con_dir = if d_apex_sph > 1e-15 {
            (sph_center - con_apex).normalize_or_zero()
        } else {
            con_axis.normalize_or_zero()
        };
        let tga = con.semi_angle().tan();
        // OCCT L1947-2033: math_DirectPolynomialRoots → quadratic
        let tgatga = tga * tga;
        let a_coeff = 1.0 + tgatga;
        let b_coeff = 2.0 * tgatga * d_apex_sph;
        let c_coeff = -sph_radius * sph_radius + d_apex_sph * d_apex_sph * tgatga;
        let disc = b_coeff * b_coeff - 4.0 * a_coeff * c_coeff;
        if disc < 0.0 {
            return;
        }
        let sqrt_disc = disc.sqrt();
        let nbsol = if disc.abs() < 1e-15 { 1 } else { 2 };
        self.typeres = AnaResultType::Circle;
        if nbsol >= 1 {
            let x1 = (-b_coeff + sqrt_disc) / (2.0 * a_coeff);
            let d1 = d_apex_sph + x1;
            self.pt1 = con_apex + d1 * con_dir;
            self.param1 = (tga * d1).abs();
            self.dir1 = con_axis.normalize_or_zero();
            self.nbint = 1;
        }
        if nbsol >= 2 {
            let x2 = (-b_coeff - sqrt_disc) / (2.0 * a_coeff);
            let d2 = d_apex_sph + x2;
            self.pt2 = con_apex + d2 * con_dir;
            self.param2 = (tga * d2).abs();
            self.dir2 = self.dir1;
            self.nbint = 2;
        }
    }

    /// OCCT L2034: Sphere/Sphere
    pub fn perform_sphere_sphere(&mut self, s1: &Quadric, s2: &Quadric, tol: f64) {
        self.init_tolerances();
        self.done = true;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;
        let o1 = s1.axis_loc();
        let o2 = s2.axis_loc();
        let d = o1.distance(o2);
        let r1 = s1.radius();
        let r2 = s2.radius();
        let (rmin, rmax) = if r1 > r2 { (r2, r1) } else { (r1, r2) };
        if d <= tol && (r1 - r2).abs() <= tol {
            self.typeres = AnaResultType::Same;
            return;
        }
        if d <= tol {
            return;
        }
        let dir = (o2 - o1).normalize_or_zero();
        let t = rmax - d - rmin;
        if t >= 0.0 && t <= tol {
            self.typeres = AnaResultType::Point;
            self.nbint = 1;
            let t2 = if r1 == rmax {
                (r1 + r2 + d) * 0.5
            } else {
                (-r1 + d - r2) * 0.5
            };
            self.pt1 = o1 + t2 * dir;
        } else if d > r1 + r2 + tol || rmax > d + rmin + tol {
            return;
        } else {
            let alpha = 0.5 * (r1 * r1 - r2 * r2 + d * d) / d;
            let beta_sq = r1 * r1 - alpha * alpha;
            let beta = if beta_sq > 0.0 { beta_sq.sqrt() } else { 0.0 };
            if beta <= TOLERANCE_CLAMP_MIN {
                self.typeres = AnaResultType::Point;
                self.nbint = 1;
                self.pt1 = o1 + ((r1 + d - r2) * 0.5) * dir;
            } else {
                self.typeres = AnaResultType::Circle;
                self.nbint = 1;
                self.dir1 = dir;
                self.param1 = beta;
                self.pt1 = o1 + alpha * dir;
            }
        }
    }

    /// OCCT L2163: Plane/Torus
    /// OCCT IntAna_QuadQuadGeo.cxx L2163-2249: Plane/Torus
    pub fn perform_plane_torus(&mut self, p: &Quadric, tor: &Quadric, tol: f64) {
        self.init_tolerances();
        self.done = true;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;
        // OCCT L2169-2170: aRMin, aRMaj
        let a_rmaj = tor.major_radius();
        let a_rmin = tor.minor_radius();
        // OCCT L2171-2175: if aRMin >= aRMaj -> NoGeometricSolution
        if a_rmin >= a_rmaj {
            self.typeres = AnaResultType::NoGeometricSolution;
            return;
        }

        let p_axis_dir = p.axis_dir();
        let t_axis_dir = tor.axis_dir();
        // OCCT L2182-2183: bParallel, bNormal
        let is_parallel = t_axis_dir.cross(p_axis_dir).length() <= self.my_epsilon_axes_para;
        let is_normal = if !is_parallel {
            t_axis_dir.dot(p_axis_dir).abs() <= self.my_epsilon_axes_para
        } else {
            false
        };
        // OCCT L2184-2188: if neither ∥ nor ⟂ -> NoGeometricSolution
        if !is_normal && !is_parallel {
            self.typeres = AnaResultType::NoGeometricSolution;
            return;
        }

        let t_loc = tor.axis_loc();
        let (a, b, c_n, d) = p.plane_coeffs();

        if is_parallel {
            // OCCT L2193-2228: parallel case
            let a_dist = a * t_loc.x + b * t_loc.y + c_n * t_loc.z + d;
            let a_dr = a_dist.abs() - a_rmin;
            if a_dr > self.my_epsilon_cylinder_delta_radius {
                self.typeres = AnaResultType::Empty;
                return;
            }
            let adj_dist = if a_dr.abs() < self.my_epsilon_cylinder_delta_radius {
                if a_dist < 0.0 { -a_rmin } else { a_rmin }
            } else {
                a_dist
            };
            self.typeres = AnaResultType::Circle;
            self.pt1 = DVec3::new(
                t_loc.x - adj_dist * a,
                t_loc.y - adj_dist * b,
                t_loc.z - adj_dist * c_n,
            );
            let a_dt = (a_rmin * a_rmin - a_dist * a_dist).abs().sqrt();
            self.param1 = a_rmaj + a_dt;
            self.dir1 = t_axis_dir.normalize_or_zero();
            self.nbint = 1;
            if a_dr < -self.my_epsilon_cylinder_delta_radius && a_dt > tol {
                self.pt2 = self.pt1;
                self.param2 = a_rmaj - a_dt;
                self.dir2 = self.dir1;
                self.nbint = 2;
            }
        } else {
            // OCCT L2231-2248: normal case
            let a_dist = (a * t_loc.x + b * t_loc.y + c_n * t_loc.z + d).abs()
                / (a * a + b * b + c_n * c_n).sqrt();
            if a_dist > self.my_epsilon_distance {
                self.typeres = AnaResultType::NoGeometricSolution;
                return;
            }
            self.typeres = AnaResultType::Circle;
            self.param1 = a_rmin;
            self.param2 = a_rmin;
            self.dir1 = p_axis_dir.normalize_or_zero();
            self.dir2 = self.dir1;
            self.nbint = 2;
            let a_dir = t_axis_dir.cross(self.dir1).normalize_or_zero();
            self.pt1 = t_loc + a_rmaj * a_dir;
            self.pt2 = t_loc - a_rmaj * a_dir;
        }
    }

    /// OCCT IntAna_QuadQuadGeo.cxx L2278-2330: Cylinder/Torus (coaxial only)
    pub fn perform_cylinder_torus(&mut self, cyl: &Quadric, tor: &Quadric, tol: f64) {
        self.init_tolerances();
        self.done = true;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;
        // OCCT L2284-2290: minor >= major -> NoGeometricSolution
        let a_rmin = tor.minor_radius();
        let a_rmaj = tor.major_radius();
        if a_rmin >= a_rmaj {
            self.typeres = AnaResultType::NoGeometricSolution;
            return;
        }

        let cyl_axis = cyl.axis_dir();
        let tor_axis = tor.axis_dir();
        // OCCT L2298-2303: check coaxial + distance
        let para = tor_axis.cross(cyl_axis).length() <= self.my_epsilon_axes_para;
        let cyl_loc = cyl.axis_loc();
        let tor_loc = tor.axis_loc();
        let perp = (tor_loc - cyl_loc).cross(tor_axis).length();
        if !para || perp > self.my_epsilon_distance {
            self.typeres = AnaResultType::NoGeometricSolution;
            return;
        }
        // OCCT L2305-2312: cylinder radius vs torus radii
        let r_cyl = cyl.radius();
        if (r_cyl + tol) < (a_rmaj - a_rmin) || (r_cyl - tol) > (a_rmaj + a_rmin) {
            self.typeres = AnaResultType::Empty;
            return;
        }
        // OCCT L2314-2329: 1 or 2 circles
        self.typeres = AnaResultType::Circle;
        let a_dist = (a_rmin * a_rmin - (r_cyl - a_rmaj) * (r_cyl - a_rmaj))
            .abs()
            .sqrt();
        let tor_loc_xyz = tor_loc;
        self.dir1 = tor_axis.normalize_or_zero();
        self.pt1 = tor_loc_xyz + a_dist * self.dir1;
        self.param1 = r_cyl;
        self.nbint = 1;
        if a_dist > tol && r_cyl > (a_rmaj - a_rmin) && r_cyl < (a_rmaj + a_rmin) {
            self.dir2 = self.dir1;
            self.pt2 = tor_loc_xyz - a_dist * self.dir2;
            self.param2 = self.param1;
            self.nbint = 2;
        }
    }

    /// OCCT IntAna_QuadQuadGeo.cxx L2357-2469: Cone/Torus (coaxial only)
    pub fn perform_cone_torus(&mut self, con: &Quadric, tor: &Quadric, tol: f64) {
        self.init_tolerances();
        self.done = true;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;
        let a_rmin = tor.minor_radius();
        let a_rmaj = tor.major_radius();
        if a_rmin >= a_rmaj {
            self.typeres = AnaResultType::NoGeometricSolution;
            return;
        }
        let con_axis = con.axis_dir();
        let tor_axis = tor.axis_dir();
        let is_parallel = tor_axis.cross(con_axis).length() <= self.my_epsilon_axes_para;
        let con_apex = con.axis_loc();
        let tor_loc = tor.axis_loc();
        let perp = (tor_loc - con_apex).cross(tor_axis).length();
        if !is_parallel || perp > self.my_epsilon_distance {
            self.typeres = AnaResultType::NoGeometricSolution;
            return;
        }
        // OCCT L2389-2468: rotate cone generatrix around torus, find circles
        // GAP: rcad does not implement the full cone-torus coaxial solution.
        // Falls back to Empty.
    }

    /// OCCT IntAna_QuadQuadGeo.cxx L2496-2561: Sphere/Torus (sphere center on torus axis only)
    pub fn perform_sphere_torus(&mut self, sph: &Quadric, tor: &Quadric, tol: f64) {
        self.init_tolerances();
        self.done = true;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;
        // OCCT L2502-2508: minor >= major -> NoGeometricSolution
        let a_rmin = tor.minor_radius();
        let a_rmaj = tor.major_radius();
        if a_rmin >= a_rmaj {
            self.typeres = AnaResultType::NoGeometricSolution;
            return;
        }
        // OCCT L2514-2518: sphere center on torus axis?
        let tor_axis = tor.axis_dir();
        let tor_loc = tor.axis_loc();
        let sph_loc = sph.axis_loc();
        let perp = (sph_loc - tor_loc).cross(tor_axis).length();
        if perp > self.my_epsilon_distance {
            self.typeres = AnaResultType::NoGeometricSolution;
            return;
        }
        // OCCT L2523-2533: distance check
        let r_sph = sph.radius();
        let a_x_dir = {
            let up = if tor_axis.x.abs() > 0.1 || tor_axis.y.abs() > 0.1 {
                DVec3::Z
            } else {
                DVec3::X
            };
            tor_axis.cross(up).cross(tor_axis).normalize_or_zero()
        };
        let a_tor_loc = tor_loc + a_rmaj * a_x_dir;
        let a_dist = (sph_loc - a_tor_loc).length();
        if (a_dist - tol) > (a_rmin + r_sph) || (a_dist + tol) < (a_rmin - r_sph).abs() {
            self.typeres = AnaResultType::Empty;
            return;
        }
        // OCCT L2535-2560: circle
        self.typeres = AnaResultType::Circle;
        let an_alpha =
            0.5 * (a_rmin * a_rmin - r_sph * r_sph + a_dist * a_dist) / a_dist.max(1e-15);
        let a_beta = (a_rmin * a_rmin - an_alpha * an_alpha).abs().sqrt();
        let a_dir12 = (sph_loc - a_tor_loc).normalize_or_zero();
        let a_ph = a_tor_loc + an_alpha * a_dir12;
        let a_dc = tor_axis.cross(a_dir12).normalize_or_zero();
        let a_dval = a_beta * a_dc;
        let a_p = a_ph + a_dval;
        let a_lin = tor_axis.normalize_or_zero();
        let param1 = (sph_loc - a_p).dot(a_lin).abs();
        self.pt1 = a_p - param1 * a_x_dir;
        self.dir1 = a_lin;
        self.param1 = a_lin.dot(self.pt1 - tor_loc).abs();
        self.nbint = 1;
        if a_dist < (r_sph + a_rmin) && a_dist > (r_sph - a_rmin).abs() && a_dval.length() > tol {
            let a_p2 = a_ph - a_dval;
            let param2 = (sph_loc - a_p2).dot(a_lin).abs();
            self.pt2 = a_p2 - param2 * a_x_dir;
            self.dir2 = self.dir1;
            self.nbint = 2;
        }
    }

    /// OCCT IntAna_QuadQuadGeo.cxx L2588-2650: Torus/Torus (coaxial only)
    pub fn perform_torus_torus(&mut self, t1: &Quadric, t2: &Quadric, tol: f64) {
        self.init_tolerances();
        self.done = true;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;
        let a1 = t1.axis_dir();
        let loc1 = t1.axis_loc();
        let a2 = t2.axis_dir();
        let loc2 = t2.axis_loc();
        let rmj1 = t1.major_radius();
        let rmn1 = t1.minor_radius();
        let rmj2 = t2.major_radius();
        let rmn2 = t2.minor_radius();
        // OCCT L2605-2610: coaxial check
        let is_parallel = a1.cross(a2).length() <= self.my_epsilon_axes_para;
        let perp = (loc2 - loc1).cross(a1).length();
        if !is_parallel || perp > self.my_epsilon_distance {
            self.typeres = AnaResultType::NoGeometricSolution;
            return;
        }
        // OCCT L2612-2617: same torus check
        if loc1.distance(loc2) <= tol && (rmn1 - rmn2).abs() <= tol && (rmj1 - rmj2).abs() <= tol {
            self.typeres = AnaResultType::Same;
            return;
        }
        // OCCT L2619-2623: invalid geometry
        if rmn1 >= rmj1 || rmn2 >= rmj2 {
            self.typeres = AnaResultType::NoGeometricSolution;
            return;
        }
        // OCCT L2625-2650: distance-based circles
        let a_x_dir = {
            let up = if a1.x.abs() > 0.1 || a1.y.abs() > 0.1 {
                DVec3::Z
            } else {
                DVec3::X
            };
            a1.cross(up).cross(a1).normalize_or_zero()
        };
        let a_p1 = loc1 + rmj1 * a_x_dir;
        let a_p2 = loc2 + rmj2 * a_x_dir;
        let a_dist = (a_p1 - a_p2).length();
        if (a_dist - tol) > (rmn1 + rmn2) || (a_dist + tol) < (rmn1 - rmn2).abs() {
            self.typeres = AnaResultType::Empty;
            return;
        }
        self.typeres = AnaResultType::Circle;
        let an_alpha = 0.5 * (rmn1 * rmn1 - rmn2 * rmn2 + a_dist * a_dist) / a_dist.max(1e-15);
        let a_beta = (rmn1 * rmn1 - an_alpha * an_alpha).abs().sqrt();
        let a_dir12 = (a_p2 - a_p1).normalize_or_zero();
        let a_ph = a_p1 + an_alpha * a_dir12;
        self.pt1 = a_ph;
        self.dir1 = a1.normalize_or_zero();
        self.param1 = a_beta;
        self.nbint = 1;
        if a_dist < (rmn1 + rmn2) && a_dist > (rmn1 - rmn2).abs() && a_beta > tol {
            self.pt2 = a_ph;
            self.param2 = a_beta;
            self.dir2 = self.dir1;
            self.nbint = 2;
        }
    }
}

/// Helper: compute a perpendicular pair to a given direction ().
fn any_perpendicular_pair(dir: DVec3) -> (DVec3, DVec3) {
    let x_dir = rcad_kernel::geom::any_perpendicular(dir);
    let y_dir = dir.cross(x_dir).normalize_or_zero();
    (x_dir, y_dir)
}

// =========================================================================
// Helper functions used by Perform methods
// =========================================================================

/// Compute a point on line at parameter t
fn line_point(origin: DVec3, dir: DVec3, t: f64) -> DVec3 {
    origin + t * dir
}

/// Project point onto line, return parameter t
fn project_on_line(point: DVec3, origin: DVec3, dir: DVec3) -> f64 {
    (point - origin).dot(dir)
}
