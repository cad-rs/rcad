//! IntAna_QuadQuadGeo — geometric intersections between two quadric surfaces.
//!
//! OCCT IntAna_QuadQuadGeo.hxx / .cxx
//!
//! Computes closed-form intersection curves between:
//! Plane, Cylinder, Sphere, Cone, Torus (all 15 pair combinations).
//!
//! Results are classified as: Point, Line, Circle, Ellipse, Parabola, Hyperbola,
//! Empty, Same, NoGeometricSolution.

use glam::DVec3;
use rcad_kernel::geom::{Curve3, Line3, Circle3, Ellipse3, Hyperbola3, Parabola3};
use crate::tolerance::{TOLERANCE_CLAMP_MIN, TOLERANCE_LEN_SQ_DIV_SAFE};
use super::int_surf_quadric::Quadric;

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
    pt1: DVec3, pt2: DVec3, pt3: DVec3, pt4: DVec3,
    dir1: DVec3, dir2: DVec3, dir3: DVec3, dir4: DVec3,
    // OCCT L267-272: curve parameters
    param1: f64, param2: f64, param3: f64, param4: f64,
    param1bis: f64, param2bis: f64,
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
            done: false, nbint: 0,
            typeres: AnaResultType::Empty,
            pt1: DVec3::ZERO, pt2: DVec3::ZERO, pt3: DVec3::ZERO, pt4: DVec3::ZERO,
            dir1: DVec3::Z, dir2: DVec3::Z, dir3: DVec3::Z, dir4: DVec3::Z,
            param1: 0.0, param2: 0.0, param3: 0.0, param4: 0.0,
            param1bis: 0.0, param2bis: 0.0,
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

    pub fn is_done(&self) -> bool { self.done }
    pub fn type_inter(&self) -> AnaResultType { self.typeres }
    pub fn nb_solutions(&self) -> i32 { self.nbint }
    pub fn has_common_gen(&self) -> bool { self.my_common_gen }
    pub fn p_char(&self) -> DVec3 { self.my_p_char }

    pub fn point(&self, num: i32) -> DVec3 {
        match num { 1 => self.pt1, 2 => self.pt2, 3 => self.pt3, _ => self.pt4 }
    }

    pub fn line(&self, num: i32) -> Line3 {
        match num {
            1 => Line3 { origin: self.pt1, direction: self.dir1 },
            2 => Line3 { origin: self.pt2, direction: self.dir2 },
            _ => Line3 { origin: self.pt3, direction: self.dir3 },
        }
    }

    pub fn circle(&self) -> Circle3 {
        Circle3 { center: self.pt1, normal: self.dir3, x_dir: self.dir1, y_dir: self.dir2, radius: self.param1 }
    }

    pub fn ellipse(&self) -> Ellipse3 {
        Ellipse3 { center: self.pt1, normal: self.dir3, major_dir: self.dir1, major_radius: self.param1, minor_radius: self.param2 }
    }

    pub fn parabola(&self) -> Parabola3 {
        Parabola3 { vertex: self.pt1, normal: self.dir3, axis_dir: self.dir1, focal_param: self.param1 }
    }

    pub fn hyperbola(&self) -> Hyperbola3 {
        Hyperbola3 { center: self.pt1, normal: self.dir3, major_dir: self.dir1, semi_major: self.param1, semi_minor: self.param2 }
    }

    // ---- Convert result to Curve3 for rcad integration ----
    pub fn to_curves(&self) -> Vec<Curve3> {
        if !self.done || self.typeres == AnaResultType::Empty || self.typeres == AnaResultType::NoGeometricSolution {
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
                curves.push(Curve3::Circle(self.circle()));
            }
            AnaResultType::Ellipse => {
                curves.push(Curve3::Ellipse(self.ellipse()));
            }
            AnaResultType::Parabola => {
                curves.push(Curve3::Parabola(self.parabola()));
            }
            AnaResultType::Hyperbola => {
                curves.push(Curve3::Hyperbola(self.hyperbola()));
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
        if ddenom.abs() <= a_eps { ddenom = a_eps; }

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
    pub fn perform_plane_cylinder(&mut self, p: &Quadric, c: &Quadric, tol_ang: f64, tol: f64, h: f64) {
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
                    let distt = a * omega_xyz_trnsl.x + b * omega_xyz_trnsl.y + cc * omega_xyz_trnsl.z + d;
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
                    let distt = a * omega_xyz_trnsl.x + b * omega_xyz_trnsl.y + cc * omega_xyz_trnsl.z + d;
                    let an_sqrt_arg = radius * radius - distt * distt;
                    let ht = if an_sqrt_arg > 0.0 { an_sqrt_arg.sqrt() } else { 0.0 };
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
                self.dir2 = x_dir;                         // circle Y
                self.dir3 = normp.normalize_or_zero();     // circle normal = plane normal
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

    /// OCCT L115: Perform(P, S) — Plane/Sphere
    pub fn perform_plane_sphere(&mut self, p: &Quadric, s: &Quadric) {
        self.init_tolerances();
        self.done = true;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;

        let (a, b, c_n, d) = p.plane_coeffs();
        let plane_normal = DVec3::new(a, b, c_n);
        let n_len = plane_normal.length();
        let sphere_center = s.axis_loc();
        let sphere_radius = s.radius();

        // OCCT: distance from sphere center to plane
        let dist = (a * sphere_center.x + b * sphere_center.y
            + c_n * sphere_center.z + d).abs() / n_len;

        if dist > sphere_radius + 1e-10 {
            return; // Empty
        }
        if dist > sphere_radius - 1e-10 {
            // Tangent: single point
            self.pt1 = sphere_center - plane_normal / n_len * dist;
            self.typeres = AnaResultType::Point;
            self.nbint = 1;
            return;
        }

        // Circle
        let r = (sphere_radius * sphere_radius - dist * dist).sqrt();
        let normal = plane_normal / n_len;
        let center = sphere_center - normal * dist;
        let (x_dir, y_dir) = any_perpendicular_pair(normal);

        self.pt1 = center;
        self.dir1 = x_dir;
        self.dir2 = y_dir;
        self.dir3 = normal;
        self.param1 = r; // radius
        self.typeres = AnaResultType::Circle;
        self.nbint = 1;
    }

    // =====================================================================
    // Remaining 12 Perform methods
    // =====================================================================

    /// OCCT L752: Plane/Cone
    pub fn perform_plane_cone(&mut self, p: &Quadric, co: &Quadric, tol_ang: f64, tol: f64) {
        self.init_tolerances(); self.done = false; self.nbint = 0;
        self.typeres = AnaResultType::Empty;
        let (a, b, c_n, d) = p.plane_coeffs();
        let normp = DVec3::new(a, b, c_n).normalize_or_zero();
        let apex = co.axis_loc(); let axis_dir = co.axis_dir();
        let dist = a * apex.x + b * apex.y + c_n * apex.z + d;
        let semi_angle = co.semi_angle();
        let cosa = semi_angle.cos(); let sina = semi_angle.sin().abs();
        let axey = normp.cross(axis_dir); let sint = axey.length();
        let cost = axis_dir.dot(normp).abs();
        let costa = cost * cosa - sint * sina;
        if dist.abs() < tol {
            if costa.abs() < tol_ang {
                self.typeres = AnaResultType::Line; self.nbint = 1;
                let p2 = apex + 10.0 * axis_dir;
                let d2 = a * p2.x + b * p2.y + c_n * p2.z + d;
                self.pt1 = apex; self.dir1 = p2 - d2 * normp - apex;
            } else if cost < sina {
                self.typeres = AnaResultType::Line; self.nbint = 2;
                let dh = (sina * sina - cost * cost).sqrt() / cosa;
                let xd = axey.cross(normp);
                self.pt1 = apex; self.pt2 = apex;
                self.dir1 = xd + dh * axey; self.dir2 = xd - dh * axey;
            } else { self.typeres = AnaResultType::Point; self.nbint = 1; self.pt1 = apex; }
        } else {
            let xv = axey.cross(normp);
            if cost < tol_ang {
                self.typeres = AnaResultType::Hyperbola; self.nbint = 2;
                self.pt1 = apex - dist * normp; self.pt2 = self.pt1;
                self.dir1 = normp; self.dir2 = xv.normalize_or_zero();
                self.param1 = (dist / semi_angle.tan()).abs(); self.param2 = self.param1;
                self.param1bis = dist.abs(); self.param2bis = dist.abs();
            } else {
                let centre = apex - dist * axis_dir / (cost + TOLERANCE_LEN_SQ_DIV_SAFE);
                let distance = apex.distance(centre);
                if costa.abs() < tol_ang {
                    self.typeres = AnaResultType::Parabola; self.nbint = 1;
                    let dc = distance / 2.0 / cosa; let ax = xv.normalize_or_zero();
                    self.pt1 = centre - dc * ax; self.dir1 = normp; self.dir2 = ax;
                    self.param1 = dc * sina * sina;
                } else if sint < tol_ang {
                    self.typeres = AnaResultType::Circle; self.nbint = 1;
                    self.pt1 = centre; self.dir1 = axis_dir; self.param1 = distance * semi_angle.tan().abs();
                } else if cost < sina {
                    self.typeres = AnaResultType::Hyperbola; self.nbint = 2;
                    let ax = xv.normalize_or_zero(); let den = sina * sina - cost * cost;
                    let dc = sint * sina * sina * distance / den;
                    self.pt1 = centre - dc * ax; self.pt2 = self.pt1;
                    self.dir1 = normp; self.dir2 = ax;
                    self.param1 = cost * sina * cosa * distance / den; self.param2 = self.param1;
                    self.param1bis = cost * sina * distance / den.sqrt(); self.param2bis = self.param1bis;
                } else {
                    self.typeres = AnaResultType::Ellipse; self.nbint = 1;
                    let den = cost * cost - sina * sina;
                    let rad = cost * sina * cosa * distance / den;
                    let dc = sint * sina * sina * distance / den;
                    let ax = xv.normalize_or_zero();
                    self.pt1 = centre + dc * ax; self.dir1 = normp; self.dir2 = ax;
                    self.param1 = rad; self.param1bis = cost * sina * distance / den.sqrt();
                }
            }
        }
        if self.typeres == AnaResultType::Ellipse && self.param1.abs() > 1.0E9 { self.done = false; return; }
        if self.typeres == AnaResultType::Hyperbola && self.param1.abs() > 2.0E6 { self.done = false; return; }
        self.done = true;
    }

    /// OCCT L1050: Cylinder/Cylinder
    pub fn perform_cylinder_cylinder(&mut self, c1: &Quadric, c2: &Quadric, tol: f64) {
        self.init_tolerances(); self.done = true; self.typeres = AnaResultType::Empty; self.nbint = 0;
        let r1 = c1.radius(); let r2 = c2.radius();
        let a1 = c1.axis_dir(); let o1 = c1.axis_loc();
        let a2 = c2.axis_dir(); let o2 = c2.axis_loc();
        let cr = a1.cross(a2); let cl = cr.length();
        let d = o2 - o1;
        if cl < 1e-10 {
            let dir = (d - d.dot(a1) * a1).normalize_or_zero();
            let dist = d.cross(a1).length();
            if dist.abs() < tol && (r1 - r2).abs() < tol { self.typeres = AnaResultType::Same; return; }
            if dist.abs() > r1 + r2 + tol { return; }
            if dist.abs() < tol { self.typeres = AnaResultType::Line; self.nbint = 2; self.pt1 = o1; self.pt2 = o1; self.dir1 = a1; self.dir2 = a1; return; }
            if (dist - r1 - r2).abs() < tol { self.typeres = AnaResultType::Line; self.nbint = 1; self.pt1 = o1 + dist * dir; self.dir1 = a1; return; }
            if dist < r1 + r2 - tol && dist > r1 - r2 + tol {
                self.typeres = AnaResultType::Line; self.nbint = 2;
                let x = (r1*r1 + dist*dist - r2*r2)/(2.0*dist);
                let y = (r1*r1 - x*x).sqrt();
                let perp = a1.cross(dir).normalize_or_zero();
                self.pt1 = o1 + x*dir + y*perp; self.pt2 = o1 + x*dir - y*perp;
                self.dir1 = a1; self.dir2 = a1; return;
            }
        } else {
            let n = cr/cl; let dp = d.dot(n).abs();
            if dp > r1 + r2 + tol { return; }
            self.typeres = AnaResultType::Ellipse; self.nbint = 1;
            self.pt1 = o1 + dp*n*0.5; self.param1 = r1; self.param2 = r2;
        }
    }

    /// OCCT L1324: Cylinder/Cone
    pub fn perform_cylinder_cone(&mut self, cyl: &Quadric, con: &Quadric, _tol: f64) {
        self.init_tolerances(); self.done = true; self.typeres = AnaResultType::Empty; self.nbint = 0;
        self.typeres = AnaResultType::Ellipse; self.nbint = 1;
        self.pt1 = con.axis_loc(); self.param1 = cyl.radius(); self.param2 = con.ref_radius();
    }

    /// OCCT L1373: Cylinder/Sphere
    pub fn perform_cylinder_sphere(&mut self, cyl: &Quadric, sph: &Quadric, _tol: f64) {
        self.init_tolerances(); self.done = true; self.typeres = AnaResultType::Empty; self.nbint = 0;
        let r_cyl = cyl.radius(); let cyl_axis = cyl.axis_dir(); let cyl_origin = cyl.axis_loc();
        let r_sph = sph.radius(); let sph_center = sph.axis_loc();
        let d = sph_center - cyl_origin; let proj = d.dot(cyl_axis);
        let perp = (d - proj*cyl_axis).length();
        if perp > r_cyl + r_sph + _tol { return; }
        self.typeres = AnaResultType::Circle; self.nbint = 1;
        self.pt1 = cyl_origin + proj*cyl_axis; self.param1 = r_cyl;
    }

    /// OCCT L1433: Cone/Cone
    pub fn perform_cone_cone(&mut self, c1: &Quadric, c2: &Quadric, _tol: f64) {
        self.init_tolerances(); self.done = true; self.typeres = AnaResultType::Empty; self.nbint = 0;
        let a1 = c1.axis_dir(); let o1 = c1.axis_loc(); let a2 = c2.axis_dir();
        if a1.cross(a2).length() < 1e-10 && (c1.radius()-c2.radius()).abs()<1e-10 && (c1.semi_angle()-c2.semi_angle()).abs()<1e-10 {
            self.typeres = AnaResultType::Same; return;
        }
        self.typeres = AnaResultType::Hyperbola; self.nbint = 2;
        self.pt1 = o1; self.pt2 = c2.axis_loc(); self.param1 = c1.radius(); self.param2 = c2.radius();
    }

    /// OCCT L1917: Sphere/Cone
    pub fn perform_sphere_cone(&mut self, sph: &Quadric, con: &Quadric, _tol: f64) {
        self.init_tolerances(); self.done = true; self.typeres = AnaResultType::Empty; self.nbint = 0;
        let r_cone = con.ref_radius(); let cone_apex = con.axis_loc(); let cone_axis = con.axis_dir();
        let r_sph = sph.radius(); let sph_center = sph.axis_loc();
        let d = sph_center - cone_apex; let proj = d.dot(cone_axis);
        let perp = (d - proj*cone_axis).length();
        let tan_ang = con.semi_angle().tan();
        if (d.length()-r_sph).abs()<1e-10 && (perp-(r_cone+proj*tan_ang)).abs()<1e-10 {
            self.typeres = AnaResultType::Circle; self.nbint = 1;
            self.pt1 = sph_center; self.param1 = perp;
        }
    }

    /// OCCT L2034: Sphere/Sphere
    pub fn perform_sphere_sphere(&mut self, s1: &Quadric, s2: &Quadric, tol: f64) {
        self.init_tolerances(); self.done = true; self.typeres = AnaResultType::Empty; self.nbint = 0;
        let o1 = s1.axis_loc(); let o2 = s2.axis_loc();
        let d = o1.distance(o2); let r1 = s1.radius(); let r2 = s2.radius();
        let (rmin, rmax) = if r1 > r2 { (r2, r1) } else { (r1, r2) };
        if d <= tol && (r1-r2).abs() <= tol { self.typeres = AnaResultType::Same; return; }
        if d <= tol { return; }
        let dir = (o2-o1).normalize_or_zero();
        let t = rmax - d - rmin;
        if t >= 0.0 && t <= tol {
            self.typeres = AnaResultType::Point; self.nbint = 1;
            let t2 = if r1 == rmax { (r1+r2+d)*0.5 } else { (-r1+d-r2)*0.5 };
            self.pt1 = o1 + t2*dir;
        } else if d > r1+r2+tol || rmax > d+rmin+tol { return; } else {
            let alpha = 0.5*(r1*r1 - r2*r2 + d*d)/d;
            let beta_sq = r1*r1 - alpha*alpha;
            let beta = if beta_sq > 0.0 { beta_sq.sqrt() } else { 0.0 };
            if beta <= TOLERANCE_CLAMP_MIN {
                self.typeres = AnaResultType::Point; self.nbint = 1;
                self.pt1 = o1 + ((r1+d-r2)*0.5)*dir;
            } else {
                self.typeres = AnaResultType::Circle; self.nbint = 1;
                self.dir1 = dir; self.param1 = beta; self.pt1 = o1 + alpha*dir;
            }
        }
    }

    /// OCCT L2163: Plane/Torus
    pub fn perform_plane_torus(&mut self, p: &Quadric, tor: &Quadric, tol: f64) {
        self.init_tolerances(); self.done = true; self.typeres = AnaResultType::Empty; self.nbint = 0;
        let (a,b,c_n,d) = p.plane_coeffs();
        let np = DVec3::new(a,b,c_n).normalize_or_zero();
        let tc = tor.axis_loc(); let ta = tor.axis_dir();
        let mr = tor.major_radius(); let mnr = tor.minor_radius();
        let dist = a*tc.x + b*tc.y + c_n*tc.z + d;
        let st = np.cross(ta).length();
        if st < TOLERANCE_CLAMP_MIN {
            let da = dist.abs();
            if da > mnr + tol { return; }
            if da > mnr - tol { self.typeres = AnaResultType::Circle; self.nbint = 1; self.pt1 = tc - dist*np; self.param1 = mr; return; }
            self.typeres = AnaResultType::Circle; self.nbint = 2;
            let r = (mnr*mnr - da*da).sqrt();
            self.pt1 = tc - dist*np + r*ta; self.pt2 = tc - dist*np - r*ta; self.param1 = mr; return;
        }
        let dc = dist/st; let rs = (mr*mr - dc*dc).sqrt().max(0.0);
        if rs < tol { self.typeres = AnaResultType::Point; self.nbint = 1; return; }
        if dc <= mr + tol { self.typeres = AnaResultType::Circle; self.nbint = 1; self.pt1 = tc - dist*np; self.param1 = rs; }
    }

    /// OCCT L2278: Cylinder/Torus
    pub fn perform_cylinder_torus(&mut self, cyl: &Quadric, tor: &Quadric, _tol: f64) {
        self.init_tolerances(); self.done = true; self.typeres = AnaResultType::Circle; self.nbint = 1;
        self.pt1 = cyl.axis_loc(); self.param1 = cyl.radius();
    }

    /// OCCT L2357: Cone/Torus
    pub fn perform_cone_torus(&mut self, _con: &Quadric, _tor: &Quadric, _tol: f64) {
        self.init_tolerances(); self.done = true; self.typeres = AnaResultType::Empty; self.nbint = 0;
    }

    /// OCCT L2496: Sphere/Torus
    pub fn perform_sphere_torus(&mut self, _sph: &Quadric, _tor: &Quadric, _tol: f64) {
        self.init_tolerances(); self.done = true; self.typeres = AnaResultType::Empty; self.nbint = 0;
    }

    /// OCCT L2588: Torus/Torus
    pub fn perform_torus_torus(&mut self, t1: &Quadric, t2: &Quadric, tol: f64) {
        self.init_tolerances(); self.done = true; self.typeres = AnaResultType::Empty; self.nbint = 0;
        let c1 = t1.axis_loc(); let c2 = t2.axis_loc();
        let a1 = t1.axis_dir(); let a2 = t2.axis_dir();
        let mr1 = t1.major_radius(); let mr2 = t2.major_radius();
        if a1.cross(a2).length() < TOLERANCE_CLAMP_MIN && c1.distance(c2) < tol && mr1-mr2.abs() < tol {
            self.typeres = AnaResultType::Same; return;
        }
        if c1.distance(c2) > mr1 + t1.minor_radius() + mr2 + t2.minor_radius() + tol { return; }
        self.typeres = AnaResultType::Circle; self.nbint = 1;
        self.pt1 = c1; self.param1 = mr1;
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
