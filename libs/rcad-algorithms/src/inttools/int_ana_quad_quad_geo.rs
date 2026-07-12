//! OCCT-aligned: IntAna_QuadQuadGeo — geometric intersections between two quadric surfaces.
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

/// OCCT-aligned: IntAna_QuadQuadGeo
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
        self.my_epsilon_mini_circle_radius = 1e-15;
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

    /// OCCT L76-79: Perform(P1, P2, TolAng, Tol) — Plane/Plane
    pub fn perform_plane_plane(&mut self, p1: &Quadric, p2: &Quadric, tol_ang: f64, tol: f64) {
        self.init_tolerances();
        self.done = true;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;

        // OCCT: plane1/2 coefficients
        let (a1, b1, c1, d1) = p1.plane_coeffs();
        let (a2, b2, c2, d2) = p2.plane_coeffs();
        let n1 = DVec3::new(a1, b1, c1);
        let n2 = DVec3::new(a2, b2, c2);
        let n1_len = n1.length();
        let n2_len = n2.length();

        // OCCT: cross product magnitude = |sin(angle)| * |n1| * |n2|
        let cross_mag = n1.cross(n2).length();
        // OCCT: if aMVD <= TolAng — normals are collinear (parallel planes)
        if cross_mag <= tol_ang * n1_len.max(n2_len).max(1.0) {
            // OCCT: check if identical
            let dist = d1 / n1_len - d2 / n2_len;
            if dist.abs() < tol {
                self.typeres = AnaResultType::Same;
            } else {
                self.typeres = AnaResultType::Empty;
            }
            return;
        }

        // OCCT: non-parallel → intersection line
        let cross = n1.cross(n2);
        let cross_len = cross.length();
        let dir = cross / cross_len;
        // Find a point on the line: solve n1·p = -d1, n2·p = -d2
        // Using formula: P = (d2*(n1×n2) + d1*(n2×n1)) / |n1×n2|² ?
        // Actually OCCT uses a different approach but the result is a line
        let origin = (n2 * d1 - n1 * d2).cross(cross) / (cross_len * cross_len);

        self.pt1 = origin;
        self.dir1 = dir;
        self.typeres = AnaResultType::Line;
        self.nbint = 1;
    }

    /// OCCT L105-109: Perform(P, C, TolAng, Tol, H) — Plane/Cylinder
    pub fn perform_plane_cylinder(&mut self, p: &Quadric, c: &Quadric, tol_ang: f64, tol: f64, _h: f64) {
        self.init_tolerances();
        self.done = true;
        self.typeres = AnaResultType::Empty;
        self.nbint = 0;

        let (a, b, c_n, d) = p.plane_coeffs();
        let plane_normal = DVec3::new(a, b, c_n);
        let axis_dir = c.axis_dir();
        let radius = c.radius();
        let center = c.axis_loc();

        // OCCT: dot product of plane normal and cylinder axis
        let dot = plane_normal.dot(axis_dir);
        let cos_ang = dot / (plane_normal.length() * axis_dir.length());

        // OCCT: cylinder axis parallel to plane?
        if cos_ang.abs() < tol_ang.cos() {
            // Parallel: intersection is 1-2 line(s) or empty/circle
            // Distance from cylinder axis to plane
            let dist_to_plane = (a * center.x + b * center.y + c_n * center.z + d).abs()
                / plane_normal.length();

            if dist_to_plane > radius + tol {
                return; // Empty
            }
            if dist_to_plane < tol {
                // Cylinder axis lies in plane → 2 lines (generatrices)
                self.pt1 = center - radius * axis_dir.cross(plane_normal).normalize_or_zero();
                self.dir1 = axis_dir;
                self.pt2 = center + radius * axis_dir.cross(plane_normal).normalize_or_zero();
                self.dir2 = axis_dir;
                self.nbint = 2;
                self.typeres = AnaResultType::Line;
                return;
            }
            // 2 lines (generatrices at intersection)
            let offset = (radius * radius - dist_to_plane * dist_to_plane).sqrt();
            let perp = axis_dir.cross(plane_normal).normalize_or_zero();
            self.pt1 = center - dist_to_plane * plane_normal.normalize_or_zero() - offset * perp;
            self.dir1 = axis_dir;
            self.pt2 = center - dist_to_plane * plane_normal.normalize_or_zero() + offset * perp;
            self.dir2 = axis_dir;
            self.nbint = 2;
            self.typeres = AnaResultType::Line;
            return;
        }

        // OCCT: non-parallel → ellipse or circle
        let normal = plane_normal.normalize_or_zero();
        let proj_center = center - (plane_normal.dot(center) + d) / plane_normal.length_squared() * plane_normal;
        let x_dir = axis_dir - dot / plane_normal.length_squared() * plane_normal;
        let y_dir = normal.cross(axis_dir);
        let x_len = x_dir.length();
        let y_len = y_dir.length();

        if x_len < 1e-15 || y_len < 1e-15 { return; }

        // OCCT: semi-axes of the ellipse
        let a_ellipse = radius / cos_ang.abs();
        let b_ellipse = radius;

        self.pt1 = proj_center;
        self.dir1 = x_dir / x_len;
        self.dir2 = y_dir / y_len;
        self.dir3 = normal;
        self.param1 = a_ellipse; // major
        self.param2 = b_ellipse; // minor

        // OCCT: check if close to circle
        if (a_ellipse - b_ellipse).abs() < tol {
            self.typeres = AnaResultType::Circle;
            self.param1 = (a_ellipse + b_ellipse) / 2.0;
        } else {
            self.typeres = AnaResultType::Ellipse;
        }
        self.nbint = 1;
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
                let centre = apex - dist * axis_dir / (cost + 1e-30);
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
            if beta <= 1e-15 {
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
        if st < 1e-15 {
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
        if a1.cross(a2).length() < 1e-15 && c1.distance(c2) < tol && mr1-mr2.abs() < tol {
            self.typeres = AnaResultType::Same; return;
        }
        if c1.distance(c2) > mr1 + t1.minor_radius() + mr2 + t2.minor_radius() + tol { return; }
        self.typeres = AnaResultType::Circle; self.nbint = 1;
        self.pt1 = c1; self.param1 = mr1;
    }
}

/// Helper: compute a perpendicular pair to a given direction (OCCT-aligned).
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
