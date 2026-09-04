// OCCT Contap_ContAna (TKHLR) — contours of quadric surfaces (sphere /
// cylinder / cone vs direction, draft angle or eye point).
//
// Contap_ContAna.hxx L36-86 + Contap_ContAna.cxx L29-483.
// NOTE the file constant: `static const double Tolpetit = 1.0e-8;` here,
// versus Contap_Contour.cxx L46 `Tolpetit = 1.e-10` — both are kept
// verbatim (the "conflicting Tolpetit" quirk of the port plan §5).

use glam::DVec3;
use rcad_kernel::geom::{Circle3, ConicalSurface, CylindricalSurface, Point3, SphericalSurface, Vec3};
use rcad_kernel::topo::topods as _topods_shape;

use crate::geomalgo::geom2d_int::Curve2dType;

/// OCCT Contap_ContAna.cxx L29 — `static const double Tolpetit = 1.0e-8;`
const TOLPETIT: f64 = 1.0e-8;

/// OCCT Contap_ContAna.
#[derive(Debug, Clone)]
pub struct ContAna {
    done: bool,
    nb_sol: i32,
    typ_l: Curve2dType,
    pt1: Point3,
    pt2: Point3,
    pt3: Point3,
    pt4: Point3,
    dir1: Vec3,
    dir2: Vec3,
    dir3: Vec3,
    dir4: Vec3,
    prm: f64,
}

impl ContAna {
    /// OCCT Contap_ContAna() (cxx L31-37).
    pub fn new() -> Self {
        ContAna {
            done: false,
            nb_sol: 0,
            typ_l: Curve2dType::OtherCurve,
            pt1: Point3::ZERO,
            pt2: Point3::ZERO,
            pt3: Point3::ZERO,
            pt4: Point3::ZERO,
            dir1: Vec3::ZERO,
            dir2: Vec3::ZERO,
            dir3: Vec3::ZERO,
            dir4: Vec3::ZERO,
            prm: 0.0,
        }
    }

    /// OCCT Perform(const gp_Sphere& S, const gp_Dir& D) (cxx L39-56) — the
    /// great circle of the plane through the center normal to D.
    pub fn perform_sphere_dir(&mut self, s: &SphericalSurface, d: Vec3) {
        self.done = false;
        self.typ_l = Curve2dType::Circle;
        self.pt1 = s.center;
        self.dir1 = d;
        let xax = s.ref_dir.normalize();
        let yax = s.axis.cross(s.ref_dir).normalize();
        if d.dot(xax).abs() < 0.9999999999999 {
            self.dir2 = d.cross(xax);
        } else {
            self.dir2 = d.cross(yax);
        }
        self.prm = s.radius;
        self.nb_sol = 1;
        self.done = true;
    }

    /// OCCT Perform(const gp_Sphere& S, const gp_Dir& D, const double Angle)
    /// (cxx L58-77) — the draft circle.
    pub fn perform_sphere_dir_angle(&mut self, s: &SphericalSurface, d: Vec3, angle: f64) {
        self.done = false;
        self.typ_l = Curve2dType::Circle;

        self.dir1 = d;
        let xax = s.ref_dir.normalize();
        let yax = s.axis.cross(s.ref_dir).normalize();
        if d.dot(xax).abs() < 0.9999999999999 {
            self.dir2 = d.cross(xax);
        } else {
            self.dir2 = d.cross(yax);
        }
        let direct = ax3_direct(s.ref_dir, s.axis.cross(s.ref_dir), s.axis);
        let alpha = if direct { angle } else { -angle };
        self.pt1 = s.center - s.radius * alpha.sin() * d;
        self.prm = s.radius * alpha.cos();
        self.nb_sol = 1;
        self.done = true;
    }

    /// OCCT Perform(const gp_Sphere& S, const gp_Pnt& Eye) (cxx L79-114) —
    /// the viewing circle.
    pub fn perform_sphere_eye(&mut self, s: &SphericalSurface, eye: Point3) {
        self.done = false;

        let radius = s.radius;
        let dist = eye.distance(s.center);
        if dist <= radius {
            self.nb_sol = 0;
        } else {
            self.prm = radius * (1.0 - radius * radius / (dist * dist)).sqrt();
            if self.prm < TOLPETIT {
                self.nb_sol = 0;
            } else {
                let locxyz = s.center;
                self.dir1 = eye - locxyz;
                self.pt1 = locxyz + (radius * radius / dist) * self.dir1;
                let xax = s.ref_dir.normalize();
                let yax = s.axis.cross(s.ref_dir).normalize();
                if self.dir1.dot(xax).abs() < 0.9999999999999 {
                    self.dir2 = self.dir1.cross(xax);
                } else {
                    self.dir2 = self.dir1.cross(yax);
                }
                self.nb_sol = 1;
                self.typ_l = Curve2dType::Circle;
            }
        }
        self.done = true;
    }

    /// OCCT Perform(const gp_Cylinder& C, const gp_Dir& D) (cxx L116-138) —
    /// the two straight contour generatrices.
    pub fn perform_cylinder_dir(&mut self, c: &CylindricalSurface, d: Vec3) {
        self.done = false;

        let mut normale = c.axis;
        normale = normale.cross(d);
        if normale.length() <= 1e-15 {
            self.nb_sol = 0;
        } else {
            normale = normale.normalize();
            self.typ_l = Curve2dType::Line;
            self.dir1 = c.axis;
            self.dir2 = self.dir1;
            self.pt1 = c.origin + c.radius * normale;
            self.pt2 = c.origin - c.radius * normale;
            self.nb_sol = 2;
        }

        self.done = true;
    }

    /// OCCT Perform(const gp_Cylinder& C, const gp_Dir& D, const double
    /// Angle) (cxx L140-198) — the two draft generatrices.
    pub fn perform_cylinder_dir_angle(&mut self, c: &CylindricalSurface, d: Vec3, angle: f64) {
        self.done = false;

        let mut coefcos = d.dot(c.ref_dir.normalize_or_zero());
        let mut coefsin = d.dot(c.y_axis());
        let coefcst = (std::f64::consts::PI * 0.5 + angle).cos();

        let norm1 = coefcos * coefcos + coefsin * coefsin;
        let norm2 = norm1.sqrt();

        if coefcst.abs() < norm2 {
            self.typ_l = Curve2dType::Line;
            self.nb_sol = 2;
            self.dir1 = c.axis;
            self.dir2 = self.dir1;

            let direct = ax3_direct(c.ref_dir, c.y_axis(), c.axis);
            if !direct {
                // The normal is inverted.
                coefcos = -coefcos;
                coefsin = -coefsin;
            }

            self.prm = (norm1 - coefcst * coefcst).sqrt();
            let cost0 = (coefcos * coefcst - coefsin * self.prm) / norm1;
            let cost1 = (coefcos * coefcst + coefsin * self.prm) / norm1;
            let sint0 = (coefcos * self.prm + coefsin * coefcst) / norm1;
            let sint1 = (-coefcos * self.prm + coefsin * coefcst) / norm1;

            let xdir = c.ref_dir.normalize_or_zero();
            let ydir = c.y_axis();

            let dirxyz = cost0 * xdir + sint0 * ydir;
            self.pt1 = c.origin + c.radius * dirxyz;

            let dirxyz = cost1 * xdir + sint1 * ydir;
            self.pt2 = c.origin + c.radius * dirxyz;
        } else {
            self.nb_sol = 0;
        }

        self.done = true;
    }

    /// OCCT Perform(const gp_Cylinder& C, const gp_Pnt& Eye) (cxx L200-226).
    pub fn perform_cylinder_eye(&mut self, c: &CylindricalSurface, eye: Point3) {
        self.done = false;

        let radius = c.radius;
        // gp_Lin theaxis(C.Axis()); dist = theaxis.Distance(Eye).
        let axis = c.axis.normalize();
        let w = eye - c.origin;
        let dist = (w - w.dot(axis) * axis).length();
        if dist <= radius {
            self.nb_sol = 0;
        } else {
            self.typ_l = Curve2dType::Line;
            self.prm = radius * (1.0 - radius * radius / (dist * dist)).sqrt();
            self.dir1 = c.axis;
            self.dir2 = self.dir1;
            let axeye = (w - w.dot(axis) * axis).normalize(); // orientate the axis to the outside
            let normale = axis.cross(axeye);
            self.pt1 = c.origin + (radius * radius / dist) * axeye;
            self.pt2 = self.pt1 - self.prm * normale;
            self.pt1 = self.pt1 + self.prm * normale;
            self.nb_sol = 2;
        }
        self.done = true;
    }

    /// OCCT Perform(const gp_Cone& C, const gp_Dir& D) (cxx L228-287).
    pub fn perform_cone_dir(&mut self, c: &ConicalSurface, d: Vec3) {
        self.done = false;

        let tgtalpha = c.half_angle_rad.tan();
        let xdir = c.ref_dir.normalize_or_zero();
        let ydir = c.axis.cross(c.ref_dir).normalize_or_zero();

        let coefcos = d.dot(xdir);
        let coefsin = d.dot(ydir);
        let coefcst = d.dot(c.axis) * tgtalpha;

        let norm1 = coefcos * coefcos + coefsin * coefsin;
        let norm2 = norm1.sqrt();

        if coefcst.abs() < norm2 {
            self.typ_l = Curve2dType::Line;
            self.nb_sol = 2;
            self.pt1 = c.apex;
            self.pt2 = self.pt1;

            self.prm = (norm1 - coefcst * coefcst).sqrt();
            let cost0 = (coefcos * coefcst - coefsin * self.prm) / norm1;
            let cost1 = (coefcos * coefcst + coefsin * self.prm) / norm1;
            let sint0 = (coefcos * self.prm + coefsin * coefcst) / norm1;
            let sint1 = (-coefcos * self.prm + coefsin * coefcst) / norm1;

            let zdir = c.axis;
            let dirxyz = cost0 * xdir + sint0 * ydir + (1.0 / tgtalpha) * zdir;
            self.dir1 = dirxyz;
            self.pt1 = self.pt1 + dirxyz;
            let dirxyz = cost1 * xdir + sint1 * ydir + (1.0 / tgtalpha) * zdir;
            self.dir2 = dirxyz;
            self.pt2 = self.pt2 + dirxyz;
        } else {
            self.nb_sol = 0;
        }
        self.done = true;
    }

    /// OCCT Perform(const gp_Cone& C, const gp_Dir& D, const double Angle)
    /// (cxx L289-391).
    pub fn perform_cone_dir_angle(&mut self, c: &ConicalSurface, d: Vec3, angle: f64) {
        self.done = false;
        self.nb_sol = 0;

        let ang = c.half_angle_rad;
        let cosa = ang.cos();
        let sina = ang.sin();
        let xdir = c.ref_dir.normalize_or_zero();
        let ydir = c.axis.cross(c.ref_dir).normalize_or_zero();

        let coefcos = d.dot(xdir);
        let coefsin = d.dot(ydir);

        let coefcst1 = (std::f64::consts::PI * 0.5 + angle).cos();

        let norm1 = coefcos * coefcos + coefsin * coefsin;
        let norm2 = norm1.sqrt();

        let coefnz = d.dot(c.axis) * sina;
        let mut coefcst = (coefcst1 + coefnz) / cosa;

        let direct = ax3_direct(c.ref_dir, ydir, c.axis);

        if coefcst.abs() < norm2 {
            self.typ_l = Curve2dType::Line;
            self.nb_sol += 2;
            self.pt1 = c.apex;
            self.pt2 = self.pt1;

            self.prm = (norm1 - coefcst * coefcst).sqrt();
            let cost0 = (coefcos * coefcst - coefsin * self.prm) / norm1;
            let cost1 = (coefcos * coefcst + coefsin * self.prm) / norm1;
            let sint0 = (coefcos * self.prm + coefsin * coefcst) / norm1;
            let sint1 = (-coefcos * self.prm + coefsin * coefcst) / norm1;

            let mut zdir = c.axis;
            if !direct {
                zdir = -zdir;
            }
            let dirxyz = cost0 * xdir + sint0 * ydir + (cosa / sina) * zdir;
            self.dir1 = dirxyz;
            self.pt1 = self.pt1 + dirxyz;
            let dirxyz = cost1 * xdir + sint1 * ydir + (cosa / sina) * zdir;
            self.dir2 = dirxyz;
            self.pt2 = self.pt2 + dirxyz;
        }

        coefcst = (coefcst1 - coefnz) / cosa;

        if coefcst.abs() < norm2 {
            self.typ_l = Curve2dType::Line;
            self.nb_sol += 2;
            self.pt3 = c.apex;
            self.pt4 = self.pt3;

            self.prm = (norm1 - coefcst * coefcst).sqrt();
            let cost0 = (coefcos * coefcst - coefsin * self.prm) / norm1;
            let cost1 = (coefcos * coefcst + coefsin * self.prm) / norm1;
            let sint0 = (coefcos * self.prm + coefsin * coefcst) / norm1;
            let sint1 = (-coefcos * self.prm + coefsin * coefcst) / norm1;

            let mut zdir = c.axis;
            if !direct {
                zdir = -zdir;
            }
            let dirxyz = cost0 * xdir + sint0 * ydir + (-cosa / sina) * zdir;
            self.dir3 = dirxyz;
            self.pt3 = self.pt3 + dirxyz;
            let dirxyz = cost1 * xdir + sint1 * ydir + (-cosa / sina) * zdir;
            self.dir4 = dirxyz;
            self.pt4 = self.pt4 + dirxyz;
            if self.nb_sol == 2 {
                self.pt1 = self.pt3;
                self.pt2 = self.pt4;
                self.dir1 = self.dir3;
                self.dir2 = self.dir4;
            }
        }

        self.done = true;
    }

    /// OCCT Perform(const gp_Cone& C, const gp_Pnt& Eye) (cxx L393-455).
    pub fn perform_cone_eye(&mut self, c: &ConicalSurface, eye: Point3) {
        self.done = false;

        let tgtalpha = c.half_angle_rad.tan();
        let xdir = c.ref_dir.normalize_or_zero();
        let ydir = c.axis.cross(c.ref_dir).normalize_or_zero();

        let apexeye = eye - c.apex;

        let coefcos = apexeye.dot(xdir);
        let coefsin = apexeye.dot(ydir);
        let coefcst = apexeye.dot(c.axis) * tgtalpha;

        let norm1 = coefcos * coefcos + coefsin * coefsin;
        let norm2 = norm1.sqrt();

        if coefcst.abs() < norm2 {
            self.typ_l = Curve2dType::Line;
            self.nb_sol = 2;
            self.pt1 = c.apex;
            self.pt2 = self.pt1;

            self.prm = (norm1 - coefcst * coefcst).sqrt();
            let cost0 = (coefcos * coefcst - coefsin * self.prm) / norm1;
            let cost1 = (coefcos * coefcst + coefsin * self.prm) / norm1;
            let sint0 = (coefcos * self.prm + coefsin * coefcst) / norm1;
            let sint1 = (-coefcos * self.prm + coefsin * coefcst) / norm1;

            let zdir = c.axis;
            let dirxyz = cost0 * xdir + sint0 * ydir + (1.0 / tgtalpha) * zdir;
            self.dir1 = dirxyz;
            self.pt1 = self.pt1 + dirxyz;
            let dirxyz = cost1 * xdir + sint1 * ydir + (1.0 / tgtalpha) * zdir;
            self.dir2 = dirxyz;
            self.pt2 = self.pt2 + dirxyz;
        } else {
            self.nb_sol = 0;
        }
        self.done = true;
    }

    /// OCCT IsDone (Contap_ContAna.lxx) — done flag.
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT NbContours.
    pub fn nb_contours(&self) -> i32 {
        self.nb_sol
    }

    /// OCCT TypeContour — GeomAbs_Line or GeomAbs_Circle when done; rcad
    /// carries the Curve2dType (Line / Circle / OtherCurve).
    pub fn type_contour(&self) -> Curve2dType {
        self.typ_l
    }

    /// OCCT Circle (lxx) — gp_Circ(gp_Ax2(pt1, dir1, dir2), prm).
    pub fn circle(&self) -> Circle3 {
        Circle3 {
            center: self.pt1,
            normal: self.dir1,
            x_dir: self.dir2,
            y_dir: self.dir1.cross(self.dir2),
            radius: self.prm,
        }
    }

    /// OCCT Line(Index) (cxx L457-483) — 1-based.
    pub fn line(&self, index: i32) -> rcad_kernel::geom::Line3 {
        if !self.done {
            panic!("StdFail_NotDone: Contap_ContAna::Line");
        }
        if self.typ_l != Curve2dType::Line || self.nb_sol == 0 {
            panic!("Standard_DomainError: Contap_ContAna::Line");
        }
        if index <= 0 || index > self.nb_sol {
            panic!("Standard_OutOfRange: Contap_ContAna::Line");
        }
        match index {
            1 => rcad_kernel::geom::Line3 {
                origin: self.pt1,
                direction: self.dir1,
            },
            2 => rcad_kernel::geom::Line3 {
                origin: self.pt2,
                direction: self.dir2,
            },
            3 => rcad_kernel::geom::Line3 {
                origin: self.pt3,
                direction: self.dir3,
            },
            4 => rcad_kernel::geom::Line3 {
                origin: self.pt4,
                direction: self.dir4,
            },
            _ => panic!("Standard_OutOfRange: Program error in Contap_ContAna"),
        }
    }
}

impl Default for ContAna {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT gp_Ax3::Direct — X ^ Y . Z > 0.
fn ax3_direct(x: Vec3, y: Vec3, z: Vec3) -> bool {
    x.cross(y).dot(z) > 0.0
}

#[allow(unused)]
fn _shape_parity(p: Point3, d: DVec3) -> (Point3, DVec3) {
    (p, d)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_cylinder() -> CylindricalSurface {
        CylindricalSurface {
            origin: DVec3::new(1.0, 2.0, 3.0),
            axis: DVec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
            ref_dir: DVec3::new(1.0, 0.0, 0.0),
            y_dir: None,
        }
    }

    /// OCCT anchor: cylinder viewed along -Z (perpendicular to the axis)
    /// — the two silhouette lines at x = center +/- R (cxx L116-138).
    #[test]
    fn cont_ana_cylinder_dir() {
        let mut ana = ContAna::new();
        ana.perform_cylinder_dir(&unit_cylinder(), DVec3::new(0.0, -1.0, 0.0));
        assert!(ana.is_done());
        assert_eq!(ana.nb_contours(), 2);
        assert_eq!(ana.type_contour(), Curve2dType::Line);
        let l1 = ana.line(1);
        assert!((l1.origin.x - 3.0).abs() < 1e-12 && (l1.origin.y - 2.0).abs() < 1e-12);
        let l2 = ana.line(2);
        assert!((l2.origin.x + 1.0).abs() < 1e-12 && (l2.origin.y - 2.0).abs() < 1e-12);
        // Both direction along +Z.
        assert!((l1.direction.z - 1.0).abs() < 1e-12);
        assert!((l2.direction.z - 1.0).abs() < 1e-12);
    }

    /// OCCT anchor: cylinder viewed along the axis — no silhouette
    /// (normale.Modulus() <= 1e-15).
    #[test]
    fn cont_ana_cylinder_axis_view() {
        let mut ana = ContAna::new();
        ana.perform_cylinder_dir(&unit_cylinder(), DVec3::new(0.0, 0.0, 1.0));
        assert!(ana.is_done());
        assert_eq!(ana.nb_contours(), 0);
    }

    /// OCCT anchor: sphere viewed along Z — the great circle in the
    /// equatorial plane, radius = sphere radius (cxx L39-56).
    #[test]
    fn cont_ana_sphere_dir() {
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::new(0.0, 0.0, 1.0),
            radius: 3.0,
            ref_dir: DVec3::new(1.0, 0.0, 0.0),
        };
        let mut ana = ContAna::new();
        ana.perform_sphere_dir(&sphere, DVec3::new(0.0, 0.0, 1.0));
        assert!(ana.is_done());
        assert_eq!(ana.nb_contours(), 1);
        assert_eq!(ana.type_contour(), Curve2dType::Circle);
        let c = ana.circle();
        assert!((c.radius - 3.0).abs() < 1e-12);
        assert!(c.center.distance(DVec3::ZERO) < 1e-12);
    }

    /// OCCT anchor: cone draft at 0 angle along axis — two silhouette
    /// generatrices through the apex (cxx L228-287).
    #[test]
    fn cont_ana_cone_dir() {
        let cone = ConicalSurface {
            apex: DVec3::new(0.0, 0.0, 0.0),
            axis: DVec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
            half_angle_rad: (30f64).to_radians(),
            ref_dir: DVec3::new(1.0, 0.0, 0.0),
        };
        let mut ana = ContAna::new();
        ana.perform_cone_dir(&cone, DVec3::new(0.0, -1.0, 0.0));
        assert!(ana.is_done());
        assert_eq!(ana.nb_contours(), 2);
        assert_eq!(ana.type_contour(), Curve2dType::Line);
        // Solution direction z-component = 1 / tan(sida) = cot(30 deg).
        let l1 = ana.line(1);
        let expected = 1.0 / (30f64).to_radians().tan();
        assert!((l1.direction.z - expected).abs() < 1e-12);
    }

    // The gp_TopAbs shape import stays referenced (rcad topods parity).
    #[allow(unused)]
    fn _t() {
        let _ = _topods_shape::Orientation::Forward;
    }
}
