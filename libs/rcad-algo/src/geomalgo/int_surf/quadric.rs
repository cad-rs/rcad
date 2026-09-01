// OCCT IntSurf_Quadric — an implicit (analytic) surface in a unified form.
//
// OCCT IntSurf_Quadric.hxx / .cxx / .lxx. Stores a quadric (Plane, Cylinder,
// Sphere, Cone, Torus) as a gp_Ax3 frame + parameter list, and provides the
// algebraic distance / gradient / value / parameters used by the intersection
// algorithms (IntPatch_ImpImpIntersection, IntWalk).
//
// rcad data-model notes:
// - OCCT `gp_Ax3 ax3` maps to (loc, x_dir, y_dir, z_dir) + ax3_direct.
// - OCCT `gp_Lin lin` (the axis line) maps to (lin_origin, lin_dir); only
//   set for Cylinder/Sphere/Cone/Torus (SetPosition(ax3.Axis())).
// - Constructors take the rcad-kernel surface types; the frame is built from
//   the surface's own axes (Plane::u_dir/v_dir, Cylinder/Sphere::ref_dir,
//   any_perpendicular(axis) for Cone/Torus) so that Value(u,v) coincides
//   with the rcad surface's point_at(u,v).

use glam::{DVec2, DVec3};

use rcad_kernel::geom::{
    any_perpendicular, ConicalSurface, CylindricalSurface, Plane, SphericalSurface,
    ToroidalSurface,
};

/// OCCT GeomAbs_SurfaceType — the analytic surface types the quadric handles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QuadricType {
    Plane,
    Cylinder,
    Sphere,
    Cone,
    Torus,
    Other,
}

/// OCCT IntSurf_Quadric — unified quadric surface representation.
///
/// Fields (IntSurf_Quadric.hxx L102-109):
///   ax3: gp_Ax3            -> loc, x_dir, y_dir, z_dir (+ ax3_direct)
///   lin: gp_Lin            -> lin_origin, lin_dir
///   typ: GeomAbs_SurfaceType
///   prm1..prm4: double
///   ax3direc: bool
#[derive(Debug, Clone)]
pub struct Quadric {
    // gp_Ax3
    loc: DVec3,
    x_dir: DVec3,
    y_dir: DVec3,
    z_dir: DVec3,
    // gp_Lin (axis line)
    lin_origin: DVec3,
    lin_dir: DVec3,
    // typ
    typ: QuadricType,
    // prm1..prm4
    prm1: f64,
    prm2: f64,
    prm3: f64,
    prm4: f64,
    // ax3direc
    ax3_direct: bool,
}

impl Quadric {
    /// OCCT IntSurf_Quadric() — default: typ = GeomAbs_OtherSurface.
    pub fn new() -> Self {
        Quadric {
            loc: DVec3::ZERO,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
            z_dir: DVec3::Z,
            lin_origin: DVec3::ZERO,
            lin_dir: DVec3::Z,
            typ: QuadricType::Other,
            prm1: 0.0,
            prm2: 0.0,
            prm3: 0.0,
            prm4: 0.0,
            ax3_direct: false,
        }
    }

    /// OCCT IntSurf_Quadric(const gp_Pln&) (IntSurf_Quadric.cxx L49-55).
    pub fn from_plane(p: &Plane) -> Self {
        // P.Coefficients(prm1, prm2, prm3, prm4)
        let d = -p.normal.dot(p.origin);
        Quadric {
            loc: p.origin,
            x_dir: p.u_dir,
            y_dir: p.v_dir,
            z_dir: p.normal,
            lin_origin: DVec3::ZERO,
            lin_dir: DVec3::Z,
            typ: QuadricType::Plane,
            prm1: p.normal.x,
            prm2: p.normal.y,
            prm3: p.normal.z,
            prm4: d,
            ax3_direct: p.normal.cross(p.u_dir).dot(p.v_dir) > 0.0,
        }
    }

    /// OCCT IntSurf_Quadric(const gp_Cylinder&) (L58-68).
    pub fn from_cylinder(c: &CylindricalSurface) -> Self {
        let z = c.axis.normalize_or_zero();
        let x = c.ref_dir.normalize_or_zero();
        let y = z.cross(x).normalize_or_zero();
        Quadric {
            loc: c.origin,
            x_dir: x,
            y_dir: y,
            z_dir: z,
            lin_origin: c.origin,
            lin_dir: z,
            typ: QuadricType::Cylinder,
            prm1: c.radius,
            prm2: 0.0,
            prm3: 0.0,
            prm4: 0.0,
            ax3_direct: z.cross(x).dot(y) > 0.0,
        }
    }

    /// OCCT IntSurf_Quadric(const gp_Sphere&) (L71-81).
    pub fn from_sphere(s: &SphericalSurface) -> Self {
        let z = s.axis.normalize_or_zero();
        let x = s.ref_dir.normalize_or_zero();
        let y = z.cross(x).normalize_or_zero();
        Quadric {
            loc: s.center,
            x_dir: x,
            y_dir: y,
            z_dir: z,
            lin_origin: s.center,
            lin_dir: z,
            typ: QuadricType::Sphere,
            prm1: s.radius,
            prm2: 0.0,
            prm3: 0.0,
            prm4: 0.0,
            ax3_direct: z.cross(x).dot(y) > 0.0,
        }
    }

    /// OCCT IntSurf_Quadric(const gp_Cone&) (L84-96) — the cone's gp_Ax3
    /// (including its XDirection) is preserved, defining the u=0 generatrix.
    pub fn from_cone(c: &ConicalSurface) -> Self {
        let z = c.axis.normalize_or_zero();
        let x = c.ref_dir.normalize_or_zero();
        let y = z.cross(x).normalize_or_zero();
        let prm2 = c.half_angle_rad;
        Quadric {
            loc: c.apex,
            x_dir: x,
            y_dir: y,
            z_dir: z,
            lin_origin: c.apex,
            lin_dir: z,
            typ: QuadricType::Cone,
            prm1: c.radius,
            prm2,
            prm3: prm2.cos(),
            prm4: 0.0,
            ax3_direct: z.cross(x).dot(y) > 0.0,
        }
    }

    /// OCCT IntSurf_Quadric(const gp_Torus&) (L99-111).
    pub fn from_torus(t: &ToroidalSurface) -> Self {
        let z = t.axis.normalize_or_zero();
        let x = any_perpendicular(z).normalize_or_zero();
        let y = z.cross(x).normalize_or_zero();
        Quadric {
            loc: t.center,
            x_dir: x,
            y_dir: y,
            z_dir: z,
            lin_origin: t.center,
            lin_dir: z,
            typ: QuadricType::Torus,
            prm1: t.major_radius,
            prm2: t.minor_radius,
            prm3: 0.0,
            prm4: 0.0,
            ax3_direct: z.cross(x).dot(y) > 0.0,
        }
    }

    /// Build from any rcad Surface3. Returns None for non-quadric surfaces.
    pub fn from_surface3(surf: &rcad_kernel::geom::Surface3) -> Option<Self> {
        match surf {
            rcad_kernel::geom::Surface3::Plane(p) => Some(Self::from_plane(p)),
            rcad_kernel::geom::Surface3::Cylinder(c) => Some(Self::from_cylinder(c)),
            rcad_kernel::geom::Surface3::Sphere(s) => Some(Self::from_sphere(s)),
            rcad_kernel::geom::Surface3::Cone(c) => Some(Self::from_cone(c)),
            rcad_kernel::geom::Surface3::Torus(t) => Some(Self::from_torus(t)),
            _ => None,
        }
    }

    /// OCCT TypeQuadric().
    pub fn type_quadric(&self) -> QuadricType {
        self.typ
    }

    // === Accessors used by IntAna_QuadQuadGeo / IntPatch_ImpImpIntersection ===

    /// OCCT IntSurf_Quadric::TypeQuadric — alias of type_quadric.
    pub fn surface_type(&self) -> QuadricType {
        self.typ
    }

    /// OCCT ax3.Location().
    pub fn location(&self) -> DVec3 {
        self.loc
    }

    /// OCCT ax3.XDirection().
    pub fn x_dir(&self) -> DVec3 {
        self.x_dir
    }

    /// OCCT ax3.YDirection().
    pub fn y_dir(&self) -> DVec3 {
        self.y_dir
    }

    /// OCCT ax3.Direction().
    pub fn z_dir(&self) -> DVec3 {
        self.z_dir
    }

    /// OCCT lin.Direction() — axis line direction (quadrics); plane uses normal.
    pub fn axis_dir(&self) -> DVec3 {
        match self.typ {
            QuadricType::Plane => self.z_dir,
            _ => self.lin_dir,
        }
    }

    /// OCCT lin.Location() — axis line origin (quadrics); plane uses location.
    pub fn axis_loc(&self) -> DVec3 {
        match self.typ {
            QuadricType::Plane => self.loc,
            _ => self.lin_origin,
        }
    }

    /// OCCT ax3.Direct().
    pub fn ax3_direct(&self) -> bool {
        self.ax3_direct
    }

    /// OCCT prm1 — radius for Cylinder/Sphere, RefRadius for Cone, MajorRadius for Torus.
    pub fn radius(&self) -> f64 {
        self.prm1
    }
    pub fn ref_radius(&self) -> f64 {
        self.prm1
    }
    pub fn major_radius(&self) -> f64 {
        self.prm1
    }

    /// OCCT prm2 — SemiAngle for Cone, MinorRadius for Torus.
    pub fn semi_angle(&self) -> f64 {
        self.prm2
    }
    pub fn minor_radius(&self) -> f64 {
        self.prm2
    }

    /// OCCT P.Coefficients(A, B, C, D) — plane coefficients.
    pub fn plane_coeffs(&self) -> (f64, f64, f64, f64) {
        (self.prm1, self.prm2, self.prm3, self.prm4)
    }

    /// OCCT Plane() (IntSurf_Quadric.lxx L27-30) — reconstruct the plane.
    pub fn plane(&self) -> Plane {
        Plane {
            origin: self.loc,
            normal: self.z_dir,
            u_dir: self.x_dir,
            v_dir: self.y_dir,
        }
    }

    /// OCCT Sphere() (L35-36) — gp_Sphere(ax3, prm1).
    pub fn sphere(&self) -> SphericalSurface {
        SphericalSurface {
            center: self.loc,
            axis: self.z_dir,
            radius: self.prm1,
            ref_dir: self.x_dir,
        }
    }

    /// OCCT Cylinder() (L39-40) — gp_Cylinder(ax3, prm1).
    pub fn cylinder(&self) -> CylindricalSurface {
        CylindricalSurface {
            origin: self.loc,
            axis: self.z_dir,
            radius: self.prm1,
            ref_dir: self.x_dir,
            y_dir: None,
        }
    }

    /// OCCT Cone() (L45-47) — gp_Cone(ax3, prm2, prm1).
    pub fn cone(&self) -> ConicalSurface {
        ConicalSurface {
            apex: self.loc,
            axis: self.z_dir,
            radius: self.prm1,
            half_angle_rad: self.prm2,
            ref_dir: self.x_dir,
        }
    }

    /// OCCT Torus() (L51-53) — gp_Torus(ax3, prm1, prm2).
    pub fn torus(&self) -> ToroidalSurface {
        ToroidalSurface {
            center: self.loc,
            axis: self.z_dir,
            ref_dir: self.x_dir,
            major_radius: self.prm1,
            minor_radius: self.prm2,
        }
    }

    // === Distance / Gradient / ValAndGrad (IntSurf_Quadric.cxx L171-403) ===

    /// OCCT Distance(P) — algebraic distance from the point to the quadric.
    pub fn distance(&self, p: DVec3) -> f64 {
        match self.typ {
            QuadricType::Plane => {
                self.prm1 * p.x + self.prm2 * p.y + self.prm3 * p.z + self.prm4
            }
            QuadricType::Cylinder => line_distance(&self.lin_origin, &self.lin_dir, p) - self.prm1,
            QuadricType::Sphere => (self.lin_origin - p).length() - self.prm1,
            QuadricType::Cone => {
                let dist = line_distance(&self.lin_origin, &self.lin_dir, p);
                let (u, v) = self.parameters(p);
                let pp = self.value(u, v);
                let distp = line_distance(&self.lin_origin, &self.lin_dir, pp);
                (dist - distp) / self.prm3
            }
            QuadricType::Torus => {
                // O = ax3.Location(); Pp = P translated by -(O->P . OZ) * OZ
                let o = self.loc;
                let oz = self.z_dir;
                let pp = p - oz * ((p - o).dot(oz));
                // DOPp = (O.SquareDistance(Pp) < 1e-14) ? ax3.XDirection() : Dir(O, Pp)
                let dop_p = if (o - pp).length_squared() < 1e-14 {
                    self.x_dir
                } else {
                    (pp - o).normalize_or_zero()
                };
                let pt = o + dop_p * self.prm1;
                (p - pt).length() - self.prm2
            }
            QuadricType::Other => 0.0,
        }
    }

    /// OCCT Gradient(P) — gradient of the algebraic distance.
    pub fn gradient(&self, p: DVec3) -> DVec3 {
        match self.typ {
            QuadricType::Plane => DVec3::new(self.prm1, self.prm2, self.prm3),
            QuadricType::Cylinder => {
                // PP = lin.Location() + Parameter(lin, P) * lin.Direction()
                let pp = self.lin_origin + elclib_parameter(&self.lin_origin, &self.lin_dir, p) * self.lin_dir;
                let mut grad = p - pp;
                let n = grad.length();
                if n > 1e-14 {
                    grad / n
                } else {
                    DVec3::ZERO
                }
            }
            QuadricType::Sphere => {
                let mut grad = p - self.lin_origin;
                let n = grad.length();
                if n > 1e-14 {
                    grad / n
                } else {
                    DVec3::ZERO
                }
            }
            QuadricType::Cone => {
                let (u, v) = self.parameters(p);
                let (_, d1u, d1v) = self.cone_d1(u, v);
                let mut grad = d1u.cross(d1v);
                if !self.ax3_direct {
                    grad = -grad;
                }
                if has_magnitude_for_normalization(&grad) {
                    grad.normalize()
                } else {
                    DVec3::ZERO
                }
            }
            QuadricType::Torus => {
                let o = self.loc;
                let oz = self.z_dir;
                let pp = p - oz * ((p - o).dot(oz));
                let dop_p = if (o - pp).length_squared() < 1e-14 {
                    self.x_dir
                } else {
                    (pp - o).normalize_or_zero()
                };
                let pt = o + dop_p * self.prm1;
                let mut grad = p - pt;
                let n = grad.length();
                if n > 1e-14 {
                    grad / n
                } else {
                    DVec3::ZERO
                }
            }
            QuadricType::Other => DVec3::ZERO,
        }
    }

    /// OCCT ValAndGrad(P, Dist, Grad).
    pub fn val_and_grad(&self, p: DVec3) -> (f64, DVec3) {
        match self.typ {
            QuadricType::Plane => (
                self.prm1 * p.x + self.prm2 * p.y + self.prm3 * p.z + self.prm4,
                DVec3::new(self.prm1, self.prm2, self.prm3),
            ),
            QuadricType::Cylinder => {
                let dist = line_distance(&self.lin_origin, &self.lin_dir, p) - self.prm1;
                let pp = self.lin_origin + elclib_parameter(&self.lin_origin, &self.lin_dir, p) * self.lin_dir;
                let mut grad = p - pp;
                let n = grad.length();
                if n > 1e-14 {
                    grad /= n;
                } else {
                    grad = DVec3::ZERO;
                }
                (dist, grad)
            }
            QuadricType::Sphere => {
                let dist = (self.lin_origin - p).length() - self.prm1;
                let mut grad = p - self.lin_origin;
                let n = grad.length();
                if n > 1e-14 {
                    grad /= n;
                } else {
                    grad = DVec3::ZERO;
                }
                (dist, grad)
            }
            QuadricType::Cone => {
                let dist = line_distance(&self.lin_origin, &self.lin_dir, p);
                let (u, v) = self.parameters(p);
                let (_, d1u, d1v) = self.cone_d1(u, v);
                let distp = line_distance(&self.lin_origin, &self.lin_dir, self.value(u, v));
                let dist = (dist - distp) / self.prm3;
                let mut grad = d1u.cross(d1v);
                if !self.ax3_direct {
                    grad = -grad;
                }
                if has_magnitude_for_normalization(&grad) {
                    grad = grad.normalize();
                } else {
                    grad = DVec3::ZERO;
                }
                (dist, grad)
            }
            QuadricType::Torus => {
                let o = self.loc;
                let oz = self.z_dir;
                let pp = p - oz * ((p - o).dot(oz));
                let dop_p = if (o - pp).length_squared() < 1e-14 {
                    self.x_dir
                } else {
                    (pp - o).normalize_or_zero()
                };
                let pt = o + dop_p * self.prm1;
                let dist = (p - pt).length() - self.prm2;
                let mut grad = p - pt;
                let n = grad.length();
                if n > 1e-14 {
                    grad /= n;
                } else {
                    grad = DVec3::ZERO;
                }
                (dist, grad)
            }
            QuadricType::Other => (0.0, DVec3::ZERO),
        }
    }

    // === Value / D1 / DN (IntSurf_Quadric.cxx L406-480) ===

    /// OCCT Value(U, V) — ElSLib::*Value.
    pub fn value(&self, u: f64, v: f64) -> DVec3 {
        match self.typ {
            QuadricType::Plane => self.loc + u * self.x_dir + v * self.y_dir,
            QuadricType::Cylinder => {
                self.loc + self.prm1 * (u.cos() * self.x_dir + u.sin() * self.y_dir) + v * self.z_dir
            }
            QuadricType::Sphere => {
                // ElSLib::SphereValue — latitude convention (V=0 at equator).
                let r = self.prm1 * v.cos();
                let a3 = self.prm1 * v.sin();
                self.loc + r * (u.cos() * self.x_dir + u.sin() * self.y_dir) + a3 * self.z_dir
            }
            QuadricType::Cone => {
                let r = self.prm1 + v * self.prm2.sin();
                let a3 = v * self.prm2.cos();
                self.loc + r * (u.cos() * self.x_dir + u.sin() * self.y_dir) + a3 * self.z_dir
            }
            QuadricType::Torus => {
                let r = self.prm1 + self.prm2 * v.cos();
                let a3 = self.prm2 * v.sin();
                self.loc + r * (u.cos() * self.x_dir + u.sin() * self.y_dir) + a3 * self.z_dir
            }
            QuadricType::Other => DVec3::ZERO,
        }
    }

    /// OCCT D1(U, V, P, D1U, D1V).
    pub fn d1(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        match self.typ {
            QuadricType::Plane => (
                self.loc + u * self.x_dir + v * self.y_dir,
                self.x_dir,
                self.y_dir,
            ),
            QuadricType::Cylinder => {
                let (su, cu) = u.sin_cos();
                let p = self.loc + self.prm1 * (cu * self.x_dir + su * self.y_dir) + v * self.z_dir;
                let d1u = self.prm1 * (-su * self.x_dir + cu * self.y_dir);
                (p, d1u, self.z_dir)
            }
            QuadricType::Sphere => {
                // ElSLib::SphereD1 (latitude convention).
                let (su, cu) = u.sin_cos();
                let (sv, cv) = v.sin_cos();
                let r = self.prm1;
                let radial = cu * self.x_dir + su * self.y_dir;
                let d1u = r * cv * (-su * self.x_dir + cu * self.y_dir);
                let d1v = r * (-sv * radial + cv * self.z_dir);
                let p = self.loc + r * (cv * radial + sv * self.z_dir);
                (p, d1u, d1v)
            }
            QuadricType::Cone => self.cone_d1(u, v),
            QuadricType::Torus => {
                let (su, cu) = u.sin_cos();
                let (sv, cv) = v.sin_cos();
                let r_maj = self.prm1;
                let r_min = self.prm2;
                let r = r_maj + r_min * cv;
                let radial = cu * self.x_dir + su * self.y_dir;
                let r_perp = -su * self.x_dir + cu * self.y_dir;
                let p = self.loc + r * radial + r_min * sv * self.z_dir;
                let d1u = r * r_perp;
                let d1v = -r_min * sv * radial + r_min * cv * self.z_dir;
                (p, d1u, d1v)
            }
            QuadricType::Other => (DVec3::ZERO, DVec3::ZERO, DVec3::ZERO),
        }
    }

    /// OCCT ElSLib::ConeD1 — used by Gradient/ValAndGrad and Normale.
    fn cone_d1(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let (su, cu) = u.sin_cos();
        let sang = self.prm2;
        let (ss, cs) = sang.sin_cos();
        let r = self.prm1 + v * ss;
        let radial = cu * self.x_dir + su * self.y_dir;
        let r_perp = -su * self.x_dir + cu * self.y_dir;
        let p = self.loc + r * radial + v * cs * self.z_dir;
        let d1u = r * r_perp;
        let d1v = ss * radial + cs * self.z_dir;
        (p, d1u, d1v)
    }

    /// OCCT DN(U, V, Nu, Nv) — the Nu+Nv-th partial derivative (ElSLib::*DN).
    pub fn dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> DVec3 {
        match self.typ {
            QuadricType::Plane => {
                if nu == 0 && nv == 1 {
                    self.y_dir
                } else if nu == 1 && nv == 0 {
                    self.x_dir
                } else {
                    DVec3::ZERO
                }
            }
            QuadricType::Cylinder => {
                // ElSLib::CylinderDN: Nu==1 -> -r sin(u) X + r cos(u) Y, Nu==0 -> Z.
                let r = self.prm1;
                match (nu, nv) {
                    (1, 0) => r * (-u.sin() * self.x_dir + u.cos() * self.y_dir),
                    (0, 1) => self.z_dir,
                    _ => DVec3::ZERO,
                }
            }
            QuadricType::Sphere => {
                // ElSLib::SphereDN (latitude): Nu==1,Nv==0 ->
                // r cos(v)(-sin(u)X + cos(u)Y); Nu==0,Nv==1 ->
                // r(-sin(v)(cos(u)X + sin(u)Y) + cos(v)Z).
                let r = self.prm1;
                let (su, cu) = u.sin_cos();
                let (sv, cv) = v.sin_cos();
                match (nu, nv) {
                    (1, 0) => r * cv * (-su * self.x_dir + cu * self.y_dir),
                    (0, 1) => r * (-sv * (cu * self.x_dir + su * self.y_dir) + cv * self.z_dir),
                    _ => DVec3::ZERO,
                }
            }
            QuadricType::Cone => {
                let sang = self.prm2;
                let (ss, cs) = sang.sin_cos();
                match (nu, nv) {
                    (1, 0) => {
                        let r = self.prm1 + v * ss;
                        r * (-u.sin() * self.x_dir + u.cos() * self.y_dir)
                    }
                    (0, 1) => ss * (u.cos() * self.x_dir + u.sin() * self.y_dir) + cs * self.z_dir,
                    _ => DVec3::ZERO,
                }
            }
            QuadricType::Torus => {
                let r_maj = self.prm1;
                let r_min = self.prm2;
                let (su, cu) = u.sin_cos();
                let (sv, cv) = v.sin_cos();
                match (nu, nv) {
                    (1, 0) => {
                        let r = r_maj + r_min * cv;
                        r * (-su * self.x_dir + cu * self.y_dir)
                    }
                    (0, 1) => {
                        -r_min * sv * (cu * self.x_dir + su * self.y_dir)
                            + r_min * cv * self.z_dir
                    }
                    _ => DVec3::ZERO,
                }
            }
            QuadricType::Other => DVec3::ZERO,
        }
    }

    // === Normale (IntSurf_Quadric.cxx L483-587) ===

    /// OCCT Normale(U, V).
    pub fn normale_uv(&self, u: f64, v: f64) -> DVec3 {
        match self.typ {
            QuadricType::Plane => {
                if self.ax3_direct {
                    self.z_dir
                } else {
                    -self.z_dir
                }
            }
            QuadricType::Cylinder => self.normale(self.value(u, v)),
            QuadricType::Sphere => self.normale(self.value(u, v)),
            QuadricType::Cone => {
                let (_, d1u, d1v) = self.cone_d1(u, v);
                if d1u.length() < 0.0000001 {
                    DVec3::ZERO
                } else {
                    d1u.cross(d1v)
                }
            }
            QuadricType::Torus => self.normale(self.value(u, v)),
            QuadricType::Other => DVec3::ZERO,
        }
    }

    /// OCCT Normale(P).
    pub fn normale(&self, p: DVec3) -> DVec3 {
        match self.typ {
            QuadricType::Plane => {
                if self.ax3_direct {
                    self.z_dir
                } else {
                    -self.z_dir
                }
            }
            QuadricType::Cylinder => {
                // lin.Normal(P).Direction(), reversed if not direct.
                let d = line_normal(&self.lin_origin, &self.lin_dir, p);
                if self.ax3_direct {
                    d
                } else {
                    -d
                }
            }
            QuadricType::Sphere => {
                if self.ax3_direct {
                    (p - self.loc).normalize_or_zero()
                } else {
                    (self.loc - p).normalize_or_zero()
                }
            }
            QuadricType::Cone => {
                let (u, v) = self.parameters(p);
                self.normale_uv(u, v)
            }
            QuadricType::Torus => {
                let o = self.loc;
                let oz = self.z_dir;
                let pp = p - oz * ((p - o).dot(oz));
                let dop_p = if (o - pp).length_squared() < 1e-14 {
                    self.x_dir
                } else {
                    (pp - o).normalize_or_zero()
                };
                let pt = o + dop_p * self.prm1;
                if (pt - p).length_squared() < 1e-14 {
                    return oz;
                }
                let d = if self.ax3_direct {
                    (p - pt).normalize_or_zero()
                } else {
                    (pt - p).normalize_or_zero()
                };
                d
            }
            QuadricType::Other => DVec3::ZERO,
        }
    }

    /// OCCT Parameters(P, U, V) — ElSLib::*Parameters.
    pub fn parameters(&self, p: DVec3) -> (f64, f64) {
        match self.typ {
            QuadricType::Plane => ((p - self.loc).dot(self.x_dir), (p - self.loc).dot(self.y_dir)),
            QuadricType::Cylinder => {
                let ploc = p - self.loc;
                let mut u = ploc.y.atan2(ploc.x);
                normalize_angle(&mut u);
                (u, ploc.z)
            }
            QuadricType::Sphere => {
                let ploc = p - self.loc;
                let l = (ploc.x * ploc.x + ploc.y * ploc.y).sqrt();
                let (mut u, v);
                if l < f64::MIN_POSITIVE {
                    // point on axis Z of the sphere
                    v = if ploc.z > 0.0 { std::f64::consts::FRAC_PI_2 } else { -std::f64::consts::FRAC_PI_2 };
                    u = 0.0;
                } else {
                    v = (ploc.z / l).atan();
                    u = ploc.y.atan2(ploc.x);
                    normalize_angle(&mut u);
                }
                (u, v)
            }
            QuadricType::Cone => {
                // OCCT ElSLib::ConeParameters (L1574-1611): the point is
                // transformed into the cone's local frame (gp_Ax3) before
                // U/V are computed.
                let ploc = p - self.loc;
                let px = ploc.dot(self.x_dir);
                let py = ploc.dot(self.y_dir);
                let pz = ploc.dot(self.z_dir);
                let sang = self.prm2;
                let (ss, cs) = sang.sin_cos();
                let mut u;
                if px.abs() < f64::MIN_POSITIVE && py.abs() < f64::MIN_POSITIVE {
                    u = 0.0;
                } else if -self.prm1 > pz * sang.tan() {
                    // the point is at the wrong side of the apex
                    u = (-py).atan2(-px);
                } else {
                    u = py.atan2(px);
                }
                normalize_angle(&mut u);
                // V = sin(Sang) * (x cosU + y sinU - R) + z * cos(Sang)
                let v = ss * (px * u.cos() + py * u.sin() - self.prm1) + cs * pz;
                (u, v)
            }
            QuadricType::Torus => {
                let ploc = p - self.loc;
                let (x, y, z) = (ploc.x, ploc.y, ploc.z);
                let mut u = y.atan2(x);
                let r_maj = self.prm1;
                let r_min = self.prm2;
                if r_maj < r_min {
                    let cosu = u.cos();
                    let sinu = u.sin();
                    let z2 = z * z;
                    let min_r2 = r_min * r_min;
                    let r_cosu = r_maj * cosu;
                    let r_sinu = r_maj * sinu;
                    let xm = x - r_cosu;
                    let ym = y - r_sinu;
                    let xp = x + r_cosu;
                    let yp = y + r_sinu;
                    let d1 = xm * xm + ym * ym + z2 - min_r2;
                    let d2 = xp * xp + yp * yp + z2 - min_r2;
                    if d2.abs() < d1.abs() {
                        u += std::f64::consts::PI;
                    }
                }
                normalize_angle(&mut u);
                let cosu = u.cos();
                let sinu = u.sin();
                let dx = DVec3::new(cosu, sinu, 0.0);
                let dpv = DVec3::new(x - r_maj * cosu, y - r_maj * sinu, z);
                let mut v;
                let mag = dpv.length();
                if mag <= f64::MIN_POSITIVE {
                    v = 0.0;
                } else {
                    let dp = dpv / mag;
                    v = angle_with_ref(dx, dp, dx.cross(DVec3::Z));
                }
                normalize_angle(&mut v);
                (u, v)
            }
            QuadricType::Other => (0.0, 0.0),
        }
    }
}

impl Default for Quadric {
    fn default() -> Self {
        Self::new()
    }
}

// === OCCT IntSurf_QuadricTool (IntSurf_QuadricTool.hxx / .cxx) ===

impl Quadric {
    /// OCCT IntSurf_QuadricTool::Value(Quad, X, Y, Z).
    pub fn tool_value(&self, x: f64, y: f64, z: f64) -> f64 {
        self.distance(DVec3::new(x, y, z))
    }

    /// OCCT IntSurf_QuadricTool::Gradient(Quad, X, Y, Z).
    pub fn tool_gradient(&self, x: f64, y: f64, z: f64) -> DVec3 {
        self.gradient(DVec3::new(x, y, z))
    }

    /// OCCT IntSurf_QuadricTool::ValueAndGradient.
    pub fn tool_value_and_gradient(&self, x: f64, y: f64, z: f64) -> (f64, DVec3) {
        self.val_and_grad(DVec3::new(x, y, z))
    }

    /// OCCT IntSurf_QuadricTool::Tolerance (IntSurf_QuadricTool.cxx L17-29).
    pub fn tool_tolerance(&self) -> f64 {
        match self.typ {
            QuadricType::Sphere => 2.0e-6 * self.prm1,
            QuadricType::Cylinder => 2.0e-6 * self.prm1,
            _ => 1.0e-6,
        }
    }
}

// === OCCT support helpers ===

/// OCCT gp::Resolution() == RealSmall() == DBL_MIN.
const RESOLUTION: f64 = f64::MIN_POSITIVE;

/// OCCT Precision::Computational() == RealEpsilon() == DBL_EPSILON.
const COMPUTATIONAL: f64 = f64::EPSILON;

/// OCCT IntSurf_Quadric.cxx L30-34 — hasMagnitudeForNormalization.
fn has_magnitude_for_normalization(v: &DVec3) -> bool {
    v.length_squared() > RESOLUTION * RESOLUTION
}

/// OCCT gp_Lin::Distance(P) — distance from a point to the axis line.
fn line_distance(origin: &DVec3, dir: &DVec3, p: DVec3) -> f64 {
    let d = p - *origin;
    let along = d.dot(*dir);
    let perp = d - *dir * along;
    perp.length()
}

/// OCCT gp_Lin::Normal(P).Direction() — unit vector from the axis to P.
fn line_normal(origin: &DVec3, dir: &DVec3, p: DVec3) -> DVec3 {
    let d = p - *origin;
    let along = d.dot(*dir);
    let perp = d - *dir * along;
    perp.normalize_or_zero()
}

/// OCCT ElCLib::Parameter(gp_Lin, gp_Pnt) — parameter of a point on the line.
fn elclib_parameter(origin: &DVec3, dir: &DVec3, p: DVec3) -> f64 {
    (p - *origin).dot(*dir)
}

/// OCCT ElSLib normalizeAngle (ElSLib.cxx L33-58) — normalize to [0, 2*PI].
fn normalize_angle(a: &mut f64) {
    let pipi = std::f64::consts::TAU;
    let negative_resolution = -COMPUTATIONAL;
    while *a < negative_resolution {
        *a += pipi;
    }
    while *a > pipi * (1.0 + RESOLUTION) {
        *a -= pipi;
    }
    if *a < 0.0 {
        *a = 0.0;
    }
}

/// OCCT gp_Dir::AngleWithRef(V, VRef) — signed angle from V to Ref about their
/// cross product (used by ElSLib::TorusParameters).
fn angle_with_ref(v: DVec3, v_ref: DVec3, cross_axis: DVec3) -> f64 {
    let mut ang = v_ref.angle_between(v);
    if cross_axis.dot(v.cross(v_ref)) < 0.0 {
        ang = -ang;
    }
    ang
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::SurfaceEval;

    #[test]
    fn plane_value_matches_surface() {
        let p = Plane::new(DVec3::new(1.0, 2.0, 3.0), DVec3::new(0.0, 1.0, 0.0));
        let q = Quadric::from_plane(&p);
        for u in [-0.5, 0.0, 0.7] {
            for v in [-1.0, 0.0, 1.5] {
                let a = q.value(u, v);
                let b = p.point_at(u, v);
                assert!((a - b).length() < 1e-12, "plane value mismatch: {a:?} vs {b:?}");
            }
        }
        // distance of a point on the plane is ~0
        let on = p.point_at(0.3, -0.2);
        assert!(q.distance(on).abs() < 1e-12, "plane distance: {}", q.distance(on));
    }

    #[test]
    fn cylinder_value_matches_surface() {
        let c = CylindricalSurface {
            origin: DVec3::new(0.0, 0.0, 1.0),
            axis: DVec3::new(0.0, 0.0, 1.0),
            radius: 2.0,
            ref_dir: DVec3::X,
            y_dir: None,
        };
        let q = Quadric::from_cylinder(&c);
        for u in [0.0, 1.0, 2.5] {
            for v in [-2.0, 0.0, 3.0] {
                let a = q.value(u, v);
                let b = c.point_at(u, v);
                assert!((a - b).length() < 1e-12, "cylinder value mismatch: {a:?} vs {b:?}");
            }
        }
        // distance from a point on the surface to the quadric == 0
        let on = c.point_at(0.5, 1.0);
        assert!(q.distance(on).abs() < 1e-12, "cylinder distance: {}", q.distance(on));
    }

    #[test]
    fn cone_value_matches_surface() {
        let co = ConicalSurface::new(DVec3::new(0.0, 0.0, 0.0), DVec3::Z, 1.0, 0.4);
        let q = Quadric::from_cone(&co);
        for u in [0.0, 1.0, 3.0] {
            for v in [-1.0, 0.0, 2.0] {
                let a = q.value(u, v);
                let b = co.point_at(u, v);
                assert!((a - b).length() < 1e-12, "cone value mismatch: {a:?} vs {b:?}");
            }
        }
        let on = co.point_at(0.5, 1.0);
        assert!(q.distance(on).abs() < 1e-9, "cone distance: {}", q.distance(on));
    }

    #[test]
    fn torus_value_matches_surface() {
        let t = ToroidalSurface {
            center: DVec3::new(1.0, 0.0, 0.0),
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            major_radius: 3.0,
            minor_radius: 1.0,
        };
        let q = Quadric::from_torus(&t);
        for u in [0.0, 1.0, 4.0] {
            for v in [0.0, 1.5, 3.0] {
                let a = q.value(u, v);
                let b = t.point_at(u, v);
                assert!((a - b).length() < 1e-12, "torus value mismatch: {a:?} vs {b:?}");
            }
        }
        let on = t.point_at(0.5, 1.0);
        assert!(q.distance(on).abs() < 1e-9, "torus distance: {}", q.distance(on));
    }

    #[test]
    fn parameters_round_trip() {
        // cylinder
        let c = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.5,
            ref_dir: DVec3::X,
            y_dir: None,
        };
        let q = Quadric::from_cylinder(&c);
        let on = c.point_at(0.7, 2.0);
        let (u, v) = q.parameters(on);
        assert!((u - 0.7).abs() < 1e-9, "cylinder u: {u} vs 0.7");
        assert!((v - 2.0).abs() < 1e-9, "cylinder v: {v} vs 2.0");
    }
}
