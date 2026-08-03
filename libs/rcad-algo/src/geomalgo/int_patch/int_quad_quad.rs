// 1:1 translation of OCCT IntAna_IntQuadQuad (IntAna_IntQuadQuad.cxx),
// IntAna_Curve (IntAna_Curve.cxx), IntAna_Quadric (IntAna_Quadric.cxx) and
// math_TrigonometricFunctionRoots (math_TrigonometricFunctionRoots.cxx).
//
// This is the general quadric-quadric intersection used when
// IntAna_QuadQuadGeo returns IntAna_NoGeometricSolution (e.g. non-coaxial
// cylinder-sphere).  It parameterizes the explicit cylinder/cone by its
// (theta, axial) coordinates, forms the discriminant of the quadratic-in-Z,
// finds its trig-polynomial roots, and builds one IntAnaCurve per positive
// interval.  Each IntAnaCurve evaluates the curve point at a parameter by
// solving the quadratic in Z and mapping (U,V) through the surface frame.

use glam::DVec3;
use rcad_kernel::precision::PCONFUSION;

const TWO_PI: f64 = std::f64::consts::TAU;

// ============================================================================
// IntAna_Quadric — general quadric coefficients (IntAna_Quadric.hxx)
// ============================================================================

/// OCCT IntAna_Quadric: f(x,y,z) =
///   CXX x^2 + CYY y^2 + CZZ z^2
///   + 2 (CXY x y + CXZ x z + CYZ y z)
///   + 2 (CX x + CY y + CZ z)
///   + CCte
#[derive(Debug, Clone, Default)]
pub struct IntAnaQuadric {
    pub cxx: f64,
    pub cyy: f64,
    pub czz: f64,
    pub cxy: f64,
    pub cxz: f64,
    pub cyz: f64,
    pub cx: f64,
    pub cy: f64,
    pub cz: f64,
    pub ccte: f64,
    pub special_points: Vec<DVec3>,
}

impl IntAnaQuadric {
    pub fn new() -> Self {
        IntAnaQuadric {
            cxx: 0.0, cyy: 0.0, czz: 0.0, cxy: 0.0, cxz: 0.0, cyz: 0.0,
            cx: 0.0, cy: 0.0, cz: 0.0, ccte: 1.0, special_points: Vec::new(),
        }
    }

    /// OCCT IntAna_Quadric::SetQuadric(gp_Pln): P.Coefficients(A,B,C,D) with
    /// f = A x + B y + C z + D = 2(CX x + CY y + CZ z) + CCte.
    pub fn set_plane(&mut self, p: &rcad_kernel::geom::Plane) {
        let d = -p.normal.dot(p.origin);
        self.cx = p.normal.x * 0.5;
        self.cy = p.normal.y * 0.5;
        self.cz = p.normal.z * 0.5;
        self.ccte = d;
        self.cxx = 0.0; self.cyy = 0.0; self.czz = 0.0;
        self.cxy = 0.0; self.cxz = 0.0; self.cyz = 0.0;
    }

    /// OCCT IntAna_Quadric::SetQuadric(gp_Cylinder): Cyl.Coefficients.
    pub fn set_cylinder(&mut self, c: &rcad_kernel::geom::CylindricalSurface) {
        let coeffs = gp_cylinder_coefficients(c);
        self.cxx = coeffs.0; self.cyy = coeffs.1; self.czz = coeffs.2;
        self.cxy = coeffs.3; self.cxz = coeffs.4; self.cyz = coeffs.5;
        self.cx = coeffs.6; self.cy = coeffs.7; self.cz = coeffs.8;
        self.ccte = coeffs.9;
    }

    /// OCCT IntAna_Quadric::SetQuadric(gp_Sphere).
    pub fn set_sphere(&mut self, s: &rcad_kernel::geom::SphericalSurface) {
        let coeffs = gp_sphere_coefficients(s);
        self.cxx = coeffs.0; self.cyy = coeffs.1; self.czz = coeffs.2;
        self.cxy = coeffs.3; self.cxz = coeffs.4; self.cyz = coeffs.5;
        self.cx = coeffs.6; self.cy = coeffs.7; self.cz = coeffs.8;
        self.ccte = coeffs.9;
        // OCCT SetQuadric(gp_Sphere): special points = the two poles.
        let axis = s.axis.normalize_or_zero();
        self.special_points.push(s.center + axis * s.radius);
        self.special_points.push(s.center - axis * s.radius);
    }

    /// OCCT IntAna_Quadric::SetQuadric(gp_Cone): Cone.Coefficients + the
    /// apex as a special point (ElSLib::Value(0.0, -RefRadius/sin(Angle))).
    pub fn set_cone(&mut self, c: &rcad_kernel::geom::ConicalSurface) {
        let coeffs = gp_cone_coefficients(c);
        self.cxx = coeffs.0; self.cyy = coeffs.1; self.czz = coeffs.2;
        self.cxy = coeffs.3; self.cxz = coeffs.4; self.cyz = coeffs.5;
        self.cx = coeffs.6; self.cy = coeffs.7; self.cz = coeffs.8;
        self.ccte = coeffs.9;
        self.special_points.push(c.apex_point());
    }

    /// Build from a rcad Surface3 (OCCT constructors).  Returns None for
    /// non-quadric surfaces (BSpline etc.).
    pub fn from_surface3(surf: &rcad_kernel::geom::Surface3) -> Option<IntAnaQuadric> {
        use rcad_kernel::geom::Surface3;
        let mut q = IntAnaQuadric::new();
        match surf {
            Surface3::Plane(p) => q.set_plane(p),
            Surface3::Cylinder(c) => q.set_cylinder(c),
            Surface3::Sphere(s) => q.set_sphere(s),
            Surface3::Cone(c) => q.set_cone(c),
            _ => return None,
        }
        Some(q)
    }

    /// OCCT IntAna_Quadric::NewCoefficients(..., Axis): rewrite the quadric in
    /// the local frame `axis` (an orthonormal frame given by x/y/z dirs + loc).
    pub fn new_coefficients(&self, x_dir: DVec3, y_dir: DVec3, z_dir: DVec3, loc: DVec3) -> IntAnaQuadric {
        // The inverse of the local->global transform (global = loc + X*x_dir +
        // Y*y_dir + Z*z_dir) is the global->local one, whose matrix is the
        // transpose with translation.  OCCT builds gp_Trsf::SetTransformation
        // (absolute -> Axis) then Inverts it, giving the local->global map:
        //   x = t11 X + t12 Y + t13 Z + t14  (t14 = loc.x) etc.
        let t11 = x_dir.x; let t12 = y_dir.x; let t13 = z_dir.x; let t14 = loc.x;
        let t21 = x_dir.y; let t22 = y_dir.y; let t23 = z_dir.y; let t24 = loc.y;
        let t31 = x_dir.z; let t32 = y_dir.z; let t33 = z_dir.z; let t34 = loc.z;

        let cxx = self.cxx; let cyy = self.cyy; let czz = self.czz;
        let cxy = self.cxy; let cxz = self.cxz; let cyz = self.cyz;
        let cx = self.cx; let cy = self.cy; let cz = self.cz;
        let ccte = self.ccte;

        let t11_p2 = t11 * t11;
        let t21_p2 = t21 * t21;
        let t31_p2 = t31 * t31;
        let t12_p2 = t12 * t12;
        let t22_p2 = t22 * t22;
        let t32_p2 = t32 * t32;
        let t13_p2 = t13 * t13;
        let t23_p2 = t23 * t23;
        let t33_p2 = t33 * t33;
        let t14_p2 = t14 * t14;
        let t24_p2 = t24 * t24;
        let t34_p2 = t34 * t34;

        let n_ccte = ccte + t14_p2 * cxx + t24_p2 * cyy + t34_p2 * czz
            + 2.0 * (t14 * (cx + t24 * cxy + t34 * cxz)
                     + t24 * (cy + t34 * cyz) + t34 * cz);

        let n_cxx = t11_p2 * cxx + t21_p2 * cyy + t31_p2 * czz
            + 2.0 * (t11 * (t21 * cxy + t31 * cxz) + t21 * t31 * cyz);

        let n_cyy = t12_p2 * cxx + t22_p2 * cyy + t32_p2 * czz
            + 2.0 * (t12 * (t22 * cxy + t32 * cxz) + t22 * t32 * cyz);

        let n_czz = t13_p2 * cxx + t33_p2 * czz + t23_p2 * cyy
            + 2.0 * (t13 * (t23 * cxy + t33 * cxz) + t23 * t33 * cyz);

        let n_cz = t13 * cx + t13 * (t14 * cxx + t24 * cxy + t34 * cxz)
            + t14 * (t23 * cxy + t33 * cxz)
            + t23 * (cy + t24 * cyy + t34 * cyz)
            + t33 * (t24 * cyz + cz + t34 * czz);

        let n_cx = t11 * (cx + t14 * cxx + t24 * cxy + t34 * cxz)
            + t14 * (t21 * cxy + t31 * cxz)
            + t21 * (cy + t24 * cyy + t34 * cyz)
            + t31 * (t24 * cyz + cz + t34 * czz);

        let n_cxy = t11 * (t12 * cxx + t22 * cxy + t32 * cxz)
            + t12 * (t21 * cxy + t31 * cxz)
            + t21 * (t22 * cyy + t32 * cyz)
            + t31 * (t22 * cyz + t32 * czz);

        let n_cxz = t11 * (t13 * cxx + t23 * cxy + t33 * cxz)
            + t13 * (t21 * cxy + t31 * cxz)
            + t21 * (t23 * cyy + t33 * cyz)
            + t31 * (t23 * cyz + t33 * czz);

        let n_cy = t12 * (cx + t14 * cxx + t24 * cxy + t34 * cxz)
            + t14 * (t22 * cxy + t32 * cxz)
            + t22 * (cy + t24 * cyy + t34 * cyz)
            + t32 * (cz + t24 * cyz + t34 * czz);

        let n_cyz = t12 * (t13 * cxx + t23 * cxy + t33 * cxz)
            + t13 * (t22 * cxy + t32 * cxz)
            + t22 * (t23 * cyy + t33 * cyz)
            + t32 * (t23 * cyz + t33 * czz);

        IntAnaQuadric {
            cxx: n_cxx, cyy: n_cyy, czz: n_czz,
            cxy: n_cxy, cxz: n_cxz, cyz: n_cyz,
            cx: n_cx, cy: n_cy, cz: n_cz, ccte: n_ccte,
            special_points: Vec::new(),
        }
    }
}

/// OCCT gp_Cylinder::Coefficients (gp_Cylinder.cxx L27-61).  The local frame
/// transform maps global -> local; rows 1-2 (Xdir, Ydir) with translation.
fn gp_cylinder_coefficients(c: &rcad_kernel::geom::CylindricalSurface) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    let origin = c.origin;
    let z = c.axis.normalize_or_zero();
    let x = c.ref_dir.normalize_or_zero();
    let y = z.cross(x).normalize_or_zero();
    let r = c.radius;
    let t11 = x.x; let t12 = x.y; let t13 = x.z; let t14 = -(origin.dot(x));
    let t21 = y.x; let t22 = y.y; let t23 = y.z; let t24 = -(origin.dot(y));
    let a1 = t11 * t11 + t21 * t21;
    let a2 = t12 * t12 + t22 * t22;
    let a3 = t13 * t13 + t23 * t23;
    let b1 = t11 * t12 + t21 * t22;
    let b2 = t11 * t13 + t21 * t23;
    let b3 = t12 * t13 + t22 * t23;
    let c1 = t11 * t14 + t21 * t24;
    let c2 = t12 * t14 + t22 * t24;
    let c3 = t13 * t14 + t23 * t24;
    let d = t14 * t14 + t24 * t24 - r * r;
    (a1, a2, a3, b1, b2, b3, c1, c2, c3, d)
}

/// OCCT gp_Sphere::Coefficients (gp_Sphere.cxx L23-60).
fn gp_sphere_coefficients(s: &rcad_kernel::geom::SphericalSurface) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    let origin = s.center;
    let z = s.axis.normalize_or_zero();
    let x = s.ref_dir.normalize_or_zero();
    let y = z.cross(x).normalize_or_zero();
    let r = s.radius;
    let t11 = x.x; let t12 = x.y; let t13 = x.z; let t14 = -(origin.dot(x));
    let t21 = y.x; let t22 = y.y; let t23 = y.z; let t24 = -(origin.dot(y));
    let t31 = z.x; let t32 = z.y; let t33 = z.z; let t34 = -(origin.dot(z));
    let a1 = t11 * t11 + t21 * t21 + t31 * t31;
    let a2 = t12 * t12 + t22 * t22 + t32 * t32;
    let a3 = t13 * t13 + t23 * t23 + t33 * t33;
    let b1 = t11 * t12 + t21 * t22 + t31 * t32;
    let b2 = t11 * t13 + t21 * t23 + t31 * t33;
    let b3 = t12 * t13 + t22 * t23 + t32 * t33;
    let c1 = t11 * t14 + t21 * t24 + t31 * t34;
    let c2 = t12 * t14 + t22 * t24 + t32 * t34;
    let c3 = t13 * t14 + t23 * t24 + t33 * t34;
    let d = t14 * t14 + t24 * t24 + t34 * t34 - r * r;
    (a1, a2, a3, b1, b2, b3, c1, c2, c3, d)
}

/// OCCT gp_Cone::Coefficients (gp_Cone.cxx L26-64).  Local equation:
/// X^2 + Y^2 - (radius + Z*tan(semiAngle))^2 = 0.
fn gp_cone_coefficients(c: &rcad_kernel::geom::ConicalSurface) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    // rcad ConicalSurface: the reference circle of radius `radius` is at the
    // `apex` point; the geometric apex is apex - axis*(radius/tan).  The frame
    // for the Coefficients is the natural cone frame at the reference point.
    let origin = c.apex;
    let z = c.axis.normalize_or_zero();
    let x = any_perpendicular_axis(z);
    let y = z.cross(x).normalize_or_zero();
    let r = c.radius;
    let k_ang = c.half_angle_rad.tan();
    let t11 = x.x; let t12 = x.y; let t13 = x.z; let t14 = -(origin.dot(x));
    let t21 = y.x; let t22 = y.y; let t23 = y.z; let t24 = -(origin.dot(y));
    let t31 = z.x * k_ang; let t32 = z.y * k_ang; let t33 = z.z * k_ang; let t34 = -(origin.dot(z)) * k_ang;
    let a1 = t11 * t11 + t21 * t21 - t31 * t31;
    let a2 = t12 * t12 + t22 * t22 - t32 * t32;
    let a3 = t13 * t13 + t23 * t23 - t33 * t33;
    let b1 = t11 * t12 + t21 * t22 - t31 * t32;
    let b2 = t11 * t13 + t21 * t23 - t31 * t33;
    let b3 = t12 * t13 + t22 * t23 - t32 * t33;
    let c1 = t11 * t14 + t21 * t24 - t31 * (r + t34);
    let c2 = t12 * t14 + t22 * t24 - t32 * (r + t34);
    let c3 = t13 * t14 + t23 * t24 - t33 * (r + t34);
    let d = t14 * t14 + t24 * t24 - r * r - t34 * t34 - 2.0 * r * t34;
    (a1, a2, a3, b1, b2, b3, c1, c2, c3, d)
}

/// A unit vector perpendicular to `v` (any).  Mirrors OCCT's any_perpendicular
/// where the canonical choice is the coordinate axis with the smallest
/// component — but any orthonormal X works because the quadric coefficients and
/// the curve evaluation use the same frame.
pub fn any_perpendicular_axis(v: DVec3) -> DVec3 {
    let avx = v.x.abs();
    let avy = v.y.abs();
    let avz = v.z.abs();
    let basis = if avx <= avy && avx <= avz {
        DVec3::X
    } else if avy <= avz {
        DVec3::Y
    } else {
        DVec3::Z
    };
    basis.cross(v).normalize_or_zero()
}

// ============================================================================
// math_TrigonometricFunctionRoots — roots of
//   a cos^2(x) + 2 b cos(x) sin(x) + c cos(x) + d sin(x) + e
// ============================================================================

/// OCCT math_TrigonometricFunctionRoots.  `Perform` translates the trig
/// polynomial into a quartic in t = tan(x/2) and solves it with the direct
/// polynomial solver, then filters to [InfBound, SupBound].
pub struct MathTrigFunctionRoots {
    nb_sol: i32,
    sol: [f64; 4],
    infinite_status: bool,
    done: bool,
}

impl MathTrigFunctionRoots {
    pub fn new(a: f64, b: f64, c: f64, d: f64, e: f64, inf_bound: f64, sup_bound: f64) -> Self {
        let mut r = MathTrigFunctionRoots { nb_sol: -1, sol: [0.0; 4], infinite_status: false, done: false };
        r.perform(a, b, c, d, e, inf_bound, sup_bound);
        r
    }

    pub fn is_done(&self) -> bool { self.done }
    pub fn infinite_roots(&self) -> bool { self.infinite_status }
    pub fn nb_solutions(&self) -> usize { self.nb_sol.max(0) as usize }
    pub fn value(&self, index: usize) -> f64 { self.sol[index] }

    fn perform(&mut self, a: f64, b: f64, c: f64, d: f64, e: f64, inf_bound: f64, sup_bound: f64) {
        let eps = 1.5e-12;
        let depi = TWO_PI;
        let (my_borne_inf, delta, modv) = if inf_bound <= f64::NEG_INFINITY && sup_bound >= f64::INFINITY {
            (0.0, depi, 0.0)
        } else if sup_bound >= f64::INFINITY {
            (inf_bound, depi, inf_bound / depi)
        } else if inf_bound <= f64::NEG_INFINITY {
            (sup_bound - depi, depi, (sup_bound - depi) / depi)
        } else {
            let mut delta = sup_bound - inf_bound;
            if delta > depi { delta = depi; }
            (inf_bound, delta, inf_bound / depi)
        };
        let modv = modv;

        self.infinite_status = false;
        self.done = true;
        let mut nzer = 0usize;
        let mut zer = [0.0f64; 4];

        if a.abs() <= eps && b.abs() <= eps {
            if c.abs() <= eps {
                if d.abs() <= eps {
                    if e.abs() <= eps {
                        self.infinite_status = true;
                        return;
                    } else {
                        self.nb_sol = 0;
                        return;
                    }
                } else {
                    // d*sin(x) + e = 0
                    self.nb_sol = 0;
                    let aa = -e / d;
                    if aa.abs() > 1.0 { return; }
                    zer[0] = aa.asin();
                    zer[1] = std::f64::consts::PI - zer[0];
                    nzer = 2;
                    for i in 0..nzer {
                        if zer[i] <= -eps {
                            zer[i] = depi - zer[i].abs();
                        }
                        zer[i] += modv.trunc() * depi;
                        let x = zer[i] - my_borne_inf;
                        if x > -delta_eps(delta) && x < delta + delta_eps(delta) {
                            self.nb_sol += 1;
                            self.sol[(self.nb_sol - 1) as usize] = zer[i];
                        }
                    }
                    return;
                }
            } else if d.abs() <= eps {
                // c*cos(x) + e = 0
                self.nb_sol = 0;
                let aa = -e / c;
                if aa.abs() > 1.0 { return; }
                zer[0] = aa.acos();
                zer[1] = -zer[0];
                nzer = 2;
                for i in 0..nzer {
                    if zer[i] <= -eps {
                        zer[i] = depi - zer[i].abs();
                    }
                    zer[i] += modv.trunc() * TWO_PI;
                    let x = zer[i] - my_borne_inf;
                    if x >= -delta_eps(delta) && x <= delta + delta_eps(delta) {
                        self.nb_sol += 1;
                        self.sol[(self.nb_sol - 1) as usize] = zer[i];
                    }
                }
                return;
            } else {
                // quadratic: (E-C) + 2D t + (E+C) t^2 for t=tan(x/2)  — the OCCT
                // path solves the cos/sin quadratic via the direct polynomial
                // in t.  Equivalent to solving c*cos + d*sin + e = 0 directly.
                let aa = e - c;
                let bb = 2.0 * d;
                let cc = e + c;
                let roots = solve_poly2(aa, bb, cc);
                if roots.0 == -2 {
                    self.done = false;
                    return;
                } else if roots.0 == -1 {
                    self.infinite_status = true;
                    return;
                } else {
                    nzer = roots.0 as usize;
                    for i in 0..nzer { zer[i] = roots.1[i]; }
                }
            }
        } else {
            // Two additional analytical cases.
            if a.abs() <= eps && e.abs() <= eps {
                if c.abs() <= eps {
                    // 2 B sin cos + D sin = 0
                    nzer = 2;
                    zer[0] = 0.0;
                    zer[1] = std::f64::consts::PI;
                    let aa = -d / (b * 2.0);
                    if aa.abs() <= 1.0 + PCONFUSION {
                        nzer = 4;
                        if aa >= 1.0 {
                            zer[2] = 0.0; zer[3] = 0.0;
                        } else if aa <= -1.0 {
                            zer[2] = std::f64::consts::PI; zer[3] = std::f64::consts::PI;
                        } else {
                            zer[2] = aa.acos();
                            zer[3] = depi - zer[2];
                        }
                    }
                    self.nb_sol = 0;
                    for i in 0..nzer {
                        if zer[i] <= my_borne_inf - eps { zer[i] += depi; }
                        zer[i] += modv.trunc() * TWO_PI;
                        let x = zer[i] - my_borne_inf;
                        if x >= -PCONFUSION && x <= delta + PCONFUSION {
                            if zer[i] < inf_bound { zer[i] = inf_bound; }
                            if zer[i] > sup_bound { zer[i] = sup_bound; }
                            self.nb_sol += 1;
                            self.sol[(self.nb_sol - 1) as usize] = zer[i];
                        }
                    }
                    return;
                }
                if d.abs() <= eps {
                    // 2 B sin cos + C cos = 0
                    nzer = 2;
                    zer[0] = std::f64::consts::FRAC_PI_2;
                    zer[1] = std::f64::consts::PI * 3.0 / 2.0;
                    let aa = -c / (b * 2.0);
                    if aa.abs() <= 1.0 + PCONFUSION {
                        nzer = 4;
                        if aa >= 1.0 {
                            zer[2] = std::f64::consts::FRAC_PI_2; zer[3] = std::f64::consts::FRAC_PI_2;
                        } else if aa <= -1.0 {
                            zer[2] = std::f64::consts::PI * 3.0 / 2.0; zer[3] = std::f64::consts::PI * 3.0 / 2.0;
                        } else {
                            zer[2] = aa.asin();
                            zer[3] = std::f64::consts::PI - zer[2];
                        }
                    }
                    self.nb_sol = 0;
                    for i in 0..nzer {
                        if zer[i] <= my_borne_inf - eps { zer[i] += depi; }
                        zer[i] += modv.trunc() * TWO_PI;
                        let x = zer[i] - my_borne_inf;
                        if x >= -PCONFUSION && x <= delta + PCONFUSION {
                            if zer[i] < inf_bound { zer[i] = inf_bound; }
                            if zer[i] > sup_bound { zer[i] = sup_bound; }
                            self.nb_sol += 1;
                            self.sol[(self.nb_sol - 1) as usize] = zer[i];
                        }
                    }
                    return;
                }
            }

            // Equation of the 4th degree.  OCCT maps the trig polynomial to a
            // quartic in t = tan(x/2):
            //   ko(1) = A - C + E, ko(2) = 2D - 4B, ko(3) = 2E - 2A,
            //   ko(4) = 4B + 2D,    ko(5) = A + C + E
            let mut ko = [a - c + e, 2.0 * d - 4.0 * b, 2.0 * e - 2.0 * a, 4.0 * b + 2.0 * d, a + c + e];
            let mut bko;
            let mut nzer4 = 0usize;
            loop {
                bko = false;
                let roots = solve_poly4(ko[0], ko[1], ko[2], ko[3], ko[4]);
                if roots.0 == -2 {
                    self.done = false;
                    return;
                } else if roots.0 == -1 {
                    self.infinite_status = true;
                    return;
                } else {
                    nzer4 = roots.0 as usize;
                    for i in 0..nzer4 { zer[i] = roots.1[i]; }
                }
                // sort
                for i in 0..nzer4 {
                    for j in (i + 1)..nzer4 {
                        if zer[j] < zer[i] {
                            zer.swap(i, j);
                        }
                    }
                }
                // detect numerically-double roots that are not true doubles
                for i in 0..nzer4.saturating_sub(1) {
                    if (zer[i + 1] - zer[i]).abs() < eps {
                        let qw = zer[i + 1];
                        let va = ko[3] + qw * (2.0 * ko[2] + qw * (3.0 * ko[1] + qw * (4.0 * ko[0])));
                        if va.abs() > eps {
                            bko = true;
                            break;
                        }
                    }
                }
                if bko {
                    for k in ko.iter_mut() { *k *= 0.0001; }
                } else {
                    break;
                }
            }
            nzer = nzer4;
        }

        // Verification of the solutions w.r.t. the bounds.
        let sup_minus_inf_100 = (sup_bound - inf_bound) * 0.01;
        self.nb_sol = 0;
        let tol1 = 1e-15;
        for i in 0..nzer {
            let mut teta = zer[i].atan();
            teta += teta;
            if zer[i] <= -eps {
                teta = depi - teta.abs();
            }
            teta += modv.trunc() * depi;
            if teta - my_borne_inf < 0.0 {
                teta += depi;
            }
            let x = teta - my_borne_inf;
            if x >= -delta_eps(delta) && x <= delta + delta_eps(delta) {
                let mut t = teta;
                // Newton refinement of the trig equation.
                let mut nit = 0;
                while nit < 10 {
                    let (f, df) = trig_value_deriv(a, b, c, d, e, t);
                    if df.abs() < 1e-15 { break; }
                    let dt = f / df;
                    t -= dt;
                    if dt.abs() < tol1 { break; }
                    nit += 1;
                }
                let delta_newton = t - teta;
                if delta_newton <= sup_minus_inf_100 && delta_newton >= -sup_minus_inf_100 {
                    teta = t;
                }
                // insert teta in increasing order
                let mut flag4 = false;
                for k in 0..self.nb_sol.max(0) as usize {
                    if teta < self.sol[k] {
                        for l in (k..self.nb_sol.max(0) as usize).rev() {
                            self.sol[l + 1] = self.sol[l];
                        }
                        self.sol[k] = teta;
                        self.nb_sol += 1;
                        flag4 = true;
                        break;
                    }
                }
                if !flag4 {
                    self.nb_sol += 1;
                    self.sol[(self.nb_sol - 1) as usize] = teta;
                }
            }
        }
        // Special case of PI.
        if self.nb_sol < 4 {
            let start_index = self.nb_sol.max(0) + 1;
            for sol_it in start_index..=4 {
                let mut teta = std::f64::consts::PI + modv.trunc() * TWO_PI;
                let x = teta - my_borne_inf;
                if x >= -delta_eps(delta) && x <= delta + delta_eps(delta) {
                    if (a - c + e).abs() <= eps {
                        let mut flag4 = false;
                        let mut j = 0usize;
                        for k in 0..self.nb_sol.max(0) as usize {
                            j = k;
                            if teta < self.sol[k] {
                                flag4 = true;
                                break;
                            }
                            if sol_it == start_index && (teta - self.sol[k]).abs() <= eps {
                                return;
                            }
                        }
                        if !flag4 {
                            self.nb_sol += 1;
                            self.sol[(self.nb_sol - 1) as usize] = teta;
                        } else {
                            for k in (j..self.nb_sol.max(0) as usize).rev() {
                                self.sol[k + 1] = self.sol[k];
                            }
                            self.sol[j] = teta;
                            self.nb_sol += 1;
                        }
                    }
                }
            }
        }
    }
}

fn delta_eps(delta: f64) -> f64 {
    // Epsilon(Delta) in OCCT — the relative machine epsilon scaled by |Delta|.
    delta.abs() * f64::EPSILON
}

fn trig_value_deriv(a: f64, b: f64, c: f64, d: f64, e: f64, u: f64) -> (f64, f64) {
    let (s, co) = u.sin_cos();
    let f = a * co * co + 2.0 * b * co * s + c * co + d * s + e;
    let df = 2.0 * ((s * co) * (a - a) + d * co - c * s + b * (co * co - s * s));
    // Correct derivative of a cos^2 + 2b cos sin + c cos + d sin:
    let df = -2.0 * a * co * s + 2.0 * b * (co * co - s * s) - c * s + d * co;
    (f, df)
}

/// OCCT math_DirectPolynomialRoots for a quadratic.  Returns (-2, _) = not
/// done, (-1, _) = infinite roots, else (n, roots).
fn solve_poly2(a: f64, b: f64, c: f64) -> (i32, [f64; 4]) {
    let mut out = [0.0f64; 4];
    if a.abs() < 1e-30 && b.abs() < 1e-30 && c.abs() < 1e-30 {
        return (-1, out);
    }
    if a.abs() < 1e-30 {
        // linear b t + c = 0
        if b.abs() < 1e-30 {
            return (0, out);
        }
        out[0] = -c / b;
        return (1, out);
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return (0, out);
    }
    let sq = disc.sqrt();
    out[0] = (-b + sq) / (2.0 * a);
    out[1] = (-b - sq) / (2.0 * a);
    let n = if disc.abs() < 1e-30 { 1 } else { 2 };
    (n, out)
}

/// OCCT math_DirectPolynomialRoots for a quartic.  Returns (-2, _) = not done,
/// (-1, _) = infinite, else (n, roots).
fn solve_poly4(a: f64, b: f64, c: f64, d: f64, e: f64) -> (i32, [f64; 4]) {
    let mut out = [0.0f64; 4];
    let eps = 1e-30;
    if a.abs() < eps && b.abs() < eps && c.abs() < eps && d.abs() < eps && e.abs() < eps {
        return (-1, out);
    }
    if a.abs() < eps {
        return solve_poly3(b, c, d, e);
    }
    let roots = rcad_kernel::math::math_poly::solve_quartic(a, b, c, d, e);
    let n = roots.len().min(4) as i32;
    for (i, r) in roots.iter().enumerate().take(n as usize) {
        out[i] = *r;
    }
    (n, out)
}

fn solve_poly3(a: f64, b: f64, c: f64, d: f64) -> (i32, [f64; 4]) {
    let mut out = [0.0f64; 4];
    if a.abs() < 1e-30 {
        return solve_poly2(b, c, d);
    }
    let roots = rcad_kernel::math::math_poly::solve_cubic(a, b, c, d);
    let n = roots.len().min(4) as i32;
    for (i, r) in roots.iter().enumerate().take(n as usize) {
        out[i] = *r;
    }
    (n, out)
}

// ============================================================================
// TrigonometricRoots — internal OCCT class wrapping MathTrigFunctionRoots with
// a residual check and sorting.
// ============================================================================

struct TrigRoots {
    roots: [f64; 4],
    done: bool,
    nb_roots: usize,
    infinite_roots: bool,
}

impl TrigRoots {
    /// OCCT TrigonometricRoots constructor (L199-285).  F = AA*CN^2 +
    /// 2*BB*CN*SN + CC*CN + DD*SN + EE.
    fn new(cc: f64, sc: f64, c: f64, s: f64, cte: f64, binf: f64, bsup: f64) -> Self {
        let pippi = TWO_PI;
        let mut tr = TrigRoots { roots: [0.0; 4], done: false, nb_roots: 0, infinite_roots: false };
        let mtfr = MathTrigFunctionRoots::new(cc, sc, c, s, cte, binf, bsup);
        if !mtfr.is_done() {
            return tr;
        }
        tr.done = true;
        if mtfr.infinite_roots() {
            tr.infinite_roots = true;
            return tr;
        }
        tr.nb_roots = mtfr.nb_solutions();
        for i in 0..tr.nb_roots {
            tr.roots[i] = mtfr.value(i);
            if tr.roots[i] < 0.0 {
                tr.roots[i] += pippi;
            }
            if tr.roots[i] > pippi {
                tr.roots[i] -= pippi;
            }
        }
        // The direct search gives unreliable results — check each root against
        // the original polynomial.
        let sv_nb_roots = tr.nb_roots;
        for i in 0..sv_nb_roots {
            let co = tr.roots[i].cos();
            let si = tr.roots[i].sin();
            let y = co * (cc * co + (sc + sc) * si + c) + s * si + cte;
            if y.abs() > 1e-8 {
                tr.done = false;
                return tr;
            }
        }
        // bubble sort
        for i in 0..sv_nb_roots {
            for j in (i + 1)..sv_nb_roots {
                if tr.roots[j] < tr.roots[i] {
                    tr.roots.swap(i, j);
                }
            }
        }
        tr.infinite_roots = false;
        if tr.nb_roots == 0 {
            if (cc.abs() + sc.abs() + c.abs() + s.abs()) < 1e-10 {
                if cte.abs() < 1e-10 {
                    tr.infinite_roots = true;
                }
            }
        }
        tr
    }

    fn is_done(&self) -> bool { self.done }
    fn infinite_roots(&self) -> bool { self.infinite_roots }
    fn nb_solutions(&self) -> usize { self.nb_roots }
    fn value(&self, i: usize) -> f64 { self.roots[i] }
}

/// OCCT MyTrigonometricFunction (L292-339).
struct MyTrigFunction {
    cc: f64,
    ss: f64,
    sc: f64,
    s: f64,
    c: f64,
    cte: f64,
}

impl MyTrigFunction {
    fn new(cc: f64, ss: f64, sc: f64, c: f64, s: f64, cte: f64) -> Self {
        MyTrigFunction { cc, ss, sc, s, c, cte }
    }
    fn value(&self, u: f64) -> f64 {
        let (s, co) = u.sin_cos();
        self.cc * co * co + self.ss * s * s + 2.0 * (s * (self.sc * co + self.s) + co * self.c) + self.cte
    }
}

// ============================================================================
// IntAna_Curve — parametric curve result of quadric-quadric intersection.
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum IntAnaCurveType {
    Cylinder,
    Cone,
}

/// OCCT IntAna_Curve: a curve defined on [DomainInf, DomainSup] whose point at
/// parameter Theta solves A(Theta) V^2 + B(Theta) V + C(Theta) = 0, where
/// A/B/C are the quadric's coefficients restricted to the surface parametrization.
#[derive(Debug, Clone)]
pub struct IntAnaCurve {
    // Z0/Z1/Z2 trig-polynomial coefficients (A = Z2, B = Z1, C = Z0 in V).
    z0_cte: f64, z0_sin: f64, z0_cos: f64, z0_sin_sin: f64, z0_cos_cos: f64, z0_cos_sin: f64,
    z1_cte: f64, z1_sin: f64, z1_cos: f64, z1_sin_sin: f64, z1_cos_cos: f64, z1_cos_sin: f64,
    z2_cte: f64, z2_sin: f64, z2_cos: f64, z2_sin_sin: f64, z2_cos_cos: f64, z2_cos_sin: f64,
    two_curves: bool,
    take_z_positive: bool,
    tolerance: f64,
    domain_inf: f64,
    domain_sup: f64,
    my_first_parameter: f64,
    my_last_parameter: f64,
    // Surface definition.
    surf_type: IntAnaCurveType,
    // cylinder: origin/axis/x/y, radius
    cyl_origin: DVec3,
    cyl_axis: DVec3,
    cyl_x: DVec3,
    cyl_y: DVec3,
    cyl_r: f64,
    // cone: apex (reference point), axis, x, y, radius, angle
    con_apex: DVec3,
    con_axis: DVec3,
    con_x: DVec3,
    con_y: DVec3,
    con_r: f64,
    con_angle: f64,
    // OCCT IntPatch_ALine vertices (IntPatch_Point): boundary points on the
    // line, added by ProcessBounds (IntCySp/IntCyCo/IntCoCo/IntCoSp).  They are
    // consumed by IntPatch_ALineToWLine::MakeWLine.
    pub vertices: Vec<crate::geomalgo::int_patch::special_points::PatchPoint>,
}

impl IntAnaCurve {
    /// OCCT IntAna_Curve::SetCylinderQuadValues (L161-217).
    fn set_cylinder_quad_values(
        cyl: &rcad_kernel::geom::CylindricalSurface,
        qxx: f64, qyy: f64, qzz: f64, qxy: f64, qxz: f64, qyz: f64,
        qx: f64, qy: f64, qz: f64, q1: f64,
        tol: f64, dom_inf: f64, dom_sup: f64,
        two_curves: bool, take_z_positive: bool,
    ) -> IntAnaCurve {
        let r = cyl.radius;
        let origin = cyl.origin;
        let axis = cyl.axis.normalize_or_zero();
        let x = cyl.ref_dir.normalize_or_zero();
        let y = axis.cross(x).normalize_or_zero();
        let r_mul2 = r + r;
        IntAnaCurve {
            z0_cte: q1,
            z0_sin: r_mul2 * qy,
            z0_cos: r_mul2 * qx,
            z0_cos_cos: qxx * r * r,
            z0_sin_sin: qyy * r * r,
            z0_cos_sin: r * r * qxy,
            z1_cte: qz + qz,
            z1_sin: r_mul2 * qyz,
            z1_cos: r_mul2 * qxz,
            z1_sin_sin: 0.0, z1_cos_cos: 0.0, z1_cos_sin: 0.0,
            z2_cte: qzz,
            z2_sin: 0.0, z2_cos: 0.0, z2_sin_sin: 0.0, z2_cos_cos: 0.0, z2_cos_sin: 0.0,
            two_curves, take_z_positive, tolerance: tol,
            domain_inf: dom_inf, domain_sup: dom_sup,
            my_first_parameter: dom_inf,
            my_last_parameter: if two_curves { dom_sup + dom_sup - dom_inf } else { dom_sup },
            surf_type: IntAnaCurveType::Cylinder,
            cyl_origin: origin, cyl_axis: axis, cyl_x: x, cyl_y: y, cyl_r: r,
            con_apex: DVec3::ZERO, con_axis: DVec3::Z, con_x: DVec3::X, con_y: DVec3::Y,
            con_r: 0.0, con_angle: 0.0,
            vertices: Vec::new(),
        }
    }

    /// OCCT IntAna_Curve::SetConeQuadValues (L95-155).
    fn set_cone_quad_values(
        con: &rcad_kernel::geom::ConicalSurface,
        qxx: f64, qyy: f64, qzz: f64, qxy: f64, qxz: f64, qyz: f64,
        qx: f64, qy: f64, qz: f64, q1: f64,
        tol: f64, dom_inf: f64, dom_sup: f64,
        two_curves: bool, take_z_positive: bool,
    ) -> IntAnaCurve {
        let r = con.radius;
        let apex = con.apex;
        let axis = con.axis.normalize_or_zero();
        let x = any_perpendicular_axis(axis);
        let y = axis.cross(x).normalize_or_zero();
        let angle = con.half_angle_rad;
        let un_sur_tg_angle = 1.0 / angle.tan();
        IntAnaCurve {
            z0_cte: q1,
            z0_sin: 0.0, z0_cos: 0.0, z0_cos_cos: 0.0, z0_sin_sin: 0.0, z0_cos_sin: 0.0,
            z1_cte: 2.0 * un_sur_tg_angle * qz,
            z1_sin: qy + qy,
            z1_cos: qx + qx,
            z1_sin_sin: 0.0, z1_cos_cos: 0.0, z1_cos_sin: 0.0,
            z2_cte: qzz * un_sur_tg_angle * un_sur_tg_angle,
            z2_sin: (un_sur_tg_angle + un_sur_tg_angle) * qyz,
            z2_cos: (un_sur_tg_angle + un_sur_tg_angle) * qxz,
            z2_cos_cos: qxx,
            z2_sin_sin: qyy,
            z2_cos_sin: qxy,
            two_curves, take_z_positive, tolerance: tol,
            domain_inf: dom_inf, domain_sup: dom_sup,
            my_first_parameter: dom_inf,
            my_last_parameter: if two_curves { dom_sup + dom_sup - dom_inf } else { dom_sup },
            surf_type: IntAnaCurveType::Cone,
            cyl_origin: DVec3::ZERO, cyl_axis: DVec3::Z, cyl_x: DVec3::X, cyl_y: DVec3::Y, cyl_r: 0.0,
            con_apex: apex, con_axis: axis, con_x: x, con_y: y, con_r: r, con_angle: angle,
            vertices: Vec::new(),
        }
    }

    pub fn domain(&self) -> [f64; 2] {
        [self.my_first_parameter, self.my_last_parameter]
    }

    /// OCCT IntAna_Curve::SetDomain (IntAna_Curve.cxx) — restrict the curve to
    /// [theFirst, theLast].
    pub fn set_domain(&mut self, the_first: f64, the_last: f64) {
        self.my_first_parameter = the_first;
        self.my_last_parameter = the_last;
    }

    /// OCCT IntAna_Curve::IsFirstOpen (IntAna_Curve.hxx L89) — the domain is
    /// bounded (firstbounded=true) by default in the IntXX flow.
    pub fn is_first_open(&self) -> bool {
        !self.my_first_parameter.is_finite()
    }

    /// OCCT IntAna_Curve::IsLastOpen (IntAna_Curve.hxx L92).
    pub fn is_last_open(&self) -> bool {
        !self.my_last_parameter.is_finite()
    }

    // OCCT IntPatch_ALine vertex access (the ALine wraps this IntAna_Curve and
    // carries IntPatch_Point vertices added by ProcessBounds).
    pub fn has_vertices(&self) -> bool {
        !self.vertices.is_empty()
    }
    pub fn vertex_params(&self) -> Vec<f64> {
        self.vertices.iter().map(|v| v.param_on_line).collect()
    }
    pub fn vertex_at(&self, i: usize) -> crate::geomalgo::int_patch::special_points::PatchPoint {
        self.vertices[i].clone()
    }

    /// OCCT IntAna_Curve::InternalUVValue (L279-374).
    fn internal_uv_value(&self, mut theta: f64) -> Option<(f64, f64, f64, f64, f64, f64)> {
        let rel_tolp = 1.0 + f64::EPSILON;
        let rel_tolm = 1.0 - f64::EPSILON;
        let a_dt = 100.0 * f64::EPSILON * (self.domain_sup + self.domain_sup - self.domain_inf).abs();

        let mut second_solution = false;
        if (theta < self.domain_inf * rel_tolm)
            || ((theta > self.domain_sup * rel_tolp) && !self.two_curves)
            || (theta > (self.domain_sup + self.domain_sup - self.domain_inf) * rel_tolp)
        {
            return None;
        }
        if (theta - self.domain_sup).abs() < a_dt {
            theta = self.domain_sup;
        } else if theta > self.domain_sup {
            theta = self.domain_sup + self.domain_sup - theta;
            second_solution = true;
        }
        let u = theta;
        if !self.two_curves {
            second_solution = self.take_z_positive;
        }
        let co = theta.cos();
        let si = theta.sin();
        let a_sin2t = (theta + theta).sin();
        let a_cos2t = (theta + theta).cos();

        let a = self.z2_cte + si * (self.z2_sin + si * self.z2_sin_sin)
            + co * (self.z2_cos + co * self.z2_cos_cos) + self.z2_cos_sin * a_sin2t;
        let a_da = co * self.z2_sin - si * self.z2_cos
            + a_sin2t * (self.z2_sin_sin - self.z2_cos_cos)
            + a_cos2t * (self.z2_cos_sin * self.z2_cos_sin);

        let b = self.z1_cte + si * (self.z1_sin + si * self.z1_sin_sin)
            + co * (self.z1_cos + co * self.z1_cos_cos) + self.z1_cos_sin * a_sin2t;
        let a_db = self.z1_sin * co - self.z1_cos * si
            + a_sin2t * (self.z1_sin_sin - self.z1_cos_cos)
            + a_cos2t * (self.z1_cos_sin + self.z1_cos_sin);

        let c = self.z0_cte + si * (self.z0_sin + si * self.z0_sin_sin)
            + co * (self.z0_cos + co * self.z0_cos_cos) + self.z0_cos_sin * a_sin2t;
        let a_dc = self.z0_sin * co - self.z0_cos * si
            + a_sin2t * (self.z0_sin_sin - self.z0_cos_cos)
            + a_cos2t * (self.z0_cos_sin + self.z0_cos_sin);

        let mut a_discriminant = b * b - 4.0 * a * c;
        let a_tol_d = 2.0 * a_dt * (b * a_db - 2.0 * (a * a_dc + c * a_da)).abs();
        if a_discriminant < a_tol_d {
            a_discriminant = 0.0;
        }

        let mut param2;
        let mut signe_sqrt_dis = 0.0;
        if a.abs() <= PCONFUSION {
            if b.abs() <= PCONFUSION {
                param2 = 0.0;
            } else {
                param2 = -c / b;
            }
        } else {
            signe_sqrt_dis = if second_solution { a_discriminant.sqrt() } else { -a_discriminant.sqrt() };
            param2 = (-b + signe_sqrt_dis) / (a + a);
        }
        Some((u, param2, a, b, c, signe_sqrt_dis))
    }

    /// OCCT IntAna_Curve::Value (L378-393).
    pub fn value(&self, theta: f64) -> Option<DVec3> {
        let (u, v, _, _, _, _) = self.internal_uv_value(theta)?;
        Some(self.internal_value(u, v))
    }

    /// Sample the curve into `n` 3D points across its parametric domain.  Used
    /// to convert the analytic IntAnaCurve into an IntPatch WLine (the OCCT
    /// IntPatch_ALineToWLine step happens inside IntPatch_Intersection before
    /// MakeCurve).
    pub fn sample(&self, n: usize) -> Vec<DVec3> {
        let d = self.domain();
        let (first, last) = (d[0], d[1]);
        let span = last - first;
        let mut out = Vec::with_capacity(n);
        for k in 0..n {
            let t = first + span * (k as f64 / (n - 1) as f64);
            if let Some(p) = self.value(t) {
                out.push(p);
            }
        }
        out
    }

    /// OCCT IntAna_Curve::FindParameter (IntAna_Curve.cxx L434-531).
    /// Projects P to the ALine, returning the list of parameters as a result
    /// of the projection.
    pub fn find_parameter(&self, p: DVec3) -> Vec<f64> {
        const A_PI_PI: f64 = std::f64::consts::TAU;
        const AN_EPS_ANG: f64 = 1.0e-8;
        const INTERNAL_PRECISION: f64 = 1.0e-8;
        let a_sq_tol_precision = rcad_kernel::precision::SQUARE_CONFUSION;

        let mut a_theta = 0.0;
        match self.surf_type {
            IntAnaCurveType::Cylinder => {
                // ElSLib::CylinderParameters(Ax3, RCyl, theP, aTheta, aZ).
                let v = p - self.cyl_origin;
                a_theta = v.y.atan2(v.x);
                let _ = self.cyl_y;
            }
            IntAnaCurveType::Cone => {
                let v = p - self.con_apex;
                a_theta = v.y.atan2(v.x);
            }
            _ => return Vec::new(),
        }
        // OCCT L468-475: the domain-boundary snap uses firstbounded/lastbounded
        // (SetCylinderQuadValues sets them false, so both conditions hold; the
        // take_z_positive branch flag is NOT part of this test).  rcad's
        // is_first_open/is_last_open (domain-parameter finiteness) mirrors the
        // bounded flags for the IntXX flow.
        if !self.is_first_open() && (self.domain_inf > a_theta)
            && ((self.domain_inf - a_theta) <= AN_EPS_ANG)
        {
            a_theta = self.domain_inf;
        } else if !self.is_last_open() && (a_theta > self.domain_sup)
            && ((a_theta - self.domain_sup) <= AN_EPS_ANG)
        {
            a_theta = self.domain_sup;
        }

        if a_theta < self.domain_inf {
            a_theta += A_PI_PI;
        } else if a_theta > self.domain_sup {
            a_theta -= A_PI_PI;
        }

        const A_MAX_PAR: usize = 5;
        let mut a_params = [self.domain_inf, self.domain_sup, a_theta,
                            if self.two_curves { self.domain_sup + self.domain_sup - a_theta } else { f64::MAX },
                            if self.two_curves { self.domain_sup + self.domain_sup - self.domain_inf } else { f64::MAX }];
        a_params[0..A_MAX_PAR - 1].sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let mut result = Vec::new();
        for i in 0..A_MAX_PAR {
            if a_params[i] > self.my_last_parameter {
                break;
            }
            if a_params[i] < self.my_first_parameter {
                continue;
            }
            if i > 0 && (a_params[i] - a_params[i - 1]) < rcad_kernel::precision::PCONFUSION {
                continue;
            }
            let Some((u, v, _, _, _, _)) = self.internal_uv_value(a_params[i]) else { continue };
            let a_p = self.internal_value(u, v);
            let mut a_sq_tol;
            if a_params[i] == a_theta
                || (self.two_curves && a_params[i] == self.domain_sup + self.domain_sup - a_theta)
            {
                a_sq_tol = INTERNAL_PRECISION;
            } else {
                a_sq_tol = a_sq_tol_precision;
            }
            if a_p.distance_squared(p) < a_sq_tol {
                result.push(a_params[i]);
            }
        }
        result
    }

    /// OCCT IntAna_Curve::D1u (L397-431).
    pub fn d1u(&self, theta: f64) -> Option<(DVec3, DVec3)> {
        let (_, _, a, _, _, signe_sqrt_dis) = self.internal_uv_value(theta)?;
        let pt = self.value(theta)?;
        if a.abs() < 1e-7 || signe_sqrt_dis.abs() < 1e-10 {
            return None;
        }
        let mut dtheta = (self.domain_sup - self.domain_inf) * 1e-6;
        let mut theta2 = theta + dtheta;
        if (theta2 < self.domain_inf)
            || ((theta2 > self.domain_sup) && !self.two_curves)
            || (theta2 > (self.domain_sup + self.domain_sup - self.domain_inf + 1e-14))
        {
            dtheta = -dtheta;
            theta2 = theta + dtheta;
        }
        let p2 = self.value(theta2)?;
        let inv = 1.0 / dtheta;
        Some((pt, (p2 - pt) * inv))
    }

    /// OCCT IntAna_Curve::InternalValue (L535-568).
    fn internal_value(&self, u: f64, v: f64) -> DVec3 {
        let v = v.clamp(-100000.0, 100000.0);
        match self.surf_type {
            IntAnaCurveType::Cylinder => {
                let co = u.cos();
                let si = u.sin();
                self.cyl_origin + self.cyl_x * (self.cyl_r * co) + self.cyl_y * (self.cyl_r * si) + self.cyl_axis * v
            }
            IntAnaCurveType::Cone => {
                // OCCT: ConeValue(U, (V-RCyl)/sin(Angle), Ax3, RCyl, Angle).
                let slant = (v - self.con_r) / self.con_angle.sin();
                let radial = self.con_r + slant * self.con_angle.sin();
                let axial = slant * self.con_angle.cos();
                let co = u.cos();
                let si = u.sin();
                self.con_apex + self.con_axis * axial + self.con_x * (radial * co) + self.con_y * (radial * si)
            }
        }
    }
}

// ============================================================================
// IntAna_IntQuadQuad
// ============================================================================

/// OCCT IntAna_IntQuadQuad — intersection of a cylinder or cone with a general
/// quadric.  Used when IntAna_QuadQuadGeo returns IntAna_NoGeometricSolution.
pub struct IntQuadQuad {
    done: bool,
    identical: bool,
    curves: [Option<IntAnaCurve>; 12],
    nb_curves: usize,
    points: Vec<DVec3>,
    my_epsilon: f64,
    my_epsilon_coeff_poly_null: f64,
}

impl IntQuadQuad {
    pub fn new() -> Self {
        IntQuadQuad {
            done: false,
            identical: false,
            curves: [None, None, None, None, None, None, None, None, None, None, None, None],
            nb_curves: 0,
            points: Vec::new(),
            my_epsilon: 1e-8,
            my_epsilon_coeff_poly_null: 1e-8,
        }
    }

    pub fn is_done(&self) -> bool { self.done }
    pub fn identical_elements(&self) -> bool { self.identical }
    pub fn nb_curves(&self) -> usize { self.nb_curves }
    pub fn curve(&self, i: usize) -> Option<&IntAnaCurve> { self.curves[i].as_ref() }
    pub fn nb_points(&self) -> usize { self.points.len() }
    pub fn point(&self, i: usize) -> DVec3 { self.points[i] }

    /// OCCT IntAna_IntQuadQuad::Perform(gp_Cylinder, IntAna_Quadric) (L375-825).
    pub fn perform_cylinder(&mut self, cyl: &rcad_kernel::geom::CylindricalSurface, quad: &IntAnaQuadric) {
        self.done = true;
        self.identical = false;
        self.nb_curves = 0;
        self.points.clear();
        let un_seul_z_par_theta = false;
        let deux_z_par_theta = true;
        let z_positif = true;
        let z_indifferent = true;
        let z_negatif = false;

        let qxx = quad.cxx; let qyy = quad.cyy; let qzz = quad.czz;
        let qxy = quad.cxy; let qxz = quad.cxz; let qyz = quad.cyz;
        let qx = quad.cx; let qy = quad.cy; let qz = quad.cz; let q1 = quad.ccte;

        let rcyl = cyl.radius;
        let a_real_epsilon = f64::EPSILON;

        // Transform the quadric into the cylinder's frame.
        let q = quad.new_coefficients(
            cyl.ref_dir.normalize_or_zero(),
            cyl.axis.normalize_or_zero().cross(cyl.ref_dir.normalize_or_zero()).normalize_or_zero(),
            cyl.axis.normalize_or_zero(),
            cyl.origin,
        );
        let qxx = q.cxx; let qyy = q.cyy; let qzz = q.czz;
        let qxy = q.cxy; let qxz = q.cxz; let qyz = q.cyz;
        let qx = q.cx; let qy = q.cy; let qz = q.cz; let q1 = q.ccte;

        if qzz.abs() < self.my_epsilon_coeff_poly_null {
            self.done = false;
            return;
        }
        let r2 = rcyl * rcyl;
        let c_1 = qz * qz - qzz * q1;
        let c_ss = r2 * (qyz * qyz - qyy * qzz);
        let c_cc = r2 * (qxz * qxz - qxx * qzz);
        let c_s = rcyl * (qyz * qz - qy * qzz);
        let c_c = rcyl * (qxz * qz - qx * qzz);
        let c_sc = r2 * (qxz * qyz - qxy * qzz);
        let mtf = MyTrigFunction::new(c_cc, c_ss, c_sc, c_c, c_s, c_1);
        let pol_dis = TrigRoots::new(c_cc - c_ss, c_sc, c_c + c_c, c_s + c_s, c_1 + c_ss, 0.0, TWO_PI);
        if !pol_dis.is_done() {
            self.done = false;
            return;
        }
        if pol_dis.infinite_roots() {
            let c0 = IntAnaCurve::set_cylinder_quad_values(
                cyl, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                self.my_epsilon, 0.0, TWO_PI, un_seul_z_par_theta, z_positif);
            let c1 = IntAnaCurve::set_cylinder_quad_values(
                cyl, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                self.my_epsilon, 0.0, TWO_PI, un_seul_z_par_theta, z_negatif);
            self.curves[0] = Some(c0);
            self.curves[1] = Some(c1);
            self.nb_curves = 2;
        } else {
            let nbsol_dis = pol_dis.nb_solutions();
            if nbsol_dis == 0 {
                if mtf.value(std::f64::consts::PI) >= -a_real_epsilon {
                    let c0 = IntAnaCurve::set_cylinder_quad_values(
                        cyl, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                        self.my_epsilon, 0.0, TWO_PI, un_seul_z_par_theta, z_positif);
                    let c1 = IntAnaCurve::set_cylinder_quad_values(
                        cyl, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                        self.my_epsilon, 0.0, TWO_PI, un_seul_z_par_theta, z_negatif);
                    self.curves[0] = Some(c0);
                    self.curves[1] = Some(c1);
                    self.nb_curves = 2;
                } else {
                    self.nb_curves = 0;
                }
            } else if nbsol_dis == 1 {
                if mtf.value(pol_dis.value(0) + std::f64::consts::PI) >= -a_real_epsilon {
                    let c0 = IntAnaCurve::set_cylinder_quad_values(
                        cyl, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                        self.my_epsilon, 0.0, TWO_PI, un_seul_z_par_theta, z_positif);
                    let c1 = IntAnaCurve::set_cylinder_quad_values(
                        cyl, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                        self.my_epsilon, 0.0, TWO_PI, un_seul_z_par_theta, z_negatif);
                    self.curves[0] = Some(c0);
                    self.curves[1] = Some(c1);
                    self.nb_curves = 2;
                } else {
                    self.nb_curves = 0;
                }
            } else {
                let mut un_pt_tg = false;
                self.nb_curves = 0;
                if nbsol_dis == 2 {
                    for i in 0..nbsol_dis {
                        let theta1 = pol_dis.value(i);
                        let theta2 = if i + 1 < nbsol_dis { pol_dis.value(i + 1) } else { pol_dis.value(0) + TWO_PI };
                        if (theta2 - theta1).abs() <= a_real_epsilon {
                            un_pt_tg = true;
                            let mut autrepar = theta1 - 0.1;
                            if autrepar < 0.0 { autrepar = theta1 + 0.1; }
                            let qwet = mtf.value(autrepar);
                            if qwet >= 0.0 {
                                let a_param = theta1 + TWO_PI;
                                let t1 = theta1;
                                let c0 = IntAnaCurve::set_cylinder_quad_values(
                                    cyl, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                                    self.my_epsilon, t1, a_param, un_seul_z_par_theta, z_positif);
                                self.curves[self.nb_curves] = Some(c0);
                                self.nb_curves += 1;
                                let c1 = IntAnaCurve::set_cylinder_quad_values(
                                    cyl, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                                    self.my_epsilon, t1, a_param, un_seul_z_par_theta, z_negatif);
                                self.curves[self.nb_curves] = Some(c1);
                                self.nb_curves += 1;
                            }
                        }
                    }
                }
                if !un_pt_tg {
                    for i in 0..nbsol_dis {
                        let theta1 = pol_dis.value(i);
                        let theta2 = if i + 1 < nbsol_dis { pol_dis.value(i + 1) } else { pol_dis.value(0) + TWO_PI };
                        if (theta2 - theta1).abs() <= 1e-12 {
                            // tangent point — skip
                        } else {
                            let qwet = mtf.value(0.5 * (theta1 + theta2))
                                + mtf.value(0.4 * theta1 + 0.6 * theta2)
                                + mtf.value(0.6 * theta1 + 0.4 * theta2);
                            if qwet >= 0.0 {
                                let theta3 = if i + 2 < nbsol_dis { pol_dis.value(i + 2) } else { pol_dis.value(0) + TWO_PI };
                                if (theta3 - theta2) < 5e-8 {
                                    let c0 = IntAnaCurve::set_cylinder_quad_values(
                                        cyl, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                                        self.my_epsilon, theta1, theta2, un_seul_z_par_theta, z_positif);
                                    self.curves[self.nb_curves] = Some(c0);
                                    self.nb_curves += 1;
                                    let c1 = IntAnaCurve::set_cylinder_quad_values(
                                        cyl, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                                        self.my_epsilon, theta1, theta2, un_seul_z_par_theta, z_negatif);
                                    self.curves[self.nb_curves] = Some(c1);
                                    self.nb_curves += 1;
                                } else {
                                    let c0 = IntAnaCurve::set_cylinder_quad_values(
                                        cyl, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                                        self.my_epsilon, theta1, theta2, deux_z_par_theta, z_indifferent);
                                    self.curves[self.nb_curves] = Some(c0);
                                    self.nb_curves += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// OCCT IntAna_IntQuadQuad::Perform(gp_Cone, IntAna_Quadric) (L841-...).
    pub fn perform_cone(&mut self, con: &rcad_kernel::geom::ConicalSurface, quad: &IntAnaQuadric) {
        let un_seul_z_par_theta = false;
        let z_positif = true;
        let z_negatif = false;

        self.done = true;
        self.identical = false;
        self.nb_curves = 0;
        self.points.clear();
        for c in self.curves.iter_mut() { *c = None; }

        let qxx = quad.cxx; let qyy = quad.cyy; let qzz = quad.czz;
        let qxy = quad.cxy; let qxz = quad.cxz; let qyz = quad.cyz;
        let qx = quad.cx; let qy = quad.cy; let qz = quad.cz; let q1 = quad.ccte;

        // Transform the quadric into the cone's APEX frame.
        let apex = con.apex_point();
        let axis = con.axis.normalize_or_zero();
        let x_dir = any_perpendicular_axis(axis);
        let y_dir = axis.cross(x_dir).normalize_or_zero();
        let q = quad.new_coefficients(x_dir, y_dir, axis, apex);
        let qxx = q.cxx; let qyy = q.cyy; let qzz = q.czz;
        let qxy = q.cxy; let qxz = q.cxz; let qyz = q.cyz;
        let qx = q.cx; let qy = q.cy; let qz = q.cz; let q1 = q.ccte;

        let tg_angle = 1.0 / con.half_angle_rad.tan();

        // A(t) for the quadratic in the "radius" parameter (the cone
        // parameterization X = W cos(t), Y = W sin(t), Z = W/tan(beta)).
        let z2_cc = qxx;
        let z2_ss = qyy;
        let z2_cte = qzz * tg_angle * tg_angle;
        let z2_sc = qxy;
        let z2_c = tg_angle * qxz;
        let z2_s = tg_angle * qyz;
        let pol_z2 = TrigRoots::new(z2_cc - z2_ss, z2_sc, z2_c + z2_c, z2_s + z2_s, z2_cte + z2_ss, 0.0, TWO_PI);
        if !pol_z2.is_done() {
            self.done = false;
            return;
        }
        let nbsol_z2 = pol_z2.nb_solutions();

        let z1_cte = 2.0 * tg_angle * qz;
        let z1_s = qy;
        let z1_c = qx;
        let pol_z1 = TrigRoots::new(0.0, 0.0, z1_c + z1_c, z1_s + z1_s, z1_cte, 0.0, TWO_PI);
        if !pol_z1.is_done() {
            self.done = false;
            return;
        }
        let mtf_z1 = MyTrigFunction::new(0.0, 0.0, 0.0, z1_c, z1_s, z1_cte);
        let nbsol_1 = pol_z1.nb_solutions();

        if pol_z2.infinite_roots() {
            if !pol_z1.infinite_roots() {
                if nbsol_1 == 0 {
                    // B(t)*z + C(t) = 0 with C(t) != 0
                    let c0 = IntAnaCurve::set_cone_quad_values(
                        con, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                        self.my_epsilon, 0.0, TWO_PI, un_seul_z_par_theta, z_positif);
                    self.curves[0] = Some(c0);
                    self.nb_curves = 1;
                }
            } else {
                if q1.abs() <= self.my_epsilon {
                    self.done = false;
                    return;
                }
                // else no solutions
            }
            return;
        }

        // Discriminant D(t) = B(t)^2 - 4 A(t) C(t).
        let c_1 = tg_angle * tg_angle * (qz * qz - qzz * q1);
        let c_ss = qy * qy - qyy * q1;
        let c_cc = qx * qx - qxx * q1;
        let c_s = tg_angle * (qy * qz - qyz * q1);
        let c_c = tg_angle * (qx * qz - qxz * q1);
        let c_sc = qx * qy - qxy * q1;
        let pol = TrigRoots::new(c_cc - c_ss, c_sc, c_c + c_c, c_s + c_s, c_1 + c_ss, 0.0, TWO_PI);
        if !pol.is_done() {
            self.done = false;
            return;
        }
        let nbsol = pol.nb_solutions();
        let mtf = MyTrigFunction::new(c_cc, c_ss, c_sc, c_c, c_s, c_1);

        if pol.infinite_roots() {
            let c0 = IntAnaCurve::set_cone_quad_values(
                con, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                self.my_epsilon, 0.0, TWO_PI, un_seul_z_par_theta, z_positif);
            let c1 = IntAnaCurve::set_cone_quad_values(
                con, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                self.my_epsilon, 0.0, TWO_PI, un_seul_z_par_theta, z_negatif);
            self.curves[0] = Some(c0);
            self.curves[1] = Some(c1);
            self.nb_curves = 2;
        } else if nbsol == 0 {
            // Discriminant has a constant sign.
            if mtf.value(std::f64::consts::PI) >= 0.0 {
                let c0 = IntAnaCurve::set_cone_quad_values(
                    con, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                    self.my_epsilon, 0.0, TWO_PI, un_seul_z_par_theta, z_positif);
                let c1 = IntAnaCurve::set_cone_quad_values(
                    con, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                    self.my_epsilon, 0.0, TWO_PI, un_seul_z_par_theta, z_negatif);
                self.curves[0] = Some(c0);
                self.curves[1] = Some(c1);
                self.nb_curves = 2;
            } else {
                self.nb_curves = 0;
            }
        } else {
            // Intervals where the discriminant is >= 0.
            let mut un_pt_tg = false;
            self.nb_curves = 0;
            if nbsol == 2 {
                for i in 0..nbsol {
                    let theta1 = pol.value(i);
                    let theta2 = if i + 1 < nbsol { pol.value(i + 1) } else { pol.value(0) + TWO_PI };
                    if (theta2 - theta1).abs() <= f64::EPSILON {
                        un_pt_tg = true;
                        let mut autrepar = theta1 - 0.1;
                        if autrepar < 0.0 { autrepar = theta1 + 0.1; }
                        let qwet = mtf.value(autrepar);
                        if qwet >= 0.0 {
                            let a_param = theta1 + TWO_PI;
                            let c0 = IntAnaCurve::set_cone_quad_values(
                                con, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                                self.my_epsilon, theta1, a_param, un_seul_z_par_theta, z_positif);
                            self.curves[self.nb_curves] = Some(c0);
                            self.nb_curves += 1;
                            let c1 = IntAnaCurve::set_cone_quad_values(
                                con, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                                self.my_epsilon, theta1, a_param, un_seul_z_par_theta, z_negatif);
                            self.curves[self.nb_curves] = Some(c1);
                            self.nb_curves += 1;
                        }
                    }
                }
            }
            if !un_pt_tg {
                for i in 0..nbsol {
                    let theta1 = pol.value(i);
                    let theta2 = if i + 1 < nbsol { pol.value(i + 1) } else { pol.value(0) + TWO_PI };
                    if (theta2 - theta1).abs() > 1e-12 {
                        let qwet = mtf.value(0.5 * (theta1 + theta2))
                            + mtf.value(0.4 * theta1 + 0.6 * theta2)
                            + mtf.value(0.6 * theta1 + 0.4 * theta2);
                        if qwet >= 0.0 {
                            let theta3 = if i + 2 < nbsol { pol.value(i + 2) } else { pol.value(0) + TWO_PI };
                            if (theta3 - theta2) < 5e-8 {
                                let c0 = IntAnaCurve::set_cone_quad_values(
                                    con, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                                    self.my_epsilon, theta1, theta2, un_seul_z_par_theta, z_positif);
                                self.curves[self.nb_curves] = Some(c0);
                                self.nb_curves += 1;
                                let c1 = IntAnaCurve::set_cone_quad_values(
                                    con, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                                    self.my_epsilon, theta1, theta2, un_seul_z_par_theta, z_negatif);
                                self.curves[self.nb_curves] = Some(c1);
                                self.nb_curves += 1;
                            } else {
                                let c0 = IntAnaCurve::set_cone_quad_values(
                                    con, qxx, qyy, qzz, qxy, qxz, qyz, qx, qy, qz, q1,
                                    self.my_epsilon, theta1, theta2, true, true);
                                self.curves[self.nb_curves] = Some(c0);
                                self.nb_curves += 1;
                            }
                        }
                    }
                }
            }
        }
    }
}

impl Default for IntQuadQuad {
    fn default() -> Self {
        Self::new()
    }
}

/// Build an IntAnaQuadric from a rcad Surface3.  OCCT's IntAna_Quadric
/// constructors take the gp_* surface; here we route through the rcad surface
/// primitives.
pub fn quadric_from_surface3(surf: &rcad_kernel::geom::Surface3) -> Option<IntAnaQuadric> {
    IntAnaQuadric::from_surface3(surf)
}
