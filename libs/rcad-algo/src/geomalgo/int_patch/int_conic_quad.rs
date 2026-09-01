//! IntAna_IntConicQuad — analytic conic-quadric intersection (Circle branch).
//!
//! 1:1 translation of OCCT `IntAna_IntConicQuad::Perform(const gp_Circ&,
//! const IntAna_Quadric&)` (IntAna_IntConicQuad.cxx L132-195) together with
//! `math_TrigonometricFunctionRoots` (math_TrigonometricFunctionRoots.cxx
//! L75-end) and the quadric coefficient extraction
//! (`IntAna_Quadric::Coefficients` / `NewCoefficients`, IntAna_Quadric.cxx
//! L122-248).
//!
//! The circle is parameterized in its own frame (x = R cos t, y = R sin t,
//! z = 0); substituting into the quadric implicit equation yields a
//! trigonometric polynomial  a*cos² + 2*b*cos*sin + c*cos + d*sin + e = 0,
//! solved analytically by reducing to a polynomial in u = tan(t/2).

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Circle3, CurveEval, Hyperbola3, Parabola3, Plane};
use rcad_kernel::math::direct_polynomial_roots::{epsilon, DirectPolynomialRoots};

use super::quad_quad_geo::{AnaResultType, QuadQuadGeo};

/// OCCT math_TrigonometricEquationFunction.hxx — the trigonometric equation
/// a·cos²(x) + 2·b·cos(x)·sin(x) + c·cos(x) + d·sin(x) + e = 0, used as the
/// Newton target in math_TrigonometricFunctionRoots L460-478.
struct TrigEquation {
    aa: f64,
    bb: f64,
    cc: f64,
    dd: f64,
    ee: f64,
}

impl rcad_kernel::math::root::FunctionValue for TrigEquation {
    fn value(&mut self, x: f64) -> Option<f64> {
        let cn = x.cos();
        let sn = x.sin();
        Some(cn * (self.aa * cn + (self.bb + self.bb) * sn + self.cc) + self.dd * sn + self.ee)
    }
}

impl rcad_kernel::math::root::FunctionWithDerivative for TrigEquation {
    fn derivative(&mut self, x: f64) -> Option<f64> {
        let cn = x.cos();
        let sn = x.sin();
        let mut d = -self.aa * cn * sn + self.bb * (cn * cn - sn * sn);
        d += d;
        d += -self.cc * sn + self.dd * cn;
        Some(d)
    }
    fn values(&mut self, x: f64) -> Option<(f64, f64)> {
        let cn = x.cos();
        let sn = x.sin();
        let aacn = self.aa * cn;
        let bbsn = self.bb * sn;
        let f = aacn * cn + bbsn * (cn + cn) + self.cc * cn + self.dd * sn + self.ee;
        let mut d = -aacn * sn + self.bb * (cn * cn - sn * sn);
        d += d;
        d += -self.cc * sn + self.dd * cn;
        Some((f, d))
    }
}

/// OCCT math_TrigonometricFunctionRoots::Perform (math_TrigonometricFunctionRoots.cxx
/// L75-551) — solutions of
///   a·cos²(x) + 2·b·cos(x)·sin(x) + c·cos(x) + d·sin(x) + e = 0
/// within [InfBound, SupBound].  Returns None when the solver reports
/// NotDone, Some((infinite, roots)) otherwise (infinite = InfiniteStatus).
fn trig_function_roots(
    a: f64, b: f64, c: f64, d: f64, e: f64,
    inf_bound: f64, sup_bound: f64,
) -> Option<(bool, Vec<f64>)> {
    use rcad_kernel::math::direct_polynomial_roots::{epsilon, DirectPolynomialRoots};
    use rcad_kernel::math::newton_function_root::NewtonFunctionRoot;
    use rcad_kernel::precision::PCONFUSION;

    const EPS: f64 = 1.5e-12;
    const TOL1: f64 = 1.0e-15;
    const NIT: i32 = 10;
    const DEPI: f64 = 2.0 * std::f64::consts::PI;
    const REAL_FIRST: f64 = -f64::MAX;
    const REAL_LAST: f64 = f64::MAX;

    // OCCT L95-123: bound setup (MyBorneInf / Delta / Mod).
    let (my_borne_inf, delta, mod_) = if inf_bound <= REAL_FIRST && sup_bound >= REAL_LAST {
        (0.0, DEPI, 0.0)
    } else if sup_bound >= REAL_LAST {
        (inf_bound, DEPI, inf_bound / DEPI)
    } else if inf_bound <= REAL_FIRST {
        (sup_bound - DEPI, DEPI, (sup_bound - DEPI) / DEPI)
    } else {
        let mut delta = sup_bound - inf_bound;
        let mod_ = inf_bound / DEPI;
        if (sup_bound - inf_bound) > DEPI {
            delta = DEPI;
        }
        (inf_bound, delta, mod_)
    };

    let mut zer: [f64; 4] = [0.0; 4];
    let mut n_zer: usize = 0;
    let mut sol: Vec<f64> = Vec::new();

    if a.abs() <= EPS && b.abs() <= EPS {
        if c.abs() <= EPS {
            if d.abs() <= EPS {
                if e.abs() <= EPS {
                    // OCCT L133-135: infinite number of solutions.
                    return Some((true, Vec::new()));
                } else {
                    // OCCT L137-140: no solution.
                    return Some((false, Vec::new()));
                }
            } else {
                // OCCT L142-173: d·sin(x) + e = 0.
                let aa = -e / d;
                if aa.abs() > 1.0 {
                    return Some((false, Vec::new()));
                }
                zer[0] = aa.asin();
                zer[1] = std::f64::consts::PI - zer[0];
                n_zer = 2;
                for i in 0..n_zer {
                    if zer[i] <= -EPS {
                        zer[i] = DEPI - zer[i].abs();
                    }
                    zer[i] += mod_.trunc() * DEPI;
                    let x = zer[i] - my_borne_inf;
                    if x > -epsilon(delta) && x < delta + epsilon(delta) {
                        sol.push(zer[i]);
                    }
                }
            }
            return Some((false, sol));
        } else if d.abs() <= EPS {
            // OCCT L175-207: c·cos(x) + e = 0.
            let aa = -e / c;
            if aa.abs() > 1.0 {
                return Some((false, Vec::new()));
            }
            zer[0] = aa.acos();
            zer[1] = -zer[0];
            n_zer = 2;
            for i in 0..n_zer {
                if zer[i] <= -EPS {
                    zer[i] = DEPI - zer[i].abs();
                }
                zer[i] += mod_.trunc() * DEPI;
                let x = zer[i] - my_borne_inf;
                if x >= -epsilon(delta) && x <= delta + epsilon(delta) {
                    sol.push(zer[i]);
                }
            }
            return Some((false, sol));
        } else {
            // OCCT L208-236: quadratic in u = tan(x/2).
            let aa = e - c;
            let bb = 2.0 * d;
            let cc = e + c;
            let resol = DirectPolynomialRoots::new_quadratic(aa, bb, cc);
            if !resol.is_done() {
                return None;
            } else if !resol.infinite_roots() {
                n_zer = resol.nb_solutions();
                for i in 0..n_zer {
                    zer[i] = resol.value(i + 1);
                }
            } else {
                return Some((true, Vec::new()));
            }
        }
    } else {
        // OCCT L240-354: two additional analytical cases.
        if a.abs() <= EPS && e.abs() <= EPS {
            if c.abs() <= EPS {
                // OCCT L243-296: 2·B·sin·cos + D·sin = 0.
                n_zer = 2;
                zer[0] = 0.0;
                zer[1] = std::f64::consts::PI;
                let aa = -d / (2.0 * b);
                if aa.abs() <= 1.0 + PCONFUSION {
                    n_zer = 4;
                    if aa >= 1.0 {
                        zer[2] = 0.0;
                        zer[3] = 0.0;
                    } else if aa <= -1.0 {
                        zer[2] = std::f64::consts::PI;
                        zer[3] = std::f64::consts::PI;
                    } else {
                        zer[2] = aa.acos();
                        zer[3] = DEPI - zer[2];
                    }
                }
                for i in 0..n_zer {
                    if zer[i] <= my_borne_inf - EPS {
                        zer[i] += DEPI;
                    }
                    zer[i] += mod_.trunc() * DEPI;
                    let x = zer[i] - my_borne_inf;
                    if x >= -PCONFUSION && x <= delta + PCONFUSION {
                        if zer[i] < inf_bound {
                            zer[i] = inf_bound;
                        }
                        if zer[i] > sup_bound {
                            zer[i] = sup_bound;
                        }
                        sol.push(zer[i]);
                    }
                }
                return Some((false, sol));
            }
            if d.abs() <= EPS {
                // OCCT L298-352: 2·B·sin·cos + C·cos = 0.
                n_zer = 2;
                zer[0] = std::f64::consts::FRAC_PI_2;
                zer[1] = std::f64::consts::PI * 1.5;
                let aa = -c / (2.0 * b);
                if aa.abs() <= 1.0 + PCONFUSION {
                    n_zer = 4;
                    if aa >= 1.0 {
                        zer[2] = std::f64::consts::FRAC_PI_2;
                        zer[3] = std::f64::consts::FRAC_PI_2;
                    } else if aa <= -1.0 {
                        zer[2] = std::f64::consts::PI * 1.5;
                        zer[3] = std::f64::consts::PI * 1.5;
                    } else {
                        zer[2] = aa.asin();
                        zer[3] = std::f64::consts::PI - zer[2];
                    }
                }
                for i in 0..n_zer {
                    if zer[i] <= my_borne_inf - EPS {
                        zer[i] += DEPI;
                    }
                    zer[i] += mod_.trunc() * DEPI;
                    let x = zer[i] - my_borne_inf;
                    if x >= -PCONFUSION && x <= delta + PCONFUSION {
                        if zer[i] < inf_bound {
                            zer[i] = inf_bound;
                        }
                        if zer[i] > sup_bound {
                            zer[i] = sup_bound;
                        }
                        sol.push(zer[i]);
                    }
                }
                return Some((false, sol));
            }
        }

        // OCCT L356-434: the general quartic.
        let mut ko = [
            a - c + e,
            2.0 * d - 4.0 * b,
            2.0 * e - 2.0 * a,
            4.0 * b + 2.0 * d,
            a + c + e,
        ];
        loop {
            let mut bko = false;
            let resol4 = DirectPolynomialRoots::new_quartic(ko[0], ko[1], ko[2], ko[3], ko[4]);
            if !resol4.is_done() {
                return None;
            } else if !resol4.infinite_roots() {
                n_zer = resol4.nb_solutions();
                for i in 0..n_zer {
                    zer[i] = resol4.value(i + 1);
                }
            } else {
                return Some((true, Vec::new()));
            }

            // OCCT L386-400: bubble sort Zer.
            let mut triok;
            loop {
                triok = true;
                for i in 0..n_zer.saturating_sub(1) {
                    if zer[i] > zer[i + 1] {
                        zer.swap(i, i + 1);
                        triok = false;
                    }
                }
                if triok {
                    break;
                }
            }

            // OCCT L402-422: double-root check; on a numerical double root
            // scale ko by 1e-4 and re-solve.
            for i in 0..n_zer.saturating_sub(1) {
                if (zer[i + 1] - zer[i]).abs() < EPS {
                    let qw = zer[i + 1];
                    let va = ko[3] + qw * (2.0 * ko[2] + qw * (3.0 * ko[1] + qw * (4.0 * ko[0])));
                    if va.abs() > EPS {
                        bko = true;
                        break;
                    }
                }
            }
            if bko {
                for k in 0..5 {
                    ko[k] *= 0.0001;
                }
            } else {
                break;
            }
        }
    }

    // OCCT L437-504: verification of the solutions against the bounds.
    let supm_infs100 = (sup_bound - inf_bound) * 0.01;
    for i in 0..n_zer {
        let mut teta = zer[i].atan() * 2.0;
        if zer[i] <= -EPS {
            teta = DEPI - teta.abs();
        }
        teta += mod_.trunc() * DEPI;
        if teta - my_borne_inf < 0.0 {
            teta += DEPI;
        }
        let x = teta - my_borne_inf;
        if x >= -epsilon(delta) && x <= delta + epsilon(delta) {
            // OCCT L460-478: Newton refinement.
            let mut teta_newton = teta;
            let mut my_f = TrigEquation {
                aa: a,
                bb: b,
                cc: c,
                dd: d,
                ee: e,
            };
            let resol = NewtonFunctionRoot::new_full_range(&mut my_f, x, TOL1, EPS, NIT);
            if resol.is_done() {
                teta_newton = resol.root();
            }
            let delta_newton = teta_newton - teta;
            if delta_newton.abs() <= supm_infs100 {
                teta = teta_newton;
            }

            // OCCT L480-502: insert Teta into the sorted Sol.
            let mut flag4 = false;
            for k in 0..sol.len() {
                if teta < sol[k] {
                    sol.insert(k, teta);
                    flag4 = true;
                    break;
                }
            }
            if !flag4 {
                sol.push(teta);
            }
        }
    }

    // OCCT L505-550: "Cas particulier de PI" — t = pi is a solution exactly
    // when the trig polynomial vanishes there, i.e. ko(1) = A - C + E = 0
    // (at t=pi: cos=-1, sin=0).  The substitution u = tan(t/2) maps t=pi to
    // u=infinity, so the quartic (and the quadratic when e-c=0) misses it.
    if sol.len() < 4 {
        let start_index = sol.len() + 1;
        for sol_it in start_index..=4 {
            let mut teta = std::f64::consts::PI + mod_.trunc() * DEPI;
            let x = teta - my_borne_inf;
            if x >= -epsilon(delta) && x <= delta + epsilon(delta) {
                if (a - c + e).abs() <= EPS {
                    let mut flag4 = false;
                    let mut j = 0usize;
                    for k in 0..sol.len() {
                        j = k;
                        if teta < sol[k] {
                            flag4 = true;
                            break;
                        }
                        if sol_it == start_index && (teta - sol[k]).abs() <= EPS {
                            return Some((false, sol));
                        }
                    }
                    if !flag4 {
                        sol.push(teta);
                    } else {
                        sol.insert(j, teta);
                    }
                }
            }
        }
    }

    Some((false, sol))
}

/// OCCT IntAna_Quadric::NewCoefficients — transform the quadric coefficients
/// into the given frame (the circle's position).  The implicit equation
/// Cxx·x² + Cyy·y² + Czz·z² + 2·Cxy·xy + 2·Cxz·xz + 2·Cyz·yz
///   + Cx·x + Cy·y + Cz·z + Ccte = 0
/// is expressed in coordinates of `frame` (origin + orthonormal axes).
/// The linear terms follow the OCCT "2·C1·x" convention (stored halved).
pub(crate) struct FrameCoefs {
    pub(crate) xx: f64, pub(crate) yy: f64, pub(crate) zz: f64,
    pub(crate) xy: f64, pub(crate) xz: f64, pub(crate) yz: f64,
    pub(crate) x: f64, pub(crate) y: f64, pub(crate) z: f64,
    pub(crate) cte: f64,
}

fn new_coefficients(
    co: &FrameCoefs,
    origin: DVec3, x_dir: DVec3, y_dir: DVec3, z_dir: DVec3,
) -> FrameCoefs {
    // Inverse transform: world = origin + X·x_dir + Y·y_dir + Z·z_dir, so
    // X = (world - origin)·x_dir etc.  The change-of-coordinates matrix t:
    //   x = t11·X + t12·Y + t13·Z + t14  (x is the first world axis).
    let t11 = x_dir.x; let t12 = y_dir.x; let t13 = z_dir.x; let t14 = origin.x;
    let t21 = x_dir.y; let t22 = y_dir.y; let t23 = z_dir.y; let t24 = origin.y;
    let t31 = x_dir.z; let t32 = y_dir.z; let t33 = z_dir.z; let t34 = origin.z;

    let (cxx0, cyy0, czz0, cxy0, cxz0, cyz0, cx0, cy0, cz0, ccte0) =
        (co.xx, co.yy, co.zz, co.xy, co.xz, co.yz, co.x, co.y, co.z, co.cte);

    let t11p2 = t11 * t11; let t21p2 = t21 * t21; let t31p2 = t31 * t31;
    let t12p2 = t12 * t12; let t22p2 = t22 * t22; let t32p2 = t32 * t32;
    let t13p2 = t13 * t13; let t23p2 = t23 * t23; let t33p2 = t33 * t33;
    let t14p2 = t14 * t14; let t24p2 = t24 * t24; let t34p2 = t34 * t34;

    let ccte = ccte0 + t14p2 * cxx0 + t24p2 * cyy0 + t34p2 * czz0
        + 2.0 * (t14 * (cx0 + t24 * cxy0 + t34 * cxz0)
            + t24 * (cy0 + t34 * cyz0) + t34 * cz0);

    let cxx = t11p2 * cxx0 + t21p2 * cyy0 + t31p2 * czz0
        + 2.0 * (t11 * (t21 * cxy0 + t31 * cxz0) + t21 * t31 * cyz0);
    let cyy = t12p2 * cxx0 + t22p2 * cyy0 + t32p2 * czz0
        + 2.0 * (t12 * (t22 * cxy0 + t32 * cxz0) + t22 * t32 * cyz0);
    let czz = t13p2 * cxx0 + t33p2 * czz0 + t23p2 * cyy0
        + 2.0 * (t13 * (t23 * cxy0 + t33 * cxz0) + t23 * t33 * cyz0);

    let cx = t11 * (cx0 + t14 * cxx0 + t24 * cxy0 + t34 * cxz0)
        + t14 * (t21 * cxy0 + t31 * cxz0)
        + t21 * (cy0 + t24 * cyy0 + t34 * cyz0)
        + t31 * (t24 * cyz0 + cz0 + t34 * czz0);
    let cy = t12 * (cx0 + t14 * cxx0 + t24 * cxy0 + t34 * cxz0)
        + t14 * (t22 * cxy0 + t32 * cxz0)
        + t22 * (cy0 + t24 * cyy0 + t34 * cyz0)
        + t32 * (cz0 + t24 * cyz0 + t34 * czz0);
    let cz = t13 * (cx0 + t14 * cxx0 + t24 * cxy0 + t34 * cxz0)
        + t14 * (t23 * cxy0 + t33 * cxz0)
        + t23 * (cy0 + t24 * cyy0 + t34 * cyz0)
        + t33 * (cz0 + t24 * cyz0 + t34 * czz0);

    let cxy = t11 * (t12 * cxx0 + t22 * cxy0 + t32 * cxz0)
        + t12 * (t21 * cxy0 + t31 * cxz0)
        + t21 * (t22 * cyy0 + t32 * cyz0)
        + t31 * (t22 * cyz0 + t32 * czz0);
    let cxz = t11 * (t13 * cxx0 + t23 * cxy0 + t33 * cxz0)
        + t13 * (t21 * cxy0 + t31 * cxz0)
        + t21 * (t23 * cyy0 + t33 * cyz0)
        + t31 * (t23 * cyz0 + t33 * czz0);
    let cyz = t12 * (t13 * cxx0 + t23 * cxy0 + t33 * cxz0)
        + t13 * (t22 * cxy0 + t32 * cxz0)
        + t22 * (t23 * cyy0 + t33 * cyz0)
        + t32 * (t23 * cyz0 + t33 * czz0);

    FrameCoefs { xx: cxx, yy: cyy, zz: czz, xy: cxy, xz: cxz, yz: cyz, x: cx, y: cy, z: cz, cte: ccte }
}

/// OCCT IntAna_Quadric::Coefficients (IntAna_Quadric.cxx L122-143) — the
/// absolute-frame coefficients of the quadric, in the "2·C1·x" convention.
/// Matches gp_Cylinder/Cone/Sphere::Coefficients (verified by implicit
/// expansion).  None for a non-canonic quadric.
pub(crate) fn quadric_frame_coefs(
    quad: &crate::geomalgo::int_surf::quadric::Quadric,
) -> Option<FrameCoefs> {
    use crate::geomalgo::int_surf::quadric::QuadricType;
    let co = match quad.type_quadric() {
        QuadricType::Plane => {
            // Plane: pn·x + d = 0  ->  Qx = pn.x/2, Qy = pn.y/2, Qz = pn.z/2.
            let pl = quad.plane();
            let d = -pl.normal.dot(pl.origin);
            FrameCoefs {
                xx: 0.0, yy: 0.0, zz: 0.0, xy: 0.0, xz: 0.0, yz: 0.0,
                x: 0.5 * pl.normal.x,
                y: 0.5 * pl.normal.y,
                z: 0.5 * pl.normal.z,
                cte: d,
            }
        }
        QuadricType::Cylinder => {
            // gp_Cylinder implicit: X² + Y² - R² = 0 in the local frame, i.e.
            // |P-O|² - ((P-O)·a)² - R² = 0.  Coeffs follow the OCCT
            // "2·C1·x" convention (linear terms stored halved).
            let o = quad.axis_loc();
            let a = quad.axis_dir();
            let r = quad.radius();
            let oa = o.dot(a);
            FrameCoefs {
                xx: 1.0 - a.x * a.x,
                yy: 1.0 - a.y * a.y,
                zz: 1.0 - a.z * a.z,
                xy: -a.x * a.y,
                xz: -a.x * a.z,
                yz: -a.y * a.z,
                x: -(o.x - a.x * oa),
                y: -(o.y - a.y * oa),
                z: -(o.z - a.z * oa),
                cte: o.length_squared() - oa * oa - r * r,
            }
        }
        QuadricType::Cone => {
            // gp_Cone implicit: X² + Y² - (R + Z·tan)² = 0 in the local frame
            // (O = reference point where the radius is R, Z = axis).
            //   |P-O|² - sec²·((P-O)·a)² - 2·R·tan·((P-O)·a) - R² = 0.
            let o = quad.axis_loc();
            let a = quad.axis_dir();
            let r = quad.radius();
            let tg = quad.semi_angle().tan();
            let sec2 = 1.0 + tg * tg;
            let oa = o.dot(a);
            FrameCoefs {
                xx: 1.0 - sec2 * a.x * a.x,
                yy: 1.0 - sec2 * a.y * a.y,
                zz: 1.0 - sec2 * a.z * a.z,
                xy: -sec2 * a.x * a.y,
                xz: -sec2 * a.x * a.z,
                yz: -sec2 * a.y * a.z,
                x: -o.x + sec2 * a.x * oa - r * tg * a.x,
                y: -o.y + sec2 * a.y * oa - r * tg * a.y,
                z: -o.z + sec2 * a.z * oa - r * tg * a.z,
                cte: o.length_squared() - sec2 * oa * oa + 2.0 * r * tg * oa - r * r,
            }
        }
        QuadricType::Sphere => {
            let sph = quad.sphere();
            let center = sph.center;
            let rad = sph.radius;
            FrameCoefs {
                xx: 1.0, yy: 1.0, zz: 1.0, xy: 0.0, xz: 0.0, yz: 0.0,
                x: -center.x, y: -center.y, z: -center.z,
                cte: center.length_squared() - rad * rad,
            }
        }
        _ => return None,
    };
    Some(co)
}

/// OCCT IntAna_IntConicQuad::Perform(const gp_Circ&, const IntAna_Quadric&)
/// — returns (in_quadric, points) where each point is (3D point, t param).
/// A circle lying entirely on the quadric reports `in_quadric = true`.
pub fn intersect_circle_quadric(
    circle: &Circle3,
    quad: &crate::geomalgo::int_surf::quadric::Quadric,
) -> Option<(bool, Vec<(DVec3, f64)>)> {
    let r = circle.radius;
    let rr = r * r;

    // OCCT Quad.Coefficients(...) — absolute-frame coefficients.  The trig
    // solver's Eps is internal (1.5e-12).
    let co = quadric_frame_coefs(quad)?;
    let inf = 0.0;
    let sup = 2.0 * std::f64::consts::PI;

    // NewCoefficients into the circle's frame (the circle's position).
    let nco = new_coefficients(&co, circle.center, circle.x_dir, circle.y_dir, circle.normal);

    // Polynomial coefficients after substituting x=R cos t, y=R sin t, z=0.
    let p_cos_cos = rr * nco.xx;
    let p_sin_sin = rr * nco.yy;
    let p_sin = r * nco.y;
    let p_cos = r * nco.x;
    let p_cos_sin = rr * nco.xy;
    let p_cte = nco.cte;

    // a·cos² + 2b·cos·sin + c·cos + d·sin + e = 0
    // OCCT IntAna_IntConicQuad.cxx L176-194: IsDone -> inquadric (infinite) or
    // the NbSolutions points.
    let (in_quadric, ts) = match trig_function_roots(
        p_cos_cos - p_sin_sin,
        p_cos_sin,
        2.0 * p_cos,
        2.0 * p_sin,
        p_cte + p_sin_sin,
        inf,
        sup,
    ) {
        None => return None,
        Some(r) => r,
    };

    if in_quadric {
        return Some((true, Vec::new()));
    }

    let pts: Vec<(DVec3, f64)> = ts
        .into_iter()
        .map(|t| {
            let p = circle.center
                + circle.x_dir * (r * t.cos())
                + circle.y_dir * (r * t.sin());
            (p, t)
        })
        .collect();
    Some((false, pts))
}

/// OCCT IntAna_IntConicQuad::Perform(const gp_Circ&, const gp_Pln&) — the
/// plane case via IntAna_QuadQuadGeo (line ∩ circle in the circle's plane).
/// Returns (is_parallel, in_quadric, points).
pub fn intersect_circle_plane(
    circle: &Circle3,
    plane: &Plane,
    tol_ang: f64,
    tol: f64,
) -> (bool, bool, Vec<(DVec3, f64)>) {
    use crate::geomalgo::int_surf::quadric::Quadric;
    let plconic = Plane {
        origin: circle.center,
        normal: circle.normal,
        u_dir: circle.x_dir,
        v_dir: circle.y_dir,
    };
    let mut gg = QuadQuadGeo::new();
    gg.perform_plane_plane(
        &Quadric::from_plane(&plconic),
        &Quadric::from_plane(plane),
        tol_ang,
        tol,
    );
    match gg.type_inter() {
        AnaResultType::Empty => {
            // parallel planes: in the quadric when close enough.
            let distmax = plane
                .normal
                .dot(circle.center - plane.origin)
                .abs()
                .max(0.0)
                + circle.radius * tol_ang;
            (true, distmax < tol, Vec::new())
        }
        AnaResultType::Same => (false, true, Vec::new()),
        AnaResultType::Line => {
            // OCCT IntAna_IntConicQuad::Perform(Circ, Pln) L522-559: project
            // the intersection line into the circle's plane and intersect the
            // 2D line with the 2D circle via IntAna2d_AnaIntersection
            // (IntAna2d_AnaIntersection_3.cxx L24-108).
            let line_origin = gg.line(1).origin;
            let line_dir = gg.line(1).direction;
            // 2D line in the circle's frame: X = (P-O)·x_dir, Y = (P-O)·y_dir.
            let o2 = DVec2::new(
                (line_origin - circle.center).dot(circle.x_dir),
                (line_origin - circle.center).dot(circle.y_dir),
            );
            let d2 = DVec2::new(line_dir.dot(circle.x_dir), line_dir.dot(circle.y_dir));
            // OCCT gp_Lin2d::Coefficients(A, B, C0): the unit normal of the
            // line, so the signed distance d = A·cx + B·cy + C0.
            let len = d2.length();
            let (a, b) = (d2.y / len, -d2.x / len);
            let c0 = -(a * o2.x + b * o2.y);
            let (cx, cy) = (0.0, 0.0); // the circle's center in its own frame
            let r = circle.radius;
            let eps_r = rcad_kernel::math::direct_polynomial_roots::epsilon(r);
            let d = a * cx + b * cy + c0;

            let mut pts = Vec::new();
            if d.abs() - r > eps_r {
                // OCCT L38-42: no solution.
            } else if (d.abs() - r).abs() <= eps_r {
                // OCCT L52-78: tangency — one point.
                let xs = cx - d * a;
                let ys = cy - d * b;
                let p3 = circle.center + circle.x_dir * xs + circle.y_dir * ys;
                pts.push((p3, ys.atan2(xs)));
            } else {
                // OCCT L80-105: two intersection points.
                let h = (r * r - d * d).sqrt();
                let xs1 = cx - d * a - h * b;
                let ys1 = cy - d * b + h * a;
                let xs2 = cx - d * a + h * b;
                let ys2 = cy - d * b - h * a;
                pts.push((
                    circle.center + circle.x_dir * xs1 + circle.y_dir * ys1,
                    ys1.atan2(xs1),
                ));
                pts.push((
                    circle.center + circle.x_dir * xs2 + circle.y_dir * ys2,
                    ys2.atan2(xs2),
                ));
            }
            (false, false, pts)
        }
        _ => (false, false, Vec::new()),
    }
}

/// OCCT IntAna_IntConicQuad::Perform(const gp_Parab&, const IntAna_Quadric&)
/// (IntAna_IntConicQuad.cxx L273-325) — the conic-quadric intersection for a
/// parabola.  In the parabola's frame x = y²/(2p), y = y, z = 0, so the
/// substituted quadric becomes a quartic in y.  Used by Intf_Tool::Inters3d
/// with the plane quadric.  Returns None when the solver is not done,
/// Some((in_quadric, (point, param) list)) otherwise.
pub fn intersect_parabola_quadric(
    parabola: &Parabola3,
    quad: &crate::geomalgo::int_surf::quadric::Quadric,
) -> Option<(bool, Vec<(DVec3, f64)>)> {
    let co = quadric_frame_coefs(quad)?;
    // The parabola's frame: X = axis (symmetry, toward focus), Y = N × X,
    // Z = N (gp_Parab::Position() — gp_Ax2 whose Z is the normal).
    let y_dir = parabola.normal.cross(parabola.axis_dir).normalize_or_zero();
    let nco = new_coefficients(&co, parabola.vertex, parabola.axis_dir, y_dir, parabola.normal);

    // OCCT f = P.Focal() (the focal LENGTH); Un_Sur_2p = 0.25 / f.
    let f = parabola.focal_param * 0.5;
    let un_sur_2p = 0.25 / f;

    let a4 = nco.xx * un_sur_2p * un_sur_2p;
    let a3 = (nco.xy + nco.xy) * un_sur_2p;
    let a2 = nco.yy + (nco.x + nco.x) * un_sur_2p;
    let a1 = nco.y + nco.y;
    let a0 = nco.cte;

    let roots = DirectPolynomialRoots::new_quartic(a4, a3, a2, a1, a0);
    if !roots.is_done() {
        return None;
    }
    if roots.infinite_roots() {
        return Some((true, Vec::new()));
    }
    let pts: Vec<(DVec3, f64)> = (1..=roots.nb_solutions())
        .map(|i| {
            let t = roots.value(i);
            (parabola.point_at(t), t)
        })
        .collect();
    Some((false, pts))
}

/// OCCT IntAna_IntConicQuad::Perform(const gp_Hypr&, const IntAna_Quadric&)
/// (IntAna_IntConicQuad.cxx L335-397) — the conic-quadric intersection for a
/// hyperbola.  In the hyperbola's frame x = R·Ch(t), y = r·Sh(t), z = 0,
/// rewritten as a quartic in S = Exp(t); valid solutions have S >= RealEpsilon
/// and the conic parameter is t = Ln(S).  Returns None when the solver is not
/// done, Some((in_quadric, (point, param) list)) otherwise.
pub fn intersect_hyperbola_quadric(
    hyperbola: &Hyperbola3,
    quad: &crate::geomalgo::int_surf::quadric::Quadric,
) -> Option<(bool, Vec<(DVec3, f64)>)> {
    let co = quadric_frame_coefs(quad)?;
    // The hyperbola's frame: X = major, Y = N × X (minor), Z = N.
    let y_dir = hyperbola.normal.cross(hyperbola.major_dir).normalize_or_zero();
    let nco = new_coefficients(
        &co,
        hyperbola.center,
        hyperbola.major_dir,
        y_dir,
        hyperbola.normal,
    );

    let r_big = hyperbola.semi_major;
    let r_small = hyperbola.semi_minor;
    let rr_p2 = r_big * r_big;
    let rr_small_p2 = r_small * r_small;
    let rr_prod = r_big * r_small;

    let a4 = rr_p2 * nco.xx + rr_prod * (nco.xy + nco.xy) + rr_small_p2 * nco.yy;
    let a3 = 4.0 * (r_big * nco.x + r_small * nco.y);
    let a2 = 2.0 * ((nco.cte + nco.cte) + nco.xx * rr_p2 - nco.yy * rr_small_p2);
    let a1 = 4.0 * (r_big * nco.x - r_small * nco.y);
    let a0 = nco.xx * rr_p2 - rr_prod * (nco.xy + nco.xy) + nco.yy * rr_small_p2;

    let roots = DirectPolynomialRoots::new_quartic(a4, a3, a2, a1, a0);
    if !roots.is_done() {
        return None;
    }
    if roots.infinite_roots() {
        return Some((true, Vec::new()));
    }
    let mut pts = Vec::new();
    for i in 1..=roots.nb_solutions() {
        let t = roots.value(i);
        if t >= f64::EPSILON {
            let param = t.ln();
            pts.push((hyperbola.point_at(param), param));
        }
    }
    Some((false, pts))
}
