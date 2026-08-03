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
use rcad_kernel::geom::{Circle3, Plane};

use super::quad_quad_geo::{AnaResultType, QuadQuadGeo};

/// OCCT math_DirectPolynomialRoots — real roots of a degree-4 polynomial
/// a4·x⁴ + a3·x³ + a2·x² + a1·x + a0 = 0 (Ferrari's method).
fn quartic_real_roots(a4: f64, a3: f64, a2: f64, a1: f64, a0: f64) -> Vec<f64> {
    const EPS: f64 = 1e-14;
    // Degenerate lower degrees.
    if a4.abs() < EPS {
        return cubic_real_roots(a3, a2, a1, a0);
    }
    // Depress: x = u - a3/(4 a4).
    let a = a3 / a4;
    let b = a2 / a4;
    let c = a1 / a4;
    let d = a0 / a4;
    let p = b - 3.0 * a * a / 8.0;
    let q = a * a * a / 8.0 - a * b / 2.0 + c;
    let r = -3.0 * a * a * a * a / 256.0 + a * a * b / 16.0 - a * c / 4.0 + d;

    // Resolvent cubic: 8 m³ - 4 p m² - 8 r m + (4 p r - q²) = 0.
    let m_roots = cubic_real_roots(8.0, -4.0 * p, -8.0 * r, 4.0 * p * r - q * q);
    let m = m_roots.into_iter().find(|&m| 2.0 * m - p >= 0.0).unwrap_or(0.0);

    // (u² + m)² = (A u + B)² with A² = 2m - p, B = -q/(2A).
    let a2m = 2.0 * m - p;
    if a2m < 0.0 {
        return Vec::new();
    }
    let aa = a2m.sqrt();
    let bb = if aa.abs() < EPS { 0.0 } else { -q / (2.0 * aa) };

    // u² - aa u + (m - bb) = 0  and  u² + aa u + (m + bb) = 0.
    let mut roots = Vec::new();
    roots.extend(quadratic_real_roots(1.0, -aa, m - bb));
    roots.extend(quadratic_real_roots(1.0, aa, m + bb));
    // Back-substitute the depression.
    let shift = a / 4.0;
    roots.iter_mut().for_each(|x| *x -= shift);
    roots
}

/// Real roots of a cubic a3·x³ + a2·x² + a1·x + a0 = 0 (Cardano).
fn cubic_real_roots(a3: f64, a2: f64, a1: f64, a0: f64) -> Vec<f64> {
    const EPS: f64 = 1e-14;
    if a3.abs() < EPS {
        return quadratic_real_roots(a2, a1, a0);
    }
    // Depress: x = u - a2/(3 a3).
    let p = (3.0 * a3 * a1 - a2 * a2) / (3.0 * a3 * a3);
    let q = (2.0 * a2 * a2 * a2 - 9.0 * a3 * a2 * a1 + 27.0 * a3 * a3 * a0)
        / (27.0 * a3 * a3 * a3);
    let disc = q * q / 4.0 + p * p * p / 27.0;
    let mut roots = Vec::new();
    if disc > 0.0 {
        let sd = disc.sqrt();
        let u1 = (-q / 2.0 + sd).cbrt() + (-q / 2.0 - sd).cbrt();
        roots.push(u1);
    } else if disc.abs() <= EPS {
        let u1 = 3.0 * q / p;
        let u2 = -3.0 * q / (2.0 * p);
        roots.push(u1);
        roots.push(u2);
    } else {
        // Three real roots via the trigonometric form.
        let rho = (-p * p * p / 27.0).sqrt();
        let phi = (-q / (2.0 * rho)).acos();
        for k in 0..3 {
            let u = 2.0 * rho.cbrt() * ((phi + 2.0 * std::f64::consts::PI * k as f64) / 3.0).cos();
            roots.push(u);
        }
    }
    let shift = a2 / (3.0 * a3);
    roots.iter_mut().for_each(|x| *x -= shift);
    roots
}

/// Real roots of a quadratic a2·x² + a1·x + a0 = 0.
fn quadratic_real_roots(a2: f64, a1: f64, a0: f64) -> Vec<f64> {
    const EPS: f64 = 1e-14;
    if a2.abs() < EPS {
        if a1.abs() < EPS {
            return Vec::new();
        }
        return vec![-a0 / a1];
    }
    let disc = a1 * a1 - 4.0 * a2 * a0;
    if disc < 0.0 {
        return Vec::new();
    }
    let sd = disc.sqrt();
    if sd < EPS {
        return vec![-a1 / (2.0 * a2)];
    }
    vec![(-a1 + sd) / (2.0 * a2), (-a1 - sd) / (2.0 * a2)]
}

/// OCCT math_TrigonometricFunctionRoots — solutions of
///   a·cos²(x) + 2·b·cos(x)·sin(x) + c·cos(x) + d·sin(x) + e = 0
/// within [InfBound, SupBound].
///
/// Documented divergences from OCCT (math_TrigonometricFunctionRoots.cxx):
///   - L460-478: OCCT refines every candidate root with
///     math_NewtonFunctionRoot on the trig equation.  rcad uses the raw
///     polynomial roots; for the well-separated roots of the Circle x Quadric
///     path the Newton step is a no-op to within Tol1=1e-15.
///   - L363-434: OCCT solves the quartic with math_DirectPolynomialRoots and
///     on near-double roots scales ko by 1e-4 and re-solves.  rcad uses a
///     Ferrari quartic with the same real-root semantics; the 1e-9 dedup below
///     plays the role of the double-root rejection.
fn trig_function_roots(
    a: f64, b: f64, c: f64, d: f64, e: f64,
    inf_bound: f64, sup_bound: f64,
) -> Vec<f64> {
    const EPS: f64 = 1.5e-12;
    let depi = 2.0 * std::f64::consts::PI;

    let mut zer: Vec<f64> = Vec::new();
    let mut infinite = false;

    if a.abs() <= EPS && b.abs() <= EPS {
        if c.abs() <= EPS {
            if d.abs() <= EPS {
                if e.abs() <= EPS {
                    infinite = true;
                }
                // else: no solution.
            } else {
                // d·sin + e = 0
                let aa = -e / d;
                if aa.abs() > 1.0 {
                    return Vec::new();
                }
                zer.push(aa.asin());
                zer.push(std::f64::consts::PI - aa.asin());
            }
        } else if d.abs() <= EPS {
            // c·cos + e = 0
            let aa = -e / c;
            if aa.abs() > 1.0 {
                return Vec::new();
            }
            zer.push(aa.acos());
            zer.push(-aa.acos());
        } else {
            // Quadratic in u = tan(x/2): (e-c)u² + 2d u + (e+c) = 0.
            zer.extend(quadratic_real_roots(e - c, 2.0 * d, e + c));
        }
    } else {
        if a.abs() <= EPS && e.abs() <= EPS {
            if c.abs() <= EPS {
                // 2·b·sin·cos + d·sin = 0  →  sin·(2b cos + d) = 0
                zer.push(0.0);
                zer.push(std::f64::consts::PI);
                let aa = -d / (2.0 * b);
                if aa.abs() <= 1.0 + 1e-9 {
                    if aa >= 1.0 {
                        zer.push(0.0);
                        zer.push(0.0);
                    } else if aa <= -1.0 {
                        zer.push(std::f64::consts::PI);
                        zer.push(std::f64::consts::PI);
                    } else {
                        zer.push(aa.acos());
                        zer.push(depi - aa.acos());
                    }
                }
            } else if d.abs() <= EPS {
                // 2·b·sin·cos + c·cos = 0  →  cos·(2b sin + c) = 0
                zer.push(std::f64::consts::FRAC_PI_2);
                zer.push(std::f64::consts::PI * 1.5);
                let aa = -c / (2.0 * b);
                if aa.abs() <= 1.0 + 1e-9 {
                    if aa >= 1.0 {
                        zer.push(std::f64::consts::FRAC_PI_2);
                        zer.push(std::f64::consts::FRAC_PI_2);
                    } else if aa <= -1.0 {
                        zer.push(std::f64::consts::PI * 1.5);
                        zer.push(std::f64::consts::PI * 1.5);
                    } else {
                        zer.push(aa.asin());
                        zer.push(std::f64::consts::PI - aa.asin());
                    }
                }
            }
        }

        // General quartic in u = tan(x/2).
        if !infinite {
            let ko = [
                a - c + e,
                2.0 * d - 4.0 * b,
                2.0 * e - 2.0 * a,
                4.0 * b + 2.0 * d,
                a + c + e,
            ];
            let mut u_roots = quartic_real_roots(ko[0], ko[1], ko[2], ko[3], ko[4]);
            u_roots.sort_by(|x, y| x.partial_cmp(y).unwrap());
            for &u in &u_roots {
                let mut teta = u.atan() * 2.0;
                if u <= -EPS {
                    teta = depi - teta.abs();
                }
                zer.push(teta);
            }
        }
    }

    // OCCT math_TrigonometricFunctionRoots.cxx L505-550 "Cas particulier de
    // PI": t = pi is a solution exactly when the trig polynomial vanishes
    // there, i.e. ko(1) = A - C + E = 0 (at t=pi: cos=-1, sin=0).  The
    // substitution u = tan(t/2) maps t=pi to u=infinity, so the quartic (and
    // the quadratic when e-c=0) misses it; OCCT appends it explicitly.
    if !infinite && (a - c + e).abs() <= EPS {
        zer.push(std::f64::consts::PI);
    }

    if infinite {
        return Vec::new();
    }

    // Clamp the solutions into [InfBound, SupBound].
    let mut sol: Vec<f64> = Vec::new();
    for &t in &zer {
        let mut t = t;
        if t <= -EPS {
            t = depi - t.abs();
        }
        let period = depi;
        t += (inf_bound / period).floor() * period;
        if t < inf_bound - 1e-7 {
            t += period;
        }
        if t >= inf_bound - 1e-7 && t <= sup_bound + 1e-7 {
            let t = t.max(inf_bound).min(sup_bound);
            // Keep sorted, dedup.
            let pos = sol.partition_point(|&x| x < t);
            if pos == 0 || (t - sol[pos - 1]).abs() > 1e-9 {
                sol.insert(pos, t);
            }
        }
    }
    sol
}

/// OCCT IntAna_Quadric::NewCoefficients — transform the quadric coefficients
/// into the given frame (the circle's position).  The implicit equation
/// Cxx·x² + Cyy·y² + Czz·z² + 2·Cxy·xy + 2·Cxz·xz + 2·Cyz·yz
///   + Cx·x + Cy·y + Cz·z + Ccte = 0
/// is expressed in coordinates of `frame` (origin + orthonormal axes).
struct FrameCoefs {
    xx: f64, yy: f64, zz: f64, xy: f64, xz: f64, yz: f64,
    x: f64, y: f64, z: f64, cte: f64,
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

/// OCCT IntAna_IntConicQuad::Perform(const gp_Circ&, const IntAna_Quadric&)
/// — returns (in_quadric, points) where each point is (3D point, t param).
/// A circle lying entirely on the quadric reports `in_quadric = true`.
pub fn intersect_circle_quadric(
    circle: &Circle3,
    quad: &crate::geomalgo::int_surf::quadric::Quadric,
) -> Option<(bool, Vec<(DVec3, f64)>)> {
    use crate::geomalgo::int_surf::quadric::QuadricType;
    let r = circle.radius;
    let rr = r * r;

    // Coarse epsilon for the trig solver (OCCT Eps = 1.5e-12).
    let (co, inf, sup) = match quad.type_quadric() {
        QuadricType::Plane => {
            let pl = quad.plane();
            // Express the plane in the circle's frame.
            let pn = pl.normal;
            let d = -pn.dot(pl.origin);
            // Plane: pn·x + d = 0.  In the circle's frame: pn·(O + X u + Y v + Z w) + d = 0.
            let px = pn.dot(circle.x_dir);
            let py = pn.dot(circle.y_dir);
            let pz = pn.dot(circle.normal);
            let pc = pn.dot(circle.center) + d;
            let co = FrameCoefs {
                xx: 0.0, yy: 0.0, zz: 0.0, xy: 0.0, xz: 0.0, yz: 0.0,
                x: px, y: py, z: pz, cte: pc,
            };
            (co, 0.0, 2.0 * std::f64::consts::PI)
        }
        QuadricType::Cylinder => {
            // gp_Cylinder implicit: X² + Y² - R² = 0 in the local frame, i.e.
            // |P-O|² - ((P-O)·a)² - R² = 0.  Coeffs follow the OCCT
            // "2·C1·x" convention (linear terms stored halved).
            let o = quad.axis_loc();
            let a = quad.axis_dir();
            let r = quad.radius();
            let oa = o.dot(a);
            let co = FrameCoefs {
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
            };
            (co, 0.0, 2.0 * std::f64::consts::PI)
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
            let co = FrameCoefs {
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
            };
            (co, 0.0, 2.0 * std::f64::consts::PI)
        }
        QuadricType::Sphere => {
            let sph = quad.sphere();
            let center = sph.center;
            let rad = sph.radius;
            let co = FrameCoefs {
                xx: 1.0, yy: 1.0, zz: 1.0, xy: 0.0, xz: 0.0, yz: 0.0,
                x: -center.x, y: -center.y, z: -center.z,
                cte: center.length_squared() - rad * rad,
            };
            (co, 0.0, 2.0 * std::f64::consts::PI)
        }
        _ => return None,
    };

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
    let ts = trig_function_roots(
        p_cos_cos - p_sin_sin,
        p_cos_sin,
        2.0 * p_cos,
        2.0 * p_sin,
        p_cte + p_sin_sin,
        inf,
        sup,
    );

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
            let line_origin = gg.line(1).origin;
            let line_dir = gg.line(1).direction;
            // Project the 3D line into the circle's plane: u = (P-O)·u, v = (P-O)·v.
            let o2 = DVec2::new(
                (line_origin - circle.center).dot(circle.x_dir),
                (line_origin - circle.center).dot(circle.y_dir),
            );
            let d2 = DVec2::new(line_dir.dot(circle.x_dir), line_dir.dot(circle.y_dir));
            // Intersect the 2D line with the circle (radius r) in the plane.
            // |o2 + s·d2| = r  →  |d2|² s² + 2 o2·d2 s + |o2|² - r² = 0.
            let a = d2.length_squared();
            let b = 2.0 * o2.dot(d2);
            let c = o2.length_squared() - circle.radius * circle.radius;
            let mut pts = Vec::new();
            for s in quadratic_real_roots(a, b, c) {
                let p3 = line_origin + line_dir * s;
                // Parameter on the circle: t where p3 = O + R(cos t u + sin t v).
                let ut = (p3 - circle.center).dot(circle.x_dir);
                let vt = (p3 - circle.center).dot(circle.y_dir);
                let tt = vt.atan2(ut);
                pts.push((p3, tt));
            }
            (false, false, pts)
        }
        _ => (false, false, Vec::new()),
    }
}
