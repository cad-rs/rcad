//! Adaptor-layer methods (OCCT `GeomAdaptor_Curve` / `GeomAdaptor_Surface`) —
//! kernel gap completion for the TKFillet blend pipeline (Stage 1e):
//! continuity / resolution / intervals on `Curve3`, the ElSLib DN forms and
//! the BSpline / Bezier surface DN marshalling, and `u_resolution` /
//! `v_resolution` / `dn` on `Surface3`.
//!
//! Extracted verbatim from `geom/eval.rs` (project Rule 5 file-size split).
use crate::geom::*;
use crate::base::extrema_ext_elc::epsilon_of as standard_epsilon;

// =============================================================================
// Adaptor-layer methods (OCCT GeomAdaptor_Curve / GeomAdaptor_Surface) —
// kernel gap completion for the TKFillet blend pipeline (Stage 1e).
// =============================================================================

/// OCCT GeomAdaptor_Curve::LocalContinuity(U1, U2)
/// (GeomAdaptor_Curve.cxx L135-207) — computes the continuity of a BSpline
/// curve between the parameters U1 and U2: C(d - m) with d = degree,
/// m = max multiplicity of the knots between U1 and U2.
fn bspline_local_continuity(bs: &BSplineCurve3, u1: f64, u2: f64) -> crate::math::GeomAbsShape {
    use crate::math::bspl_lib::locate_parameter_main;
    use crate::math::GeomAbsShape;

    let (tk, tm) = bs.knots_mults();
    let nb = tk.len() as i32;
    let mut index1 = 0i32;
    let mut index2 = 0i32;
    let mut new_first = 0.0f64;
    let mut new_last = 0.0f64;
    let p_confusion = 1e-9; // OCCT Precision::PConfusion()
    locate_parameter_main(
        &tk,
        u1,
        bs.is_periodic,
        1,
        nb,
        &mut index1,
        &mut new_first,
        tk[0],
        tk[tk.len() - 1],
    );
    locate_parameter_main(
        &tk,
        u2,
        bs.is_periodic,
        1,
        nb,
        &mut index2,
        &mut new_last,
        tk[0],
        tk[tk.len() - 1],
    );
    if (new_first - crate::math::bspl_lib::at(&tk, index1 + 1)).abs() < p_confusion {
        if index1 < nb {
            index1 += 1;
        }
    }
    if (new_last - crate::math::bspl_lib::at(&tk, index2)).abs() < p_confusion {
        index2 -= 1;
    }
    let mult_max;
    // Handle periodic curves.
    if bs.is_periodic && index1 == nb {
        index1 = 1;
    }

    if (index2 - index1 <= 0) && !bs.is_periodic {
        mult_max = 100; // CN between 2 consecutive nodes
    } else {
        let mut m = crate::math::bspl_lib::ati(&tm, index1 + 1);
        for i in (index1 + 1)..=index2 {
            let ti = crate::math::bspl_lib::ati(&tm, i);
            if ti > m {
                m = ti;
            }
        }
        mult_max = bs.degree as i32 - m;
    }
    if mult_max <= 0 {
        GeomAbsShape::C0
    } else if mult_max == 1 {
        GeomAbsShape::C1
    } else if mult_max == 2 {
        GeomAbsShape::C2
    } else if mult_max == 3 {
        GeomAbsShape::C3
    } else {
        GeomAbsShape::CN
    }
}

impl Curve3 {
    /// OCCT GeomAdaptor_Curve::Continuity() (GeomAdaptor_Curve.cxx L330-355).
    pub fn continuity(&self) -> crate::math::GeomAbsShape {
        use crate::math::GeomAbsShape;
        match self {
            Curve3::BSpline(bs) => {
                bspline_local_continuity(bs, bs.first_parameter(), bs.last_parameter())
            }
            Curve3::Offset(off) => {
                // OCCT: GetBasisCurveContinuity() shifted one level down.
                let s = off.basis.continuity();
                match s {
                    GeomAbsShape::CN => GeomAbsShape::CN,
                    GeomAbsShape::C3 => GeomAbsShape::C2,
                    GeomAbsShape::C2 => GeomAbsShape::C1,
                    GeomAbsShape::C1 => GeomAbsShape::C0,
                    _ => GeomAbsShape::C0,
                }
            }
            _ => GeomAbsShape::CN,
        }
    }

    /// OCCT GeomAdaptor_Curve::Resolution(R3D) (GeomAdaptor_Curve.cxx
    /// L1116-1147).
    pub fn resolution(&self, r3d: f64) -> f64 {
        const PI_2: f64 = 2.0 * std::f64::consts::PI;
        match self {
            Curve3::Line(_) => r3d,
            Curve3::Circle(c) => {
                let r = c.radius;
                if r > r3d / 2.0 {
                    2.0 * (r3d / (2.0 * r)).asin()
                } else {
                    PI_2
                }
            }
            Curve3::Ellipse(e) => r3d / e.major_radius,
            Curve3::BSpline(bs) => bs.bsplclib_resolution(r3d),
            Curve3::Bezier(bz) => {
                // OCCT Geom_BezierCurve::Resolution delegates to
                // BSplCLib::Resolution with flat Bezier knots (0..0..1..1):
                // every knot span has length 1.
                let degree = bz.control_points.len().saturating_sub(1).max(1);
                let mut max_der = 0.0f64;
                for ii in 1..bz.control_points.len() {
                    let mut value = 0.0f64;
                    for kk in 0..3 {
                        let mut factor = bz.control_points[ii][kk] - bz.control_points[ii - 1][kk];
                        if factor < 0.0 {
                            factor = -factor;
                        }
                        value += factor;
                    }
                    value *= 1.0; // inverse = 1 / (FK[ii+Degree] - FK[ii]) == 1
                    if max_der < value {
                        max_der = value;
                    }
                }
                let max_derivative = max_der * degree as f64;
                if max_derivative > f64::MIN_POSITIVE {
                    r3d / max_derivative
                } else {
                    r3d / f64::MIN_POSITIVE
                }
            }
            _ => r3d * 0.01, // OCCT default: Precision::Parametric(R3D)
        }
    }

    /// OCCT GeomAdaptor_Curve::NbIntervals(S) (GeomAdaptor_Curve.cxx L371-460).
    pub fn nb_intervals(&self, s: crate::math::GeomAbsShape) -> usize {
        match self {
            Curve3::BSpline(bs) => {
                if (!bs.is_periodic && s <= self.continuity())
                    || s == crate::math::GeomAbsShape::C0
                {
                    return 1;
                }
                let a_degree = bs.degree as i32;
                let a_cont = match s {
                    crate::math::GeomAbsShape::C1 => 1,
                    crate::math::GeomAbsShape::C2 => 2,
                    crate::math::GeomAbsShape::C3 => 3,
                    crate::math::GeomAbsShape::CN => a_degree,
                    _ => panic!("Standard_DomainError: GeomAdaptor_Curve::NbIntervals()"),
                };
                let an_eps = (self.resolution(1e-7)).min(1e-9);
                bs.bsplclib_intervals(a_cont, bs.first_parameter(), bs.last_parameter(), an_eps, None)
            }
            Curve3::Offset(off) => {
                // OCCT offset branch (GeomAdaptor_Curve.cxx L413-455).
                use crate::math::GeomAbsShape;
                let base_s = match s {
                    GeomAbsShape::C0 => GeomAbsShape::C1,
                    GeomAbsShape::C1 => GeomAbsShape::C2,
                    GeomAbsShape::C2 => GeomAbsShape::C3,
                    _ => GeomAbsShape::CN,
                };
                let mut my_nb_intervals = 1usize;
                let i_nb_basis_int = off.basis.nb_intervals(base_s);
                if i_nb_basis_int > 1 {
                    let mut rdf_inter = Vec::with_capacity(1 + i_nb_basis_int);
                    off.basis.intervals(&mut rdf_inter, base_s);
                    let domain = off.basis.default_domain();
                    let (my_first, my_last) = (domain[0], domain[1]);
                    for i_int in 0..i_nb_basis_int {
                        if rdf_inter[i_int] > my_first && rdf_inter[i_int] < my_last {
                            my_nb_intervals += 1;
                        }
                    }
                }
                my_nb_intervals
            }
            _ => 1,
        }
    }

    /// OCCT GeomAdaptor_Curve::Intervals(T, S) (GeomAdaptor_Curve.cxx
    /// L466-561).  Fills `t` with the interval bounds.
    pub fn intervals(&self, t: &mut Vec<f64>, s: crate::math::GeomAbsShape) {
        match self {
            Curve3::BSpline(bs) => {
                if (!bs.is_periodic && s <= self.continuity())
                    || s == crate::math::GeomAbsShape::C0
                {
                    let domain = bs.default_domain();
                    let (first, last) = (domain[0], domain[1]);
                    t.clear();
                    t.push(first);
                    t.push(last);
                    return;
                }
                let a_degree = bs.degree as i32;
                let a_cont = match s {
                    crate::math::GeomAbsShape::C1 => 1,
                    crate::math::GeomAbsShape::C2 => 2,
                    crate::math::GeomAbsShape::C3 => 3,
                    crate::math::GeomAbsShape::CN => a_degree,
                    _ => panic!("Standard_DomainError: GeomAdaptor_Curve::Intervals()"),
                };
                let an_eps = (self.resolution(1e-7)).min(1e-9);
                bs.bsplclib_intervals(
                    a_cont,
                    bs.first_parameter(),
                    bs.last_parameter(),
                    an_eps,
                    Some(t),
                );
            }
            Curve3::Offset(off) => {
                use crate::math::GeomAbsShape;
                let base_s = match s {
                    GeomAbsShape::C0 => GeomAbsShape::C1,
                    GeomAbsShape::C1 => GeomAbsShape::C2,
                    GeomAbsShape::C2 => GeomAbsShape::C3,
                    _ => GeomAbsShape::CN,
                };
                let mut my_nb_intervals = 1usize;
                let domain = off.basis.default_domain();
                let (my_first, my_last) = (domain[0], domain[1]);
                let i_nb_basis_int = off.basis.nb_intervals(base_s);
                t.clear();
                t.push(my_first);
                if i_nb_basis_int > 1 {
                    let mut rdf_inter = Vec::with_capacity(1 + i_nb_basis_int);
                    off.basis.intervals(&mut rdf_inter, base_s);
                    for i_int in 0..i_nb_basis_int {
                        if rdf_inter[i_int] > my_first && rdf_inter[i_int] < my_last {
                            t.push(rdf_inter[i_int]);
                            my_nb_intervals += 1;
                        }
                    }
                }
                t.push(my_last);
            }
            _ => {
                let domain = self.default_domain();
                let (first, last) = (domain[0], domain[1]);
                t.clear();
                t.push(first);
                t.push(last);
            }
        }
    }
}

// --- ElSLib DN forms (OCCT ElSLib.cxx L169-500) for Surface3::dn ---

/// OCCT ElSLib::PlaneDN (ElSLib.cxx L169-180).
fn plane_dn(u_dir: DVec3, v_dir: DVec3, nu: i32, nv: i32) -> DVec3 {
    if nu == 0 && nv == 1 {
        v_dir
    } else if nu == 1 && nv == 0 {
        u_dir
    } else {
        DVec3::ZERO
    }
}

/// OCCT ElSLib::ConeDN (ElSLib.cxx L182-215).
#[allow(clippy::too_many_arguments)]
fn cone_dn(
    x_dir: DVec3,
    y_dir: DVec3,
    z_dir: DVec3,
    location: DVec3,
    radius: f64,
    s_angle: f64,
    u: f64,
    v: f64,
    nu: i32,
    nv: i32,
) -> DVec3 {
    let um = u + nu as f64 * (std::f64::consts::PI * 0.5);
    let mut xdir = x_dir * um.cos() + y_dir * um.sin();
    if nv == 0 {
        xdir *= radius + v * s_angle.sin();
        if nu == 0 {
            xdir += location;
        }
        xdir
    } else if nv == 1 {
        xdir *= s_angle.sin();
        if nu == 0 {
            xdir += z_dir * s_angle.cos();
        }
        xdir
    } else {
        DVec3::ZERO
    }
}

/// OCCT ElSLib::CylinderDN (ElSLib.cxx L217-265).
fn cylinder_dn(x_dir: DVec3, y_dir: DVec3, z_dir: DVec3, radius: f64, u: f64, nu: i32, nv: i32) -> DVec3 {
    if nu + nv < 1 || nu < 0 || nv < 0 {
        return DVec3::ZERO;
    }
    if nv == 0 {
        let r_cos_u = radius * u.cos();
        let r_sin_u = radius * u.sin();
        if (nu + 6) % 4 == 0 {
            x_dir * -r_cos_u + y_dir * -r_sin_u
        } else if (nu + 5) % 4 == 0 {
            x_dir * r_sin_u + y_dir * -r_cos_u
        } else if (nu + 3) % 4 == 0 {
            x_dir * -r_sin_u + y_dir * r_cos_u
        } else if nu % 4 == 0 {
            x_dir * r_cos_u + y_dir * r_sin_u
        } else {
            DVec3::ZERO
        }
    } else if nv == 1 && nu == 0 {
        z_dir
    } else {
        DVec3::ZERO
    }
}

/// OCCT ElSLib::SphereDN (ElSLib.cxx L267-365).
#[allow(clippy::too_many_arguments)]
fn sphere_dn(
    x_dir: DVec3,
    y_dir: DVec3,
    z_dir: DVec3,
    radius: f64,
    u: f64,
    v: f64,
    nu: i32,
    nv: i32,
) -> DVec3 {
    let is_odd = |i: i32| i % 2 != 0;
    if nu + nv < 1 || nu < 0 || nv < 0 {
        return DVec3::ZERO;
    }
    let cos_u = u.cos();
    let sin_u = u.sin();
    let r_cos_v = radius * v.cos();
    if nu == 0 {
        let r_sin_v = radius * v.sin();
        let (a1, a2, a3) = if is_odd(nv) {
            (-r_sin_v * cos_u, -r_sin_v * sin_u, r_cos_v)
        } else {
            (-r_cos_v * cos_u, -r_cos_v * sin_u, -r_sin_v)
        };
        let mut x = a1 * x_dir.x + a2 * y_dir.x + a3 * z_dir.x;
        let mut y = a1 * x_dir.y + a2 * y_dir.y + a3 * z_dir.y;
        let mut z = a1 * x_dir.z + a2 * y_dir.z + a3 * z_dir.z;
        if (nv + 2) % 4 != 0 && (nv + 3) % 4 != 0 {
            x = -x;
            y = -y;
            z = -z;
        }
        DVec3::new(x, y, z)
    } else if nv == 0 {
        let (a1, a2) = if is_odd(nu) {
            (-r_cos_v * sin_u, r_cos_v * cos_u)
        } else {
            (r_cos_v * cos_u, r_cos_v * sin_u)
        };
        let mut x = a1 * x_dir.x + a2 * y_dir.x;
        let mut y = a1 * x_dir.y + a2 * y_dir.y;
        let mut z = a1 * x_dir.z + a2 * y_dir.z;
        if (nu + 2) % 4 == 0 || (nu + 1) % 4 == 0 {
            x = -x;
            y = -y;
            z = -z;
        }
        DVec3::new(x, y, z)
    } else {
        let r_sin_v = radius * v.sin();
        let (a1, a2) = if is_odd(nu) {
            (-sin_u, cos_u)
        } else {
            (-cos_u, -sin_u)
        };
        let a3 = if is_odd(nv) { -r_sin_v } else { -r_cos_v };
        let mut x = (a1 * x_dir.x + a2 * y_dir.x) * a3;
        let mut y = (a1 * x_dir.y + a2 * y_dir.y) * a3;
        let mut z = (a1 * x_dir.z + a2 * y_dir.z) * a3;
        if ((nu + 2) % 4 != 0 && (nu + 3) % 4 != 0 && ((nv + 2) % 4 == 0 || (nv + 3) % 4 == 0))
            || (((nu + 2) % 4 == 0 || (nu + 3) % 4 == 0) && (nv + 2) % 4 != 0 && (nv + 3) % 4 != 0)
        {
            x = -x;
            y = -y;
            z = -z;
        }
        DVec3::new(x, y, z)
    }
}

/// OCCT ElSLib::TorusDN (ElSLib.cxx L367-500).
#[allow(clippy::too_many_arguments)]
fn torus_dn(
    x_dir: DVec3,
    y_dir: DVec3,
    z_dir: DVec3,
    major_radius: f64,
    minor_radius: f64,
    u: f64,
    v: f64,
    nu: i32,
    nv: i32,
) -> DVec3 {
    let is_odd = |i: i32| i % 2 != 0;
    if nu + nv < 1 || nu < 0 || nv < 0 {
        return DVec3::ZERO;
    }
    let cos_u = u.cos();
    let sin_u = u.sin();
    // OCCT OCC620: eps = 10 * (MinorRadius + MajorRadius) * RealEpsilon().
    let eps = 10.0 * (minor_radius + major_radius) * f64::EPSILON;
    if nv == 0 {
        let r = major_radius + minor_radius * v.cos();
        let (mut a1, mut a2) = if is_odd(nu) {
            (-r * sin_u, r * cos_u)
        } else {
            (-r * cos_u, -r * sin_u)
        };
        if a1.abs() <= eps {
            a1 = 0.0;
        }
        if a2.abs() <= eps {
            a2 = 0.0;
        }
        let mut x = a1 * x_dir.x + a2 * y_dir.x;
        let mut y = a1 * x_dir.y + a2 * y_dir.y;
        let mut z = a1 * x_dir.z + a2 * y_dir.z;
        if (nu + 2) % 4 != 0 && (nu + 3) % 4 != 0 {
            x = -x;
            y = -y;
            z = -z;
        }
        DVec3::new(x, y, z)
    } else if nu == 0 {
        let r_cos_v = minor_radius * v.cos();
        let r_sin_v = minor_radius * v.sin();
        let (mut a1, mut a2, mut a3) = if is_odd(nv) {
            (-r_sin_v * cos_u, -r_sin_v * sin_u, r_cos_v)
        } else {
            (-r_cos_v * cos_u, -r_cos_v * sin_u, -r_sin_v)
        };
        if a1.abs() <= eps {
            a1 = 0.0;
        }
        if a2.abs() <= eps {
            a2 = 0.0;
        }
        if a3.abs() <= eps {
            a3 = 0.0;
        }
        let mut x = a1 * x_dir.x + a2 * y_dir.x + a3 * z_dir.x;
        let mut y = a1 * x_dir.y + a2 * y_dir.y + a3 * z_dir.y;
        let mut z = a1 * x_dir.z + a2 * y_dir.z + a3 * z_dir.z;
        if (nv + 2) % 4 != 0 && (nv + 3) % 4 != 0 {
            x = -x;
            y = -y;
            z = -z;
        }
        DVec3::new(x, y, z)
    } else if is_odd(nu) && is_odd(nv) {
        let rsin_v = minor_radius * v.sin();
        let (mut a1, mut a2) = (rsin_v * sin_u, -rsin_v * cos_u);
        if a1.abs() <= eps {
            a1 = 0.0;
        }
        if a2.abs() <= eps {
            a2 = 0.0;
        }
        DVec3::new(
            a1 * x_dir.x + a2 * y_dir.x,
            a1 * x_dir.y + a2 * y_dir.y,
            a1 * x_dir.z + a2 * y_dir.z,
        )
    } else if !is_odd(nu) && !is_odd(nv) {
        let rcos_v = minor_radius * v.cos();
        let (mut a1, mut a2) = (rcos_v * cos_u, rcos_v * sin_u);
        if a1.abs() <= eps {
            a1 = 0.0;
        }
        if a2.abs() <= eps {
            a2 = 0.0;
        }
        DVec3::new(
            a1 * x_dir.x + a2 * y_dir.x,
            a1 * x_dir.y + a2 * y_dir.y,
            a1 * x_dir.z + a2 * y_dir.z,
        )
    } else {
        DVec3::ZERO
    }
}

/// Radius of the V-iso circle of a cone at parameter v
/// (OCCT Geom_ConicalSurface::VIso).
fn cone_iso_radius_at(c: &crate::geom::ConicalSurface, v: f64) -> f64 {
    c.radius + v * c.half_angle_rad.tan()
}

/// Direction-wise pole bound for the U resolution of a BSpline surface
/// (OCCT Geom_BSplineSurface::Resolution / BSplSLib::Resolution — the
/// two-direction pole-difference bound, applied per direction).
fn bspline_surface_u_resolution(bs: &BSplineSurface, tolerance_3d: f64) -> f64 {
    let degree = bs.degree_u;
    let mut max_der = 0.0f64;
    let nu = bs.control_points.len();
    let nv = if nu > 0 { bs.control_points[0].len() } else { 0 };
    for j in 0..nv {
        for i in 1..nu {
            let inverse = 1.0 / (bs.knots_u[i + degree] - bs.knots_u[i]);
            let mut value = 0.0f64;
            for kk in 0..3 {
                let mut factor = bs.control_points[i][j][kk] - bs.control_points[i - 1][j][kk];
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
            }
            value *= inverse;
            if max_der < value {
                max_der = value;
            }
        }
    }
    let max_derivative = max_der * degree as f64;
    if max_derivative > f64::MIN_POSITIVE {
        tolerance_3d / max_derivative
    } else {
        tolerance_3d / f64::MIN_POSITIVE
    }
}

/// Direction-wise pole bound for the V resolution of a BSpline surface
/// (see bspline_surface_u_resolution).
fn bspline_surface_v_resolution(bs: &BSplineSurface, tolerance_3d: f64) -> f64 {
    let degree = bs.degree_v;
    let mut max_der = 0.0f64;
    let nu = bs.control_points.len();
    let nv = if nu > 0 { bs.control_points[0].len() } else { 0 };
    for i in 0..nu {
        for j in 1..nv {
            let inverse = 1.0 / (bs.knots_v[j + degree] - bs.knots_v[j]);
            let mut value = 0.0f64;
            for kk in 0..3 {
                let mut factor = bs.control_points[i][j][kk] - bs.control_points[i][j - 1][kk];
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
            }
            value *= inverse;
            if max_der < value {
                max_der = value;
            }
        }
    }
    let max_derivative = max_der * degree as f64;
    if max_derivative > f64::MIN_POSITIVE {
        tolerance_3d / max_derivative
    } else {
        tolerance_3d / f64::MIN_POSITIVE
    }
}

/// OCCT `Geom_BSplineSurface::EvalDN` (Geom_BSplineSurface_1.cxx L279-313) —
/// `BSplSLib::DN(U, V, Nu, Nv, 0, 0, myPoles, Weights(), myUFlatKnots,
/// myVFlatKnots, NoMults, NoMults, myUDeg, myVDeg, myURational, myVRational,
/// myUPeriodic, myVPeriodic, Vn)`.
///
/// `Weights()` (Geom_BSplineSurface_1.cxx L914-921) is the OCCT null pointer
/// exactly when the surface is non-rational — mirrored by the `Option` below,
/// which [`crate::geom::eval_b::bspl_slib_dn`] forwards to `PrepareEval`.
pub(crate) fn bspline_surface_dn(bs: &BSplineSurface, u: f64, v: f64, nu: i32, nv: i32) -> DVec3 {
    let rational = bs.is_rational_u() || bs.is_rational_v();
    let weights: Option<&[Vec<f64>]> = if rational { Some(&bs.weights) } else { None };
    crate::geom::eval_b::bspl_slib_dn(
        u,
        v,
        nu,
        nv,
        // OCCT L296-297: UIndex = VIndex = 0 (the LocateParameter arms).
        0,
        0,
        bs.degree_u as i32,
        bs.degree_v as i32,
        bs.is_rational_u(),
        bs.is_rational_v(),
        bs.is_periodic_u,
        bs.is_periodic_v,
        &bs.control_points,
        weights,
        &bs.knots_u,
        &bs.knots_v,
    )
}

/// The rcad marshalling of the OCCT `Geom_BezierSurface` evaluation frame
/// shared by `EvalD0` / `EvalD1` / `EvalD2` / `EvalDN` (Geom_BezierSurface.cxx
/// L1416-1724):
///   - `UDegree = myPoles.ColLength() - 1` (`NbUPoles() - 1`),
///     `VDegree = myPoles.RowLength() - 1` (`NbVPoles() - 1`);
///   - `UIndex = VIndex = 0` and `UPer = VPer = false`;
///   - `UKnots() == VKnots()` are the compact `{0, 1}` form with
///     `UMultiplicities() == VMultiplicities() == {Degree+1, Degree+1}`, which
///     the rcad re-host carries as the equivalent clamped flat sequence (see
///     [`crate::geom::eval_b`]);
///   - `Weights()` is `&myWeights` for a rational Bezier and
///     `BSplSLib::NoWeights()` otherwise.
struct BezierCallFrame {
    u_degree: i32,
    v_degree: i32,
    rational_u: bool,
    rational_v: bool,
    knots_u: Vec<f64>,
    knots_v: Vec<f64>,
}

impl BezierCallFrame {
    fn new(bez: &BezierSurface) -> Self {
        let u_degree = bez.control_points.len().saturating_sub(1) as i32;
        let v_degree = bez
            .control_points
            .first()
            .map(|row| row.len())
            .unwrap_or(0)
            .saturating_sub(1) as i32;
        BezierCallFrame {
            u_degree,
            v_degree,
            rational_u: bezier_is_rational_u(bez),
            rational_v: bezier_is_rational_v(bez),
            knots_u: crate::geom::eval_b::bezier_flat_knots(u_degree as usize),
            knots_v: crate::geom::eval_b::bezier_flat_knots(v_degree as usize),
        }
    }

    fn weights<'a>(&self, bez: &'a BezierSurface) -> Option<&'a [Vec<f64>]> {
        if self.rational_u || self.rational_v {
            Some(&bez.weights)
        } else {
            None
        }
    }
}

/// OCCT `Geom_BezierSurface::EvalD0` (Geom_BezierSurface.cxx L1416-1466) —
/// `BSplSLib::D0(U, V, 1, 1, myPoles, Weights(), UKnots(), UKnots(),
/// &UMultiplicities(), &VMultiplicities(), ColLength - 1, RowLength - 1,
/// myURational, myVRational, false, false, P)`.
pub(crate) fn bezier_surface_d0(bez: &BezierSurface, u: f64, v: f64) -> DVec3 {
    let f = BezierCallFrame::new(bez);
    crate::geom::eval_b::bspl_slib_d0(
        u,
        v,
        0,
        0,
        f.u_degree,
        f.v_degree,
        f.rational_u,
        f.rational_v,
        false,
        false,
        &bez.control_points,
        f.weights(bez),
        &f.knots_u,
        &f.knots_v,
    )
}

/// OCCT `Geom_BezierSurface::EvalD1` (Geom_BezierSurface.cxx L1470-1524) —
/// `BSplSLib::D1` with the [`BezierCallFrame`] arguments.
pub(crate) fn bezier_surface_d1(bez: &BezierSurface, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
    let f = BezierCallFrame::new(bez);
    crate::geom::eval_b::bspl_slib_d1(
        u,
        v,
        0,
        0,
        f.u_degree,
        f.v_degree,
        f.rational_u,
        f.rational_v,
        false,
        false,
        &bez.control_points,
        f.weights(bez),
        &f.knots_u,
        &f.knots_v,
    )
}

/// OCCT `Geom_BezierSurface::EvalD2` (Geom_BezierSurface.cxx L1528-1590) —
/// `BSplSLib::D2` with the [`BezierCallFrame`] arguments, returned in the
/// rcad `SurfaceEval::derivatives2` order `(P, dP/du, dP/dv, d2P/du2,
/// d2P/dudv, d2P/dv2)`.
pub(crate) fn bezier_surface_d2(
    bez: &BezierSurface,
    u: f64,
    v: f64,
) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
    let f = BezierCallFrame::new(bez);
    crate::geom::eval_b::bspl_slib_d2(
        u,
        v,
        0,
        0,
        f.u_degree,
        f.v_degree,
        f.rational_u,
        f.rational_v,
        false,
        false,
        &bez.control_points,
        f.weights(bez),
        &f.knots_u,
        &f.knots_v,
    )
}

/// OCCT `Geom_BezierSurface::EvalDN` (Geom_BezierSurface.cxx L1666-1724) —
/// `BSplSLib::DN(U, V, Nu, Nv, 0, 0, myPoles, Weights(), UKnots(), UKnots(),
/// &UMultiplicities(), &VMultiplicities(), ColLength - 1, RowLength - 1,
/// myURational, myVRational, false, false, Derivative)`.
pub(crate) fn bezier_surface_dn(bez: &BezierSurface, u: f64, v: f64, nu: i32, nv: i32) -> DVec3 {
    let f = BezierCallFrame::new(bez);
    crate::geom::eval_b::bspl_slib_dn(
        u,
        v,
        nu,
        nv,
        0,
        0,
        f.u_degree,
        f.v_degree,
        f.rational_u,
        f.rational_v,
        // OCCT L1697-1698: the Bezier is never periodic.
        false,
        false,
        &bez.control_points,
        f.weights(bez),
        &f.knots_u,
        &f.knots_v,
    )
}

/// OCCT `Geom_BezierSurface::EvalD0` derivation of `myURational`
/// (Geom_BezierSurface.cxx L458-469): `myURational` is set when two weights of
/// the same V-row differ by more than `Epsilon(abs(w))` of the left one.
fn bezier_is_rational_u(bez: &BezierSurface) -> bool {
    for row in &bez.weights {
        for j in 0..row.len().saturating_sub(1) {
            if (row[j] - row[j + 1]).abs() > standard_epsilon(row[j].abs()) {
                return true;
            }
        }
    }
    false
}

/// OCCT `Geom_BezierSurface::EvalD1` derivation of `myVRational`
/// (Geom_BezierSurface.cxx L446-457): the V direction counterpart over the
/// weights of each U-column.
fn bezier_is_rational_v(bez: &BezierSurface) -> bool {
    let n_u = bez.weights.len();
    let n_v = bez.weights.first().map(|r| r.len()).unwrap_or(0);
    for j in 0..n_v {
        for i in 0..n_u.saturating_sub(1) {
            let w = bez.weights[i][j];
            if (w - bez.weights[i + 1][j]).abs() > standard_epsilon(w.abs()) {
                return true;
            }
        }
    }
    false
}

impl Surface3 {
    /// OCCT GeomAdaptor_Surface::UResolution(R3d)
    /// (GeomAdaptor_Surface.cxx L1818-1898).  Promoted from the local
    /// GeomAdaptor translation in bop/algo/wire_splitter.rs.
    pub fn u_resolution(&self, r3d: f64) -> f64 {
        let mut res = 0.0f64;
        match self {
            Surface3::LinearExtrusion(ext) => {
                return ext.profile.resolution(r3d);
            }
            Surface3::Torus(t) => {
                let r = t.major_radius + t.minor_radius;
                if r > crate::CONFUSION {
                    res = r3d / (2.0 * r);
                }
            }
            Surface3::Sphere(s) => {
                let r = s.radius;
                if r > crate::CONFUSION {
                    res = r3d / (2.0 * r);
                }
            }
            Surface3::Cylinder(c) => {
                let r = c.radius;
                if r > crate::CONFUSION {
                    res = r3d / (2.0 * r);
                }
            }
            Surface3::Cone(c) => {
                let domain = self.default_domain();
                let (v_first, v_last) = (domain[2], domain[3]);
                if v_last - v_first > 1e10 {
                    // Not truly bounded => unknown resolution.
                    return r3d * 0.01; // Precision::Parametric(R3d)
                }
                // OCCT: R = max(radius of VIso(myVLast), radius of VIso(myVFirst)).
                let r1 = cone_iso_radius_at(c, v_last);
                let r2 = cone_iso_radius_at(c, v_first);
                let r = r1.max(r2);
                return if r > crate::CONFUSION { r3d / r } else { 0.0 };
            }
            Surface3::Plane(_) => {
                return r3d;
            }
            Surface3::BSpline(bs) => {
                return bspline_surface_u_resolution(bs, r3d);
            }
            Surface3::Offset(off) => {
                return off.basis.u_resolution(r3d);
            }
            _ => {
                return r3d * 0.01; // Precision::Parametric(R3d)
            }
        }
        if res <= 1.0 {
            2.0 * res.asin()
        } else {
            2.0 * std::f64::consts::PI
        }
    }

    /// OCCT GeomAdaptor_Surface::VResolution(R3d)
    /// (GeomAdaptor_Surface.cxx L1900-1958).  Promoted from the local
    /// GeomAdaptor translation in bop/algo/wire_splitter.rs.
    pub fn v_resolution(&self, r3d: f64) -> f64 {
        let mut res = 0.0f64;
        match self {
            Surface3::Revolution(rev) => {
                return rev.profile.resolution(r3d);
            }
            Surface3::Torus(t) => {
                let r = t.minor_radius;
                if r > crate::CONFUSION {
                    res = r3d / (2.0 * r);
                }
            }
            Surface3::Sphere(s) => {
                let r = s.radius;
                if r > crate::CONFUSION {
                    res = r3d / (2.0 * r);
                }
            }
            Surface3::LinearExtrusion(_)
            | Surface3::Cylinder(_)
            | Surface3::Cone(_)
            | Surface3::Plane(_) => {
                return r3d;
            }
            Surface3::BSpline(bs) => {
                return bspline_surface_v_resolution(bs, r3d);
            }
            Surface3::Offset(off) => {
                return off.basis.v_resolution(r3d);
            }
            _ => {
                return r3d * 0.01; // Precision::Parametric(R3d)
            }
        }
        if res <= 1.0 {
            2.0 * res.asin()
        } else {
            2.0 * std::f64::consts::PI
        }
    }

    /// OCCT GeomAdaptor_Surface::DN(U, V, Nu, Nv)
    /// (GeomAdaptor_Surface.cxx L1697-1814) — dispatch to the per-type
    /// derivative: ElSLib::DN forms for the quadrics (L1796-1805),
    /// Geom_BSplineSurface::EvalDN / Geom_BezierSurface::EvalDN (both
    /// `BSplSLib::DN`) for the polynomial kinds (L1731-1752 /
    /// L1807-1813), the offset adaptor's DN (L1782-1794), the swept surfaces
    /// (L1754 / L1768) and — through the default arm at L1807-1813 — the
    /// `GeomEval` ellipsoid / circular helicoid.  The `Trimmed` wrapper is the
    /// `Geom_RectangularTrimmedSurface::EvalDN` delegation to the basis
    /// surface (Geom_RectangularTrimmedSurface.cxx L419-429).  The rcad-only
    /// Ruled / Coons / Pipe / TriBezier kinds have no OCCT
    /// `Geom_Surface::EvalDN` counterpart and raise the explicit gap below.
    pub fn dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> DVec3 {
        match self {
            Surface3::Plane(p) => plane_dn(p.u_dir, p.v_dir, nu, nv),
            Surface3::Cylinder(c) => {
                let x_dir = c.ref_dir;
                let y_dir = c.y_dir.unwrap_or_else(|| c.axis.cross(c.ref_dir));
                cylinder_dn(x_dir, y_dir, c.axis, c.radius, u, nu, nv)
            }
            Surface3::Cone(c) => {
                let x_dir = c.ref_dir;
                let y_dir = c.axis.cross(c.ref_dir);
                cone_dn(
                    x_dir, y_dir, c.axis, c.apex, c.radius, c.half_angle_rad, u, v, nu, nv,
                )
            }
            Surface3::Sphere(s) => {
                let x_dir = s.ref_dir;
                let y_dir = s.axis.cross(s.ref_dir);
                sphere_dn(x_dir, y_dir, s.axis, s.radius, u, v, nu, nv)
            }
            Surface3::Torus(t) => {
                let x_dir = t.ref_dir;
                let y_dir = t.axis.cross(t.ref_dir);
                torus_dn(
                    x_dir,
                    y_dir,
                    t.axis,
                    t.major_radius,
                    t.minor_radius,
                    u,
                    v,
                    nu,
                    nv,
                )
            }
            Surface3::BSpline(bs) => bspline_surface_dn(bs, u, v, nu, nv),
            Surface3::Bezier(bez) => bezier_surface_dn(bez, u, v, nu, nv),
            Surface3::Offset(of) => {
                crate::geom::offset_surface_utils::offset_payload_eval_dn(of, u, v, nu, nv)
            }
            Surface3::LinearExtrusion(le) => {
                crate::geom::extrusion_utils::linear_extrusion_eval_dn(le, u, v, nu, nv)
            }
            Surface3::Ellipsoid(el) => crate::geom::eval_c::ellipsoid_eval_dn(el, u, v, nu, nv),
            Surface3::Helicoid(h) => crate::geom::eval_c::helicoid_eval_dn(h, u, v, nu, nv),
            Surface3::Revolution(rev) => {
                crate::geom::revolution_utils::revolution_eval_dn(rev, u, v, nu, nv)
            }
            // OCCT Geom_RectangularTrimmedSurface::EvalDN
            // (Geom_RectangularTrimmedSurface.cxx L419-429).
            Surface3::Trimmed(t) => {
                // OCCT L424-427: if (Nu + Nv < 1 || Nu < 0 || Nv < 0) throw
                // Geom_UndefinedDerivative.
                assert!(
                    nu + nv >= 1 && nu >= 0 && nv >= 0,
                    "Geom_UndefinedDerivative: Geom_RectangularTrimmedSurface::EvalDN"
                );
                // OCCT L428: return basisSurf->EvalDN(U, V, Nu, Nv).
                t.basis.dn(u, v, nu, nv)
            }
            _ => {
                panic!(
                    "GAP: GeomAdaptor_Surface::DN (GeomAdaptor_Surface.cxx L1697-1814): the \
                     rcad GeomAdaptor_Surface DN engine covers ElSLib surfaces \
                     (L1796-1805), Geom_BSplineSurface (L1731-1752), Geom_BezierSurface \
                     (L1807-1813; Geom_BezierSurface::EvalDN), the offset adaptor DN \
                     (L1782-1794), the extrusion (L1754), the revolution (L1768), the \
                     GeomEval ellipsoid / circular helicoid (the L1807-1813 default arm) \
                     and Geom_RectangularTrimmedSurface (the Trimmed arm above); the \
                     rcad-only Ruled / Coons / Pipe / TriBezier kinds have no OCCT \
                     Geom_Surface::EvalDN counterpart"
                );
            }
        }
    }
}
