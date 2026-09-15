//! Parametric evaluation traits mirroring the OCCT `Geom_Curve` /
//! `Geom_Surface` / `Geom2d_Curve` abstract layers and their OCCT intermediate
//! classes (`Geom_Conic`, `Geom_BoundedCurve`, `Geom_ElementarySurface`,
//! `Geom_BoundedSurface`, `Geom_SweptSurface`, `Geom2d_Conic`,
//! `Geom2d_BoundedCurve`).
//!
//! Extracted verbatim from `geom/mod.rs` (project Rule 5 file-size split);
//! all public paths stay stable through the `pub use` re-exports in `mod.rs`.
use glam::{DVec2, DVec3};

use super::Curve3;

/// Parametric evaluation of a 3D curve: `t -> Point3`.
///
/// Mirrors OCCT `Geom_Curve::Value(t)` / `D1(t)`.
pub trait CurveEval {
    /// Point on the curve at parameter `t`.
    fn point_at(&self, t: f64) -> DVec3;
    /// Unit tangent vector at parameter `t`.
    fn tangent_at(&self, t: f64) -> DVec3;
    /// First derivative (non-unit velocity vector) at parameter `t`.
    /// OCCT-aligned: Extrema_CurveTool::D1 — used by Extrema_LocateExtPC
    /// for the Newton method solving g(u) = (P-C)·C' = 0.
    /// Default: 6-point central difference (accurate to O(h⁴)).
    fn derivative_at(&self, t: f64) -> DVec3 {
        // 6-point stencil: f'(t) ≈ [f(t-2h)-8f(t-h)+8f(t+h)-f(t+2h)] / (12h)
        let h = 1e-6;
        let fp2 = self.point_at(t + 2.0 * h);
        let fp1 = self.point_at(t + h);
        let fm1 = self.point_at(t - h);
        let fm2 = self.point_at(t - 2.0 * h);
        (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * h)
    }
    /// Natural parameter domain `[t_min, t_max]`.
    /// Unbounded producers (line / parabola / hyperbola) return
    /// `±Precision::Infinite()` (±2e100, OCCT Geom_Line.cxx L137/L144,
    /// Geom_Parabola.cxx L114/L128, Geom_Hyperbola.cxx L94/L101); circles and
    /// ellipses use `[0, 2π]`.
    fn default_domain(&self) -> [f64; 2];

    /// OCCT-aligned: IsClosed — true for periodic curves where start == end (circle, ellipse).
    fn is_closed(&self) -> bool {
        false
    }

    /// OCCT-aligned: IsPeriodic — true for curves with cyclic parameter (circle, ellipse).
    fn is_periodic(&self) -> bool {
        false
    }

    /// OCCT-aligned: ReversedParameter(t) — parameter of the same geometric point
    /// when traversed in the opposite direction.
    /// For periodic curves (circle): `period - t` (mod period).
    /// For non-periodic curves: `-t` (line).
    fn reversed_parameter(&self, t: f64) -> f64 {
        t
    }

    /// OCCT-aligned: D2(t) — second derivative d²P/dt² at parameter `t`.
    ///
    /// Default: 5-point central difference of `point_at`:
    /// ```text
    /// f''(t) = [-f(t+2h) + 16f(t+h) - 30f(t) + 16f(t-h) - f(t-2h)] / (12h²)
    /// ```
    /// Override with analytic formula when available (Line3 → 0, Circle3 → -R·N, etc.).
    fn derivative2_at(&self, t: f64) -> DVec3 {
        let h = 1e-4;
        let fp2 = self.point_at(t + 2.0 * h);
        let fp1 = self.point_at(t + h);
        let f = self.point_at(t);
        let fm1 = self.point_at(t - h);
        let fm2 = self.point_at(t - 2.0 * h);
        (-fp2 + 16.0 * fp1 - 30.0 * f + 16.0 * fm1 - fm2) / (12.0 * h * h)
    }

    /// OCCT-aligned: D3(t) — third derivative d³P/dt³ at parameter `t`.
    ///
    /// Default: central difference of `derivative2_at`.
    fn derivative3_at(&self, t: f64) -> DVec3 {
        let h = 1e-4;
        (self.derivative2_at(t + h) - self.derivative2_at(t - h)) / (2.0 * h)
    }

    /// Signed curvature at parameter `t`.
    ///
    /// For 3D curves: `k = |r' × r''| / |r'|³`.
    /// Returns 0 when the velocity is zero (degenerate point).
    /// OCCT-aligned: computed via `D1 × D2 / |D1|³`.
    fn curvature_at(&self, t: f64) -> f64 {
        let d1 = self.derivative_at(t);
        let d2 = self.derivative2_at(t);
        let speed = d1.length();
        if speed < 1e-15 {
            return 0.0;
        }
        d1.cross(d2).length() / (speed * speed * speed)
    }

    /// OCCT-aligned: `Geom_Curve::TransformedParameter(t, T)`.
    ///
    /// Returns the parameter on the transformed curve corresponding to parameter
    /// `t` on the original curve after transformation `T`. Used for curve-on-surface
    /// evaluation after a `TopLoc_Location` transform.
    ///
    /// Default: identity (`t` unchanged — correct for isometric transformations).
    /// For scaling transformations, override with `t / scale_factor`.
    fn transformed_parameter(&self, t: f64) -> f64 {
        t
    }

    /// OCCT-aligned: `Geom_Curve::ParametricTransformation(T)`.
    ///
    /// Returns the scale factor for parametric transformation. Used when computing
    /// parametric tolerance after a `TopLoc_Location` transformation.
    ///
    /// Default: 1.0 (correct for isometric / uniform scaling).
    fn parametric_transformation(&self) -> f64 {
        1.0
    }
}

/// Parametric evaluation of a 3D surface: `(u, v) -> Point3`.
///
/// Mirrors OCCT `Geom_Surface::Value(u, v)`.
pub trait SurfaceEval {
    /// Point on the surface at parameter `(u, v)`.
    fn point_at(&self, u: f64, v: f64) -> DVec3;
    /// Outward unit normal at parameter `(u, v)`.
    fn normal_at(&self, u: f64, v: f64) -> DVec3;
    /// Natural parameter domain `[u_min, u_max, v_min, v_max]`.
    fn default_domain(&self) -> [f64; 4];
    /// First partial derivatives `(point, dP/du, dP/dv)` at `(u, v)`.
    /// Default: finite-difference approximation (2-point, 1e-6 step).
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let eps = 1e-6;
        let p = self.point_at(u, v);
        let pu = self.point_at(u + eps, v);
        let pv = self.point_at(u, v + eps);
        (p, (pu - p) / eps, (pv - p) / eps)
    }

    /// OCCT-aligned: IsUClosed / IsVClosed — true if surface is closed in that direction.
    fn is_u_closed(&self) -> bool {
        false
    }
    fn is_v_closed(&self) -> bool {
        false
    }

    /// OCCT-aligned: IsUPeriodic / IsVPeriodic — true if parameter is cyclic.
    fn is_u_periodic(&self) -> bool {
        false
    }
    fn is_v_periodic(&self) -> bool {
        false
    }

    /// OCCT-aligned: UReversedParameter / VReversedParameter — parameter in reverse direction.
    fn u_reversed_parameter(&self, t: f64) -> f64 {
        t
    }
    fn v_reversed_parameter(&self, t: f64) -> f64 {
        t
    }

    /// OCCT-aligned: D2(u,v) — second-order partial derivatives.
    ///
    /// Returns `(P, dP/du, dP/dv, d²P/du², d²P/dudv, d²P/dv²)`.
    ///
    /// Default: 3-point central difference for each second-order term:
    /// ```text
    /// Puu  = [P(u+h,v) - 2P(u,v) + P(u-h,v)] / h²
    /// Pvv  = [P(u,v+h) - 2P(u,v) + P(u,v-h)] / h²
    /// Puv  = [P(u+h,v+h) - P(u+h,v-h) - P(u-h,v+h) + P(u-h,v-h)] / (4h²)
    /// ```
    /// Override with analytic formula for known surface types.
    fn derivatives2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        let h = 1e-5;
        let (p, pu, pv) = self.derivatives(u, v);
        let p_up = self.point_at(u + h, v);
        let p_um = self.point_at(u - h, v);
        let p_vp = self.point_at(u, v + h);
        let p_vm = self.point_at(u, v - h);
        let p_pp = self.point_at(u + h, v + h);
        let p_pm = self.point_at(u + h, v - h);
        let p_mp = self.point_at(u - h, v + h);
        let p_mm = self.point_at(u - h, v - h);
        let h2 = h * h;
        let puu = (p_up - 2.0 * p + p_um) / h2;
        let pvv = (p_vp - 2.0 * p + p_vm) / h2;
        let puv = (p_pp - p_pm - p_mp + p_mm) / (4.0 * h2);
        (p, pu, pv, puu, puv, pvv)
    }
}

/// Parametric evaluation of a 2D curve (PCurve): `t -> Point2`.
///
/// ✅ OCCT-aligned: Geom2d_Curve (provides Value/D0/D1/D2/D3 + domain).
pub trait Curve2dEval {
    /// Point on the 2D curve at parameter `t` (OCCT `Value(t)` / `D0(t)`).
    fn point_at(&self, t: f64) -> DVec2;

    /// Unit tangent vector at parameter `t` (OCCT `D1(t).normalize()`).
    /// Default: finite-difference approximation.
    fn tangent_at(&self, t: f64) -> DVec2 {
        let eps = 1e-7;
        let dp = self.point_at(t + eps) - self.point_at(t - eps);
        dp.normalize_or_zero()
    }

    /// First derivative (velocity) vector at parameter `t` (OCCT `D1(t)`).
    /// Default: central-difference approximation.
    fn derivative_at(&self, t: f64) -> DVec2 {
        let eps = 1e-7;
        (self.point_at(t + eps) - self.point_at(t - eps)) / (2.0 * eps)
    }

    /// Natural parameter domain `[t_min, t_max]` (OCCT `FirstParameter() / LastParameter()`).
    fn default_domain(&self) -> [f64; 2] {
        [f64::NEG_INFINITY, f64::INFINITY]
    }

    /// OCCT-aligned: IsClosed — true for closed curves (circle, ellipse).
    fn is_closed(&self) -> bool {
        false
    }

    /// OCCT-aligned: IsPeriodic — true for periodic curves (circle, ellipse).
    fn is_periodic(&self) -> bool {
        false
    }

    /// OCCT-aligned: ReversedParameter(t) — parameter of the same point in reverse.
    fn reversed_parameter(&self, t: f64) -> f64 {
        t
    }

    /// OCCT-aligned: D2(t) — second derivative d²P/dt² at parameter `t`.
    /// Default: 5-point central difference of `point_at`.
    fn derivative2_at(&self, t: f64) -> DVec2 {
        let h = 1e-4;
        let fp2 = self.point_at(t + 2.0 * h);
        let fp1 = self.point_at(t + h);
        let f = self.point_at(t);
        let fm1 = self.point_at(t - h);
        let fm2 = self.point_at(t - 2.0 * h);
        (-fp2 + 16.0 * fp1 - 30.0 * f + 16.0 * fm1 - fm2) / (12.0 * h * h)
    }

    /// OCCT-aligned: D3(t) — third derivative at parameter `t`.
    /// Default: central difference of `derivative2_at`.
    fn derivative3_at(&self, t: f64) -> DVec2 {
        let h = 1e-4;
        (self.derivative2_at(t + h) - self.derivative2_at(t - h)) / (2.0 * h)
    }

    /// OCCT-aligned: EvalDN(t, N) — the N-th derivative, N >= 1
    /// (Geom2d_Curve::EvalDN).  Default: the exact D1..D3 overrides for
    /// N <= 3; for N > 3 a central finite difference of `derivative3_at`
    /// (the per-type DN engines beyond D3 are not re-hosted yet; reached
    /// only by the offset `AdjustDerivative` Taylor climb at singular
    /// points).
    fn derivative_n_at(&self, t: f64, n: i32) -> DVec2 {
        match n {
            1 => self.derivative_at(t),
            2 => self.derivative2_at(t),
            3 => self.derivative3_at(t),
            _ => {
                let h = 1e-4;
                (self.derivative3_at(t + h) - self.derivative3_at(t - h)) / (2.0 * h)
            }
        }
    }

    /// Signed curvature at parameter `t` in the 2D plane.
    ///
    /// `k = (x'y'' - y'x'') / (x'² + y'²)^(3/2)`.
    /// Positive = counter-clockwise turning. Returns 0 at degenerate points.
    fn curvature_at(&self, t: f64) -> f64 {
        let d1 = self.derivative_at(t);
        let d2 = self.derivative2_at(t);
        let speed_sq = d1.length_squared();
        if speed_sq < 1e-30 {
            return 0.0;
        }
        (d1.x * d2.y - d1.y * d2.x) / (speed_sq * speed_sq.sqrt())
    }
}

/// OCCT-aligned: `Geom_Conic` intermediate abstract class.
///
/// Groups all conic curves (Circle, Ellipse, Hyperbola, Parabola) with their
/// shared geometric properties: position frame, eccentricity, and local axes.
///
/// In OCCT boolean algorithms this grouping is used for dynamic type checks
/// like `IsKind(STANDARD_TYPE(Geom_Conic))` before calling conic-specific
/// methods such as `XAxis()` / `YAxis()`.
pub trait ConicEval: CurveEval {
    /// The reference point of the conic:
    ///   Circle/Ellipse/Hyperbola → center
    ///   Parabola → vertex
    fn position(&self) -> DVec3;

    /// The normal of the conic's plane (gp_Ax2::Direction).
    fn normal(&self) -> DVec3;

    /// OCCT-aligned: eccentricity of the conic.
    ///   Circle → 0
    ///   Ellipse → sqrt(1 - (b/a)²)  (0 < e < 1)
    ///   Hyperbola → sqrt(1 + (b/a)²) (e > 1)
    ///   Parabola → 1
    fn eccentricity(&self) -> f64;

    /// OCCT-aligned: XAxis() — the local X direction (major axis for Circle/Ellipse).
    fn x_axis(&self) -> DVec3;

    /// OCCT-aligned: YAxis() — the local Y direction (minor axis for Circle/Ellipse).
    fn y_axis(&self) -> DVec3;
}

/// OCCT-aligned: `Geom_BoundedCurve` intermediate abstract class.
///
/// Groups bounded curves (BSpline, Bezier) that always have a finite domain.
/// All `BoundedCurve` types have `default_domain()` returning finite values.
pub trait BoundedCurveEval: CurveEval {
    /// The degree of the underlying polynomial/rational representation.
    fn degree(&self) -> usize;
}

// --- ConicEval implementations ---

pub trait Conic2dEval: Curve2dEval {
    /// Center or vertex position.
    fn position(&self) -> DVec2;
    /// OCCT-aligned: eccentricity.
    fn eccentricity(&self) -> f64;
    /// Local X-axis direction (major axis for Circle/Ellipse).
    fn x_axis(&self) -> DVec2;
    /// Local Y-axis direction (minor axis for Circle/Ellipse).
    fn y_axis(&self) -> DVec2;
}

/// OCCT-aligned: `Geom2d_BoundedCurve` intermediate abstract class.
pub trait BoundedCurve2dEval: Curve2dEval {
    fn degree(&self) -> usize;
}

pub trait ElementarySurfaceEval: SurfaceEval {
    /// Origin of the local coordinate system (gp_Ax3::Location).
    fn position(&self) -> DVec3;
    /// Normal / axis direction (gp_Ax3::Direction).
    fn axis_dir(&self) -> DVec3;
    /// U-direction (gp_Ax3::XDirection).
    fn x_axis(&self) -> DVec3;
    /// V-direction (gp_Ax3::YDirection).
    fn y_axis(&self) -> DVec3;
}

/// OCCT-aligned: `Geom_BoundedSurface` intermediate abstract class.
///
/// Groups bounded surfaces (BSpline, Bezier) whose domain is always finite.
pub trait BoundedSurfaceEval: SurfaceEval {
    fn degree_u(&self) -> usize;
    fn degree_v(&self) -> usize;
}

/// OCCT-aligned: `Geom_SweptSurface` intermediate abstract class.
///
/// Groups swept surfaces (LinearExtrusion, Revolution) that share a profile curve.
pub trait SweptSurfaceEval: SurfaceEval {
    fn profile(&self) -> &Curve3;
}
