//! Per-kind evaluation impls for the surface kinds: elementary quadrics
//! (plane, cylinder, sphere, cone, torus), the OCCT `Geom_SweptSurface` pair
//! (linear extrusion, revolution), the trimmed wrapper and the composite
//! kinds (pipe, Coons, ruled).
//!
//! Extracted verbatim from `geom/eval.rs` (project Rule 5 file-size split);
//! the evaluation traits stay re-exported from `geom/eval.rs`.
use crate::geom::*;
use crate::core::precision::{INFINITE_VALUE, is_infinite_value};
use crate::geom::orthonormal_frame;
use std::f64::consts::PI;

// --- SurfaceEval implementations ---

impl SurfaceEval for Plane {
    /// OCCT-aligned: P(u,v) = origin + u*u_dir + v*v_dir using stored axes.
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        self.origin + u * self.u_dir + v * self.v_dir
    }
    fn normal_at(&self, _u: f64, _v: f64) -> DVec3 {
        self.normal
    }
    fn default_domain(&self) -> [f64; 4] {
        // OCCT Geom_Plane::Bounds (Geom_Plane.cxx L181-184): +/-Precision::Infinite().
        [
            -INFINITE_VALUE,
            INFINITE_VALUE,
            -INFINITE_VALUE,
            INFINITE_VALUE,
        ]
    }
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        (
            self.origin + u * self.u_dir + v * self.v_dir,
            self.u_dir,
            self.v_dir,
        )
    }
    fn derivatives2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        (
            self.origin + u * self.u_dir + v * self.v_dir,
            self.u_dir,
            self.v_dir,
            DVec3::ZERO, // Puu
            DVec3::ZERO, // Puv
            DVec3::ZERO, // Pvv
        )
    }
}

impl SurfaceEval for CylindricalSurface {
    /// u = azimuth angle [0, 2π], v = height along axis.
    /// OCCT-aligned: uses stored ref_dir for deterministic UV mapping.
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = self.y_axis();
        self.origin + self.radius * (u.cos() * x_ax + u.sin() * y_ax) + v * self.axis
    }
    fn normal_at(&self, u: f64, _v: f64) -> DVec3 {
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = self.y_axis();
        (u.cos() * x_ax + u.sin() * y_ax).normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        // OCCT Geom_CylindricalSurface::Bounds
        // (Geom_CylindricalSurface.cxx L162-163): V = +/-Precision::Infinite().
        [0.0, 2.0 * PI, -INFINITE_VALUE, INFINITE_VALUE]
    }
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = self.y_axis();
        let (su, cu) = u.sin_cos();
        let p = self.origin + self.radius * (cu * x_ax + su * y_ax) + v * self.axis;
        let dpu = self.radius * (-su * x_ax + cu * y_ax);
        (p, dpu, self.axis)
    }
    fn is_u_closed(&self) -> bool {
        true
    }
    fn is_u_periodic(&self) -> bool {
        true
    }
    fn u_reversed_parameter(&self, t: f64) -> f64 {
        2.0 * PI - t
    }
}

impl CylindricalSurface {
    /// UV coordinates of world point `p` relative to this cylindrical surface.
    ///
    /// `u` = azimuth (−π, π], `v` = height along axis, matching
    /// [`SurfaceEval::point_at`] (which uses the stored `ref_dir` for u=0).
    /// Off-surface points are radially projected onto the cylinder.
    pub fn world_to_uv(self, p: DVec3) -> DVec2 {
        let axis = self.axis.normalize_or_zero();
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = self.y_axis();
        let d = p - self.origin;
        let v = d.dot(axis);
        let radial = d - axis * v;
        let r = radial.length();
        if r < 1e-15 {
            return DVec2::new(0.0, v);
        }
        let u = radial.dot(y_ax).atan2(radial.dot(x_ax));
        DVec2::new(u, v)
    }
}

impl SphericalSurface {
    /// Construct a sphere with `ref_dir` derived from [`any_perpendicular(axis)`](any_perpendicular).
    pub fn new(center: Point3, axis: Vec3, radius: f64) -> Self {
        Self {
            center,
            axis,
            radius,
            ref_dir: any_perpendicular(axis),
        }
    }

    /// Construct a sphere with an explicit `ref_dir` (used after mirroring / transforming).
    pub fn new_with_ref_dir(center: Point3, axis: Vec3, radius: f64, ref_dir: Vec3) -> Self {
        Self {
            center,
            axis,
            radius,
            ref_dir,
        }
    }

    /// Spherical coordinates of world point `p`: longitude `u` ∈ (−π, π],
    /// latitude `v` ∈ [−π/2, π/2] (0 = equator, +π/2 = axis pole), matching
    /// [`SurfaceEval::point_at`] / `properties` sphere helpers and OCCT
    /// `ElSLib::SphereD0` (radial projection when `p` is off the surface).
    pub fn world_to_uv(self, p: DVec3) -> DVec2 {
        let ax = self.axis.normalize_or_zero();
        let r = self.radius;
        if r < 1e-15 {
            return DVec2::ZERO;
        }
        let w = (p - self.center) / r;
        if w.length_squared() < 1e-20 {
            return DVec2::ZERO;
        }
        let w = w.normalize();
        // OCCT ElSLib::SphereD0: the axis component is R*sin(V), so
        // V = asin(w.A) (0 = equator, +pi/2 = axis pole).
        let v = w.dot(ax).clamp(-1.0, 1.0).asin();
        let x_ax = self.ref_dir.normalize();
        let y_ax = ax.cross(x_ax).normalize();
        let w_t = w - ax * w.dot(ax);
        if w_t.length_squared() < 1e-12 {
            return DVec2::new(0.0, v);
        }
        let w_t = w_t.normalize();
        let u = w_t.dot(y_ax).atan2(w_t.dot(x_ax));
        DVec2::new(u, v)
    }
}

impl SurfaceEval for SphericalSurface {
    /// OCCT ElSLib::SphereD0: P = O + R*(cos(v)*radial + sin(v)*A),
    /// radial = cos(u)X + sin(u)Y.  u = longitude [0, 2π], v = latitude
    /// [−π/2, π/2] (0 = equator, +π/2 = axis pole).
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = self.ref_dir.normalize();
        let y_ax = self.axis.cross(x_ax).normalize();
        self.center
            + self.radius * (v.cos() * (u.cos() * x_ax + u.sin() * y_ax) + v.sin() * self.axis)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let p = self.point_at(u, v);
        (p - self.center).normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, -PI / 2.0, PI / 2.0]
    }
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let x_ax = self.ref_dir.normalize();
        let y_ax = self.axis.cross(x_ax).normalize();
        let (su, cu) = u.sin_cos();
        let (sv, cv) = v.sin_cos();
        let radial = cu * x_ax + su * y_ax;
        let p = self.center + self.radius * (cv * radial + sv * self.axis);
        let dpu = self.radius * cv * (-su * x_ax + cu * y_ax);
        let dpv = self.radius * (-sv * radial + cv * self.axis);
        (p, dpu, dpv)
    }
    fn derivatives2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        // OCCT Geom_SphericalSurface::D2 (ElSLib::SphereD2). In the basis
        // (X, Y, A) with radial = cos(u)X + sin(u)Y and
        // P = O + R*(cos(v)*radial + sin(v)*A):
        //   Puu = R cos(v)(-cos(u)X - sin(u)Y)
        //   Puv = -R sin(v)(-sin(u)X + cos(u)Y)
        //   Pvv = -(P - center)
        let x_ax = self.ref_dir.normalize();
        let y_ax = self.axis.cross(x_ax).normalize();
        let (su, cu) = u.sin_cos();
        let (sv, cv) = v.sin_cos();
        let radial = cu * x_ax + su * y_ax;
        let p = self.center + self.radius * (cv * radial + sv * self.axis);
        let dpu = self.radius * cv * (-su * x_ax + cu * y_ax);
        let dpv = self.radius * (-sv * radial + cv * self.axis);
        let d2u = self.radius * cv * (-cu * x_ax - su * y_ax);
        let duv = -self.radius * sv * (-su * x_ax + cu * y_ax);
        let d2v = -(p - self.center);
        (p, dpu, dpv, d2u, duv, d2v)
    }
    fn is_u_closed(&self) -> bool {
        true
    }
    fn is_u_periodic(&self) -> bool {
        true
    }
    fn u_reversed_parameter(&self, t: f64) -> f64 {
        2.0 * PI - t
    }
    fn v_reversed_parameter(&self, t: f64) -> f64 {
        PI - t
    } // OCCT: colatitude [0, π]
}

impl SurfaceEval for ConicalSurface {
    /// u = azimuth [0, 2π], v = distance along the cone generatrix from the
    /// reference circle at `self.apex`.
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let axis = self.axis_dir();
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = axis.cross(x_ax).normalize();
        let radial = self.radius_at_slant(v);
        let axial = self.axial_from_slant(v);
        self.apex + axial * axis + radial * (u.cos() * x_ax + u.sin() * y_ax)
    }
    fn normal_at(&self, u: f64, _v: f64) -> DVec3 {
        let axis = self.axis_dir();
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = axis.cross(x_ax).normalize();
        let radial = u.cos() * x_ax + u.sin() * y_ax;
        let half = self.half_angle_rad;
        (radial * half.cos() - axis * half.sin()).normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        // OCCT Geom_ConicalSurface::Bounds (Geom_ConicalSurface.cxx L212-213):
        // V1 = -Precision::Infinite() (NOT 0), V2 = Precision::Infinite().
        [0.0, 2.0 * PI, -INFINITE_VALUE, INFINITE_VALUE]
    }
    fn is_u_closed(&self) -> bool {
        true
    }
    fn is_u_periodic(&self) -> bool {
        true
    }
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let axis = self.axis_dir();
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = axis.cross(x_ax).normalize();
        let (su, cu) = u.sin_cos();
        let radial = self.radius_at_slant(v);
        let axial = self.axial_from_slant(v);
        // d(radius)/dv = sin(half_angle), d(axial)/dv = cos(half_angle)
        let half = self.half_angle_rad;
        let dr = half.sin();
        let da = half.cos();
        let r_vec = cu * x_ax + su * y_ax;
        let p = self.apex + axial * axis + radial * r_vec;
        let dpu = radial * (-su * x_ax + cu * y_ax);
        let dpv = da * axis + dr * r_vec;
        (p, dpu, dpv)
    }
}

impl SurfaceEval for ToroidalSurface {
    /// u = major angle [0, 2π], v = minor angle [0, 2π].
    /// OCCT-aligned: uses the stored ref_dir (gp_Ax3 XDirection) for u=0, so
    /// a rotated torus keeps its seam position in the UV mapping.
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = self.axis.cross(x_ax).normalize();
        let tube_center = self.center + self.major_radius * (u.cos() * x_ax + u.sin() * y_ax);
        let radial = (u.cos() * x_ax + u.sin() * y_ax).normalize();
        tube_center + self.minor_radius * (v.cos() * radial + v.sin() * self.axis)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = self.axis.cross(x_ax).normalize();
        let radial = (u.cos() * x_ax + u.sin() * y_ax).normalize();
        (v.cos() * radial + v.sin() * self.axis).normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, 0.0, 2.0 * PI]
    }
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = self.axis.cross(x_ax).normalize();
        let (su, cu) = u.sin_cos();
        let (sv, cv) = v.sin_cos();
        let r_vec = cu * x_ax + su * y_ax;
        let r_perp = -su * x_ax + cu * y_ax;
        let r_major = self.major_radius;
        let r_minor = self.minor_radius;
        let tube = r_major + r_minor * cv;
        let p = self.center + tube * r_vec + r_minor * sv * self.axis;
        let dpu = tube * r_perp;
        let dpv = -r_minor * sv * r_vec + r_minor * cv * self.axis;
        (p, dpu, dpv)
    }
    fn is_u_closed(&self) -> bool {
        true
    }
    fn is_u_periodic(&self) -> bool {
        true
    }
    fn u_reversed_parameter(&self, t: f64) -> f64 {
        2.0 * PI - t
    }
    fn is_v_closed(&self) -> bool {
        true
    }
    fn is_v_periodic(&self) -> bool {
        true
    }
    fn v_reversed_parameter(&self, t: f64) -> f64 {
        2.0 * PI - t
    }
    /// OCCT ElSLib::D2 for a torus — exact analytic second derivatives
    /// (the trait-default finite differences have no OCCT counterpart; the
    /// Contap walker consumes D2 through Contap_SurfProps::NormAndDn's
    /// generic branch for toroidal surfaces).
    fn derivatives2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = self.axis.cross(x_ax).normalize();
        let (su, cu) = u.sin_cos();
        let (sv, cv) = v.sin_cos();
        let r_vec = cu * x_ax + su * y_ax;
        let r_perp = -su * x_ax + cu * y_ax;
        let tube = self.major_radius + self.minor_radius * cv;
        let p = self.center + tube * r_vec + self.minor_radius * sv * self.axis;
        let dpu = tube * r_perp;
        let dpv = -self.minor_radius * sv * r_vec + self.minor_radius * cv * self.axis;
        let dpuu = -tube * r_vec;
        let dpuv = self.minor_radius * sv * r_perp;
        let dpvv = -self.minor_radius * (cv * r_vec + sv * self.axis);
        (p, dpu, dpv, dpuu, dpuv, dpvv)
    }
}

impl ToroidalSurface {
    /// UV coordinates of world point `p` relative to this toroidal surface.
    ///
    /// `u` = major angle (−π, π], `v` = minor angle [0, 2π),
    /// matching [`SurfaceEval::point_at`].  When `p` is on the surface
    /// the returned `(u, v)` is exact; off-surface points project onto
    /// the tube center circle in the radial direction.
    pub fn world_to_uv(self, p: DVec3) -> DVec2 {
        use std::f64::consts::TAU;
        let axis = self.axis.normalize_or_zero();
        // OCCT ElSLib::TorusParameters uses the torus's gp_Ax3 XDirection
        // (the stored ref_dir, preserved through rotation) as the u=0 ref.
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = axis.cross(x_ax).normalize();
        let local = p - self.center;
        let axial = local.dot(axis);
        let radial_vec = local - axis * axial;
        let radial_dist = radial_vec.length();

        // u = azimuth around main axis
        let u = if radial_dist < 1e-15 {
            0.0
        } else {
            let rn = radial_vec / radial_dist;
            rn.dot(y_ax).atan2(rn.dot(x_ax))
        };

        // v = angle around tube:
        // On surface: radial_dist = R + r·cos(v), axial = r·sin(v)
        //   → v = atan2(axial, radial_dist - R)
        let v_base = axial.atan2(radial_dist - self.major_radius);
        // Convert v from [-π, π] to [0, 2π)
        let v = if v_base < 0.0 { v_base + TAU } else { v_base };

        DVec2::new(u, v)
    }
}

impl SurfaceEval for EllipsoidalSurface {
    /// See [`crate::geom::eval_c::ellipsoid_eval_d1`] for the parameterisation
    /// mapping between the rcad colatitude payload and the OCCT
    /// `GeomEval_EllipsoidSurface` latitude.
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let (axis, x_axis, y_axis) = orthonormal_frame(self.axis, self.ref_dir);
        self.center
            + self.radius_x * v.sin() * u.cos() * x_axis
            + self.radius_y * v.sin() * u.sin() * y_axis
            + self.radius_z * v.cos() * axis
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let (axis, x_axis, y_axis) = orthonormal_frame(self.axis, self.ref_dir);
        let p = self.point_at(u, v) - self.center;
        let x = p.dot(x_axis);
        let y = p.dot(y_axis);
        let z = p.dot(axis);
        let grad = (x / (self.radius_x * self.radius_x)) * x_axis
            + (y / (self.radius_y * self.radius_y)) * y_axis
            + (z / (self.radius_z * self.radius_z)) * axis;
        grad.normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 4] {
        [0.0, 2.0 * PI, 0.0, PI]
    }
    /// OCCT `GeomEval_EllipsoidSurface::EvalD1` (GeomEval_EllipsoidSurface.cxx
    /// L197-221) instead of the trait's finite-difference default.
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let d1 = crate::geom::eval_c::ellipsoid_eval_d1(self, u, v);
        (d1.point, d1.d1u, d1.d1v)
    }
    /// OCCT `GeomEval_EllipsoidSurface::EvalD2` (cxx L225-259).
    fn derivatives2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        let d2 = crate::geom::eval_c::ellipsoid_eval_d2(self, u, v);
        (d2.point, d2.d1u, d2.d1v, d2.d2u, d2.d2uv, d2.d2v)
    }
}

impl SurfaceEval for HelicoidSurface {
    /// OCCT `GeomEval_CircularHelicoidSurface::EvalD0` (cxx L143-156).
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let (axis, x_axis, y_axis) = orthonormal_frame(self.axis, self.ref_dir);
        let lead = self.pitch / (2.0 * PI);
        self.origin + v * (u.cos() * x_axis + u.sin() * y_axis) + (lead * u) * axis
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let (axis, x_axis, y_axis) = orthonormal_frame(self.axis, self.ref_dir);
        let lead = self.pitch / (2.0 * PI);
        let du = v * (-u.sin() * x_axis + u.cos() * y_axis) + lead * axis;
        let dv = u.cos() * x_axis + u.sin() * y_axis;
        du.cross(dv).normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 4] {
        [-2.0 * PI, 2.0 * PI, -10.0, 10.0]
    }
    /// OCCT `GeomEval_CircularHelicoidSurface::EvalD1` (cxx L160-183) instead
    /// of the trait's finite-difference default.
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let d1 = crate::geom::eval_c::helicoid_eval_d1(self, u, v);
        (d1.point, d1.d1u, d1.d1v)
    }
    /// OCCT `GeomEval_CircularHelicoidSurface::EvalD2` (cxx L187-217).
    fn derivatives2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        let d2 = crate::geom::eval_c::helicoid_eval_d2(self, u, v);
        (d2.point, d2.d1u, d2.d1v, d2.d2u, d2.d2uv, d2.d2v)
    }
}

impl SurfaceEval for LinearExtrusionSurface {
    /// u = profile parameter, v = extrusion distance along direction.
    /// OCCT `Geom_SurfaceOfLinearExtrusion::EvalD0` (cxx L150-163).
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        crate::geom::extrusion_utils::linear_extrusion_eval_d0(self, u, v)
    }
    fn normal_at(&self, u: f64, _v: f64) -> DVec3 {
        let tangent = self.profile.tangent_at(u);
        let n = tangent.cross(self.direction);
        if n.length_squared() < 1e-20 {
            return DVec3::Z;
        }
        n.normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        // OCCT Geom_SurfaceOfLinearExtrusion::Bounds
        // (Geom_SurfaceOfLinearExtrusion.cxx L142-145): V = +/-Infinite(),
        // U = basis curve parameters.
        let [t1, t2] = self.profile.default_domain();
        [t1, t2, -INFINITE_VALUE, INFINITE_VALUE]
    }
    /// OCCT `Geom_SurfaceOfLinearExtrusion::EvalD1` (cxx L166-184) —
    /// `Geom_ExtrusionUtils::CalculateD1` instead of the trait's
    /// finite-difference default.
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let d1 = crate::geom::extrusion_utils::linear_extrusion_eval_d1(self, u, v);
        (d1.point, d1.d1u, d1.d1v)
    }
    /// OCCT `Geom_SurfaceOfLinearExtrusion::EvalD2` (cxx L188-210) —
    /// `Geom_ExtrusionUtils::CalculateD2` instead of the trait's
    /// finite-difference default.
    fn derivatives2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        let d2 = crate::geom::extrusion_utils::linear_extrusion_eval_d2(self, u, v);
        (d2.point, d2.d1u, d2.d1v, d2.d2u, d2.d2uv, d2.d2v)
    }
}

impl SurfaceEval for RevolutionSurface {
    /// u = azimuth angle [0, 2π], v = profile parameter.
    /// OCCT `Geom_SurfaceOfRevolution::EvalD0` (cxx L236-248) —
    /// `Geom_RevolutionUtils::CalculateD0`.
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        crate::geom::revolution_utils::revolution_eval_d0(self, u, v)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-6;
        let du = (self.point_at(u + eps, v) - self.point_at(u - eps, v)) / (2.0 * eps);
        let dv = (self.point_at(u, v + eps) - self.point_at(u, v - eps)) / (2.0 * eps);
        let n = du.cross(dv);
        if n.length_squared() < 1e-20 {
            return DVec3::Z;
        }
        n.normalize()
    }
    fn default_domain(&self) -> [f64; 4] {
        let [t1, t2] = self.profile.default_domain();
        [0.0, 2.0 * PI, t1, t2]
    }
    /// OCCT `Geom_SurfaceOfRevolution::EvalD1` (cxx L252-270) —
    /// `Geom_RevolutionUtils::CalculateD1` instead of the trait's
    /// finite-difference default.
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let d1 = crate::geom::revolution_utils::revolution_eval_d1(self, u, v);
        (d1.point, d1.d1u, d1.d1v)
    }
    /// OCCT `Geom_SurfaceOfRevolution::EvalD2` (cxx L274-296) —
    /// `Geom_RevolutionUtils::CalculateD2` instead of the trait's
    /// finite-difference default.
    fn derivatives2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        let d2 = crate::geom::revolution_utils::revolution_eval_d2(self, u, v);
        (d2.point, d2.d1u, d2.d1v, d2.d2u, d2.d2uv, d2.d2v)
    }
}

impl SurfaceEval for TrimmedSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        self.basis.point_at(u, v)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        self.basis.normal_at(u, v)
    }
    fn default_domain(&self) -> [f64; 4] {
        self.trim
    }
}

// --- SweptSurfaceEval implementations ---

impl SweptSurfaceEval for LinearExtrusionSurface {
    fn profile(&self) -> &Curve3 { &self.profile }
}

impl SweptSurfaceEval for RevolutionSurface {
    fn profile(&self) -> &Curve3 { &self.profile }
}

// --- ElementarySurfaceEval implementations ---

impl ElementarySurfaceEval for Plane {
    fn position(&self) -> DVec3 { self.origin }
    fn axis_dir(&self) -> DVec3 { self.normal }
    fn x_axis(&self) -> DVec3 { self.u_dir }
    fn y_axis(&self) -> DVec3 { self.v_dir }
}

impl ElementarySurfaceEval for CylindricalSurface {
    fn position(&self) -> DVec3 { self.origin }
    fn axis_dir(&self) -> DVec3 { self.axis }
    fn x_axis(&self) -> DVec3 { self.ref_dir.normalize_or_zero() }
    fn y_axis(&self) -> DVec3 { self.y_axis() }
}

impl ElementarySurfaceEval for SphericalSurface {
    fn position(&self) -> DVec3 { self.center }
    fn axis_dir(&self) -> DVec3 { self.axis }
    fn x_axis(&self) -> DVec3 { self.ref_dir.normalize_or_zero() }
    fn y_axis(&self) -> DVec3 { self.axis.cross(self.ref_dir).normalize_or_zero() }
}

impl ElementarySurfaceEval for ConicalSurface {
    fn position(&self) -> DVec3 { self.apex_point() }
    fn axis_dir(&self) -> DVec3 { self.axis_dir() }
    fn x_axis(&self) -> DVec3 { self.ref_dir.normalize_or_zero() }
    fn y_axis(&self) -> DVec3 { self.axis_dir().cross(self.ref_dir.normalize_or_zero()).normalize_or_zero() }
}

impl ElementarySurfaceEval for ToroidalSurface {
    fn position(&self) -> DVec3 { self.center }
    fn axis_dir(&self) -> DVec3 { self.axis }
    fn x_axis(&self) -> DVec3 { any_perpendicular(self.axis) }
    fn y_axis(&self) -> DVec3 { self.axis.cross(any_perpendicular(self.axis)).normalize_or_zero() }
}

fn remap_unit_to_curve_domain(curve: &Curve3, t: f64) -> f64 {
    let [t0, t1] = curve.default_domain();
    // OCCT Precision::IsInfinite (Precision.hxx L350-353): an unbounded
    // domain leaves the unit parameter untouched.
    if is_infinite_value(t0) || is_infinite_value(t1) {
        return t;
    }
    t0 + (t1 - t0) * t
}

fn projected_frame_from_tangent(tangent: DVec3, ref_dir: DVec3) -> (DVec3, DVec3) {
    let tangent = tangent.normalize_or_zero();
    let mut x_axis = ref_dir - tangent * ref_dir.dot(tangent);
    if x_axis.length_squared() <= 1e-24 {
        x_axis = any_perpendicular(tangent);
    } else {
        x_axis = x_axis.normalize();
    }
    let y_axis = tangent.cross(x_axis).normalize_or_zero();
    (x_axis, y_axis)
}

impl SurfaceEval for PipeSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let center = self.spine.point_at(v);
        let tangent = self.spine.tangent_at(v);
        let (x_axis, y_axis) = projected_frame_from_tangent(tangent, self.ref_dir);
        center + self.radius * (u.cos() * x_axis + u.sin() * y_axis)
    }

    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let du = self.point_at(u + eps, v) - self.point_at(u - eps, v);
        let dv = self.point_at(u, v + eps) - self.point_at(u, v - eps);
        du.cross(dv).normalize_or_zero()
    }

    fn default_domain(&self) -> [f64; 4] {
        let [v0, v1] = self.spine.default_domain();
        [0.0, 2.0 * PI, v0, v1]
    }
}

impl SurfaceEval for CoonsSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let south = self
            .south
            .point_at(remap_unit_to_curve_domain(&self.south, u));
        let north = self
            .north
            .point_at(remap_unit_to_curve_domain(&self.north, u));
        let west = self
            .west
            .point_at(remap_unit_to_curve_domain(&self.west, v));
        let east = self
            .east
            .point_at(remap_unit_to_curve_domain(&self.east, v));

        let p00 = self
            .south
            .point_at(remap_unit_to_curve_domain(&self.south, 0.0));
        let p10 = self
            .south
            .point_at(remap_unit_to_curve_domain(&self.south, 1.0));
        let p01 = self
            .north
            .point_at(remap_unit_to_curve_domain(&self.north, 0.0));
        let p11 = self
            .north
            .point_at(remap_unit_to_curve_domain(&self.north, 1.0));

        let linear_u = south * (1.0 - v) + north * v;
        let linear_v = west * (1.0 - u) + east * u;
        let bilinear = p00 * ((1.0 - u) * (1.0 - v))
            + p10 * (u * (1.0 - v))
            + p01 * ((1.0 - u) * v)
            + p11 * (u * v);
        linear_u + linear_v - bilinear
    }

    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let du = self.point_at((u + eps).clamp(0.0, 1.0), v)
            - self.point_at((u - eps).clamp(0.0, 1.0), v);
        let dv = self.point_at(u, (v + eps).clamp(0.0, 1.0))
            - self.point_at(u, (v - eps).clamp(0.0, 1.0));
        du.cross(dv).normalize_or_zero()
    }

    fn default_domain(&self) -> [f64; 4] {
        [0.0, 1.0, 0.0, 1.0]
    }
}

impl SurfaceEval for RuledSurface {
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let start = self.start.point_at(u);
        let end = self.end.point_at(u);
        start.lerp(end, v)
    }

    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let eps = 1e-5;
        let du = (self.point_at(u + eps, v) - self.point_at(u - eps, v)) / (2.0 * eps);
        let dv = self.end.point_at(u) - self.start.point_at(u);
        du.cross(dv).normalize_or_zero()
    }

    fn default_domain(&self) -> [f64; 4] {
        let [u0, u1] = self.start.default_domain();
        [u0, u1, 0.0, 1.0]
    }
}
