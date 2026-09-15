//! Per-kind evaluation impls for the offset kinds: `Geom_OffsetCurve`,
//! `Geom_OffsetSurface` and `Geom2d_OffsetCurve`.  The exact derivative
//! engines live in `geom::offset_surface_utils` / `geom::offset2d_dn`.
//!
//! Extracted verbatim from `geom/eval.rs` (project Rule 5 file-size split);
//! the evaluation traits stay re-exported from `geom/eval.rs`.
use crate::geom::*;

impl CurveEval for OffsetCurve3 {
    fn point_at(&self, t: f64) -> DVec3 {
        let base_pt = self.basis.point_at(t);
        let tangent = self.basis.tangent_at(t);
        let perp = tangent.cross(self.offset_dir);
        let perp_len = perp.length();
        if perp_len < 1e-15 {
            return base_pt;
        }
        base_pt + self.offset_distance * (perp / perp_len)
    }
    fn tangent_at(&self, t: f64) -> DVec3 {
        let eps = 1e-6;
        let [t0, t1] = self.basis.default_domain();
        let t_lo = (t - eps).max(t0);
        let t_hi = (t + eps).min(t1);
        let dp = self.point_at(t_hi) - self.point_at(t_lo);
        let len = dp.length();
        if len < 1e-15 { DVec3::X } else { dp / len }
    }
    fn default_domain(&self) -> [f64; 2] {
        self.basis.default_domain()
    }
}

impl SurfaceEval for OffsetSurface {
    /// OCCT `Geom_OffsetSurface::EvalD0` (Geom_OffsetSurface.cxx L338-358).
    fn point_at(&self, u: f64, v: f64) -> DVec3 {
        crate::geom::offset_surface_utils::offset_payload_eval_d0(self, u, v)
    }
    fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        // Offset preserves the normal direction (first-order approximation)
        self.basis.normal_at(u, v)
    }
    fn default_domain(&self) -> [f64; 4] {
        self.basis.default_domain()
    }
    /// OCCT `Geom_OffsetSurface::EvalD1` (Geom_OffsetSurface.cxx L362-389) —
    /// the `Geom_OffsetSurfaceUtils::EvaluateD1` body instead of the trait's
    /// finite-difference default.
    fn derivatives(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let d1 = crate::geom::offset_surface_utils::offset_payload_eval_d1(self, u, v);
        (d1.point, d1.d1u, d1.d1v)
    }
    /// OCCT `Geom_OffsetSurface::EvalD2` (Geom_OffsetSurface.cxx L393-424) —
    /// the `Geom_OffsetSurfaceUtils::EvaluateD2` body instead of the trait's
    /// finite-difference default.
    fn derivatives2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        let d2 = crate::geom::offset_surface_utils::offset_payload_eval_d2(self, u, v);
        (d2.point, d2.d1u, d2.d1v, d2.d2u, d2.d2uv, d2.d2v)
    }
}

impl Curve2dEval for OffsetCurve2d {
    fn point_at(&self, t: f64) -> DVec2 {
        // OCCT Geom2d_OffsetCurve::EvalD0 (Geom2d_OffsetCurve.cxx L214-229)
        // — the basis point plus Geom2d_OffsetCurveUtils::CalculateD0 over
        // the basis D1; the analytic engine lives in offset2d_dn.
        crate::geom::offset2d_dn::eval_d0(self, t)
    }
    fn derivative_at(&self, t: f64) -> DVec2 {
        // OCCT Geom2d_OffsetCurve::EvalD1 (Geom2d_OffsetCurve.cxx L233-249)
        // — the basis D2 through CalculateD1 (offset2d_dn).
        crate::geom::offset2d_dn::eval_d1(self, t).d1
    }
    fn derivative2_at(&self, t: f64) -> DVec2 {
        // OCCT Geom2d_OffsetCurve::EvalD2 (Geom2d_OffsetCurve.cxx L253-285)
        // — the basis D3, the AdjustDerivative(3) singular guard and
        // CalculateD2 (offset2d_dn).
        crate::geom::offset2d_dn::eval_d2(self, t).d2
    }
    fn derivative3_at(&self, t: f64) -> DVec2 {
        // OCCT Geom2d_OffsetCurve::EvalD3 (Geom2d_OffsetCurve.cxx L289-328)
        // — the basis D3/DN(4), the AdjustDerivative(4) singular guard and
        // CalculateD3 (offset2d_dn).
        crate::geom::offset2d_dn::eval_d3(self, t).d3
    }
    fn derivative_n_at(&self, t: f64, n: i32) -> DVec2 {
        // OCCT Geom2d_OffsetCurve::EvalDN (Geom2d_OffsetCurve.cxx L332-357).
        crate::geom::offset2d_dn::eval_dn(self, t, n)
    }
    fn tangent_at(&self, t: f64) -> DVec2 {
        // The rcad tangent_at is the normalized D1 (the BSplineCurve2 impl
        // convention); no separate OCCT counterpart.
        self.derivative_at(t).normalize_or_zero()
    }
    fn default_domain(&self) -> [f64; 2] {
        // OCCT FirstParameter/LastParameter (Geom2d_OffsetCurve.cxx
        // L361-371) delegate to the basis curve.
        self.basis.default_domain()
    }
}
