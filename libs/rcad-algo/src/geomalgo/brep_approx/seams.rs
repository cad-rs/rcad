// OCCT BRepApprox_Approx_0.cxx L67-101 — the alias table of the
// ApproxInt_Approx.gxx instantiation that is BRepApprox_Approx:
//   TheWLine      = BRepApprox_ApproxLine
//   TheMultiLine  = BRepApprox_TheMultiLineOfApprox
// Binding the landed generic engines (approx_int::ApproxIntWLine /
// ApproxIntMultiLine seams) to the BRepApprox data layer; the impls live in
// this child module so they can reach the MultiLine's private fields, exactly
// as the OCCT _0.cxx aliases reach the class internals.
//
// The concrete BRepApprox_Approx driver (the Perform overloads with the
// ThePSurfaceTool quadric switch) is OCCT ApproxInt_Approx.gxx L226-397 and
// lands separately.

use glam::{DVec2, DVec3};

use crate::geomalgo::approx_int::{
    ApproxIntMultiLine, ApproxIntWLine, ApproxLineTrsf, ApproxStatus,
};

use super::{ApproxLine, TheMultiLineOfApprox};

// ---------------------------------------------------------------------------
// The generic ApproxInt_Approx engine seams (approx_int::ApproxIntWLine /
// ApproxIntMultiLine) bound to the BRepApprox data layer — the alias table of
// BRepApprox_Approx_0.cxx L67-101 (TheWLine = BRepApprox_ApproxLine,
// TheMultiLine = BRepApprox_TheMultiLineOfApprox).  With these impls the
// landed engine (approx_int::WLineApprox) runs the BRepApprox instantiation
// without any engine change.
// ---------------------------------------------------------------------------

/// OCCT TheWLine = BRepApprox_ApproxLine binding of the BRepApprox_Approx
/// chain (BRepApprox_Approx_0.cxx L68-71).
impl ApproxIntWLine for ApproxLine {
    /// OCCT BRepApprox_ApproxLine::NbPnts (hxx L43; cxx L47-62).
    fn nb_pnts(&self) -> usize {
        ApproxLine::nb_pnts(self)
    }
    /// OCCT BRepApprox_ApproxLine::Point(i).Value() (cxx L66-94) — the 3d
    /// point of the 1-based Index.
    fn point_p3d(&self, i: usize) -> DVec3 {
        ApproxLine::point(self, i).value()
    }
    /// OCCT BRepApprox_ApproxLine::Point(i).ParametersOnS1 (cxx L66-94).
    fn point_uv1(&self, i: usize) -> DVec2 {
        let (u1, v1, _u2, _v2) = ApproxLine::point(self, i).parameters();
        DVec2::new(u1, v1)
    }
    /// OCCT BRepApprox_ApproxLine::Point(i).ParametersOnS2 (cxx L66-94).
    fn point_uv2(&self, i: usize) -> DVec2 {
        let (_u1, _v1, u2, v2) = ApproxLine::point(self, i).parameters();
        DVec2::new(u2, v2)
    }
}

/// OCCT TheMultiLine = BRepApprox_TheMultiLineOfApprox binding
/// (BRepApprox_Approx_0.cxx L100-101) — forwards to the inherent
/// ApproxInt_MultiLine methods above.
impl<'a> ApproxIntMultiLine for TheMultiLineOfApprox<'a> {
    type TheWLine = ApproxLine;

    /// OCCT BRepApprox_TheMultiLineOfApprox::Line (hxx L44).
    fn the_wline(&self) -> &ApproxLine {
        self.my_line
            .as_deref()
            .expect("BRepApprox_TheMultiLineOfApprox: null line")
    }
    fn first_point(&self) -> usize {
        TheMultiLineOfApprox::first_point(self)
    }
    fn last_point(&self) -> usize {
        TheMultiLineOfApprox::last_point(self)
    }
    fn nb_p3d(&self) -> usize {
        TheMultiLineOfApprox::nb_p3d(self)
    }
    fn nb_p2d(&self) -> usize {
        TheMultiLineOfApprox::nb_p2d(self)
    }
    fn what_status(&self) -> ApproxStatus {
        TheMultiLineOfApprox::what_status(self)
    }
    /// OCCT LineTool::Value(ML, Index, tabP) (lxx L61-69) — slot 1.
    fn value_p3d(&self, index: usize) -> DVec3 {
        let mut tab_pnt = [DVec3::ZERO; 1];
        TheMultiLineOfApprox::value_3d(self, index, &mut tab_pnt);
        tab_pnt[0]
    }
    /// OCCT LineTool::Value(ML, Index, tabP2d) (lxx L70-79) — the NbP2d
    /// slots (the S1 slot is filled unconditionally like the engine seam of
    /// the GeomInt binding, so NbP2d == 0 stays safe for the caller).
    fn value_p2d(&self, index: usize) -> Vec<DVec2> {
        let mut tab_pnt2d = vec![DVec2::ZERO; self.nb_p2d().max(1)];
        TheMultiLineOfApprox::value_2d(self, index, &mut tab_pnt2d);
        tab_pnt2d
    }
    /// OCCT LineTool::Tangency(ML, Index, tabV, tabV2d) (lxx L111-118) —
    /// ML.Tangency(Index, tabV) && ML.Tangency(Index, tabV2d); a false
    /// result degrades the constraint to PassPoint (None here).
    fn tangency(&self, index: usize) -> Option<(Vec<DVec3>, Vec<DVec2>)> {
        let mut tab_v = vec![DVec3::ZERO; self.nb_p3d().max(1)];
        let mut tab_v2d = vec![DVec2::ZERO; self.nb_p2d().max(1)];
        let ok3 = TheMultiLineOfApprox::tangency_3d(self, index, &mut tab_v);
        let ok2 = TheMultiLineOfApprox::tangency_2d(self, index, &mut tab_v2d);
        if ok3 && ok2 {
            Some((tab_v, tab_v2d))
        } else {
            None
        }
    }
    /// OCCT LineTool::MakeMLBetween (lxx L120-126) — returns a (possibly
    /// empty) sub-line; the empty line reports FirstPoint == LastPoint, which
    /// the engine reads as nbpds == 0 (the OCCT failure path).
    fn make_ml_between(&self, low: usize, high: usize, n: usize) -> Option<Self> {
        Some(TheMultiLineOfApprox::make_ml_between(self, low, high, n))
    }
    /// OCCT LineTool::MakeMLOneMorePoint (lxx L128-131) — the OCCT
    /// (bool, TheMultiLine&) out-parameter pair as Option.
    fn make_ml_one_more_point(&self, low: usize, high: usize, indbad: usize) -> Option<Self> {
        let mut the_new_multi_line = TheMultiLineOfApprox::new();
        if TheMultiLineOfApprox::make_ml_one_more_point(
            self,
            low,
            high,
            indbad,
            &mut the_new_multi_line,
        ) {
            Some(the_new_multi_line)
        } else {
            None
        }
    }
    /// OCCT ApproxInt_TheMultiLine constructor as invoked by buildKnots
    /// (ApproxInt_Approx.gxx L548-564) and buildCurve (L638-654) — the same
    /// multi-line restricted to [indicemin, indicemax] carrying the engine
    /// translation; the OCCT P2DOnFirst argument is myData.ApproxU1V1.
    fn sub_line(&self, indicemin: usize, indicemax: usize, the_trsf: &ApproxLineTrsf) -> Self {
        TheMultiLineOfApprox::new_sv_surfaces(
            self.my_line
                .clone()
                .expect("BRepApprox_TheMultiLineOfApprox: null line"),
            self.ptr_on_my_sv_surfaces.clone(),
            self.nbp3d,
            self.nbp2d,
            the_trsf.approx_u1v1,
            the_trsf.approx_u2v2,
            the_trsf.xo,
            the_trsf.yo,
            the_trsf.zo,
            the_trsf.u1o,
            the_trsf.v1o,
            the_trsf.u2o,
            the_trsf.v2o,
            the_trsf.approx_u1v1,
            indicemin,
            indicemax,
        )
    }
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;
    use std::sync::Arc;

    use glam::{DVec2, DVec3};

    use crate::geomalgo::approx_int::{
        ApproxIntMultiLine, ApproxIntWLine, ApproxLineTrsf, ApproxParamType, ApproxStatus,
        WLineApprox,
    };
    use crate::geomalgo::int_surf::{LineOn2S, PntOn2S};

    use crate::geomalgo::brep_approx::{ApproxLine, SharedSvSurfaces, SvSurfaces, TheMultiLineOfApprox};

    /// Walking-line fixture: point i (1-based, i = 1..=n) sits at
    /// p(t) = (t, t*t, 0) with t = (i-1)/(n-1); UV1 = (t, 0.5),
    /// UV2 = (0.25, t) — the analytic images of the parabola on two planes.
    fn walking_line(n: usize) -> Arc<LineOn2S> {
        let mut line = LineOn2S::new();
        for i in 0..n {
            let t = i as f64 / (n - 1) as f64;
            let mut p = PntOn2S::new();
            p.set_value_all(DVec3::new(t, t * t, 0.0), t, 0.5, 0.25, t);
            line.add(&p);
        }
        Arc::new(line)
    }

    /// Analytic SvSurfaces of the fixture: the 3d tangent of the parabola
    /// (1, 2t, 0) normalized, the constant pcurve tangents (1,0) and (0,1).
    struct ParabolaSv {
        use_solver: Cell<bool>,
    }

    impl ParabolaSv {
        fn new() -> Arc<ParabolaSv> {
            Arc::new(ParabolaSv {
                use_solver: Cell::new(false),
            })
        }
    }

    impl SvSurfaces for ParabolaSv {
        fn compute(
            &self,
            u1: &mut f64,
            _v1: &mut f64,
            _u2: &mut f64,
            _v2: &mut f64,
            pt: &mut DVec3,
            tg: &mut DVec3,
            tguv1: &mut DVec2,
            tguv2: &mut DVec2,
        ) -> bool {
            *pt = DVec3::new(*u1, *u1 * *u1, 0.0);
            *tg = DVec3::new(1.0, 2.0 * *u1, 0.0).normalize();
            *tguv1 = DVec2::new(1.0, 0.0);
            *tguv2 = DVec2::new(0.0, 1.0);
            true
        }

        fn pnt(&self, u1: f64, _v1: f64, _u2: f64, _v2: f64, p: &mut DVec3) {
            *p = DVec3::new(u1, u1 * u1, 0.0);
        }

        fn seek_point(&self, _u1: f64, _v1: f64, _u2: f64, _v2: f64, _point: &mut PntOn2S) -> bool {
            false
        }

        fn tangency(&self, u1: f64, _v1: f64, _u2: f64, _v2: f64, tg: &mut DVec3) -> bool {
            *tg = DVec3::new(1.0, 2.0 * u1, 0.0).normalize();
            true
        }

        fn tangency_on_surf_1(&self, _u1: f64, _v1: f64, _u2: f64, _v2: f64, tg: &mut DVec2) -> bool {
            *tg = DVec2::new(1.0, 0.0);
            true
        }

        fn tangency_on_surf_2(&self, _u1: f64, _v1: f64, _u2: f64, _v2: f64, tg: &mut DVec2) -> bool {
            *tg = DVec2::new(0.0, 1.0);
            true
        }

        fn set_use_solver(&self, the_use_sol: bool) {
            self.use_solver.set(the_use_sol);
        }

        fn get_use_solver(&self) -> bool {
            self.use_solver.get()
        }
    }

    /// OCCT anchor: the ApproxIntWLine seam reads BRepApprox_ApproxLine
    /// exactly like ApproxInt_Approx::ComputeTrsf3d/ComputeTrsf2d do
    /// (ApproxInt_Approx.gxx L35-84 over BRepApprox_Approx_0.cxx L68-71).
    #[test]
    fn engine_seam_wline_binding() {
        let wl = ApproxLine::new_line_on_2s(Some(walking_line(5)), false);
        let w: &dyn ApproxIntWLine = &wl;
        assert_eq!(w.nb_pnts(), 5);
        let p = w.point_p3d(3);
        let t = 0.5;
        assert!((p.x - t).abs() < 1e-15 && (p.y - t * t).abs() < 1e-15 && p.z.abs() < 1e-15);
        let uv1 = w.point_uv1(3);
        assert!((uv1.x - t).abs() < 1e-15 && (uv1.y - 0.5).abs() < 1e-15);
        let uv2 = w.point_uv2(3);
        assert!((uv2.x - 0.25).abs() < 1e-15 && (uv2.y - t).abs() < 1e-15);
    }

    /// OCCT anchor: the ApproxIntMultiLine seam forwards Value/Tangency/
    /// WhatStatus through the BRepApprox_TheMultiLineOfApprox bodies
    /// (ApproxInt_MultiLine.gxx L154-298), and the null-SvSurfaces failure
    /// paths return the OCCT literals (no tangency, no point insertion, the
    /// empty MakeMLBetween line over [1, 1]).
    #[test]
    fn engine_seam_multi_line_forwarding() {
        let sv: Option<SharedSvSurfaces> = Some(ParabolaSv::new());
        let line = Arc::new(ApproxLine::new_line_on_2s(Some(walking_line(5)), false));
        let ml = TheMultiLineOfApprox::new_sv_surfaces(
            line, sv, 1, 2, true, true, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, true, 1, 5,
        );
        let m = &ml;
        assert_eq!(m.first_point(), 1);
        assert_eq!(m.last_point(), 5);
        assert_eq!(m.nb_p3d(), 1);
        assert_eq!(m.nb_p2d(), 2);
        assert_eq!(m.what_status(), ApproxStatus::PointsAdded);

        // Value: the translated 3d point and both 2d slots (gxx L178-220).
        let p3 = m.value_p3d(4);
        let t = 3.0 / 4.0;
        assert!((p3.x - t).abs() < 1e-15 && (p3.y - t * t).abs() < 1e-15);
        let p2 = m.value_p2d(4);
        assert_eq!(p2.len(), 2);
        assert!((p2[0].x - t).abs() < 1e-15 && (p2[0].y - 0.5).abs() < 1e-15);
        assert!((p2[1].x - 0.25).abs() < 1e-15 && (p2[1].y - t).abs() < 1e-15);

        // Tangency: both slots from the SvSurfaces (gxx L224-298).
        let (v3, v2) = m.tangency(4).expect("tangency with SvSurfaces");
        let v3_expect = DVec3::new(1.0, 2.0 * t, 0.0).normalize();
        assert!(v3[0].distance(v3_expect) < 1e-15);
        assert_eq!(v2.len(), 2);
        assert_eq!(v2[0], DVec2::new(1.0, 0.0));
        assert_eq!(v2[1], DVec2::new(0.0, 1.0));

        // The null-SvSurfaces line: WhatStatus = NoPointsAdded (gxx L154-160),
        // Tangency false (gxx L229), MakeMLBetween returns the empty line over
        // [1, 1] (gxx L306-328 — FirstPoint == LastPoint is the engine's
        // failure read), MakeMLOneMorePoint false (gxx L637-638).
        let line = Arc::new(ApproxLine::new_line_on_2s(Some(walking_line(5)), false));
        let ml = TheMultiLineOfApprox::new_sv_surfaces(
            line, None, 1, 2, true, true, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, true, 1, 5,
        );
        let m = &ml;
        assert_eq!(m.what_status(), ApproxStatus::NoPointsAdded);
        assert!(m.tangency(2).is_none());
        // The inherent methods take precedence for the concrete type, so the
        // seam methods are called with explicit qualification.
        assert!(
            <TheMultiLineOfApprox as ApproxIntMultiLine>::make_ml_one_more_point(&ml, 1, 5, 2)
                .is_none()
        );
        let other = <TheMultiLineOfApprox as ApproxIntMultiLine>::make_ml_between(&ml, 1, 5, 3)
            .expect("empty sub-line, not None");
        assert_eq!(other.first_point(), 1);
        assert_eq!(other.last_point(), 1);
    }

    /// OCCT anchor: sub_line reproduces the ApproxInt_TheMultiLine
    /// constructor of buildKnots (ApproxInt_Approx.gxx L548-564) — the
    /// restricted [indicemin, indicemax], the engine translation, and the
    /// literal P2DOnFirst = myData.ApproxU1V1.
    #[test]
    fn engine_seam_sub_line_occt_literal() {
        let line = Arc::new(ApproxLine::new_line_on_2s(Some(walking_line(9)), false));
        // NbP3d/NbP2d follow the buildKnots formula (gxx L553-555):
        // (ApproxXYZ?1:0) = 1 and (ApproxU1V1?1:0)+(ApproxU2V2?1:0) = 1 here.
        let ml = TheMultiLineOfApprox::new_sv_surfaces(
            line,
            Some(ParabolaSv::new()),
            1,
            1,
            true,
            true,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            true,
            1,
            9,
        );
        // OCCT passes myData.ApproxU1V1 as P2DOnFirst; with ApproxU1V1 false
        // the single-2d Value must read the S2 slot even though the outer
        // multi-line was built with P2DOnFirst = true.
        let sub = ml.sub_line(
            3,
            7,
            &ApproxLineTrsf {
                approx_u1v1: false,
                approx_u2v2: true,
                xo: 1.0,
                yo: 2.0,
                zo: 3.0,
                u1o: 0.0,
                v1o: 0.0,
                u2o: 10.0,
                v2o: 20.0,
            },
        );
        assert_eq!(sub.first_point(), 3);
        assert_eq!(sub.last_point(), 7);
        // t = (4-1)/8 = 0.375 at line index 4.
        let p3 = sub.value_p3d(4);
        assert!((p3.x - (0.375 + 1.0)).abs() < 1e-15);
        assert!((p3.y - (0.375 * 0.375 + 2.0)).abs() < 1e-15);
        let p2 = sub.value_p2d(4);
        assert_eq!(p2.len(), 1);
        // P2DOnFirst == ApproxU1V1 == false (gxx L564) -> the single 2d slot
        // is the S2 image (u2 + U2o, v2 + V2o).
        assert!((p2[0].x - (0.25 + 10.0)).abs() < 1e-15);
        assert!((p2[0].y - (0.375 + 20.0)).abs() < 1e-15);
    }

    /// OCCT anchor: the landed ApproxInt_Approx engine (approx_int::
    /// WLineApprox) driven through the BRepApprox binding approximates the
    /// walking line to tolerance — the 3d parabola and both pcurve images.
    #[test]
    fn engine_seam_approx_engine_round_trip() {
        let line = Arc::new(ApproxLine::new_line_on_2s(Some(walking_line(9)), false));
        let ml = TheMultiLineOfApprox::new_sv_surfaces(
            line,
            Some(ParabolaSv::new()),
            1,
            2,
            true,
            true,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            true,
            1,
            9,
        );
        let mut apx = WLineApprox::new();
        // The consumer pattern (bop/int_tools/face_make_curve.rs L1407):
        // SetParameters before Perform; the tolerances are divided by the
        // ratio 1.5 internally (ApproxInt_Approx.gxx L399-427).
        apx.set_parameters(
            1.0e-6, 1.0e-6, 4, 8, 40, 30, true, ApproxParamType::ChordLength,
        );
        apx.perform(&ml, true, true, true, 1, 9);
        assert!(apx.my_done, "approximation not done");
        // The single-pass ComputeLine fit of the 9-point parabola with
        // tangency constraints reports its achieved tolerance.
        assert!(
            apx.my_tol_reached3d < 1.0e-3,
            "tol3d {}",
            apx.my_tol_reached3d
        );

        // Endpoint constraints (PassPoint/Tangency): first and last poles are
        // the exact walking-line endpoints in 3d and on both pcurves.
        let poles = &apx.value().poles;
        assert!(poles.len() >= 2);
        let first = &poles[0];
        let last = &poles[poles.len() - 1];
        assert!(first.p3d[0].distance(DVec3::ZERO) < 1e-10);
        assert!(last.p3d[0].distance(DVec3::new(1.0, 1.0, 0.0)) < 1e-10);
        assert!(first.p2d[0].distance(DVec2::new(0.0, 0.5)) < 1e-10);
        assert!(last.p2d[0].distance(DVec2::new(1.0, 0.5)) < 1e-10);
        assert!(first.p2d[1].distance(DVec2::new(0.25, 0.0)) < 1e-10);
        assert!(last.p2d[1].distance(DVec2::new(0.25, 1.0)) < 1e-10);
    }
}
