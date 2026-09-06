// OCCT BRepApprox_Approx (BRepApprox_Approx.hxx L51-177) — the
// ApproxInt_Approx.gxx instantiation of BRepApprox_Approx_0.cxx L67-101:
//   ThePSurface     = BRepAdaptor_Surface
//   ThePSurfaceTool = BRepApprox_SurfaceTool
//   TheISurface     = IntSurf_Quadric
//   TheWLine        = BRepApprox_ApproxLine
//   TheMultiLine    = BRepApprox_TheMultiLineOfApprox
//   ThePrmPrmSvSurfaces = BRepApprox_ThePrmPrmSvSurfacesOfApprox
//   TheImpPrmSvSurfaces = BRepApprox_TheImpPrmSvSurfacesOfApprox
//
// The engine bodies (prepareDS / fillData / buildKnots / buildCurve and the
// tolerance bookkeeping) are the landed [`WLineApprox`] — the shared
// ApproxInt_Approx core generic over the (TheWLine, TheMultiLine) template
// parameters.  This shell is the concrete BRepApprox driver: the Perform
// overloads with the GetType quadric switch and the concrete SvSurfaces
// construction (ApproxInt_Approx.gxx L184-397).
//
// The hxx-only members without a gxx body: the no-argument Perform()
// (BRepApprox_Approx.hxx L113) is declared but has no definition in OCCT
// (there is no BRepApprox_Approx.cxx); it is never called by the boolean
// consumers and has no rcad counterpart.  The hxx default arguments
// (ApproxXYZ = ApproxU1V1 = ApproxU2V2 = true, indicemin = indicemax = 0,
// hxx L90-94/L100-101) are not modeled — the rcad Perform calls take all
// arguments explicitly.

use std::sync::Arc;

use crate::geomalgo::approx_int::{ApproxParamType, MultiBSpCurve, WLineApprox};
use crate::geomalgo::brep_approx::surface_tool;
use crate::geomalgo::brep_approx::{ApproxLine, SvSurfaces, TheMultiLineOfApprox};
use crate::geomalgo::brep_approx_imp_prm::ImpPrmSvSurfaces;
use crate::geomalgo::brep_approx_prm_prm::PrmPrmSvSurfaces;
use crate::geomalgo::int_patch::GeomAbsSurfaceType;
use crate::geomalgo::int_surf::Quadric;
use crate::topalgo::brep_adaptor::surface::BRepAdaptorSurface;

// OCCT ApproxInt_Approx.gxx L26-28: if the quantity of points is less than
// aMinNbPointsForApprox then interpolation is used.
const A_MIN_NB_POINTS_FOR_APPROX: usize = 5;

/// OCCT BRepApprox_Approx — the ApproxInt_Approx.gxx template instantiated
/// over the BRepApprox data layer (BRepApprox_Approx_0.cxx L67-101).
pub struct BRepApproxApprox {
    /// The shared ApproxInt_Approx core: myComputeLine / myComputeLineBezier
    /// / myBezToBSpl / myWithTangency / myTol3d / myTol2d / myDegMin /
    /// myDegMax / myNbIterMax / myTolReached3d / myTolReached2d / myData /
    /// myKnots of ApproxInt_Approx.gxx.
    pub apx: WLineApprox,
}

impl BRepApproxApprox {
    /// OCCT ApproxInt_Approx::ApproxInt_Approx() (L166-180) —
    /// myComputeLine(4, 8, 0.001, 0.001, 5), myComputeLineBezier(4, 8, 0.001,
    /// 0.001, 5), myWithTangency = true, myTol3d = myTol2d = 0.001,
    /// myDegMin = 4, myDegMax = 8, myNbIterMax = 5, myTolReached3d =
    /// myTolReached2d = 0 (WLineApprox::new carries the same defaults).
    /// OCCT L178: myComputeLine.SetContinuity(2) — the landed ComputeLine has
    /// no continuity field (the non-Bezier ComputeLine branch of the engine
    /// is not landed), so the statement has no counterpart.
    pub fn new() -> Self {
        BRepApproxApprox {
            apx: WLineApprox::new(),
        }
    }

    /// OCCT ApproxInt_Approx::Perform(Surf1, Surf2, theline, ApproxXYZ,
    /// ApproxU1V1, ApproxU2V2, indicemin, indicemax) (L226-343) — definition
    /// of the next steps according to the surface types (the quadric switch,
    /// i.e. coordination algorithm).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_prm_prm(
        &mut self,
        surf1: &BRepAdaptorSurface,
        surf2: &BRepAdaptorSurface,
        line: Arc<ApproxLine>,
        approx_xyz: bool,
        approx_u1v1: bool,
        approx_u2v2: bool,
        indicemin: usize,
        indicemax: usize,
    ) {
        // OCCT L236.
        self.apx.my_tol_reached3d = 0.0;
        self.apx.my_tol_reached2d = 0.0;

        // OCCT L238-244.
        let type_s1 = surface_tool::get_type(surf1);
        let type_s2 = surface_tool::get_type(surf2);

        let is_quadric = matches!(
            type_s1,
            GeomAbsSurfaceType::Plane
                | GeomAbsSurfaceType::Cylinder
                | GeomAbsSurfaceType::Sphere
                | GeomAbsSurfaceType::Cone
        ) || matches!(
            type_s2,
            GeomAbsSurfaceType::Plane
                | GeomAbsSurfaceType::Cylinder
                | GeomAbsSurfaceType::Sphere
                | GeomAbsSurfaceType::Cone
        );

        if is_quadric {
            // OCCT L248-294: the quadric switch.  The OCCT default-constructed
            // IntSurf_Quadric followed by Quad.SetValue(...) maps to the
            // direct from_* construction in every case arm.
            let mut second_is_implicit = false;
            let quad: Quadric;
            match type_s1 {
                GeomAbsSurfaceType::Plane => {
                    quad = Quadric::from_plane(&surface_tool::plane(surf1));
                }
                GeomAbsSurfaceType::Cylinder => {
                    quad = Quadric::from_cylinder(&surface_tool::cylinder(surf1));
                }
                GeomAbsSurfaceType::Sphere => {
                    quad = Quadric::from_sphere(&surface_tool::sphere(surf1));
                }
                GeomAbsSurfaceType::Cone => {
                    quad = Quadric::from_cone(&surface_tool::cone(surf1));
                }
                _ => {
                    second_is_implicit = true;
                    match type_s2 {
                        GeomAbsSurfaceType::Plane => {
                            quad = Quadric::from_plane(&surface_tool::plane(surf2));
                        }
                        GeomAbsSurfaceType::Cylinder => {
                            quad = Quadric::from_cylinder(&surface_tool::cylinder(surf2));
                        }
                        GeomAbsSurfaceType::Sphere => {
                            quad = Quadric::from_sphere(&surface_tool::sphere(surf2));
                        }
                        GeomAbsSurfaceType::Cone => {
                            quad = Quadric::from_cone(&surface_tool::cone(surf2));
                        }
                        _ => {
                            // OCCT L288-289: default — Quad stays the default
                            // IntSurf_Quadric (unreachable when isQuadric).
                            quad = Quadric::new();
                        }
                    } // switch (typeS2)
                }
            } // switch (typeS1)

            // OCCT L296-306: Perform(Quad, (SecondIsImplicit ? Surf1 : Surf2),
            // ..., !SecondIsImplicit); return.
            let psurf = if second_is_implicit { surf1 } else { surf2 };
            self.perform_implicit(
                &quad,
                psurf,
                line,
                approx_xyz,
                approx_u1v1,
                approx_u2v2,
                indicemin,
                indicemax,
                !second_is_implicit,
            );
            return;
        }

        // Here, isQuadric == FALSE (OCCT L309-343): the Param-Param perform.
        // OCCT L315: ApproxInt_ThePrmPrmSvSurfaces myPrmPrmSvSurfaces(Surf1,
        // Surf2); the rcad slot carries the lifetime instead of the OCCT
        // void* ptrsvsurf (L332).
        let my_prm_prm_sv_surfaces = PrmPrmSvSurfaces::new(surf1, surf2);

        // OCCT L317-326: the nbpntbez determination.
        let nbpntbez = indicemax - indicemin;
        if nbpntbez < A_MIN_NB_POINTS_FOR_APPROX {
            self.apx.my_bezier_approx = false;
        } else {
            self.apx.my_bezier_approx = true;
        }
        // OCCT L331: const bool cut = myData.myBezierApprox;
        let cut = self.apx.my_bezier_approx;

        // OCCT L335 buildKnots / L342 buildCurve construct the
        // ApproxInt_TheMultiLine inside the engine with the freshly computed
        // myData translation and P2DOnFirst = myData.ApproxU1V1 (L562/L652);
        // the landed engine reads the multi-line's own translation, so the
        // outer multi-line is built here with the translation zeroed (the
        // engine fillData computes the translation itself, and the knot
        // cutting is translation-invariant — landed behavior verified on the
        // GeomInt chain).
        let ml = TheMultiLineOfApprox::new_sv_surfaces(
            line,
            Some(Arc::new(my_prm_prm_sv_surfaces)),
            (approx_xyz as usize), // (ApproxXYZ ? 1 : 0)
            (approx_u1v1 as usize) + (approx_u2v2 as usize),
            approx_u1v1,
            approx_u2v2,
            0.0, // Xo
            0.0, // Yo
            0.0, // Zo
            0.0, // U1o
            0.0, // V1o
            0.0, // U2o
            0.0, // V2o
            approx_u1v1,
            indicemin,
            indicemax,
        );

        // OCCT L334-342: buildKnots(theline, ptrsvsurf); Init(..., cut, ...);
        // buildCurve(theline, ptrsvsurf) — the engine shared sequence.
        self.apx
            .perform_impl(&ml, cut, approx_xyz, approx_u1v1, approx_u2v2, indicemin, indicemax);
    }

    /// OCCT ApproxInt_Approx::Perform(ISurf, PSurf, theline, ApproxXYZ,
    /// ApproxU1V1, ApproxU2V2, indicemin, indicemax, isTheQuadFirst)
    /// (L350-395) — the Analytic-Param perform.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_implicit(
        &mut self,
        i_surf: &Quadric,
        p_surf: &BRepAdaptorSurface,
        line: Arc<ApproxLine>,
        approx_xyz: bool,
        approx_u1v1: bool,
        approx_u2v2: bool,
        indicemin: usize,
        indicemax: usize,
        is_the_quad_first: bool,
    ) {
        // OCCT L364-366: the ImpPrmSvSurfaces selects the first surface by
        // isTheQuadFirst.
        let my_imp_prm_sv_surfaces = if is_the_quad_first {
            ImpPrmSvSurfaces::new_implicit_first(i_surf, p_surf)
        } else {
            ImpPrmSvSurfaces::new_parametric_first(p_surf, i_surf)
        };
        // OCCT L368: SetUseSolver(false).
        my_imp_prm_sv_surfaces.set_use_solver(false);
        // OCCT L381: void* ptrsvsurf = &myImpPrmSvSurfaces.
        let ptrsvsurf = Arc::new(my_imp_prm_sv_surfaces);

        // OCCT L370-378: the nbpntbez determination.
        let nbpntbez = indicemax - indicemin;
        if nbpntbez < A_MIN_NB_POINTS_FOR_APPROX {
            self.apx.my_bezier_approx = false;
        } else {
            self.apx.my_bezier_approx = true;
        }
        // OCCT L380: const bool cut = myData.myBezierApprox;
        let cut = self.apx.my_bezier_approx;

        // OCCT L384 fillData / L387 buildKnots / L389-394 Init + buildCurve.
        // The SvSurfaces functor is constructed (L364-368) before fillData
        // runs inside the engine sequence — the OCCT order is preserved.
        // The outer multi-line carries the zero translation and
        // P2DOnFirst = ApproxU1V1 (L562/L652) like the PrmPrm perform.
        let ml = TheMultiLineOfApprox::new_sv_surfaces(
            line,
            Some(ptrsvsurf),
            (approx_xyz as usize),
            (approx_u1v1 as usize) + (approx_u2v2 as usize),
            approx_u1v1,
            approx_u2v2,
            0.0, // Xo
            0.0, // Yo
            0.0, // Zo
            0.0, // U1o
            0.0, // V1o
            0.0, // U2o
            0.0, // V2o
            approx_u1v1,
            indicemin,
            indicemax,
        );
        self.apx
            .perform_impl(&ml, cut, approx_xyz, approx_u1v1, approx_u2v2, indicemin, indicemax);
    }

    /// OCCT ApproxInt_Approx::Perform(theline, ApproxXYZ, ApproxU1V1,
    /// ApproxU2V2, indicemin, indicemax) (L184-218) — the WLine perform.
    /// The multi-line is built without the SvSurfaces (the OCCT NULL
    /// PtrOnmySvSurfaces), so WhatStatus reports NoPointsAdded and no extra
    /// points are added on the line.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_wline(
        &mut self,
        line: Arc<ApproxLine>,
        approx_xyz: bool,
        approx_u1v1: bool,
        approx_u2v2: bool,
        indicemin: usize,
        indicemax: usize,
    ) {
        let ml = TheMultiLineOfApprox::new_line(
            line,
            (approx_xyz as usize),
            (approx_u1v1 as usize) + (approx_u2v2 as usize),
            approx_u1v1,
            approx_u2v2,
            0.0, // Xo
            0.0, // Yo
            0.0, // Zo
            0.0, // U1o
            0.0, // V1o
            0.0, // U2o
            0.0, // V2o
            approx_u1v1,
            indicemin,
            indicemax,
        );
        // OCCT L192 prepareDS, L201 fillData, L204 buildKnots(theline, NULL),
        // L205-210 the two-knot split, L212-215 Init(..., true, ...), L217
        // buildCurve — the engine Perform(WLine) sequence.
        self.apx.perform(&ml, approx_xyz, approx_u1v1, approx_u2v2, indicemin, indicemax);
    }

    /// OCCT ApproxInt_Approx::SetParameters (L399-427) — forwarded to the
    /// engine (the RatioTol 1.5 division, the PassPoint constraints for
    /// !ApproxWithTangency and myBezierApprox = true are the engine body).
    #[allow(clippy::too_many_arguments)]
    pub fn set_parameters(
        &mut self,
        tol3d: f64,
        tol2d: f64,
        deg_min: usize,
        deg_max: usize,
        nb_iter_max: i32,
        nb_pnt_max: usize,
        approx_with_tangency: bool,
        parametrization: ApproxParamType,
    ) {
        self.apx.set_parameters(
            tol3d,
            tol2d,
            deg_min,
            deg_max,
            nb_iter_max,
            nb_pnt_max,
            approx_with_tangency,
            parametrization,
        );
    }

    /// OCCT ApproxInt_Approx::NbMultiCurves (L431-434) — returns 1.
    pub fn nb_multi_curves(&self) -> usize {
        self.apx.nb_multi_curves()
    }

    /// OCCT ApproxInt_Approx::TolReached3d (L459-462) — myTolReached3d *
    /// RatioTol (1.5).
    pub fn tol_reached3d(&self) -> f64 {
        self.apx.tol_reached3d()
    }

    /// OCCT ApproxInt_Approx::TolReached2d (L466-469) — myTolReached2d *
    /// RatioTol (1.5).
    pub fn tol_reached2d(&self) -> f64 {
        self.apx.tol_reached2d()
    }

    /// OCCT ApproxInt_Approx::IsDone (L473-483) — myComputeLineBezier.
    /// NbMultiCurves() > 0 in the Bezier regime (the engine sets my_done at
    /// the end of the shared Perform sequence).
    pub fn is_done(&self) -> bool {
        self.apx.is_done()
    }

    /// OCCT ApproxInt_Approx::Value(Index) (L487-497) — the Index argument is
    /// ignored (NbMultiCurves == 1); the Bezier regime returns the
    /// myBezToBSpl spline.
    pub fn value(&self) -> MultiBSpCurve {
        self.apx.value()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use glam::{DVec2, DVec3};
    use rcad_kernel::geom::{CylindricalSurface, Line3, Plane, Surface3};
    use rcad_kernel::topo::topods::{BRepBuilder, Shape};

    use super::*;
    use crate::geomalgo::approx_int::ApproxParamType;
    use crate::geomalgo::brep_approx::ApproxLine;
    use crate::geomalgo::int_surf::{LineOn2S, PntOn2S};

    /// A plane face over the argument surface frame and UV window
    /// [umin, umax, vmin, vmax] (the build_plane_brep pattern of
    /// brep_approx/mod.rs tests).
    fn build_plane_face(
        origin: DVec3,
        normal: DVec3,
        u_dir: DVec3,
        v_dir: DVec3,
        window: [f64; 4],
    ) -> (rcad_kernel::BRep, Shape) {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(&mut brep, origin, 1e-7);
        let v2 = b.add_vertex(&mut brep, origin + u_dir, 1e-7);
        let e = b.add_edge(
            &mut brep,
            Some(rcad_kernel::geom::Curve3::Line(Line3 {
                origin,
                direction: u_dir,
            })),
            v1,
            v2,
            [0.0, 1.0],
        );
        let wire = brep.add_twire(vec![e]);
        let face = brep.add_tface(
            Some(Surface3::Plane(Plane {
                origin,
                normal,
                u_dir,
                v_dir,
            })),
            wire,
            Vec::new(),
            None,
            Some(window),
            Vec::new(),
            true,
        );
        (brep, face)
    }

    /// A unit cylinder face around Z over the UV window [0, 2] x [0, 1]
    /// (u = angle in radians, v = z).
    fn build_cylinder_face_z() -> (rcad_kernel::BRep, Shape) {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(&mut brep, DVec3::new(1.0, 0.0, 0.0), 1e-7);
        let v2 = b.add_vertex(&mut brep, DVec3::new(1.0, 0.0, 1.0), 1e-7);
        let e = b.add_edge(
            &mut brep,
            Some(rcad_kernel::geom::Curve3::Line(Line3::new(DVec3::new(1.0, 0.0, 0.0), DVec3::new(0.0, 0.0, 1.0)))),
            v1,
            v2,
            [0.0, 1.0],
        );
        let wire = brep.add_twire(vec![e]);
        let face = brep.add_tface(
            Some(Surface3::Cylinder(CylindricalSurface {
                origin: DVec3::ZERO,
                axis: DVec3::new(0.0, 0.0, 1.0),
                radius: 1.0,
                ref_dir: DVec3::new(1.0, 0.0, 0.0),
                y_dir: None,
            })),
            wire,
            Vec::new(),
            None,
            Some([0.0, 2.0, 0.0, 1.0]),
            Vec::new(),
            true,
        );
        (brep, face)
    }

    /// Walking line of `n` points on the pair of planes z = 0 (u = X, v = Y)
    /// and x = 0 (u = Y, v = Z): point i (1-based) is p(t) = (0, w(t), 0)
    /// with t = (i-1)/(n-1) and w(t) = 0.1 + 0.6 t + 0.2 t^2 — the parabola
    /// along the intersection line (the Y axis) of the two planes.  UV1 =
    /// (0, w), UV2 = (w, 0).
    fn walking_line_two_planes(n: usize) -> Arc<LineOn2S> {
        let mut line = LineOn2S::new();
        for i in 0..n {
            let t = i as f64 / (n - 1) as f64;
            let w = 0.1 + 0.6 * t + 0.2 * t * t;
            let mut p = PntOn2S::new();
            p.set_value_all(DVec3::new(0.0, w, 0.0), 0.0, w, w, 0.0);
            line.add(&p);
        }
        Arc::new(line)
    }

    /// The endpoint assertions shared by the anchors: with the endpoint
    /// constraints the first/last poles are the exact walking-line endpoints
    /// in 3d and on both pcurves (after the engine un-translation).
    fn assert_endpoints_two_planes(value: &MultiBSpCurve) {
        let poles = &value.poles;
        assert!(poles.len() >= 2);
        let first = &poles[0];
        let last = &poles[poles.len() - 1];
        // 3D endpoints: (0, 0.1, 0) and (0, 0.9, 0).
        assert!(
            first.p3d[0].distance(DVec3::new(0.0, 0.1, 0.0)) < 1e-8,
            "first pole {:?}",
            first.p3d[0]
        );
        assert!(
            last.p3d[0].distance(DVec3::new(0.0, 0.9, 0.0)) < 1e-8,
            "last pole {:?}",
            last.p3d[0]
        );
        // UV1 on the plane z = 0 (u = X, v = Y): (0, 0.1) and (0, 0.9).
        assert!(first.p2d[0].distance(DVec2::new(0.0, 0.1)) < 1e-8);
        assert!(last.p2d[0].distance(DVec2::new(0.0, 0.9)) < 1e-8);
        // UV2 on the plane x = 0 (u = Y, v = Z): (0.1, 0) and (0.9, 0).
        assert!(first.p2d[1].distance(DVec2::new(0.1, 0.0)) < 1e-8);
        assert!(last.p2d[1].distance(DVec2::new(0.9, 0.0)) < 1e-8);
    }

    fn set_consumer_parameters(apx: &mut BRepApproxApprox) {
        // The consumer pattern (IntTools_FaceFace::MakeCurve): SetParameters
        // before Perform; the tolerances are divided by RatioTol internally
        // (ApproxInt_Approx.gxx L399-427).
        apx.set_parameters(1.0e-6, 1.0e-6, 4, 8, 40, 30, true, ApproxParamType::ChordLength);
    }

    /// Anchor: Perform(Surf1, Surf2, ...) over two orthogonal plane faces —
    /// both GetType's are GeomAbs_Plane, so the quadric switch dispatches to
    /// Perform(Quad = Plane(Surf1), Surf2, ..., isTheQuadFirst = true) and
    /// the approximation completes end-to-end with the exact endpoint poles
    /// in 3d and on both pcurves (ApproxInt_Approx.gxx L226-343).
    #[test]
    fn approx_prm_prm_two_planes_end_to_end() {
        let (brep1, face1) = build_plane_face(
            DVec3::ZERO,
            DVec3::new(0.0, 0.0, 1.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            [0.0, 1.0, 0.0, 1.0],
        );
        let (brep2, face2) = build_plane_face(
            DVec3::ZERO,
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.0, 0.0, 1.0),
            [0.0, 1.0, 0.0, 1.0],
        );
        let s1 = BRepAdaptorSurface::initialize_face(&brep1, &face1, true);
        let s2 = BRepAdaptorSurface::initialize_face(&brep2, &face2, true);
        let line = Arc::new(ApproxLine::new_line_on_2s(Some(walking_line_two_planes(9)), false));

        let mut apx = BRepApproxApprox::new();
        set_consumer_parameters(&mut apx);
        apx.perform_prm_prm(&s1, &s2, line, true, true, true, 1, 9);
        assert!(apx.is_done(), "prm-prm approximation not done");
        // TolReached3d/2d carry the RatioTol factor (L459-469).
        assert!(apx.tol_reached3d() >= 0.0 && apx.tol_reached2d() >= 0.0);
        assert_eq!(apx.nb_multi_curves(), 1); // L431-434: returns 1
        assert_endpoints_two_planes(&apx.value());
    }

    /// Anchor: Perform(ISurf, PSurf, ..., isTheQuadFirst = true) — the plane
    /// quadric z = 0 against the tilted plane face P(u, v) = u X +
    /// v (Y + Z)/sqrt(2); the intersection is the line v = 0 (3D along X)
    /// and the approximation completes end-to-end with the exact endpoint
    /// poles (ApproxInt_Approx.gxx L350-395).
    #[test]
    fn approx_implicit_plane_quadric_end_to_end() {
        // The implicit quadric: the plane z = 0 (u = X, v = Y).
        let quad = Quadric::from_plane(&Plane {
            origin: DVec3::ZERO,
            normal: DVec3::new(0.0, 0.0, 1.0),
            u_dir: DVec3::new(1.0, 0.0, 0.0),
            v_dir: DVec3::new(0.0, 1.0, 0.0),
        });
        // The parametric plane face over the UV window [0, 1]^2.
        let u_dir = DVec3::new(1.0, 0.0, 0.0);
        let v_dir = DVec3::new(0.0, 1.0, 1.0).normalize();
        let (brep, face) = build_plane_face(DVec3::ZERO, u_dir.cross(v_dir).normalize(), u_dir, v_dir, [0.0, 1.0, 0.0, 1.0]);
        let p_surf = BRepAdaptorSurface::initialize_face(&brep, &face, true);

        // Walking line: p(t) = (w(t), 0, 0) with w(t) = 0.1 + 0.6 t + 0.2 t^2
        // — on the intersection line of the two planes; UV1 (the quadric
        // plane) = (w, 0), UV2 (the parametric plane) = (w, 0).
        let n = 9usize;
        let mut line = LineOn2S::new();
        for i in 0..n {
            let t = i as f64 / (n - 1) as f64;
            let w = 0.1 + 0.6 * t + 0.2 * t * t;
            let mut p = PntOn2S::new();
            p.set_value_all(DVec3::new(w, 0.0, 0.0), w, 0.0, w, 0.0);
            line.add(&p);
        }
        let line = Arc::new(ApproxLine::new_line_on_2s(Some(Arc::new(line)), false));

        let mut apx = BRepApproxApprox::new();
        set_consumer_parameters(&mut apx);
        apx.perform_implicit(&quad, &p_surf, line, true, true, true, 1, n, true);
        assert!(apx.is_done(), "implicit-prm approximation not done");

        let poles = apx.value().poles;
        assert!(poles.len() >= 2);
        let first = &poles[0];
        let last = &poles[poles.len() - 1];
        // 3D endpoints: (0.1, 0, 0) and (0.9, 0, 0).
        assert!(
            first.p3d[0].distance(DVec3::new(0.1, 0.0, 0.0)) < 1e-8,
            "first pole {:?}",
            first.p3d[0]
        );
        assert!(
            last.p3d[0].distance(DVec3::new(0.9, 0.0, 0.0)) < 1e-8,
            "last pole {:?}",
            last.p3d[0]
        );
        // UV1 on the quadric plane: (0.1, 0) and (0.9, 0).
        assert!(first.p2d[0].distance(DVec2::new(0.1, 0.0)) < 1e-8);
        assert!(last.p2d[0].distance(DVec2::new(0.9, 0.0)) < 1e-8);
        // UV2 on the parametric plane: (0.1, 0) and (0.9, 0).
        assert!(first.p2d[1].distance(DVec2::new(0.1, 0.0)) < 1e-8);
        assert!(last.p2d[1].distance(DVec2::new(0.9, 0.0)) < 1e-8);
    }

    /// Anchor: the quadric switch dispatch of Perform(Surf1, Surf2, ...).
    /// Two planes: GetType(S1) = GeomAbs_Plane selects Quad = Plane(Surf1)
    /// with SecondIsImplicit = false, i.e. the switch must dispatch to
    /// Perform(Quad, Surf2, ..., isTheQuadFirst = true) — running that
    /// overload directly must reproduce the very same curve.  The
    /// plane x cylinder arm exercises the same dispatch with an implicit
    /// quadric first against the cylinder PSurf (no cylinder-face BRep
    /// helper existed; the face is built inline over the analytic window).
    #[test]
    fn approx_quadric_switch_dispatch() {
        // Arm 1: plane x plane — dispatch equivalence.
        let (brep1, face1) = build_plane_face(
            DVec3::ZERO,
            DVec3::new(0.0, 0.0, 1.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            [0.0, 1.0, 0.0, 1.0],
        );
        let (brep2, face2) = build_plane_face(
            DVec3::ZERO,
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.0, 0.0, 1.0),
            [0.0, 1.0, 0.0, 1.0],
        );
        let s1 = BRepAdaptorSurface::initialize_face(&brep1, &face1, true);
        let s2 = BRepAdaptorSurface::initialize_face(&brep2, &face2, true);
        let line = Arc::new(ApproxLine::new_line_on_2s(Some(walking_line_two_planes(9)), false));

        let mut apx = BRepApproxApprox::new();
        set_consumer_parameters(&mut apx);
        apx.perform_prm_prm(&s1, &s2, line.clone(), true, true, true, 1, 9);
        assert!(apx.is_done(), "dispatched approximation not done");

        // The same computation through the explicit overload: Quad =
        // Plane(Surf1), PSurf = Surf2, isTheQuadFirst = true (OCCT
        // L296-304 with SecondIsImplicit = false).
        let quad = Quadric::from_plane(&surface_tool::plane(&s1));
        let mut apx2 = BRepApproxApprox::new();
        set_consumer_parameters(&mut apx2);
        apx2.perform_implicit(&quad, &s2, line, true, true, true, 1, 9, true);
        assert!(apx2.is_done(), "explicit implicit-perform not done");

        // Both runs take the identical code path, so the resulting curves
        // are identical.
        let v1 = apx.value();
        let v2 = apx2.value();
        assert_eq!(v1.poles.len(), v2.poles.len());
        for (a, b) in v1.poles.iter().zip(v2.poles.iter()) {
            assert_eq!(a.p3d, b.p3d);
            assert_eq!(a.p2d, b.p2d);
        }
        assert_endpoints_two_planes(&v1);

        // Arm 2: plane x cylinder — typeS1 = Plane again dispatches to the
        // implicit perform, this time with the cylinder as PSurf.
        let (brep3, face3) = build_plane_face(
            DVec3::ZERO,
            DVec3::new(0.0, 0.0, 1.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            [-1.2, 1.2, -1.2, 1.2],
        );
        let (brep4, face4) = build_cylinder_face_z();
        let s3 = BRepAdaptorSurface::initialize_face(&brep3, &face3, true);
        let s4 = BRepAdaptorSurface::initialize_face(&brep4, &face4, true);

        // Walking line: the circle arc theta in [0.3, 1.3] of the unit
        // cylinder x^2 + y^2 = 1 in the plane z = 0; UV1 (the plane) =
        // (cos, sin), UV2 (the cylinder) = (theta, 0).
        let n = 9usize;
        let mut line = LineOn2S::new();
        for i in 0..n {
            let theta = 0.3 + i as f64 / (n - 1) as f64;
            let mut p = PntOn2S::new();
            p.set_value_all(
                DVec3::new(theta.cos(), theta.sin(), 0.0),
                theta.cos(),
                theta.sin(),
                theta,
                0.0,
            );
            line.add(&p);
        }
        let line = Arc::new(ApproxLine::new_line_on_2s(Some(Arc::new(line)), false));

        let mut apx3 = BRepApproxApprox::new();
        set_consumer_parameters(&mut apx3);
        apx3.perform_prm_prm(&s3, &s4, line, true, true, true, 1, n);
        assert!(apx3.is_done(), "plane x cylinder approximation not done");

        let poles = apx3.value().poles;
        assert!(poles.len() >= 2);
        let first = &poles[0];
        let last = &poles[poles.len() - 1];
        let (t0, t1) = (0.3f64, 1.3f64);
        // 3D endpoints on the circle.
        assert!(
            first.p3d[0].distance(DVec3::new(t0.cos(), t0.sin(), 0.0)) < 1e-8,
            "first pole {:?}",
            first.p3d[0]
        );
        assert!(
            last.p3d[0].distance(DVec3::new(t1.cos(), t1.sin(), 0.0)) < 1e-8,
            "last pole {:?}",
            last.p3d[0]
        );
        // UV1 on the plane: (cos, sin).
        assert!(first.p2d[0].distance(DVec2::new(t0.cos(), t0.sin())) < 1e-8);
        assert!(last.p2d[0].distance(DVec2::new(t1.cos(), t1.sin())) < 1e-8);
        // UV2 on the cylinder: (theta, 0).
        assert!(first.p2d[1].distance(DVec2::new(t0, 0.0)) < 1e-8);
        assert!(last.p2d[1].distance(DVec2::new(t1, 0.0)) < 1e-8);
    }

    /// Anchor: Perform(theline, ...) without the SvSurfaces — the multi-line
    /// reports WhatStatus = NoPointsAdded, no extra points are added, and the
    /// approximation of the 9-point line completes with the exact endpoint
    /// poles (ApproxInt_Approx.gxx L184-218 with theLine = BRepApprox_ApproxLine).
    #[test]
    fn approx_wline_no_sv_path() {
        let line = Arc::new(ApproxLine::new_line_on_2s(Some(walking_line_two_planes(9)), false));
        let mut apx = BRepApproxApprox::new();
        set_consumer_parameters(&mut apx);
        apx.perform_wline(line, true, true, true, 1, 9);
        assert!(apx.is_done(), "wline approximation not done");
        assert_eq!(apx.nb_multi_curves(), 1);
        assert_endpoints_two_planes(&apx.value());
    }
}
