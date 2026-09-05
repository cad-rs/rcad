// OCCT BRepApprox (TKTopAlgo) — the data/tool layer of the intersection-line
// approximation chain (Stage 3c-2 of the HLR port).  1:1 translation of:
//   BRepApprox_TheMultiLineOfApprox.hxx L1-152, whose bodies are the
//     ApproxInt_MultiLine.gxx template instantiation (gxx L32-778),
//   BRepApprox_TheMultiLineToolOfApprox.hxx L30-123, whose bodies are the
//     ApproxInt_MultiLineTool.lxx one-liners (lxx L29-184),
//   BRepApprox_SurfaceTool.hxx L37-164 + BRepApprox_SurfaceTool.lxx
//     L31-251 + BRepApprox_SurfaceTool.cxx L22-40,
//   BRepApprox_ApproxLine.hxx L30-54 + BRepApprox_ApproxLine.cxx L29-94.
//
// The MultiLine carries a `void* PtrOnmySvSurfaces` to an ApproxInt_SvSurfaces
// (BRepApprox_TheMultiLineOfApprox.hxx L33-34/L134); that abstract OCCT class
// is translated here as the object-safe [`SvSurfaces`] trait (the same
// encoding used for Adaptor3d_Surface in hlr::contap::surface_adaptor), so the
// pointer maps to `Option<Arc<dyn SvSurfaces>>`.
//
// OCCT handle<T> fields map to Arc<T> (Option-wrapped for nullable handles);
// NCollection_Array1 out-parameters map to slices whose index carries the
// OCCT Lower() offset; local arrays over [Low, High] map to a Vec indexed by
// `i - low`.

use std::sync::Arc;

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{BSplineCurve2, BSplineCurve3};
use rcad_kernel::math::math_matrix::Vector;
use rcad_kernel::precision::SQUARE_CONFUSION;

use super::approx_int::ApproxStatus;
use super::int_surf::{LineOn2S, PntOn2S};

// ---------------------------------------------------------------------------
// ApproxInt_SvSurfaces — the abstract base the MultiLine's void* points to.
// ---------------------------------------------------------------------------

/// OCCT ApproxInt_SvSurfaces (ApproxInt/ApproxInt_SvSurfaces.hxx L41-99) —
/// root class for classes dedicated to calculate 2d and 3d points and
/// tangents of intersection lines of two surfaces of different types for
/// given u, v parameters of intersection point on two surfaces.
///
/// The field myUseSolver manages the calculation type (hxx L35-40): when
/// true, the input u1, v1, u2, v2 are the first approximation of the exact
/// intersection point and are refined by the solver; when false they are
/// already "exact".  The OCCT base-class field becomes implementor state
/// accessed through the trait's set_use_solver / get_use_solver virtuals
/// (hxx L93-95).
pub trait SvSurfaces {
    /// OCCT ApproxInt_SvSurfaces::Compute (hxx L52-59) — returns True if Tg,
    /// Tguv1, Tguv2 can be computed; the in/out u1..v2 carry the refined
    /// intersection parameters.
    #[allow(clippy::too_many_arguments)]
    fn compute(
        &self,
        u1: &mut f64,
        v1: &mut f64,
        u2: &mut f64,
        v2: &mut f64,
        pt: &mut DVec3,
        tg: &mut DVec3,
        tguv1: &mut DVec2,
        tguv2: &mut DVec2,
    ) -> bool;

    /// OCCT ApproxInt_SvSurfaces::Pnt (hxx L61-65).
    fn pnt(&self, u1: f64, v1: f64, u2: f64, v2: f64, p: &mut DVec3);

    /// OCCT ApproxInt_SvSurfaces::SeekPoint (hxx L68-72) — computes point on
    /// curve and parameters on the surfaces.
    fn seek_point(&self, u1: f64, v1: f64, u2: f64, v2: f64, point: &mut PntOn2S) -> bool;

    /// OCCT ApproxInt_SvSurfaces::Tangency (hxx L74-78).
    fn tangency(&self, u1: f64, v1: f64, u2: f64, v2: f64, tg: &mut DVec3) -> bool;

    /// OCCT ApproxInt_SvSurfaces::TangencyOnSurf1 (hxx L80-84).
    fn tangency_on_surf_1(&self, u1: f64, v1: f64, u2: f64, v2: f64, tg: &mut DVec2) -> bool;

    /// OCCT ApproxInt_SvSurfaces::TangencyOnSurf2 (hxx L86-90).
    fn tangency_on_surf_2(&self, u1: f64, v1: f64, u2: f64, v2: f64, tg: &mut DVec2) -> bool;

    /// OCCT ApproxInt_SvSurfaces::SetUseSolver (hxx L93).
    fn set_use_solver(&self, the_use_sol: bool);

    /// OCCT ApproxInt_SvSurfaces::GetUseSolver (hxx L95).
    fn get_use_solver(&self) -> bool;
}

/// OCCT `void* PtrOnmySvSurfaces` over an ApproxInt_SvSurfaces — shared
/// access to the surfaces functor (BRepApprox_TheMultiLineOfApprox.hxx
/// L134).
pub type SharedSvSurfaces = Arc<dyn SvSurfaces>;

// ---------------------------------------------------------------------------
// BRepApprox_ApproxLine
// ---------------------------------------------------------------------------

/// OCCT BRepApprox_ApproxLine (BRepApprox_ApproxLine.hxx L30-54) — holds the
/// result of an intersection line: either three B-spline curves (the 3d
/// curve and the two pcurves) or an IntSurf_LineOn2S of points on two
/// surfaces.  The OCCT Standard_Transient handle maps to Arc.
#[derive(Debug, Clone)]
pub struct ApproxLine {
    /// OCCT myCurveXYZ (hxx L50).
    my_curve_xyz: Option<Arc<BSplineCurve3>>,
    /// OCCT myCurveUV1 (hxx L51).
    my_curve_uv1: Option<Arc<BSplineCurve2>>,
    /// OCCT myCurveUV2 (hxx L52).
    my_curve_uv2: Option<Arc<BSplineCurve2>>,
    /// OCCT myLineOn2S (hxx L53).
    my_line_on_2s: Option<Arc<LineOn2S>>,
}

impl ApproxLine {
    /// OCCT BRepApprox_ApproxLine::BRepApprox_ApproxLine(CurveXYZ, CurveUV1,
    /// CurveUV2) (hxx L34-36; cxx L29-36).
    pub fn new(
        curve_xyz: Option<Arc<BSplineCurve3>>,
        curve_uv1: Option<Arc<BSplineCurve2>>,
        curve_uv2: Option<Arc<BSplineCurve2>>,
    ) -> Self {
        ApproxLine {
            my_curve_xyz: curve_xyz,
            my_curve_uv1: curve_uv1,
            my_curve_uv2: curve_uv2,
            my_line_on_2s: None,
        }
    }

    /// OCCT BRepApprox_ApproxLine::BRepApprox_ApproxLine(lin, theTang)
    /// (hxx L40-41; cxx L40-43) — theTang has been entered only for
    /// compatibility with the alias IntPatch_WLine and is not used in this
    /// class (hxx L38-39); the OCCT .cxx leaves the parameter unnamed.
    pub fn new_line_on_2s(lin: Option<Arc<LineOn2S>>, _the_tang: bool) -> Self {
        ApproxLine {
            my_curve_xyz: None,
            my_curve_uv1: None,
            my_curve_uv2: None,
            my_line_on_2s: lin,
        }
    }

    /// OCCT BRepApprox_ApproxLine::NbPnts (hxx L43; cxx L47-62).
    pub fn nb_pnts(&self) -> usize {
        if let Some(c) = &self.my_curve_xyz {
            // OCCT return (myCurveXYZ->NbPoles());
            return c.control_points.len();
        }
        if let Some(c) = &self.my_curve_uv1 {
            // OCCT return (myCurveUV1->NbPoles());
            return c.control_points.len();
        }
        if let Some(c) = &self.my_curve_uv2 {
            // OCCT return (myCurveUV2->NbPoles());
            return c.control_points.len();
        }
        // OCCT return (myLineOn2S->NbPoints());
        self.my_line_on_2s
            .as_ref()
            .expect("BRepApprox_ApproxLine: null line")
            .nb_points()
    }

    /// OCCT BRepApprox_ApproxLine::Point(Index) (hxx L45; cxx L66-94) —
    /// Index is 1-based.
    pub fn point(&self, index: usize) -> PntOn2S {
        if let Some(line) = &self.my_line_on_2s {
            if line.nb_points() != 0 {
                // OCCT return (myLineOn2S->Value(Index)); the rcad
                // LineOn2S::value is 0-based.
                return line.value(index - 1).clone();
            }
        }
        // OCCT L75-88: gp_Pnt2d P1, P2; gp_Pnt P; the Pole(Index) reads
        // (Pole is 1-based in OCCT).
        let mut p1 = DVec2::ZERO;
        let mut p2 = DVec2::ZERO;
        let mut p = DVec3::ZERO;
        if let Some(c) = &self.my_curve_xyz {
            p = c.control_points[index - 1];
        }
        if let Some(c) = &self.my_curve_uv1 {
            p1 = c.control_points[index - 1];
        }
        if let Some(c) = &self.my_curve_uv2 {
            p2 = c.control_points[index - 1];
        }

        // OCCT L90-93: IntSurf_PntOn2S aPntOn2S;
        // aPntOn2S.SetValue(P, P1.X(), P1.Y(), P2.X(), P2.Y());
        let mut a_pnt_on_2s = PntOn2S::new();
        a_pnt_on_2s.set_value_all(p, p1.x, p1.y, p2.x, p2.y);

        a_pnt_on_2s
    }
}

// ---------------------------------------------------------------------------
// BRepApprox_TheMultiLineOfApprox
// ---------------------------------------------------------------------------

/// OCCT BRepApprox_TheMultiLineOfApprox (BRepApprox_TheMultiLineOfApprox.hxx
/// L36-150) — the ApproxInt_MultiLine instantiation over
/// BRepApprox_ApproxLine.  The class SvSurfaces is used when the
/// approximation algorithm needs some extra points on the line; a new line is
/// then created which shares the same surfaces and functions (hxx L43-48).
#[derive(Clone)]
pub struct TheMultiLineOfApprox {
    /// OCCT PtrOnmySvSurfaces (hxx L134).
    ptr_on_my_sv_surfaces: Option<SharedSvSurfaces>,
    /// OCCT myLine (hxx L135).
    my_line: Option<Arc<ApproxLine>>,
    /// OCCT indicemin (hxx L136).
    indicemin: usize,
    /// OCCT indicemax (hxx L137).
    indicemax: usize,
    /// OCCT nbp3d (hxx L138).
    nbp3d: usize,
    /// OCCT nbp2d (hxx L139).
    nbp2d: usize,
    /// OCCT myApproxU1V1 (hxx L140).
    my_approx_u1v1: bool,
    /// OCCT myApproxU2V2 (hxx L141).
    my_approx_u2v2: bool,
    /// OCCT p2donfirst (hxx L142).
    p2donfirst: bool,
    /// OCCT Xo (hxx L143).
    xo: f64,
    /// OCCT Yo (hxx L144).
    yo: f64,
    /// OCCT Zo (hxx L145).
    zo: f64,
    /// OCCT U1o (hxx L146).
    u1o: f64,
    /// OCCT V1o (hxx L147).
    v1o: f64,
    /// OCCT U2o (hxx L148).
    u2o: f64,
    /// OCCT V2o (hxx L149).
    v2o: f64,
}

impl TheMultiLineOfApprox {
    /// OCCT BRepApprox_TheMultiLineOfApprox() (hxx L41; gxx L32-50).
    pub fn new() -> Self {
        TheMultiLineOfApprox {
            ptr_on_my_sv_surfaces: None, // OCCT PtrOnmySvSurfaces = NULL
            my_line: None,               // OCCT myLine = NULL
            indicemin: 0,
            indicemax: 0,
            nbp3d: 0,
            nbp2d: 0,
            my_approx_u1v1: false,
            my_approx_u2v2: false,
            p2donfirst: true,
            xo: 0.0,
            yo: 0.0,
            zo: 0.0,
            u1o: 0.0,
            v1o: 0.0,
            u2o: 0.0,
            v2o: 0.0,
        }
    }

    /// OCCT BRepApprox_TheMultiLineOfApprox(line, PtrSvSurfaces, NbP3d,
    /// NbP2d, ApproxU1V1, ApproxU2V2, xo, yo, zo, u1o, v1o, u2o, v2o,
    /// P2DOnFirst, IndMin, IndMax) (hxx L49-64; gxx L54-94).
    #[allow(clippy::too_many_arguments)]
    pub fn new_sv_surfaces(
        line: Arc<ApproxLine>,
        sv_surfaces: Option<SharedSvSurfaces>,
        nb_p3d: usize,
        nb_p2d: usize,
        approx_u1v1: bool,
        approx_u2v2: bool,
        xo: f64,
        yo: f64,
        zo: f64,
        u1o: f64,
        v1o: f64,
        u2o: f64,
        v2o: f64,
        p2d_on_first: bool,
        ind_min: usize,
        ind_max: usize,
    ) -> Self {
        TheMultiLineOfApprox {
            ptr_on_my_sv_surfaces: sv_surfaces,
            my_line: Some(line),
            indicemin: ind_min.min(ind_max), // OCCT std::min(IndMin, IndMax)
            indicemax: ind_min.max(ind_max), // OCCT std::max(IndMin, IndMax)
            nbp3d: nb_p3d,
            nbp2d: nb_p2d,
            my_approx_u1v1: approx_u1v1,
            my_approx_u2v2: approx_u2v2,
            p2donfirst: p2d_on_first,
            xo,
            yo,
            zo,
            u1o,
            v1o,
            u2o,
            v2o,
        }
    }

    /// OCCT BRepApprox_TheMultiLineOfApprox(line, NbP3d, NbP2d, ApproxU1V1,
    /// ApproxU2V2, xo, yo, zo, u1o, v1o, u2o, v2o, P2DOnFirst, IndMin,
    /// IndMax) (hxx L67-81; gxx L98-136) — no extra points will be added on
    /// the current line.
    #[allow(clippy::too_many_arguments)]
    pub fn new_line(
        line: Arc<ApproxLine>,
        nb_p3d: usize,
        nb_p2d: usize,
        approx_u1v1: bool,
        approx_u2v2: bool,
        xo: f64,
        yo: f64,
        zo: f64,
        u1o: f64,
        v1o: f64,
        u2o: f64,
        v2o: f64,
        p2d_on_first: bool,
        ind_min: usize,
        ind_max: usize,
    ) -> Self {
        TheMultiLineOfApprox {
            ptr_on_my_sv_surfaces: None, // OCCT PtrOnmySvSurfaces(nullptr)
            my_line: Some(line),
            indicemin: ind_min.min(ind_max),
            indicemax: ind_min.max(ind_max),
            nbp3d: nb_p3d,
            nbp2d: nb_p2d,
            my_approx_u1v1: approx_u1v1,
            my_approx_u2v2: approx_u2v2,
            p2donfirst: p2d_on_first,
            xo,
            yo,
            zo,
            u1o,
            v1o,
            u2o,
            v2o,
        }
    }

    /// OCCT BRepApprox_TheMultiLineOfApprox::FirstPoint (hxx L83; gxx
    /// L140-143).
    pub fn first_point(&self) -> usize {
        self.indicemin
    }

    /// OCCT BRepApprox_TheMultiLineOfApprox::LastPoint (hxx L85; gxx
    /// L147-150).
    pub fn last_point(&self) -> usize {
        self.indicemax
    }

    /// OCCT BRepApprox_TheMultiLineOfApprox::NbP2d (hxx L88; gxx L171-174) —
    /// the number of 2d points of a TheLine.
    pub fn nb_p2d(&self) -> usize {
        self.nbp2d
    }

    /// OCCT BRepApprox_TheMultiLineOfApprox::NbP3d (hxx L91; gxx L164-167) —
    /// the number of 3d points of a TheLine.
    pub fn nb_p3d(&self) -> usize {
        self.nbp3d
    }

    /// OCCT BRepApprox_TheMultiLineOfApprox::WhatStatus (hxx L93; gxx
    /// L154-160).
    pub fn what_status(&self) -> ApproxStatus {
        if self.ptr_on_my_sv_surfaces.is_some() {
            ApproxStatus::PointsAdded
        } else {
            ApproxStatus::NoPointsAdded
        }
    }

    /// OCCT BRepApprox_TheMultiLineOfApprox::Value(MPointIndex, tabPt)
    /// (hxx L96; gxx L178-182) — the 3d points of the multipoint
    /// MPointIndex when only 3d points exist.
    pub fn value_3d(&self, index: usize, tab_pnt: &mut [DVec3]) {
        let a_p = self.my_line.as_ref().unwrap().point(index).value();
        // OCCT TabPnt(1).SetCoord(aP.X() + Xo, aP.Y() + Yo, aP.Z() + Zo);
        tab_pnt[0] = DVec3::new(a_p.x + self.xo, a_p.y + self.yo, a_p.z + self.zo);
    }

    /// OCCT BRepApprox_TheMultiLineOfApprox::Value(MPointIndex, tabPt2d)
    /// (hxx L99; gxx L186-210) — the 2d points of the multipoint
    /// MPointIndex when only 2d points exist.
    pub fn value_2d(&self, index: usize, tab_pnt2d: &mut [DVec2]) {
        let p_on_2s = self.my_line.as_ref().unwrap().point(index);
        let (u1, v1, u2, v2) = p_on_2s.parameters();
        if self.nbp2d == 1 {
            if self.p2donfirst {
                // OCCT TabPnt2d(1).SetCoord(u1 + U1o, v1 + V1o);
                tab_pnt2d[0] = DVec2::new(u1 + self.u1o, v1 + self.v1o);
            } else {
                // OCCT TabPnt2d(1).SetCoord(u2 + U2o, v2 + V2o);
                tab_pnt2d[0] = DVec2::new(u2 + self.u2o, v2 + self.v2o);
            }
        } else {
            // OCCT TabPnt2d(1).SetCoord(u1 + U1o, v1 + V1o);
            tab_pnt2d[0] = DVec2::new(u1 + self.u1o, v1 + self.v1o);
            if tab_pnt2d.len() >= 2 {
                // OCCT TabPnt2d(2).SetCoord(u2 + U2o, v2 + V2o);
                tab_pnt2d[1] = DVec2::new(u2 + self.u2o, v2 + self.v2o);
            }
        }
    }

    /// OCCT BRepApprox_TheMultiLineOfApprox::Value(MPointIndex, tabPt,
    /// tabPt2d) (hxx L102-104; gxx L214-220) — the 3d and 2d points of the
    /// multipoint MPointIndex.
    pub fn value_3d_2d(&self, index: usize, tab_pnt: &mut [DVec3], tab_pnt2d: &mut [DVec2]) {
        self.value_3d(index, tab_pnt);
        self.value_2d(index, tab_pnt2d);
    }

    /// OCCT BRepApprox_TheMultiLineOfApprox::Tangency(MPointIndex, tabV)
    /// (hxx L107; gxx L224-240) — the 3d tangency points of the multipoint
    /// MPointIndex only when 3d points exist.
    pub fn tangency_3d(&self, index: usize, tab_v: &mut [DVec3]) -> bool {
        // OCCT if (PtrOnmySvSurfaces == NULL) return false;
        let sv = match &self.ptr_on_my_sv_surfaces {
            Some(sv) => sv,
            None => return false,
        };

        let p_on_2s = self.my_line.as_ref().unwrap().point(index);
        let (u1, v1, u2, v2) = p_on_2s.parameters();

        // OCCT L233: bool ret = ((TheSvSurfaces*)PtrOnmySvSurfaces)
        //   ->Tangency(u1, v1, u2, v2, TabVec(1));
        let ret = sv.tangency(u1, v1, u2, v2, &mut tab_v[0]);
        if !ret {
            // OCCT TabVec(1).SetCoord(0.0, 0.0, 0.0);
            tab_v[0] = DVec3::new(0.0, 0.0, 0.0);
        }

        ret
    }

    /// OCCT BRepApprox_TheMultiLineOfApprox::Tangency(MPointIndex, tabV2d)
    /// (hxx L110; gxx L244-289) — the 2d tangency points of the multipoint
    /// MPointIndex only when 2d points exist.
    // The OCCT statement order is preserved verbatim (ret is assigned in
    // every branch before its first read).
    #[allow(unused_assignments)]
    pub fn tangency_2d(&self, index: usize, tab_v2d: &mut [DVec2]) -> bool {
        // OCCT if (PtrOnmySvSurfaces == NULL) return false;
        let sv = match &self.ptr_on_my_sv_surfaces {
            Some(sv) => sv,
            None => return false,
        };

        let p_on_2s = self.my_line.as_ref().unwrap().point(index);
        let (u1, v1, u2, v2) = p_on_2s.parameters();

        let mut ret = false;
        if self.nbp2d == 1 {
            if self.p2donfirst {
                // OCCT L258: TangencyOnSurf1(u1, v1, u2, v2, TabVec2d(1));
                ret = sv.tangency_on_surf_1(u1, v1, u2, v2, &mut tab_v2d[0]);
            } else {
                // OCCT L262: TangencyOnSurf2(u1, v1, u2, v2, TabVec2d(1));
                ret = sv.tangency_on_surf_2(u1, v1, u2, v2, &mut tab_v2d[0]);
            }
        } else {
            // OCCT L267: TangencyOnSurf1(u1, v1, u2, v2, TabVec2d(1));
            ret = sv.tangency_on_surf_1(u1, v1, u2, v2, &mut tab_v2d[0]);
            if ret {
                if tab_v2d.len() >= 2 {
                    // OCCT L272-274: ret = (ret && TangencyOnSurf2(...));
                    ret = ret && sv.tangency_on_surf_2(u1, v1, u2, v2, &mut tab_v2d[1]);
                }
            }
        }

        if !ret {
            // OCCT L281: TabVec2d(1) = gp_Vec2d(0.0, 0.0);
            tab_v2d[0] = DVec2::new(0.0, 0.0);
            if tab_v2d.len() >= 2 {
                // OCCT L284: TabVec2d(2) = gp_Vec2d(0.0, 0.0);
                tab_v2d[1] = DVec2::new(0.0, 0.0);
            }
        }

        ret
    }

    /// OCCT BRepApprox_TheMultiLineOfApprox::Tangency(MPointIndex, tabV,
    /// tabV2d) (hxx L113-115; gxx L293-298) — the 3d and 2d points of the
    /// multipoint MPointIndex.
    pub fn tangency_3d_2d(&self, index: usize, tab_v: &mut [DVec3], tab_v2d: &mut [DVec2]) -> bool {
        self.tangency_3d(index, tab_v) && self.tangency_2d(index, tab_v2d)
    }

    /// OCCT BRepApprox_TheMultiLineOfApprox::MakeMLBetween(Low, High,
    /// NbPointsToInsert) (hxx L119-121; gxx L302-627) — tries to make a
    /// sub-line between Low and High points of this line by adding
    /// NbPointsToInsert new points.
    // The OCCT for-loop counters and the duplicated CodeErreur assignment
    // are preserved verbatim (gxx L417, L543-544).
    #[allow(unused_variables, unused_assignments)]
    pub fn make_ml_between(
        &self,
        low: usize,
        high: usize,
        a_nb_pnts_to_insert: usize,
    ) -> TheMultiLineOfApprox {
        // OCCT L306-328: without a surfaces functor, return an empty line
        // over [1, 1].
        let sv = match &self.ptr_on_my_sv_surfaces {
            Some(sv) => sv,
            None => {
                let vide1 = Arc::new(LineOn2S::new());
                let vide = Arc::new(ApproxLine::new_line_on_2s(Some(vide1), false));
                return TheMultiLineOfApprox::new_sv_surfaces(
                    vide,
                    None,
                    self.nbp3d,
                    self.nbp2d,
                    self.my_approx_u1v1,
                    self.my_approx_u2v2,
                    self.xo,
                    self.yo,
                    self.zo,
                    self.u1o,
                    self.v1o,
                    self.u2o,
                    self.v2o,
                    self.p2donfirst,
                    1,
                    1,
                );
            }
        };
        let my_line = self.my_line.as_ref().unwrap();

        // OCCT L330-334.
        let a_save_use_solver = sv.get_use_solver();
        if !a_save_use_solver {
            sv.set_use_solver(true);
        }
        let mut nb_pnts_to_insert = a_nb_pnts_to_insert;
        if nb_pnts_to_insert < (high - low) {
            nb_pnts_to_insert = high - low;
        }
        let mut nb_pnts = nb_pnts_to_insert + high - low + 1;
        let mut nb_pntsmin = high - low;
        nb_pntsmin += nb_pntsmin;

        if nb_pnts < nb_pntsmin {
            nb_pnts = nb_pntsmin;
        }

        // OCCT L345-347: gp_Vec T; gp_Vec2d TS1, TS2; gp_Pnt P;

        // OCCT L358-362: NCollection_Array1 U1/V1/U2/V2/AC(Low, High); the
        // rcad Vecs are indexed by `i - low`.
        let mut u1a = vec![0.0f64; high - low + 1];
        let mut v1a = vec![0.0f64; high - low + 1];
        let mut u2a = vec![0.0f64; high - low + 1];
        let mut v2a = vec![0.0f64; high - low + 1];
        let mut ac = vec![0.0f64; high - low + 1];
        let mut s;
        let ds;

        // OCCT L370-375.
        let (mut u1, mut v1, mut u2, mut v2) = my_line.point(low).parameters();
        u1a[0] = u1;
        v1a[0] = v1;
        u2a[0] = u2;
        v2a[0] = v2;
        ac[0] = 0.0;

        // OCCT L392-400 (#else branch: parameterization along the XYZ arc
        // length).
        for i in (low + 1)..=high {
            let p = my_line.point(i);
            (u1, v1, u2, v2) = p.parameters();
            u1a[i - low] = u1;
            v1a[i - low] = v1;
            u2a[i - low] = u2;
            v2a[i - low] = v2;
            ac[i - low] = ac[i - low - 1]
                + my_line
                    .point(i - 1)
                    .value()
                    .distance(my_line.point(i).value());
        }

        // OCCT L404-409: the result structures.
        let mut result_pnt_on_2s_line = LineOn2S::new();

        let mut start_p_on_2s = PntOn2S::new();
        let _start_params = [0.0f64; 4]; // OCCT NCollection_Array1 StartParams(1, 4)

        ds = ac[high - low] / (nb_pnts as f64 - 1.0);
        let mut indice = low;
        let mut has_been_inserted = false;
        let dsmin = ds * 0.3;
        let smax = ac[high - low];

        // OCCT L417: for (i = 2, s = ds; (s < smax && Indice <= High - 1);
        //              i++, s += ds)
        let mut i = 2usize;
        s = ds;
        while s < smax && indice <= high - 1 {
            // OCCT L425-434: the current Indice with AC(Indice) <= s <
            // AC(Indice + 1).
            while ac[indice + 1 - low] <= s {
                if !has_been_inserted {
                    result_pnt_on_2s_line.add(&my_line.point(indice));
                }

                has_been_inserted = false;
                indice += 1;
                if indice == high {
                    break;
                }
            }

            if indice == high {
                break;
            }

            if !has_been_inserted && ac[indice - low] <= s {
                result_pnt_on_2s_line.add(&my_line.point(indice));
                has_been_inserted = true;
            }

            let a = s - ac[indice - low];
            let b = ac[indice + 1 - low] - s;
            let nab = 1.0 / (a + b);

            // OCCT L449-472: if the new point is far enough from the two
            // existing points, compute and insert it.
            if (a > dsmin) && (b > dsmin) {
                u1 = (u1a[indice - low] * b + u1a[indice + 1 - low] * a) * nab;
                v1 = (v1a[indice - low] * b + v1a[indice + 1 - low] * a) * nab;
                u2 = (u2a[indice - low] * b + u2a[indice + 1 - low] * a) * nab;
                v2 = (v2a[indice - low] * b + v2a[indice + 1 - low] * a) * nab;

                let mut p = DVec3::ZERO;
                let mut t = DVec3::ZERO;
                let mut ts1 = DVec2::ZERO;
                let mut ts2 = DVec2::ZERO;
                if sv.compute(
                    &mut u1,
                    &mut v1,
                    &mut u2,
                    &mut v2,
                    &mut p,
                    &mut t,
                    &mut ts1,
                    &mut ts2,
                ) {
                    // OCCT L463: StartPOn2S.SetValue(P, u1, v1, u2, v2);
                    start_p_on_2s.set_value_all(p, u1, v1, u2, v2);
                    result_pnt_on_2s_line.add(&start_p_on_2s);
                }
            } else if b < 0.0 {
                // OCCT L478-501: the point is not far enough from the two
                // existing points; advance the index.
                while ac[indice + 1 - low] <= s {
                    if !has_been_inserted {
                        result_pnt_on_2s_line.add(&my_line.point(indice));
                    }

                    has_been_inserted = false;
                    indice += 1;

                    if indice == high {
                        break;
                    }
                }

                if indice == high {
                    break;
                }

                if !has_been_inserted && ac[indice - low] <= s {
                    result_pnt_on_2s_line.add(&my_line.point(indice));
                    has_been_inserted = true;
                }
            } else {
                // OCCT L504: s += (dsmin - ds);
                s += dsmin - ds;
            }

            i += 1;
            s += ds;
        }

        // OCCT L509: the final point.
        result_pnt_on_2s_line.add(&my_line.point(high));
        let result_pnt_on_2s_line = Arc::new(result_pnt_on_2s_line);

        // OCCT L510: the candidate line.
        let temp = Arc::new(ApproxLine::new_line_on_2s(
            Some(result_pnt_on_2s_line.clone()),
            false,
        ));

        // OCCT L512-524: a posteriori check — no too-sharp turn in 2d.
        let (u1f, v1f, u2f, v2f) = temp.point(1).parameters();
        let mut p1a = DVec2::new(u1f, v1f);
        let mut p2a = DVec2::new(u2f, v2f);

        let (u1s, v1s, u2s, v2s) = temp.point(2).parameters();
        let mut p1b = DVec2::new(u1s, v1s);
        let mut p2b = DVec2::new(u2s, v2s);

        let mut p1c;
        let mut p2c;

        let mut code_erreur = 0i32;

        // OCCT L527: for (i = 3, NbPnts = temp->NbPnts();
        //              CodeErreur == 0 && i <= NbPnts; i++)
        let nb_pnts2 = temp.nb_pnts();
        let mut i = 3usize;
        while code_erreur == 0 && i <= nb_pnts2 {
            let (mut u1, mut v1, mut u2, mut v2) = temp.point(i).parameters();
            // Turn P1A P1B P1C (OCCT L531).
            p1c = DVec2::new(u1, v1);
            let mut du = p1b.x - p1a.x;
            let mut dv = p1b.y - p1a.y;
            let mut duv2 = 0.25 * (du * du + dv * dv);
            u1 = p1b.x + du;
            v1 = p1b.y + dv;
            du = p1c.x - u1;
            dv = p1c.y - v1;
            let mut d = du * du + dv * dv;
            if d > duv2 {
                code_erreur = 1;
                code_erreur = 1; // OCCT L543-544: the duplicated assignment
                break;
            }

            // Turn P2A P2B P2C (OCCT L548).
            p2c = DVec2::new(u2, v2);
            du = p2b.x - p2a.x;
            dv = p2b.y - p2a.y;
            duv2 = 0.25 * (du * du + dv * dv);
            u2 = p2b.x + du;
            v2 = p2b.y + dv;
            du = p2c.x - u2;
            dv = p2c.y - v2;
            d = du * du + dv * dv;
            if d > duv2 {
                code_erreur = 2;
                break;
            }

            p1a = p1b;
            p2a = p2b;
            p1b = p1c;
            p2b = p2c;

            i += 1;
        }

        // OCCT L583-602: accept the new line only when enough points were
        // produced and no sharp turn was detected.
        if temp.nb_pnts() >= nb_pnts_to_insert + high - low + 1 && code_erreur == 0 {
            sv.set_use_solver(a_save_use_solver);
            // OCCT L587: (High - Low > 10) ? PtrOnmySvSurfaces : NULL
            let sv_ptr = if high - low > 10 {
                self.ptr_on_my_sv_surfaces.clone()
            } else {
                None
            };
            TheMultiLineOfApprox::new_sv_surfaces(
                temp,
                sv_ptr,
                self.nbp3d,
                self.nbp2d,
                self.my_approx_u1v1,
                self.my_approx_u2v2,
                self.xo,
                self.yo,
                self.zo,
                self.u1o,
                self.v1o,
                self.u2o,
                self.v2o,
                self.p2donfirst,
                1,
                result_pnt_on_2s_line.nb_points(),
            )
        } else {
            // OCCT L607-626: not enough points — restore the solver flag and
            // return the empty line.
            sv.set_use_solver(a_save_use_solver);
            let vide1 = Arc::new(LineOn2S::new());
            let vide = Arc::new(ApproxLine::new_line_on_2s(Some(vide1), false));
            TheMultiLineOfApprox::new_sv_surfaces(
                vide,
                None,
                self.nbp3d,
                self.nbp2d,
                self.my_approx_u1v1,
                self.my_approx_u2v2,
                self.xo,
                self.yo,
                self.zo,
                self.u1o,
                self.v1o,
                self.u2o,
                self.v2o,
                self.p2donfirst,
                1,
                1,
            )
        }
    }

    /// OCCT BRepApprox_TheMultiLineOfApprox::MakeMLOneMorePoint(theLow,
    /// theHigh, theIndbad, theNewMultiLine) (hxx L125-128; gxx L631-753) —
    /// tries to make a sub-line by adding one more point between
    /// (indbad-1)-th and indbad-th points.
    // The OCCT statement order is preserved verbatim (gxx L667-668).
    #[allow(unused_assignments)]
    pub fn make_ml_one_more_point(
        &self,
        the_low: usize,
        the_high: usize,
        the_indbad: usize,
        the_new_multi_line: &mut TheMultiLineOfApprox,
    ) -> bool {
        let mut other_line_made = false;

        // OCCT L637-638.
        let sv = match &self.ptr_on_my_sv_surfaces {
            Some(sv) => sv,
            None => return false,
        };
        let my_line = self.my_line.as_ref().unwrap();

        // OCCT L640-644.
        let a_save_use_solver = sv.get_use_solver();
        if !a_save_use_solver {
            sv.set_use_solver(true);
        }

        // OCCT L646-648.
        let sq_tol3d = SQUARE_CONFUSION;
        let mut tolerance = Vector::new(1, 2);
        tolerance.set(1, 1.0e-8);
        tolerance.set(2, 1.0e-8);

        // OCCT L650-652.
        let mut result_pnt_on_2s_line = LineOn2S::new();
        for indice in the_low..=the_high {
            result_pnt_on_2s_line.add(&my_line.point(indice));
        }

        // Insert new point between (theIndbad-1) and theIndbad
        // Using the SvSurfaces for Rsnld: it may be ImpPrm or PrmPrm
        // (OCCT L654-665).
        let prev_pnt = my_line.point(the_indbad - 1).value();
        let cur_pnt = my_line.point(the_indbad).value();
        let (uprev1, vprev1, uprev2, vprev2) = my_line.point(the_indbad - 1).parameters();
        let (ucur1, vcur1, ucur2, vcur2) = my_line.point(the_indbad).parameters();
        let umid1 = (uprev1 + ucur1) / 2.0;
        let vmid1 = (vprev1 + vcur1) / 2.0;
        let umid2 = (uprev2 + ucur2) / 2.0;
        let vmid2 = (vprev2 + vcur2) / 2.0;
        let mut mid_point = PntOn2S::new();
        let mut is_new_point_invalid = false;
        // OCCT L668-669.
        is_new_point_invalid = self.my_approx_u1v1
            && (ucur1 - umid1).abs() <= tolerance.get(1)
            && (vcur1 - vmid1).abs() <= tolerance.get(2);
        if !is_new_point_invalid {
            // OCCT L672-673.
            is_new_point_invalid = self.my_approx_u2v2
                && (ucur2 - umid2).abs() <= tolerance.get(1)
                && (vcur2 - vmid2).abs() <= tolerance.get(2);
            // OCCT L674-675.
            if !is_new_point_invalid && sv.seek_point(umid1, vmid1, umid2, vmid2, &mut mid_point) {
                // OCCT L677-680.
                let new_pnt = mid_point.value();
                let sq_dist_new_prev = new_pnt.distance_squared(prev_pnt);
                let sq_dist_new_cur = new_pnt.distance_squared(cur_pnt);
                is_new_point_invalid = sq_dist_new_prev <= sq_tol3d || sq_dist_new_cur <= sq_tol3d;
                if !is_new_point_invalid {
                    // OCCT L683-708.
                    let (unew1, vnew1, unew2, vnew2) = mid_point.parameters();
                    if self.my_approx_u1v1 {
                        let sq_dist_cur_mid1 = (ucur1 - umid1) * (ucur1 - umid1)
                            + (vcur1 - vmid1) * (vcur1 - vmid1);
                        let sq_dist_mid_new1 = (umid1 - unew1) * (umid1 - unew1)
                            + (vmid1 - vnew1) * (vmid1 - vnew1);
                        is_new_point_invalid = sq_dist_mid_new1 > sq_dist_cur_mid1;
                    }
                    if !is_new_point_invalid {
                        if self.my_approx_u2v2 {
                            let sq_dist_cur_mid2 = (ucur2 - umid2) * (ucur2 - umid2)
                                + (vcur2 - vmid2) * (vcur2 - vmid2);
                            let sq_dist_mid_new2 = (umid2 - unew2) * (umid2 - unew2)
                                + (vmid2 - vnew2) * (vmid2 - vnew2);
                            is_new_point_invalid = sq_dist_mid_new2 > sq_dist_cur_mid2;
                        }
                        if !is_new_point_invalid {
                            // OCCT L705: InsertBefore(theIndbad - theLow + 1,
                            // MidPoint); the rcad LineOn2S::insert_before is
                            // 0-based.
                            result_pnt_on_2s_line
                                .insert_before(the_indbad - the_low, &mid_point);
                            other_line_made = true;
                        }
                    }
                }
            }
        }

        if !other_line_made {
            // OCCT L713-717.
            sv.set_use_solver(a_save_use_solver);
            return false;
        }

        // OCCT L719-733: the #ifdef DRAW debug dump is not translated.

        // OCCT L734-736.
        let result_pnt_on_2s_line = Arc::new(result_pnt_on_2s_line);
        let temp = Arc::new(ApproxLine::new_line_on_2s(
            Some(result_pnt_on_2s_line.clone()),
            false,
        ));
        sv.set_use_solver(a_save_use_solver);
        *the_new_multi_line = TheMultiLineOfApprox::new_sv_surfaces(
            temp,
            self.ptr_on_my_sv_surfaces.clone(),
            self.nbp3d,
            self.nbp2d,
            self.my_approx_u1v1,
            self.my_approx_u2v2,
            self.xo,
            self.yo,
            self.zo,
            self.u1o,
            self.v1o,
            self.u2o,
            self.v2o,
            self.p2donfirst,
            1,
            result_pnt_on_2s_line.nb_points(),
        );
        true
    }

    /// OCCT BRepApprox_TheMultiLineOfApprox::Dump (hxx L131; gxx L757-778) —
    /// dump of the current multi-line.
    pub fn dump(&self) {
        // OCCT L759-760: NCollection_Array1 anArr1(1, 1); anArr2(1, 2);
        let mut an_arr1 = vec![DVec3::ZERO; 1];
        let mut an_arr2 = vec![DVec2::ZERO; 2];

        let an_ind_f = self.first_point();
        let an_ind_l = self.last_point();

        for ind in an_ind_f..=an_ind_l {
            self.value_3d_2d(ind, &mut an_arr1, &mut an_arr2);
            println!(
                "{:4}  [{:+10.20} {:+10.20} {:+10.20}]  [{:+10.20} {:+10.20}]  [{:+10.20} {:+10.20}]",
                ind,
                an_arr1[0].x,
                an_arr1[0].y,
                an_arr1[0].z,
                an_arr2[0].x,
                an_arr2[0].y,
                an_arr2[1].x,
                an_arr2[1].y
            );
        }
    }
}

impl Default for TheMultiLineOfApprox {
    fn default() -> Self {
        TheMultiLineOfApprox::new()
    }
}

// ---------------------------------------------------------------------------
// BRepApprox_TheMultiLineToolOfApprox
// ---------------------------------------------------------------------------

/// OCCT BRepApprox_TheMultiLineToolOfApprox
/// (BRepApprox_TheMultiLineToolOfApprox.hxx L30-123) — the static LineTool
/// the approximation algorithms call; every body is the one-liner of
/// ApproxInt_MultiLineTool.lxx forwarded to the MultiLine.
pub mod the_multi_line_tool {
    use super::{ApproxStatus, DVec2, DVec3, TheMultiLineOfApprox};

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::FirstPoint (hxx L36; lxx
    /// L43-46).
    pub fn first_point(ml: &TheMultiLineOfApprox) -> usize {
        ml.first_point()
    }

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::LastPoint (hxx L39; lxx
    /// L50-53).
    pub fn last_point(ml: &TheMultiLineOfApprox) -> usize {
        ml.last_point()
    }

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::NbP2d (hxx L42; lxx L29-32).
    pub fn nb_p2d(ml: &TheMultiLineOfApprox) -> usize {
        ml.nb_p2d()
    }

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::NbP3d (hxx L45; lxx L36-39).
    pub fn nb_p3d(ml: &TheMultiLineOfApprox) -> usize {
        ml.nb_p3d()
    }

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::Value(ML, MPointIndex,
    /// tabPt) (hxx L49-51; lxx L57-62) — the 3d points of the multipoint
    /// MPointIndex when only 3d points exist.
    pub fn value_3d(ml: &TheMultiLineOfApprox, index: usize, tab_pnt: &mut [DVec3]) {
        ml.value_3d(index, tab_pnt)
    }

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::Value(ML, MPointIndex,
    /// tabPt2d) (hxx L55-57; lxx L66-71) — the 2d points of the multipoint
    /// MPointIndex when only 2d points exist.
    pub fn value_2d(ml: &TheMultiLineOfApprox, index: usize, tab_pnt2d: &mut [DVec2]) {
        ml.value_2d(index, tab_pnt2d)
    }

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::Value(ML, MPointIndex,
    /// tabPt, tabPt2d) (hxx L61-64; lxx L75-81) — the 3d and 2d points of
    /// the multipoint MPointIndex.
    pub fn value_3d_2d(
        ml: &TheMultiLineOfApprox,
        index: usize,
        tab_pnt: &mut [DVec3],
        tab_pnt2d: &mut [DVec2],
    ) {
        ml.value_3d_2d(index, tab_pnt, tab_pnt2d)
    }

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::Tangency(ML, MPointIndex,
    /// tabV) (hxx L68-70; lxx L85-90) — the 3d tangency points of the
    /// multipoint MPointIndex when only 3d points exist.
    pub fn tangency_3d(ml: &TheMultiLineOfApprox, index: usize, tab_v: &mut [DVec3]) -> bool {
        ml.tangency_3d(index, tab_v)
    }

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::Tangency(ML, MPointIndex,
    /// tabV2d) (hxx L74-76; lxx L94-99) — the 2d tangency points of the
    /// multipoint MPointIndex only when 2d points exist.
    pub fn tangency_2d(ml: &TheMultiLineOfApprox, index: usize, tab_v2d: &mut [DVec2]) -> bool {
        ml.tangency_2d(index, tab_v2d)
    }

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::Tangency(ML, MPointIndex,
    /// tabV, tabV2d) (hxx L80-83; lxx L103-109) — the 3d and 2d tangency
    /// points of the multipoint MPointIndex.
    pub fn tangency_3d_2d(
        ml: &TheMultiLineOfApprox,
        index: usize,
        tab_v: &mut [DVec3],
        tab_v2d: &mut [DVec2],
    ) -> bool {
        ml.tangency_3d_2d(index, tab_v, tab_v2d)
    }

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::Curvature(ML, MPointIndex,
    /// tabV) (hxx L87-89; lxx L113-120) — the instantiation returns False
    /// unconditionally.
    pub fn curvature_3d(_ml: &TheMultiLineOfApprox, _index: usize, _tab_v: &mut [DVec3]) -> bool {
        false
    }

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::Curvature(ML, MPointIndex,
    /// tabV2d) (hxx L93-95; lxx L124-131) — the instantiation returns False
    /// unconditionally.
    pub fn curvature_2d(
        _ml: &TheMultiLineOfApprox,
        _index: usize,
        _tab_v2d: &mut [DVec2],
    ) -> bool {
        false
    }

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::Curvature(ML, MPointIndex,
    /// tabV, tabV2d) (hxx L99-102; lxx L135-144) — the instantiation returns
    /// False unconditionally.
    pub fn curvature_3d_2d(
        _ml: &TheMultiLineOfApprox,
        _index: usize,
        _tab_v: &mut [DVec3],
        _tab_v2d: &mut [DVec2],
    ) -> bool {
        false
    }

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::MakeMLBetween(ML, I1, I2,
    /// NbPMin) (hxx L105-108; lxx L161-168) — is called if WhatStatus
    /// returned "PointsAdded".
    pub fn make_ml_between(
        ml: &TheMultiLineOfApprox,
        i1: usize,
        i2: usize,
        nb_p_min: usize,
    ) -> TheMultiLineOfApprox {
        ml.make_ml_between(i1, i2, nb_p_min)
    }

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::MakeMLOneMorePoint(ML, I1,
    /// I2, indbad, OtherLine) (hxx L111-115; lxx L172-179) — is called when
    /// the Bezier curve contains a loop.
    pub fn make_ml_one_more_point(
        ml: &TheMultiLineOfApprox,
        i1: usize,
        i2: usize,
        indbad: usize,
        other_line: &mut TheMultiLineOfApprox,
    ) -> bool {
        ml.make_ml_one_more_point(i1, i2, indbad, other_line)
    }

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::WhatStatus(ML, I1, I2)
    /// (hxx L117-119; lxx L148-157).
    pub fn what_status(ml: &TheMultiLineOfApprox, _i1: usize, _i2: usize) -> ApproxStatus {
        ml.what_status()
    }

    /// OCCT BRepApprox_TheMultiLineToolOfApprox::Dump(ML) (hxx L122; lxx
    /// L181-184) — dump of the current multi-line.
    pub fn dump(ml: &TheMultiLineOfApprox) {
        ml.dump()
    }
}

// ---------------------------------------------------------------------------
// BRepApprox_SurfaceTool
// ---------------------------------------------------------------------------

/// OCCT BRepApprox_SurfaceTool (BRepApprox_SurfaceTool.hxx L37-164) — the
/// static surface tool over the BRepAdaptor_Surface.  The OCCT one-liners of
/// BRepApprox_SurfaceTool.lxx forward to the adaptor; rcad forwards to the
/// landed BRepAdaptorSurface (topalgo/brep_adaptor/surface.rs).  Methods
/// whose adaptor machinery is not landed yet keep their OCCT signature and
/// are documented unimplemented.
pub mod surface_tool {
    use std::sync::Arc;

    use rcad_kernel::geom::{BSplineSurface, BezierSurface, Point3, Vec3};
    use rcad_kernel::math::GeomAbsShape;

    use crate::geomalgo::int_curve_surface::Adaptor3dCurveBasis;
    use crate::geomalgo::int_patch::GeomAbsSurfaceType;
    use crate::topalgo::brep_adaptor::surface::BRepAdaptorSurface;

    /// OCCT BRepApprox_SurfaceTool::FirstUParameter (hxx L42; lxx L31-34).
    pub fn first_u_parameter(s: &BRepAdaptorSurface) -> f64 {
        s.first_u_parameter()
    }

    /// OCCT BRepApprox_SurfaceTool::FirstVParameter (hxx L44; lxx L36-39).
    pub fn first_v_parameter(s: &BRepAdaptorSurface) -> f64 {
        s.first_v_parameter()
    }

    /// OCCT BRepApprox_SurfaceTool::LastUParameter (hxx L46; lxx L41-44).
    pub fn last_u_parameter(s: &BRepAdaptorSurface) -> f64 {
        s.last_u_parameter()
    }

    /// OCCT BRepApprox_SurfaceTool::LastVParameter (hxx L48; lxx L46-49).
    pub fn last_v_parameter(s: &BRepAdaptorSurface) -> f64 {
        s.last_v_parameter()
    }

    /// OCCT BRepApprox_SurfaceTool::NbUIntervals (hxx L50; lxx L51-55) —
    /// unimplemented: the continuity-interval machinery of
    /// BRepAdaptor_Surface is not landed yet.
    pub fn nb_u_intervals(_s: &BRepAdaptorSurface, _sh: GeomAbsShape) -> usize {
        unimplemented!(
            "OCCT BRepApprox_SurfaceTool::NbUIntervals (SurfaceTool.lxx L51-55) needs BRepAdaptor_Surface::NbUIntervals (not yet landed)"
        )
    }

    /// OCCT BRepApprox_SurfaceTool::NbVIntervals (hxx L52; lxx L57-61) —
    /// unimplemented: the continuity-interval machinery of
    /// BRepAdaptor_Surface is not landed yet.
    pub fn nb_v_intervals(_s: &BRepAdaptorSurface, _sh: GeomAbsShape) -> usize {
        unimplemented!(
            "OCCT BRepApprox_SurfaceTool::NbVIntervals (SurfaceTool.lxx L57-61) needs BRepAdaptor_Surface::NbVIntervals (not yet landed)"
        )
    }

    /// OCCT BRepApprox_SurfaceTool::UIntervals (hxx L54-56; lxx L63-68) —
    /// unimplemented: the continuity-interval machinery of
    /// BRepAdaptor_Surface is not landed yet.
    pub fn u_intervals(_s: &BRepAdaptorSurface, _tab: &mut [f64], _sh: GeomAbsShape) {
        unimplemented!(
            "OCCT BRepApprox_SurfaceTool::UIntervals (SurfaceTool.lxx L63-68) needs BRepAdaptor_Surface::UIntervals (not yet landed)"
        )
    }

    /// OCCT BRepApprox_SurfaceTool::VIntervals (hxx L58-60; lxx L70-75) —
    /// unimplemented: the continuity-interval machinery of
    /// BRepAdaptor_Surface is not landed yet.
    pub fn v_intervals(_s: &BRepAdaptorSurface, _tab: &mut [f64], _sh: GeomAbsShape) {
        unimplemented!(
            "OCCT BRepApprox_SurfaceTool::VIntervals (SurfaceTool.lxx L70-75) needs BRepAdaptor_Surface::VIntervals (not yet landed)"
        )
    }

    /// OCCT BRepApprox_SurfaceTool::UTrim (hxx L63-66; lxx L77-83) — if
    /// First >= Last; returns the trimmed adaptor.
    pub fn u_trim<'a>(
        s: &'a BRepAdaptorSurface<'a>,
        first: f64,
        last: f64,
        tol: f64,
    ) -> BRepAdaptorSurface<'a> {
        s.u_trim(first, last, tol)
    }

    /// OCCT BRepApprox_SurfaceTool::VTrim (hxx L69-72; lxx L85-91) — if
    /// First >= Last; returns the trimmed adaptor.
    pub fn v_trim<'a>(
        s: &'a BRepAdaptorSurface<'a>,
        first: f64,
        last: f64,
        tol: f64,
    ) -> BRepAdaptorSurface<'a> {
        s.v_trim(first, last, tol)
    }

    /// OCCT BRepApprox_SurfaceTool::IsUClosed (hxx L74; lxx L93-96).
    pub fn is_u_closed(s: &BRepAdaptorSurface) -> bool {
        s.is_u_closed()
    }

    /// OCCT BRepApprox_SurfaceTool::IsVClosed (hxx L76; lxx L98-101).
    pub fn is_v_closed(s: &BRepAdaptorSurface) -> bool {
        s.is_v_closed()
    }

    /// OCCT BRepApprox_SurfaceTool::IsUPeriodic (hxx L78; lxx L103-106).
    pub fn is_u_periodic(s: &BRepAdaptorSurface) -> bool {
        s.is_u_periodic()
    }

    /// OCCT BRepApprox_SurfaceTool::UPeriod (hxx L80; lxx L108-111).
    pub fn u_period(s: &BRepAdaptorSurface) -> f64 {
        s.u_period()
    }

    /// OCCT BRepApprox_SurfaceTool::IsVPeriodic (hxx L82; lxx L113-116).
    pub fn is_v_periodic(s: &BRepAdaptorSurface) -> bool {
        s.is_v_periodic()
    }

    /// OCCT BRepApprox_SurfaceTool::VPeriod (hxx L84; lxx L118-121).
    pub fn v_period(s: &BRepAdaptorSurface) -> f64 {
        s.v_period()
    }

    /// OCCT BRepApprox_SurfaceTool::Value (hxx L86; lxx L123-128).
    pub fn value(s: &BRepAdaptorSurface, u: f64, v: f64) -> Point3 {
        s.value(u, v)
    }

    /// OCCT BRepApprox_SurfaceTool::D0 (hxx L88; lxx L130-136) — the surface
    /// point (the OCCT D0(S, U, V, P) out-parameter is the return value).
    pub fn d0(s: &BRepAdaptorSurface, u: f64, v: f64) -> Point3 {
        s.value(u, v)
    }

    /// OCCT BRepApprox_SurfaceTool::D1 (hxx L90-95; lxx L138-146) — (P, D1U,
    /// D1V).
    pub fn d1(s: &BRepAdaptorSurface, u: f64, v: f64) -> (Point3, Vec3, Vec3) {
        s.d1(u, v)
    }

    /// OCCT BRepApprox_SurfaceTool::D2 (hxx L97-105; lxx L148-159) — (P,
    /// D1U, D1V, D2U, D2V, D2UV).
    #[allow(clippy::type_complexity)]
    pub fn d2(s: &BRepAdaptorSurface, u: f64, v: f64) -> (Point3, Vec3, Vec3, Vec3, Vec3, Vec3) {
        s.d2(u, v)
    }

    /// OCCT BRepApprox_SurfaceTool::D3 (hxx L107-119; lxx L161-176) —
    /// unimplemented: BRepAdaptor_Surface has no third-order evaluation
    /// landed yet.
    #[allow(clippy::too_many_arguments)]
    pub fn d3(
        _s: &BRepAdaptorSurface,
        _u: f64,
        _v: f64,
        _p: &mut Point3,
        _d1u: &mut Vec3,
        _d1v: &mut Vec3,
        _d2u: &mut Vec3,
        _d2v: &mut Vec3,
        _d2uv: &mut Vec3,
        _d3u: &mut Vec3,
        _d3v: &mut Vec3,
        _d3uuv: &mut Vec3,
        _d3uvv: &mut Vec3,
    ) {
        unimplemented!(
            "OCCT BRepApprox_SurfaceTool::D3 (SurfaceTool.lxx L161-176) needs BRepAdaptor_Surface::D3 (not yet landed)"
        )
    }

    /// OCCT BRepApprox_SurfaceTool::DN (hxx L121-125; lxx L178-185) —
    /// unimplemented: BRepAdaptor_Surface has no DN evaluation landed yet.
    pub fn dn(_s: &BRepAdaptorSurface, _u: f64, _v: f64, _nu: usize, _nv: usize) -> Vec3 {
        unimplemented!(
            "OCCT BRepApprox_SurfaceTool::DN (SurfaceTool.lxx L178-185) needs BRepAdaptor_Surface::DN (not yet landed)"
        )
    }

    /// OCCT BRepApprox_SurfaceTool::UResolution (hxx L127; lxx L187-190).
    pub fn u_resolution(s: &BRepAdaptorSurface, r3d: f64) -> f64 {
        s.u_resolution(r3d)
    }

    /// OCCT BRepApprox_SurfaceTool::VResolution (hxx L129; lxx L192-195).
    pub fn v_resolution(s: &BRepAdaptorSurface, r3d: f64) -> f64 {
        s.v_resolution(r3d)
    }

    /// OCCT BRepApprox_SurfaceTool::GetType (hxx L131; lxx L197-200).
    pub fn get_type(s: &BRepAdaptorSurface) -> GeomAbsSurfaceType {
        s.get_type()
    }

    /// OCCT BRepApprox_SurfaceTool::Plane (hxx L133; lxx L202-205).
    pub fn plane(
        s: &BRepAdaptorSurface,
    ) -> rcad_kernel::geom::Plane {
        s.plane()
    }

    /// OCCT BRepApprox_SurfaceTool::Cylinder (hxx L135; lxx L207-210).
    pub fn cylinder(s: &BRepAdaptorSurface) -> rcad_kernel::geom::CylindricalSurface {
        s.cylinder()
    }

    /// OCCT BRepApprox_SurfaceTool::Cone (hxx L137; lxx L212-215).
    pub fn cone(s: &BRepAdaptorSurface) -> rcad_kernel::geom::ConicalSurface {
        s.cone()
    }

    /// OCCT BRepApprox_SurfaceTool::Torus (hxx L139; lxx L222-225).
    pub fn torus(s: &BRepAdaptorSurface) -> rcad_kernel::geom::ToroidalSurface {
        s.torus()
    }

    /// OCCT BRepApprox_SurfaceTool::Sphere (hxx L141; lxx L217-220).
    pub fn sphere(s: &BRepAdaptorSurface) -> rcad_kernel::geom::SphericalSurface {
        s.sphere()
    }

    /// OCCT BRepApprox_SurfaceTool::Bezier (hxx L143; lxx L227-230) —
    /// unimplemented: the Bezier-surface extraction of BRepAdaptor_Surface
    /// is not landed yet (the OCCT handle maps to Arc).
    pub fn bezier(_s: &BRepAdaptorSurface) -> Arc<BezierSurface> {
        unimplemented!(
            "OCCT BRepApprox_SurfaceTool::Bezier (SurfaceTool.lxx L227-230) needs BRepAdaptor_Surface::Bezier (not yet landed)"
        )
    }

    /// OCCT BRepApprox_SurfaceTool::BSpline (hxx L145; lxx L232-236) —
    /// unimplemented: the BSpline-surface extraction of BRepAdaptor_Surface
    /// is not landed yet (the OCCT handle maps to Arc).
    pub fn bspline(_s: &BRepAdaptorSurface) -> Arc<BSplineSurface> {
        unimplemented!(
            "OCCT BRepApprox_SurfaceTool::BSpline (SurfaceTool.lxx L232-236) needs BRepAdaptor_Surface::BSpline (not yet landed)"
        )
    }

    /// OCCT BRepApprox_SurfaceTool::AxeOfRevolution (hxx L147; lxx L238-241)
    /// — the gp_Ax1 as (Location, Direction); unimplemented: the
    /// revolution-surface adaptor machinery is not landed yet.
    pub fn axe_of_revolution(_s: &BRepAdaptorSurface) -> (Point3, Vec3) {
        unimplemented!(
            "OCCT BRepApprox_SurfaceTool::AxeOfRevolution (SurfaceTool.lxx L238-241) needs BRepAdaptor_Surface::AxeOfRevolution (not yet landed)"
        )
    }

    /// OCCT BRepApprox_SurfaceTool::Direction (hxx L149; lxx L243-246) —
    /// unimplemented: the direction adaptor machinery is not landed yet.
    pub fn direction(_s: &BRepAdaptorSurface) -> Vec3 {
        unimplemented!(
            "OCCT BRepApprox_SurfaceTool::Direction (SurfaceTool.lxx L243-246) needs BRepAdaptor_Surface::Direction (not yet landed)"
        )
    }

    /// OCCT BRepApprox_SurfaceTool::BasisCurve (hxx L151; lxx L248-251) —
    /// the handle<Adaptor3d_Curve> (the landed Adaptor3dCurveBasis trait is
    /// the rcad encoding of that handle); unimplemented: the basis-curve
    /// adaptor machinery (Adaptor3d_CurveOnSurface) is not landed yet.
    pub fn basis_curve(_s: &BRepAdaptorSurface) -> Arc<dyn Adaptor3dCurveBasis> {
        unimplemented!(
            "OCCT BRepApprox_SurfaceTool::BasisCurve (SurfaceTool.lxx L248-251) needs Adaptor3d_CurveOnSurface (not yet landed)"
        )
    }

    /// OCCT BRepApprox_SurfaceTool::NbSamplesU (hxx L153; cxx L22-25).
    pub fn nb_samples_u(_s: &BRepAdaptorSurface) -> usize {
        10
    }

    /// OCCT BRepApprox_SurfaceTool::NbSamplesV (hxx L155; cxx L27-30).
    pub fn nb_samples_v(_s: &BRepAdaptorSurface) -> usize {
        10
    }

    /// OCCT BRepApprox_SurfaceTool::NbSamplesU (hxx L157-159; cxx L32-35).
    pub fn nb_samples_u_window(_s: &BRepAdaptorSurface, _u1: f64, _u2: f64) -> usize {
        10
    }

    /// OCCT BRepApprox_SurfaceTool::NbSamplesV (hxx L161-163; cxx L37-40).
    pub fn nb_samples_v_window(_s: &BRepAdaptorSurface, _v1: f64, _v2: f64) -> usize {
        10
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use std::cell::Cell;

    use super::*;
    use crate::geomalgo::int_patch::GeomAbsSurfaceType;
    use crate::topalgo::brep_adaptor::surface::BRepAdaptorSurface;
    use rcad_kernel::geom::{Line3, Plane, Point3, Surface3, Vec3};
    use rcad_kernel::topo::topods::{BRepBuilder, Shape};

    /// A test ApproxInt_SvSurfaces implementation with fixed outputs.
    struct StubSvSurfaces {
        use_solver: Cell<bool>,
    }

    impl StubSvSurfaces {
        fn new() -> Arc<StubSvSurfaces> {
            Arc::new(StubSvSurfaces {
                use_solver: Cell::new(false),
            })
        }
    }

    impl SvSurfaces for StubSvSurfaces {
        fn compute(
            &self,
            _u1: &mut f64,
            _v1: &mut f64,
            _u2: &mut f64,
            _v2: &mut f64,
            _pt: &mut DVec3,
            _tg: &mut DVec3,
            _tguv1: &mut DVec2,
            _tguv2: &mut DVec2,
        ) -> bool {
            false
        }

        fn pnt(&self, _u1: f64, _v1: f64, _u2: f64, _v2: f64, _p: &mut DVec3) {}

        fn seek_point(&self, _u1: f64, _v1: f64, _u2: f64, _v2: f64, _point: &mut PntOn2S) -> bool {
            false
        }

        fn tangency(&self, _u1: f64, _v1: f64, _u2: f64, _v2: f64, tg: &mut DVec3) -> bool {
            *tg = DVec3::new(1.0, 2.0, 3.0);
            true
        }

        fn tangency_on_surf_1(&self, _u1: f64, _v1: f64, _u2: f64, _v2: f64, tg: &mut DVec2) -> bool {
            *tg = DVec2::new(0.5, 0.25);
            true
        }

        fn tangency_on_surf_2(&self, _u1: f64, _v1: f64, _u2: f64, _v2: f64, tg: &mut DVec2) -> bool {
            *tg = DVec2::new(0.75, 0.125);
            true
        }

        fn set_use_solver(&self, the_use_sol: bool) {
            self.use_solver.set(the_use_sol);
        }

        fn get_use_solver(&self) -> bool {
            self.use_solver.get()
        }
    }

    /// A 3-point IntSurf_LineOn2S; point i (1-based) has value
    /// (t, -t, 2t) and parameters (t, t+1, t+2, t+3) with t = i - 1.
    fn line_on_2s_3pts() -> Arc<LineOn2S> {
        let mut line = LineOn2S::new();
        for i in 0..3 {
            let t = i as f64;
            let mut p = PntOn2S::new();
            p.set_value_all(DVec3::new(t, -t, 2.0 * t), t, t + 1.0, t + 2.0, t + 3.0);
            line.add(&p);
        }
        Arc::new(line)
    }

    /// A valid multi-line over the 3-point line with the offsets
    /// Xo=10 Yo=20 Zo=30, U1o=0.1 V1o=0.2 U2o=0.3 V2o=0.4.
    fn multi_line_3pts(sv: Option<SharedSvSurfaces>) -> TheMultiLineOfApprox {
        let line = Arc::new(ApproxLine::new_line_on_2s(Some(line_on_2s_3pts()), false));
        TheMultiLineOfApprox::new_sv_surfaces(
            line, sv, 1, 1, true, true, 10.0, 20.0, 30.0, 0.1, 0.2, 0.3, 0.4, true, 1, 3,
        )
    }

    /// OCCT anchor: BRepApprox_ApproxLine spline constructor, NbPnts and the
    /// Point(Index) pole round-trip (cxx L29-36, L47-94).
    #[test]
    fn approx_line_spline_round_trip() {
        let cxyz = Arc::new(BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 2.0, 3.0)],
            weights: vec![1.0, 1.0],
            is_periodic: false,
        });
        let cuv1 = Arc::new(BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                DVec2::new(0.1, 0.2),
                DVec2::new(0.3, 0.4),
                DVec2::new(0.5, 0.6),
            ],
            weights: vec![1.0; 3],
        });
        let cuv2 = Arc::new(BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                DVec2::new(1.1, 1.2),
                DVec2::new(1.3, 1.4),
                DVec2::new(1.5, 1.6),
            ],
            weights: vec![1.0; 3],
        });

        let line = ApproxLine::new(
            Some(cxyz.clone()),
            Some(cuv1.clone()),
            Some(cuv2.clone()),
        );
        // NbPnts reads the first non-null curve: CurveXYZ (cxx L49-52).
        assert_eq!(line.nb_pnts(), 2);
        let pt = line.point(2);
        assert_eq!(pt.value(), DVec3::new(1.0, 2.0, 3.0));
        let (u1, v1, u2, v2) = pt.parameters();
        assert_eq!((u1, v1), (0.3, 0.4));
        assert_eq!((u2, v2), (1.3, 1.4));
    }

    /// OCCT anchor: BRepApprox_ApproxLine(IntSurf_LineOn2S) constructor —
    /// NbPnts is the line's point count and Point(Index) returns the stored
    /// point (cxx L40-43, L47-62, L66-74).
    #[test]
    fn approx_line_on_2s_round_trip() {
        let src = line_on_2s_3pts();
        let line = ApproxLine::new_line_on_2s(Some(src.clone()), false);
        assert_eq!(line.nb_pnts(), 3);
        for i in 1..=3 {
            let p = line.point(i);
            let s = src.value(i - 1);
            assert_eq!(p.value(), s.value());
            assert_eq!(p.parameters(), s.parameters());
        }
    }

    /// OCCT anchor: the default constructor zeroes everything except
    /// p2donfirst (gxx L32-50).
    #[test]
    fn multi_line_default_ctor() {
        let ml = TheMultiLineOfApprox::new();
        assert_eq!(ml.first_point(), 0);
        assert_eq!(ml.last_point(), 0);
        assert_eq!(ml.nb_p3d(), 0);
        assert_eq!(ml.nb_p2d(), 0);
        assert_eq!(ml.what_status(), ApproxStatus::NoPointsAdded);
    }

    /// OCCT anchor: the constructor path, the accessors, and the Value
    /// offset arithmetic (gxx L98-136, L178-210).
    #[test]
    fn multi_line_of_approx_lifecycle() {
        let ml = TheMultiLineOfApprox::new_line(
            Arc::new(ApproxLine::new_line_on_2s(Some(line_on_2s_3pts()), false)),
            1,
            1,
            true,
            true,
            10.0,
            20.0,
            30.0,
            0.1,
            0.2,
            0.3,
            0.4,
            true,
            1,
            3,
        );
        assert_eq!(ml.first_point(), 1);
        assert_eq!(ml.last_point(), 3);
        assert_eq!(ml.nb_p3d(), 1);
        assert_eq!(ml.nb_p2d(), 1);
        assert_eq!(ml.what_status(), ApproxStatus::NoPointsAdded);

        // Value 3d: point(2) = (1, -1, 2) + (10, 20, 30).
        let mut tab3 = [DVec3::ZERO];
        ml.value_3d(2, &mut tab3);
        assert_eq!(tab3[0], DVec3::new(11.0, 19.0, 32.0));

        // Value 2d with P2DOnFirst: (u1 + U1o, v1 + V1o).
        let mut tab2 = [DVec2::ZERO];
        ml.value_2d(2, &mut tab2);
        assert_eq!(tab2[0], DVec2::new(1.1, 2.2));

        // The combined form (gxx L214-220).
        let mut tab3 = [DVec3::ZERO];
        let mut tab2 = [DVec2::ZERO];
        ml.value_3d_2d(1, &mut tab3, &mut tab2);
        assert_eq!(tab3[0], DVec3::new(10.0, 20.0, 30.0));
        assert_eq!(tab2[0], DVec2::new(0.1, 1.2));
    }

    /// OCCT anchor: the two-surface 2d form fills both slots and honours the
    /// Length() >= 2 guards (gxx L202-209); with P2DOnFirst false and
    /// nbp2d == 1 slot 1 takes the second surface (gxx L198-200).
    #[test]
    fn multi_line_value_2d_two_surfaces() {
        let ml = TheMultiLineOfApprox::new_line(
            Arc::new(ApproxLine::new_line_on_2s(Some(line_on_2s_3pts()), false)),
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
            false,
            1,
            3,
        );
        // nbp2d == 2: slots (u1, v1) and (u2, v2); point(3) parameters are
        // (2, 3, 4, 5).
        let mut tab2 = [DVec2::ZERO; 2];
        ml.value_2d(3, &mut tab2);
        assert_eq!(tab2[0], DVec2::new(2.0, 3.0));
        assert_eq!(tab2[1], DVec2::new(4.0, 5.0));

        // Length() < 2: only slot 1 is written (gxx L205-208).
        let mut tab2 = [DVec2::new(-1.0, -1.0)];
        ml.value_2d(3, &mut tab2);
        assert_eq!(tab2[0], DVec2::new(2.0, 3.0));
    }

    /// OCCT anchor: with a null SvSurfaces pointer Tangency returns False
    /// and leaves the tabs untouched (gxx L226-227, L246-247).
    #[test]
    fn multi_line_tangency_without_sv_surfaces() {
        let ml = multi_line_3pts(None);
        let mut tab3 = [DVec3::new(9.0, 9.0, 9.0)];
        assert!(!ml.tangency_3d(1, &mut tab3));
        assert_eq!(tab3[0], DVec3::new(9.0, 9.0, 9.0));
        let mut tab2 = [DVec2::new(9.0, 9.0)];
        assert!(!ml.tangency_2d(1, &mut tab2));
        assert_eq!(tab2[0], DVec2::new(9.0, 9.0));
        assert!(!ml.tangency_3d_2d(1, &mut tab3, &mut tab2));
    }

    /// OCCT anchor: with a non-null SvSurfaces the Tangency forms route to
    /// Tangency / TangencyOnSurf1 (nbp2d == 1 and P2DOnFirst) and
    /// WhatStatus reports PointsAdded (gxx L154-160, L224-289).
    #[test]
    fn multi_line_tangency_with_sv_surfaces() {
        let sv: SharedSvSurfaces = StubSvSurfaces::new();
        let ml = multi_line_3pts(Some(sv.clone()));
        assert!(!sv.get_use_solver());
        assert_eq!(ml.what_status(), ApproxStatus::PointsAdded);

        let mut tab3 = [DVec3::ZERO];
        assert!(ml.tangency_3d(1, &mut tab3));
        assert_eq!(tab3[0], DVec3::new(1.0, 2.0, 3.0));

        let mut tab2 = [DVec2::ZERO];
        assert!(ml.tangency_2d(1, &mut tab2));
        assert_eq!(tab2[0], DVec2::new(0.5, 0.25));

        let mut tab3 = [DVec3::ZERO];
        let mut tab2 = [DVec2::ZERO];
        assert!(ml.tangency_3d_2d(1, &mut tab3, &mut tab2));
        assert_eq!(tab3[0], DVec3::new(1.0, 2.0, 3.0));
        assert_eq!(tab2[0], DVec2::new(0.5, 0.25));
    }

    /// OCCT anchor: MakeMLBetween with a null SvSurfaces returns the empty
    /// 1..1 line carrying the same offsets (gxx L306-328); MakeMLOneMorePoint
    /// returns False (gxx L637-638).
    #[test]
    fn multi_line_make_ml_null_sv_surfaces() {
        let ml = multi_line_3pts(None);
        let sub = ml.make_ml_between(1, 3, 2);
        assert_eq!(sub.first_point(), 1);
        assert_eq!(sub.last_point(), 1);
        assert_eq!(sub.nb_p3d(), 1);
        assert_eq!(sub.nb_p2d(), 1);

        let mut other = TheMultiLineOfApprox::new();
        assert!(!ml.make_ml_one_more_point(1, 3, 2, &mut other));
        assert_eq!(other.first_point(), 0);
    }

    /// OCCT anchor: the MakeMLOneMorePoint validation rejects the midpoint
    /// for a failed SeekPoint and restores the solver flag (gxx L668-717).
    #[test]
    fn multi_line_make_ml_one_more_point_rejected() {
        let sv: SharedSvSurfaces = StubSvSurfaces::new();
        let ml = multi_line_3pts(Some(sv.clone()));
        sv.set_use_solver(false);
        let mut other = TheMultiLineOfApprox::new();
        // The stub SeekPoint returns False -> OtherLineMade stays False.
        assert!(!ml.make_ml_one_more_point(1, 3, 2, &mut other));
        // The solver flag is restored (gxx L715).
        assert!(!sv.get_use_solver());
    }

    /// OCCT anchor: Dump() runs the point loop over [FirstPoint, LastPoint]
    /// (gxx L757-778).
    #[test]
    fn multi_line_dump() {
        let ml = multi_line_3pts(None);
        ml.dump();
    }

    /// OCCT anchor: the LineTool statics forward to the MultiLine
    /// (ApproxInt_MultiLineTool.lxx L29-184); the Curvature statics return
    /// False unconditionally (lxx L113-144).
    #[test]
    fn the_multi_line_tool_statics() {
        let ml = multi_line_3pts(None);
        assert_eq!(the_multi_line_tool::first_point(&ml), 1);
        assert_eq!(the_multi_line_tool::last_point(&ml), 3);
        assert_eq!(the_multi_line_tool::nb_p2d(&ml), 1);
        assert_eq!(the_multi_line_tool::nb_p3d(&ml), 1);

        let mut tab3 = [DVec3::ZERO];
        the_multi_line_tool::value_3d(&ml, 2, &mut tab3);
        assert_eq!(tab3[0], DVec3::new(11.0, 19.0, 32.0));

        let mut tab2 = [DVec2::ZERO];
        the_multi_line_tool::value_2d(&ml, 2, &mut tab2);
        assert_eq!(tab2[0], DVec2::new(1.1, 2.2));

        let mut tab3 = [DVec3::ZERO];
        let mut tab2 = [DVec2::ZERO];
        the_multi_line_tool::value_3d_2d(&ml, 2, &mut tab3, &mut tab2);
        assert_eq!(tab3[0], DVec3::new(11.0, 19.0, 32.0));
        assert_eq!(tab2[0], DVec2::new(1.1, 2.2));

        assert!(!the_multi_line_tool::tangency_3d(&ml, 2, &mut tab3));
        assert!(!the_multi_line_tool::tangency_2d(&ml, 2, &mut tab2));
        assert!(!the_multi_line_tool::tangency_3d_2d(&ml, 2, &mut tab3, &mut tab2));
        assert!(!the_multi_line_tool::curvature_3d(&ml, 2, &mut tab3));
        assert!(!the_multi_line_tool::curvature_2d(&ml, 2, &mut tab2));
        assert!(!the_multi_line_tool::curvature_3d_2d(&ml, 2, &mut tab3, &mut tab2));

        assert_eq!(
            the_multi_line_tool::what_status(&ml, 1, 3),
            ApproxStatus::NoPointsAdded
        );

        let sub = the_multi_line_tool::make_ml_between(&ml, 1, 2, 1);
        assert_eq!(sub.last_point(), 1);

        let mut other = TheMultiLineOfApprox::new();
        assert!(!the_multi_line_tool::make_ml_one_more_point(&ml, 1, 2, 2, &mut other));

        the_multi_line_tool::dump(&ml);
    }

    /// A plane face over the UV window [0, 1] x [0, 1].
    fn build_plane_brep() -> (rcad_kernel::BRep, Shape) {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(&mut brep, Point3::ZERO, 1e-7);
        let v2 = b.add_vertex(&mut brep, Point3::new(1.0, 0.0, 0.0), 1e-7);
        let e = b.add_edge(
            &mut brep,
            Some(rcad_kernel::geom::Curve3::Line(Line3 {
                origin: Point3::ZERO,
                direction: Vec3::new(1.0, 0.0, 0.0),
            })),
            v1,
            v2,
            [0.0, 1.0],
        );
        let wire = brep.add_twire(vec![e]);
        let face = brep.add_tface(
            Some(Surface3::Plane(Plane {
                origin: Point3::ZERO,
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_dir: Vec3::new(1.0, 0.0, 0.0),
                v_dir: Vec3::new(0.0, 1.0, 0.0),
            })),
            wire,
            Vec::new(),
            None,
            Some([0.0, 1.0, 0.0, 1.0]),
            Vec::new(),
            true,
        );
        (brep, face)
    }

    /// OCCT anchor: the SurfaceTool one-liners forward to the adaptor
    /// (BRepApprox_SurfaceTool.lxx) and the NbSamples statics return 10
    /// (BRepApprox_SurfaceTool.cxx L22-40).
    #[test]
    fn surface_tool_plane_statics() {
        let (brep, face) = build_plane_brep();
        let s = BRepAdaptorSurface::initialize_face(&brep, &face, true);

        assert_eq!(surface_tool::first_u_parameter(&s), 0.0);
        assert_eq!(surface_tool::last_u_parameter(&s), 1.0);
        assert_eq!(surface_tool::first_v_parameter(&s), 0.0);
        assert_eq!(surface_tool::last_v_parameter(&s), 1.0);
        assert_eq!(surface_tool::get_type(&s), GeomAbsSurfaceType::Plane);
        assert!(!surface_tool::is_u_closed(&s));
        assert!(!surface_tool::is_v_closed(&s));
        assert!(!surface_tool::is_u_periodic(&s));
        assert!(!surface_tool::is_v_periodic(&s));
        assert_eq!(surface_tool::u_period(&s), 0.0);
        assert_eq!(surface_tool::v_period(&s), 0.0);

        let p = surface_tool::value(&s, 0.25, 0.75);
        assert!((p.x - 0.25).abs() < 1e-12);
        assert!((p.y - 0.75).abs() < 1e-12);
        assert!(p.z.abs() < 1e-12);

        let d0p = surface_tool::d0(&s, 0.25, 0.75);
        assert_eq!(d0p, p);

        let (p1, d1u, d1v) = surface_tool::d1(&s, 0.25, 0.75);
        assert_eq!(p1, p);
        assert_eq!(d1u, Vec3::new(1.0, 0.0, 0.0));
        assert_eq!(d1v, Vec3::new(0.0, 1.0, 0.0));

        let (p2, _d1u, _d1v, d2u, d2v, d2uv) = surface_tool::d2(&s, 0.25, 0.75);
        assert_eq!(p2, p);
        assert_eq!(d2u, Vec3::ZERO);
        assert_eq!(d2v, Vec3::ZERO);
        assert_eq!(d2uv, Vec3::ZERO);

        assert_eq!(surface_tool::u_resolution(&s, 0.01), 0.01);
        assert_eq!(surface_tool::v_resolution(&s, 0.01), 0.01);

        let pln = surface_tool::plane(&s);
        assert_eq!(pln.origin, Point3::ZERO);
        assert_eq!(pln.normal, Vec3::new(0.0, 0.0, 1.0));

        assert_eq!(surface_tool::nb_samples_u(&s), 10);
        assert_eq!(surface_tool::nb_samples_v(&s), 10);
        assert_eq!(surface_tool::nb_samples_u_window(&s, 0.0, 0.5), 10);
        assert_eq!(surface_tool::nb_samples_v_window(&s, 0.0, 0.5), 10);
    }

    /// OCCT anchor: UTrim/VTrim narrow the window on the same surface
    /// (BRepApprox_SurfaceTool.lxx L77-91).
    #[test]
    fn surface_tool_trims() {
        let (brep, face) = build_plane_brep();
        let s = BRepAdaptorSurface::initialize_face(&brep, &face, true);

        let ut = surface_tool::u_trim(&s, 0.2, 0.7, 1e-7);
        assert_eq!(surface_tool::first_u_parameter(&ut), 0.2);
        assert_eq!(surface_tool::last_u_parameter(&ut), 0.7);
        assert_eq!(surface_tool::first_v_parameter(&ut), 0.0);
        assert_eq!(surface_tool::last_v_parameter(&ut), 1.0);

        let vt = surface_tool::v_trim(&s, 0.3, 0.6, 1e-7);
        assert_eq!(surface_tool::first_u_parameter(&vt), 0.0);
        assert_eq!(surface_tool::last_u_parameter(&vt), 1.0);
        assert_eq!(surface_tool::first_v_parameter(&vt), 0.3);
        assert_eq!(surface_tool::last_v_parameter(&vt), 0.6);
    }
}
