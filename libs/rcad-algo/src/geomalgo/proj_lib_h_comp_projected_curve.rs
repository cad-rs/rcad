// OCCT ProjLib_CompProjectedCurve.hxx L37-278 + ProjLib_CompProjectedCurve.cxx
// L52-2391 and ProjLib_HCompProjectedCurve.hxx L17-25
// (ModelingData/TKGeomBase/ProjLib).
//
// 1:1 translation of the compound projected curve (normal projection of a 3D
// curve on a surface) and of the HCompProjectedCurve alias.  The computation
// machinery (Init, Perform, the split statics and the interval builder) lives
// in the sibling module [b] (`proj_lib_h_comp_projected_curve_b.rs`), which is
// the #[path]-declared continuation of this file (the pair models the single
// OCCT class; the split keeps both files under the 2000-line limit).
//
// Dependency encodings:
// - handle(Adaptor3d_Surface) / handle(Adaptor3d_Curve) /
//   handle(Adaptor2d_Curve2d) -> rcad_kernel::base::proj_lib::adaptor,
// - ProjLib_PrjResolve / ProjLib_PrjFunc ->
//   rcad_kernel::base::proj_lib::prj_resolve,
// - NCollection_HSequence<HSequence<gp_Pnt>> -> Vec<Vec<DVec3>> (the triple
//   (t, u, v) stored in the X/Y/Z slots, as in OCCT),
// - NCollection_HArray1<T> -> Option<Vec<T>> (None models the null handle),
// - Extrema_ExtCS / Extrema_ExtCC -> GAP carriers in [b] (staged);
//   Extrema_ExtPS routes to the kernel Surface3 engine through the
//   kernel_surface() bridge,
// - GeomLib::FuseIntervals -> rcad_kernel::base::geom_lib::fuse_intervals.
//
// OCCT ProjLib_HCompProjectedCurve.hxx L23:
//   typedef ProjLib_CompProjectedCurve ProjLib_HCompProjectedCurve;
// — the alias is mirrored as a Rust type alias ([HCompProjectedCurve]).

use std::sync::Arc;

use glam::{DVec2, DVec3};

use rcad_kernel::base::geom_lib::fuse_intervals;
use rcad_kernel::base::proj_lib::adaptor::{
    Adaptor2dCurve2d, Adaptor3dCurve, Adaptor3dSurface, Curve2dHandle, Geom2dCurveAdaptor,
    SurfaceHandle,
};
use rcad_kernel::base::proj_lib::prj_resolve::PrjResolve;
use rcad_kernel::core::precision::{p_confusion, CONFUSION};
use rcad_kernel::geom::{Curve2d, Curve3, Line2d};
use rcad_kernel::math::GeomAbsShape;

use crate::geomalgo::approx_curve_on_surface::ApproxCurveOnSurface;
use crate::geomalgo::int_patch::GeomAbsSurfaceType;

use b::{d1, d2, d2_curv_on_surf, dich_exact_bound, exact_bound, initial_point};

#[path = "proj_lib_h_comp_projected_curve_b.rs"]
pub(crate) mod b;

/// OCCT ProjLib_CompProjectedCurve.cxx L52: #define FuncTol 1.e-10.
pub(crate) const FUNC_TOL: f64 = 1.0e-10;

// ---------------------------------------------------------------------------
// ProjLib_CompProjectedCurve
// ---------------------------------------------------------------------------

/// OCCT ProjLib_CompProjectedCurve (hxx L37-278) — search of all solutions of
/// the normal projection of a 3D curve on a surface.
pub struct CompProjectedCurve {
    /// hxx L251: handle(Adaptor3d_Surface) mySurface.
    pub(crate) my_surface: Option<SurfaceHandle>,
    /// hxx L252: handle(Adaptor3d_Curve) myCurve.
    pub(crate) my_curve: Option<CurveHandleAlias>,
    /// hxx L253: Standard_Integer myNbCurves.
    pub(crate) my_nb_curves: i32,
    /// hxx L254: the (t, u, v) triples per continuous part.
    pub(crate) my_sequence: Option<Vec<Vec<DVec3>>>,
    /// hxx L255: handle(NCollection_HArray1(bool)) myUIso.
    pub(crate) my_u_iso: Option<Vec<bool>>,
    /// hxx L256: handle(NCollection_HArray1(bool)) myVIso.
    pub(crate) my_v_iso: Option<Vec<bool>>,
    /// hxx L257: handle(NCollection_HArray1(bool)) mySnglPnts.
    pub(crate) my_sngl_pnts: Option<Vec<bool>>,
    /// hxx L258: handle(NCollection_HArray1(double)) myMaxDistance.
    pub(crate) my_max_distance: Option<Vec<f64>>,
    /// hxx L259: handle(NCollection_HArray1(double)) myTabInt.
    pub(crate) my_tab_int: Option<Vec<f64>>,
    /// hxx L260: double myTol3d.
    pub(crate) my_tol3d: f64,
    /// hxx L261: GeomAbs_Shape myContinuity.
    pub(crate) my_continuity: GeomAbsShape,
    /// hxx L262: Standard_Integer myMaxDegree.
    pub(crate) my_max_degree: i32,
    /// hxx L263: Standard_Integer myMaxSeg.
    pub(crate) my_max_seg: i32,
    /// hxx L264: bool myProj2d.
    pub(crate) my_proj2d: bool,
    /// hxx L265: bool myProj3d.
    pub(crate) my_proj3d: bool,
    /// hxx L266: double myMaxDist.
    pub(crate) my_max_dist: f64,
    /// hxx L267: double myTolU.
    pub(crate) my_tol_u: f64,
    /// hxx L268: double myTolV.
    pub(crate) my_tol_v: f64,
    /// hxx L270: handle(NCollection_HArray1(bool)) myResultIsPoint.
    pub(crate) my_result_is_point: Option<Vec<bool>>,
    /// hxx L271: handle(NCollection_HArray1(double)) myResult2dUApproxError.
    pub(crate) my_result2d_u_approx_error: Option<Vec<f64>>,
    /// hxx L272: handle(NCollection_HArray1(double)) myResult2dVApproxError.
    pub(crate) my_result2d_v_approx_error: Option<Vec<f64>>,
    /// hxx L273: handle(NCollection_HArray1(double)) myResult3dApproxError.
    pub(crate) my_result3d_approx_error: Option<Vec<f64>>,
    /// hxx L274: handle(NCollection_HArray1(gp_Pnt)) myResult3dPoint.
    pub(crate) my_result3d_point: Option<Vec<DVec3>>,
    /// hxx L275: handle(NCollection_HArray1(gp_Pnt2d)) myResult2dPoint.
    pub(crate) my_result2d_point: Option<Vec<DVec2>>,
    /// hxx L276: handle(NCollection_HArray1(handle(Geom_Curve)))
    /// myResult3dCurve.
    pub(crate) my_result3d_curve: Option<Vec<Option<Curve3>>>,
    /// hxx L277: handle(NCollection_HArray1(handle(Geom2d_Curve)))
    /// myResult2dCurve.
    pub(crate) my_result2d_curve: Option<Vec<Option<Curve2d>>>,
}

/// OCCT `handle(Adaptor3d_Curve)`.
pub(crate) type CurveHandleAlias = Arc<dyn Adaptor3dCurve>;

impl Clone for CompProjectedCurve {
    /// OCCT ShallowCopy (cxx L631-654): the sequence handles are shared; the
    /// rcad Vec fields are cloned (value semantics over shared data — the
    /// copies are never mutated after Init in the OCCT usage).
    fn clone(&self) -> Self {
        CompProjectedCurve {
            my_surface: self.my_surface.clone(),
            my_curve: self.my_curve.clone(),
            my_nb_curves: self.my_nb_curves,
            my_sequence: self.my_sequence.clone(),
            my_u_iso: self.my_u_iso.clone(),
            my_v_iso: self.my_v_iso.clone(),
            my_sngl_pnts: self.my_sngl_pnts.clone(),
            my_max_distance: self.my_max_distance.clone(),
            my_tab_int: self.my_tab_int.clone(),
            my_tol3d: self.my_tol3d,
            my_continuity: self.my_continuity,
            my_max_degree: self.my_max_degree,
            my_max_seg: self.my_max_seg,
            my_proj2d: self.my_proj2d,
            my_proj3d: self.my_proj3d,
            my_max_dist: self.my_max_dist,
            my_tol_u: self.my_tol_u,
            my_tol_v: self.my_tol_v,
            my_result_is_point: self.my_result_is_point.clone(),
            my_result2d_u_approx_error: self.my_result2d_u_approx_error.clone(),
            my_result2d_v_approx_error: self.my_result2d_v_approx_error.clone(),
            my_result3d_approx_error: self.my_result3d_approx_error.clone(),
            my_result3d_point: self.my_result3d_point.clone(),
            my_result2d_point: self.my_result2d_point.clone(),
            my_result3d_curve: self.my_result3d_curve.clone(),
            my_result2d_curve: self.my_result2d_curve.clone(),
        }
    }
}

impl CompProjectedCurve {
    /// OCCT ProjLib_CompProjectedCurve() (cxx L547-553).  The C++ members not
    /// in the initializer list are default-initialized (indeterminate PODs);
    /// rcad uses zero/false defaults — the object is only usable after
    /// Load + Init.
    pub fn new_empty() -> Self {
        CompProjectedCurve {
            my_surface: None,
            my_curve: None,
            my_nb_curves: 0,
            my_sequence: None,
            my_u_iso: None,
            my_v_iso: None,
            my_sngl_pnts: None,
            my_max_distance: None,
            my_tab_int: None,
            my_tol3d: 0.0,
            my_continuity: GeomAbsShape::C0,
            my_max_degree: 0,
            my_max_seg: 0,
            my_proj2d: false,
            my_proj3d: false,
            my_max_dist: 0.0,
            my_tol_u: 0.0,
            my_tol_v: 0.0,
            my_result_is_point: None,
            my_result2d_u_approx_error: None,
            my_result2d_v_approx_error: None,
            my_result3d_approx_error: None,
            my_result3d_point: None,
            my_result2d_point: None,
            my_result3d_curve: None,
            my_result2d_curve: None,
        }
    }

    /// OCCT ProjLib_CompProjectedCurve(S, C, TolU, TolV) (cxx L557-577) —
    /// try to find all solutions.
    pub fn new(
        the_surface: SurfaceHandle,
        the_curve: CurveHandleAlias,
        the_tol_u: f64,
        the_tol_v: f64,
    ) -> Self {
        let mut this = CompProjectedCurve::with_defaults(the_surface, the_curve, the_tol_u, the_tol_v);
        this.init();
        this
    }

    /// OCCT ProjLib_CompProjectedCurve(S, C, TolU, TolV, MaxDist)
    /// (cxx L581-602) — the MaxDist-optimized search; MaxDist < 0 behaves as
    /// the full search.
    pub fn new_with_max_dist(
        the_surface: SurfaceHandle,
        the_curve: CurveHandleAlias,
        the_tol_u: f64,
        the_tol_v: f64,
        the_max_dist: f64,
    ) -> Self {
        let mut this = CompProjectedCurve::with_defaults(the_surface, the_curve, the_tol_u, the_tol_v);
        this.my_max_dist = the_max_dist;
        this.init();
        this
    }

    /// The common initializer of the two TolU/TolV constructors (cxx L562-574).
    fn with_defaults(
        the_surface: SurfaceHandle,
        the_curve: CurveHandleAlias,
        the_tol_u: f64,
        the_tol_v: f64,
    ) -> Self {
        CompProjectedCurve {
            my_surface: Some(the_surface),
            my_curve: Some(the_curve),
            my_nb_curves: 0,
            my_sequence: Some(Vec::new()),
            my_u_iso: None,
            my_v_iso: None,
            my_sngl_pnts: None,
            my_max_distance: None,
            my_tab_int: None,
            my_tol3d: 1.0e-6,
            my_continuity: GeomAbsShape::C2,
            my_max_degree: 14,
            my_max_seg: 16,
            my_proj2d: true,
            my_proj3d: false,
            my_max_dist: -1.0,
            my_tol_u: the_tol_u,
            my_tol_v: the_tol_v,
            my_result_is_point: None,
            my_result2d_u_approx_error: None,
            my_result2d_v_approx_error: None,
            my_result3d_approx_error: None,
            my_result3d_point: None,
            my_result2d_point: None,
            my_result3d_curve: None,
            my_result2d_curve: None,
        }
    }

    /// OCCT ProjLib_CompProjectedCurve(Tol3d, S, C, MaxDist = -1.0)
    /// (cxx L606-627) — tolerances computed from the surface resolutions.
    pub fn new_tol3d(
        the_tol3d: f64,
        the_surface: SurfaceHandle,
        the_curve: CurveHandleAlias,
        the_max_dist: f64,
    ) -> Self {
        let mut this = CompProjectedCurve {
            my_surface: Some(the_surface),
            my_curve: Some(the_curve),
            my_nb_curves: 0,
            my_sequence: Some(Vec::new()),
            my_u_iso: None,
            my_v_iso: None,
            my_sngl_pnts: None,
            my_max_distance: None,
            my_tab_int: None,
            my_tol3d: the_tol3d,
            my_continuity: GeomAbsShape::C2,
            my_max_degree: 14,
            my_max_seg: 16,
            my_proj2d: true,
            my_proj3d: false,
            my_max_dist: the_max_dist,
            my_tol_u: 0.0,
            my_tol_v: 0.0,
            my_result_is_point: None,
            my_result2d_u_approx_error: None,
            my_result2d_v_approx_error: None,
            my_result3d_approx_error: None,
            my_result3d_point: None,
            my_result2d_point: None,
            my_result3d_curve: None,
            my_result2d_curve: None,
        };
        // OCCT L623-624.
        let surf = this.my_surface.as_ref().unwrap();
        this.my_tol_u = p_confusion().max(surf.u_resolution(the_tol3d));
        this.my_tol_v = p_confusion().max(surf.v_resolution(the_tol3d));

        this.init();
        this
    }

    /// OCCT ProjLib_CompProjectedCurve::Init (cxx L658-1225) — delegates to
    /// the sibling module (the 1:1 body).
    pub(crate) fn init(&mut self) {
        b::init_body(self);
    }

    /// OCCT ProjLib_CompProjectedCurve::Perform (cxx L1229-1408) — delegates
    /// to the sibling module.
    pub(crate) fn perform(&mut self) {
        b::perform_body(self);
    }

    /// OCCT SetTol3d (cxx L1412-1415).
    pub fn set_tol3d(&mut self, the_tol3d: f64) {
        self.my_tol3d = the_tol3d;
    }

    /// OCCT SetContinuity (cxx L1419-1422).
    pub fn set_continuity(&mut self, the_continuity: GeomAbsShape) {
        self.my_continuity = the_continuity;
    }

    /// OCCT SetMaxDegree (cxx L1426-1433).
    pub fn set_max_degree(&mut self, the_max_degree: i32) {
        if the_max_degree < 1 {
            return;
        }
        self.my_max_degree = the_max_degree;
    }

    /// OCCT SetMaxSeg (cxx L1437-1444).
    pub fn set_max_seg(&mut self, the_max_seg: i32) {
        if the_max_seg < 1 {
            return;
        }
        self.my_max_seg = the_max_seg;
    }

    /// OCCT SetProj3d (cxx L1448-1451).
    pub fn set_proj3d(&mut self, the_proj3d: bool) {
        self.my_proj3d = the_proj3d;
    }

    /// OCCT SetProj2d (cxx L1455-1458).
    pub fn set_proj2d(&mut self, the_proj2d: bool) {
        self.my_proj2d = the_proj2d;
    }

    /// OCCT Load(const handle(Adaptor3d_Surface)& S) (cxx L1462-1465).
    pub fn load_surface(&mut self, s: SurfaceHandle) {
        self.my_surface = Some(s);
    }

    /// OCCT Load(const handle(Adaptor3d_Curve)& C) (cxx L1469-1472).
    pub fn load_curve(&mut self, c: CurveHandleAlias) {
        self.my_curve = Some(c);
    }

    /// OCCT GetSurface (cxx L1476-1479).
    pub fn get_surface(&self) -> &SurfaceHandle {
        self.my_surface.as_ref().expect("GetSurface()")
    }

    /// OCCT GetCurve (cxx L1483-1486).
    pub fn get_curve(&self) -> &CurveHandleAlias {
        self.my_curve.as_ref().expect("GetCurve()")
    }

    /// OCCT GetTolerance (cxx L1490-1494).
    pub fn get_tolerance(&self, tol_u: &mut f64, tol_v: &mut f64) {
        *tol_u = self.my_tol_u;
        *tol_v = self.my_tol_v;
    }

    /// OCCT NbCurves (cxx L1498-1501).
    pub fn nb_curves(&self) -> i32 {
        self.my_nb_curves
    }

    /// OCCT Bounds (cxx L1505-1513) — throws Standard_NoSuchObject when the
    /// index is out of range.
    pub fn bounds(&self, index: i32, udeb: &mut f64, ufin: &mut f64) {
        if index < 1 || index > self.my_nb_curves {
            panic!("Standard_NoSuchObject");
        }
        let seq = self.my_sequence.as_ref().unwrap();
        let part = &seq[(index - 1) as usize];
        *udeb = part[0].x;
        *ufin = part[part.len() - 1].x;
    }

    /// OCCT IsSinglePnt (cxx L1517-1525).
    pub fn is_single_pnt(&self, index: i32, p: &mut DVec2) -> bool {
        if index < 1 || index > self.my_nb_curves {
            panic!("Standard_NoSuchObject");
        }
        let seq = self.my_sequence.as_ref().unwrap();
        let first = &seq[(index - 1) as usize][0];
        *p = DVec2::new(first.y, first.z);
        self.my_sngl_pnts.as_ref().unwrap()[(index - 1) as usize]
    }

    /// OCCT IsUIso (cxx L1529-1537).
    pub fn is_u_iso(&self, index: i32, u: &mut f64) -> bool {
        if index < 1 || index > self.my_nb_curves {
            panic!("Standard_NoSuchObject");
        }
        let seq = self.my_sequence.as_ref().unwrap();
        let first = &seq[(index - 1) as usize][0];
        *u = first.y;
        self.my_u_iso.as_ref().unwrap()[(index - 1) as usize]
    }

    /// OCCT IsVIso (cxx L1541-1549).
    pub fn is_v_iso(&self, index: i32, v: &mut f64) -> bool {
        if index < 1 || index > self.my_nb_curves {
            panic!("Standard_NoSuchObject");
        }
        let seq = self.my_sequence.as_ref().unwrap();
        let first = &seq[(index - 1) as usize][0];
        *v = first.z;
        self.my_v_iso.as_ref().unwrap()[(index - 1) as usize]
    }

    /// OCCT Value (cxx L1553-1558).
    pub fn value(&self, t: f64) -> DVec2 {
        self.d0(t)
    }

    /// OCCT D0 (cxx L1562-1704).
    pub fn d0(&self, u: f64) -> DVec2 {
        let mut found = false;
        let mut i = 0usize;
        let mut udeb = 0.0f64;
        let mut ufin = 0.0f64;

        for idx in 1..=self.my_nb_curves as usize {
            self.bounds(idx as i32, &mut udeb, &mut ufin);
            if u >= udeb && u <= ufin {
                found = true;
                i = idx;
                break;
            }
        }
        if !found {
            // OCCT: throw Standard_DomainError("ProjLib_CompProjectedCurve::D0");
            panic!("Standard_DomainError: ProjLib_CompProjectedCurve::D0");
        }

        let seq = self.my_sequence.as_ref().unwrap();
        let part = &seq[i - 1];
        let end = part.len();
        let mut j = 1usize; // 1-based j as in OCCT (part[j], part[j+1])
        while j < end {
            if u >= part[j - 1].x && u <= part[j].x {
                break;
            }
            j += 1;
        }

        let mut u0;
        let mut v0;

        // Cubic Interpolation (cxx L1597-1663).
        if part.len() < 4 || (u - part[j - 1].x).abs() <= p_confusion() {
            u0 = part[j - 1].y;
            v0 = part[j - 1].z;
        } else if (u - part[j].x).abs() <= p_confusion() {
            u0 = part[j].y;
            v0 = part[j].z;
        } else {
            let mut jj = j;
            if jj == 1 {
                jj = 2;
            }
            if jj > part.len() - 2 {
                jj = part.len() - 2;
            }

            let x1 = part[jj - 2].x;
            let x2 = part[jj - 1].x;
            let x3 = part[jj].x;
            let x4 = part[jj + 1].x;

            let y1 = DVec2::new(part[jj - 2].y, part[jj - 2].z);
            let y2 = DVec2::new(part[jj - 1].y, part[jj - 1].z);
            let y3 = DVec2::new(part[jj].y, part[jj].z);
            let y4 = DVec2::new(part[jj + 1].y, part[jj + 1].z);

            let i1 = (y1 - y2) / (x1 - x2);
            let i2 = (y2 - y3) / (x2 - x3);
            let i3 = (y3 - y4) / (x3 - x4);

            let i21 = (i1 - i2) / (x1 - x3);
            let i22 = (i2 - i3) / (x2 - x4);

            let i31 = (i21 - i22) / (x1 - x4);

            let res = y1 + (u - x1) * (i1 + (u - x2) * (i21 + (u - x3) * i31));

            u0 = res.x;
            v0 = res.y;

            let surf = self.my_surface.as_ref().unwrap();
            if u0 < surf.first_u_parameter() {
                u0 = surf.first_u_parameter();
            } else if u0 > surf.last_u_parameter() {
                u0 = surf.last_u_parameter();
            }

            if v0 < surf.first_v_parameter() {
                v0 = surf.first_v_parameter();
            } else if v0 > surf.last_v_parameter() {
                v0 = surf.last_v_parameter();
            }
        }
        // End of cubic interpolation

        // OCCT L1666-1673: ProjLib_PrjResolve aPrjPS(*myCurve, *mySurface, 1)
        // + Perform.
        let curve = self.my_curve.as_ref().unwrap();
        let surf = self.my_surface.as_ref().unwrap();
        let mut p = DVec2::new(0.0, 0.0);
        {
            let mut a_prj_ps = PrjResolve::new(curve.as_ref(), surf.as_ref(), 1);
            a_prj_ps.perform(
                u,
                u0,
                v0,
                DVec2::new(self.my_tol_u, self.my_tol_v),
                DVec2::new(surf.first_u_parameter(), surf.first_v_parameter()),
                DVec2::new(surf.last_u_parameter(), surf.last_v_parameter()),
                FUNC_TOL,
                false,
            );
            if a_prj_ps.is_done() {
                p = a_prj_ps.solution();
            } else {
                // OCCT L1679-1703: the Extrema_ExtPS fallback
                // (Extrema_ExtFlag_MIN) routed to the kernel engine through
                // the kernel_surface bridge.
                let the_point = curve.value(u);
                if let Some(kernel_surf) = surf.kernel_surface() {
                    let a_ext_ps = rcad_kernel::base::extrema::ExtPS::new(
                        the_point,
                        kernel_surf,
                        self.my_tol_u,
                        self.my_tol_v,
                    );
                    if a_ext_ps.is_done() && a_ext_ps.nb_ext() > 0 {
                        let mut imin = 1usize;
                        for k in 2..=a_ext_ps.nb_ext() {
                            if a_ext_ps.square_distance(k) < a_ext_ps.square_distance(imin) {
                                imin = k;
                            }
                        }
                        let p_ons = a_ext_ps.point(imin);
                        p = DVec2::new(p_ons.u, p_ons.v);
                    } else {
                        p = DVec2::new(u0, v0);
                    }
                } else {
                    p = DVec2::new(u0, v0);
                }
            }
        }
        p
    }

    /// OCCT D1 (cxx L1708-1715).
    pub fn d1(&self, t: f64) -> (DVec2, DVec2) {
        let p = self.d0(t);
        let u = p.x;
        let v = p.y;
        let curve = self.my_curve.as_ref().unwrap();
        let surf = self.my_surface.as_ref().unwrap();
        let mut vel = DVec2::ZERO;
        d1(t, u, v, &mut vel, curve.as_ref(), surf.as_ref());
        (p, vel)
    }

    /// OCCT D2 (cxx L1719-1726).
    pub fn d2(&self, t: f64) -> (DVec2, DVec2, DVec2) {
        let p = self.d0(t);
        let u = p.x;
        let v = p.y;
        let curve = self.my_curve.as_ref().unwrap();
        let surf = self.my_surface.as_ref().unwrap();
        let mut v1 = DVec2::ZERO;
        let mut v2 = DVec2::ZERO;
        d2(t, u, v, &mut v1, &mut v2, curve.as_ref(), surf.as_ref());
        (p, v1, v2)
    }

    /// OCCT DN (cxx L1730-1755).
    pub fn dn(&self, t: f64, n: i32) -> DVec2 {
        if n < 1 {
            panic!("Standard_OutOfRange: ProjLib_CompProjectedCurve : N must be greater than 0");
        } else if n == 1 {
            let (_p, v) = self.d1(t);
            return v;
        } else if n == 2 {
            let (_p, _v1, v2) = self.d2(t);
            return v2;
        } else if n > 2 {
            panic!("Standard_NotImplemented: ProjLib_CompProjectedCurve::DN");
        }
        DVec2::ZERO
    }

    /// OCCT GetSequence (cxx L1759-1763).
    pub fn get_sequence(&self) -> &Vec<Vec<DVec3>> {
        self.my_sequence.as_ref().expect("GetSequence()")
    }

    /// OCCT FirstParameter (cxx L1767-1770).
    pub fn first_parameter(&self) -> f64 {
        self.my_curve.as_ref().unwrap().first_parameter()
    }

    /// OCCT LastParameter (cxx L1774-1777).
    pub fn last_parameter(&self) -> f64 {
        self.my_curve.as_ref().unwrap().last_parameter()
    }

    /// OCCT Continuity (cxx L1781-1796).
    pub fn continuity(&self) -> GeomAbsShape {
        let mut cont_c = self.my_curve.as_ref().unwrap().continuity();
        let cont_su = self.my_surface.as_ref().unwrap().u_continuity();
        if (cont_su as u8) < (cont_c as u8) {
            cont_c = cont_su;
        }
        let cont_sv = self.my_surface.as_ref().unwrap().v_continuity();
        if (cont_sv as u8) < (cont_c as u8) {
            cont_c = cont_sv;
        }
        cont_c
    }

    /// OCCT MaxDistance (cxx L1800-1807).
    pub fn max_distance(&self, index: i32) -> f64 {
        if index < 1 || index > self.my_nb_curves {
            panic!("Standard_NoSuchObject");
        }
        self.my_max_distance.as_ref().unwrap()[(index - 1) as usize]
    }

    /// OCCT NbIntervals (cxx L1811-1816).
    pub fn nb_intervals(&mut self, s: GeomAbsShape) -> i32 {
        self.my_tab_int = None;
        self.build_intervals(s);
        self.my_tab_int.as_ref().unwrap().len() as i32 - 1
    }

    /// OCCT Intervals (cxx L1820-1828).
    pub fn intervals(&mut self, t: &mut Vec<f64>, s: GeomAbsShape) {
        if self.my_tab_int.is_none() {
            self.build_intervals(s);
        }
        *t = self.my_tab_int.as_ref().unwrap().clone();
    }

    /// OCCT BuildIntervals (cxx L1832-2092).
    fn build_intervals(&mut self, s: GeomAbsShape) {
        // OCCT L1834-1854: SforS mapping.
        let sfor_s = match s {
            GeomAbsShape::C0 => GeomAbsShape::C1,
            GeomAbsShape::C1 => GeomAbsShape::C2,
            GeomAbsShape::C2 => GeomAbsShape::C3,
            GeomAbsShape::C3 => GeomAbsShape::CN,
            GeomAbsShape::CN => GeomAbsShape::CN,
        };
        let curve = self.my_curve.as_ref().unwrap();
        let surf = self.my_surface.as_ref().unwrap();

        let nb_int_cur = curve.nb_intervals(s);
        let nb_int_sur_u = surf.nb_u_intervals(sfor_s);
        let nb_int_sur_v = surf.nb_v_intervals(sfor_s);

        let cut_pnts_t: Vec<f64> = curve.intervals(s);
        let cut_pnts_u: Vec<f64> = surf.u_intervals(sfor_s);
        let cut_pnts_v: Vec<f64> = surf.v_intervals(sfor_s);

        // processing projection bounds (cxx L1873-1878).
        let mut b_arr = vec![0.0f64; 2 * self.my_nb_curves as usize];
        for i in 1..=self.my_nb_curves {
            self.bounds(
                i,
                &mut b_arr[(2 * i - 1) as usize - 1],
                &mut b_arr[(2 * i) as usize - 1],
            );
        }

        // processing curve discontinuities (cxx L1880-1888).
        let mut c_arr: Option<Vec<f64>> = None;
        if nb_int_cur > 1 {
            let mut arr = vec![0.0f64; nb_int_cur - 1];
            for i in 1..=(nb_int_cur - 1) {
                arr[i - 1] = cut_pnts_t[i]; // OCCT: CutPntsT(i + 1)
            }
            c_arr = Some(arr);
        }

        // processing U-surface discontinuities (cxx L1890-1964).
        let mut t_udisc: Vec<f64> = Vec::new();
        for k in 2..=nb_int_sur_u as i32 {
            for i in 1..=self.my_nb_curves {
                let seq = self.my_sequence.as_ref().unwrap();
                let part = &seq[(i - 1) as usize];
                for j in 1..part.len() as i32 {
                    let ul = part[(j - 1) as usize].y;
                    let ur = part[j as usize].y;

                    if (ul - cut_pnts_u[(k - 1) as usize]).abs() <= self.my_tol_u {
                        t_udisc.push(part[(j - 1) as usize].x);
                    } else if (ur - cut_pnts_u[(k - 1) as usize]).abs() <= self.my_tol_u {
                        t_udisc.push(part[j as usize].x);
                    } else if (ul < cut_pnts_u[(k - 1) as usize]
                        && cut_pnts_u[(k - 1) as usize] < ur)
                        || (ur < cut_pnts_u[(k - 1) as usize]
                            && cut_pnts_u[(k - 1) as usize] < ul)
                    {
                        let v = (part[(j - 1) as usize].z + part[j as usize].z) / 2.0;
                        let mut solver =
                            PrjResolve::new(curve.as_ref(), surf.as_ref(), 2);

                        let mut d = DVec2::ZERO;
                        let triple = part[(j - 1) as usize];
                        d1(triple.x, triple.y, triple.z, &mut d, curve.as_ref(), surf.as_ref());
                        let tol = if d.x.abs() < CONFUSION {
                            self.my_tol_u
                        } else {
                            self.my_tol_u.min(self.my_tol_u / d.x.abs())
                        };

                        let tl = part[(j - 1) as usize].x;
                        let tr = part[j as usize].x;

                        solver.perform(
                            (tl + tr) / 2.0,
                            cut_pnts_u[(k - 1) as usize],
                            v,
                            DVec2::new(tol, self.my_tol_v),
                            DVec2::new(tl, surf.first_v_parameter()),
                            DVec2::new(tr, surf.last_v_parameter()),
                            FUNC_TOL,
                            false,
                        );
                        if solver.is_done() {
                            t_udisc.push(solver.solution().x);
                        }
                    }
                }
            }
        }
        let mut i2 = 1usize;
        while i2 < t_udisc.len() {
            if t_udisc[i2] - t_udisc[i2 - 1] < p_confusion() {
                t_udisc.remove(i2);
            } else {
                i2 += 1;
            }
        }

        // processing V-surface discontinuities (cxx L1965-2041).
        let mut t_vdisc: Vec<f64> = Vec::new();
        for k in 2..=nb_int_sur_v as i32 {
            for i in 1..=self.my_nb_curves {
                let seq = self.my_sequence.as_ref().unwrap();
                let part = &seq[(i - 1) as usize];
                for j in 1..part.len() as i32 {
                    let vl = part[(j - 1) as usize].z;
                    let vr = part[j as usize].z;

                    if (vl - cut_pnts_v[(k - 1) as usize]).abs() <= self.my_tol_v {
                        t_vdisc.push(part[(j - 1) as usize].x);
                    } else if (vr - cut_pnts_v[(k - 1) as usize]).abs() <= self.my_tol_v {
                        t_vdisc.push(part[j as usize].x);
                    } else if (vl < cut_pnts_v[(k - 1) as usize]
                        && cut_pnts_v[(k - 1) as usize] < vr)
                        || (vr < cut_pnts_v[(k - 1) as usize]
                            && cut_pnts_v[(k - 1) as usize] < vl)
                    {
                        let u = (part[(j - 1) as usize].y + part[j as usize].y) / 2.0;
                        let mut solver =
                            PrjResolve::new(curve.as_ref(), surf.as_ref(), 3);

                        let mut d = DVec2::ZERO;
                        let triple = part[(j - 1) as usize];
                        d1(triple.x, triple.y, triple.z, &mut d, curve.as_ref(), surf.as_ref());
                        let tol = if d.y.abs() < CONFUSION {
                            self.my_tol_v
                        } else {
                            self.my_tol_v.min(self.my_tol_v / d.y.abs())
                        };

                        let tl = part[(j - 1) as usize].x;
                        let tr = part[j as usize].x;

                        solver.perform(
                            (tl + tr) / 2.0,
                            u,
                            cut_pnts_v[(k - 1) as usize],
                            DVec2::new(tol, self.my_tol_v),
                            DVec2::new(tl, surf.first_u_parameter()),
                            DVec2::new(tr, surf.last_u_parameter()),
                            FUNC_TOL,
                            false,
                        );
                        if solver.is_done() {
                            t_vdisc.push(solver.solution().x);
                        }
                    }
                }
            }
        }
        let mut i3 = 1usize;
        while i3 < t_vdisc.len() {
            if t_vdisc[i3] - t_vdisc[i3 - 1] < p_confusion() {
                t_vdisc.remove(i3);
            } else {
                i3 += 1;
            }
        }

        // fusion (cxx L2043-2084).
        let mut fusion: Vec<f64> = Vec::new();
        if let Some(c_arr) = &c_arr {
            fuse_intervals(&b_arr, c_arr, &mut fusion, p_confusion(), false);
            b_arr = fusion.clone();
            fusion.clear();
        }
        if !t_udisc.is_empty() {
            fuse_intervals(&b_arr, &t_udisc, &mut fusion, p_confusion(), false);
            b_arr = fusion.clone();
            fusion.clear();
        }
        if !t_vdisc.is_empty() {
            fuse_intervals(&b_arr, &t_vdisc, &mut fusion, p_confusion(), false);
            b_arr = fusion.clone();
        }

        self.my_tab_int = Some(b_arr);
    }

    /// OCCT Trim (cxx L2096-2104).
    fn trim_curve(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor2dCurve2d> {
        // OCCT: HCS = new ProjLib_HCompProjectedCurve(*this); (the alias is
        // ProjLib_CompProjectedCurve).
        let mut hcs = self.clone();
        hcs.load_surface(self.my_surface.as_ref().unwrap().clone());
        hcs.load_curve(self.my_curve.as_ref().unwrap().trim(first, last, tol));
        Arc::new(hcs)
    }

    /// OCCT GetType (cxx L2108-2111).
    pub fn get_type(&self) -> rcad_kernel::base::proj_lib::CurveType {
        rcad_kernel::base::proj_lib::CurveType::Other
    }

    /// OCCT ResultIsPoint (cxx L2115-2118).
    pub fn result_is_point(&self, the_index: i32) -> bool {
        self.my_result_is_point.as_ref().unwrap()[(the_index - 1) as usize]
    }

    /// OCCT GetResult2dUApproxError (cxx L2122-2125).
    pub fn get_result2d_u_approx_error(&self, the_index: i32) -> f64 {
        self.my_result2d_u_approx_error.as_ref().unwrap()[(the_index - 1) as usize]
    }

    /// OCCT GetResult2dVApproxError (cxx L2129-2132).
    pub fn get_result2d_v_approx_error(&self, the_index: i32) -> f64 {
        self.my_result2d_v_approx_error.as_ref().unwrap()[(the_index - 1) as usize]
    }

    /// OCCT GetResult3dApproxError (cxx L2136-2139).
    pub fn get_result3d_approx_error(&self, the_index: i32) -> f64 {
        self.my_result3d_approx_error.as_ref().unwrap()[(the_index - 1) as usize]
    }

    /// OCCT GetResult2dC (cxx L2143-2146).
    pub fn get_result2d_c(&self, the_index: i32) -> Option<Curve2d> {
        self.my_result2d_curve.as_ref().unwrap()[(the_index - 1) as usize].clone()
    }

    /// OCCT GetResult3dC (cxx L2150-2153).
    pub fn get_result3d_c(&self, the_index: i32) -> Option<Curve3> {
        self.my_result3d_curve.as_ref().unwrap()[(the_index - 1) as usize].clone()
    }

    /// OCCT GetResult2dP (cxx L2157-2162).
    pub fn get_result2d_p(&self, the_index: i32) -> DVec2 {
        assert!(
            self.result_is_point(the_index),
            "ProjLib_CompProjectedCurve : result is not a point 2d"
        );
        self.my_result2d_point.as_ref().unwrap()[(the_index - 1) as usize]
    }

    /// OCCT GetResult3dP (cxx L2166-2171).
    pub fn get_result3d_p(&self, the_index: i32) -> DVec3 {
        assert!(
            self.result_is_point(the_index),
            "ProjLib_CompProjectedCurve : result is not a point 3d"
        );
        self.my_result3d_point.as_ref().unwrap()[(the_index - 1) as usize]
    }

    /// OCCT GetProj2d (hxx L238).
    pub fn get_proj2d(&self) -> bool {
        self.my_proj2d
    }

    /// OCCT GetProj3d (hxx L241).
    pub fn get_proj3d(&self) -> bool {
        self.my_proj3d
    }

    /// OCCT UpdateTripleByTrapCriteria (cxx L2175-2238).
    pub(crate) fn update_triple_by_trap_criteria(&self, the_point: &mut DVec3) {
        let mut is_problems_possible = false;
        let surf = self.my_surface.as_ref().unwrap();

        // Check possible traps cases:

        // 25892 bug (cxx L2181-2192).
        if surf.get_type() == GeomAbsSurfaceType::SurfaceOfRevolution {
            let a_v_res = surf.v_resolution(CONFUSION);
            let a_max_tol = p_confusion().max(a_v_res);

            if (the_point.z - surf.first_v_parameter()).abs() < a_max_tol
                || (the_point.z - surf.last_v_parameter()).abs() < a_max_tol
            {
                is_problems_possible = true;
            }
        }

        // 27135 bug. Trap on degenerated edge (cxx L2195-2202).
        if surf.get_type() == GeomAbsSurfaceType::Sphere
            && ((the_point.z - surf.first_v_parameter()).abs() < p_confusion()
                || (the_point.z - surf.last_v_parameter()).abs() < p_confusion()
                || (the_point.y - surf.first_u_parameter()).abs() < p_confusion()
                || (the_point.y - surf.last_u_parameter()).abs() < p_confusion())
        {
            is_problems_possible = true;
        }

        if !is_problems_possible {
            return;
        }

        let curve = self.my_curve.as_ref().unwrap();
        let mut u = 0.0f64;
        let mut v = 0.0f64;
        let is_done = initial_point(
            curve.value(the_point.x),
            the_point.x,
            curve.as_ref(),
            surf.as_ref(),
            p_confusion(),
            p_confusion(),
            &mut u,
            &mut v,
            self.my_max_dist,
        );

        if !is_done {
            return;
        }

        // Restore original position in case of period jump (cxx L2226-2235).
        if surf.is_u_periodic()
            && ((u - the_point.y).abs() - surf.u_period()).abs() < p_confusion()
        {
            u = the_point.y;
        }
        if surf.is_v_periodic()
            && ((v - the_point.z).abs() - surf.v_period()).abs() < p_confusion()
        {
            v = the_point.z;
        }
        the_point.y = u;
        the_point.z = v;
    }
}

// ---------------------------------------------------------------------------
// The Adaptor2d_Curve2d inheritance
// ---------------------------------------------------------------------------

/// OCCT: class ProjLib_CompProjectedCurve : public Adaptor2d_Curve2d — the
/// adaptor interface routed to the members.
impl Adaptor2dCurve2d for CompProjectedCurve {
    /// OCCT FirstParameter.
    fn first_parameter(&self) -> f64 {
        CompProjectedCurve::first_parameter(self)
    }

    /// OCCT LastParameter.
    fn last_parameter(&self) -> f64 {
        CompProjectedCurve::last_parameter(self)
    }

    /// OCCT Value.
    fn value(&self, u: f64) -> DVec2 {
        CompProjectedCurve::value(self, u)
    }

    /// OCCT D0.
    fn d0(&self, u: f64) -> DVec2 {
        CompProjectedCurve::d0(self, u)
    }

    /// OCCT D1.
    fn d1(&self, u: f64) -> (DVec2, DVec2) {
        CompProjectedCurve::d1(self, u)
    }

    /// OCCT D2.
    fn d2(&self, u: f64) -> (DVec2, DVec2, DVec2) {
        CompProjectedCurve::d2(self, u)
    }

    /// OCCT Continuity.
    fn continuity(&self) -> GeomAbsShape {
        CompProjectedCurve::continuity(self)
    }

    /// OCCT GetType.
    fn get_type(&self) -> rcad_kernel::base::proj_lib::CurveType {
        CompProjectedCurve::get_type(self)
    }

    /// OCCT Line() — never valid (GetType() is OtherCurve).
    fn line(&self) -> Line2d {
        unreachable!("ProjLib_CompProjectedCurve::Line")
    }

    /// OCCT BSpline() — never valid (GetType() is OtherCurve).
    fn bspline(&self) -> Option<rcad_kernel::geom::BSplineCurve2> {
        None
    }

    /// OCCT Bezier() — never valid (GetType() is OtherCurve).
    fn bezier(&self) -> Option<rcad_kernel::geom::BezierCurve2> {
        None
    }

    /// OCCT Trim.
    fn trim(&self, first_param: f64, last_param: f64, tol: f64) -> Arc<dyn Adaptor2dCurve2d> {
        self.trim_curve(first_param, last_param, tol)
    }
}

// ---------------------------------------------------------------------------
// ProjLib_HCompProjectedCurve
// ---------------------------------------------------------------------------

/// OCCT ProjLib_HCompProjectedCurve.hxx L23:
/// `typedef ProjLib_CompProjectedCurve ProjLib_HCompProjectedCurve;`
pub type HCompProjectedCurve = CompProjectedCurve;

/// OCCT `new ProjLib_HCompProjectedCurve(S, C, TolU, TolV, MaxDist)` — the
/// 5-argument constructor form the BRepAlgoAPI consumers use
/// (BRepLib_MakeBRep FaceAttacher normal-projection path).
pub fn new_h_comp_projected_curve(
    the_surface: SurfaceHandle,
    the_curve: CurveHandleAlias,
    the_tol_u: f64,
    the_tol_v: f64,
    the_max_dist: f64,
) -> HCompProjectedCurve {
    CompProjectedCurve::new_with_max_dist(
        the_surface,
        the_curve,
        the_tol_u,
        the_tol_v,
        the_max_dist,
    )
}

/// Build the Geom2dAdaptor_Curve handle over a 2D curve (OCCT
/// `new Geom2dAdaptor_Curve(PCur2d)` in Perform).
pub(crate) fn geom2d_adaptor_curve(curve: Curve2d) -> Curve2dHandle {
    Arc::new(Geom2dCurveAdaptor::new(curve))
}
