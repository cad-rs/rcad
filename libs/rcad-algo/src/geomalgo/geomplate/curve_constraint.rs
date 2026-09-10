//! OCCT GeomPlate_CurveConstraint (TKGeomAlgo/GeomPlate) — 1:1 port of
//! GeomPlate_CurveConstraint.cxx (whole file) and of the hxx member block.
//!
//! Architecture mappings:
//! - `occ::handle<Adaptor3d_CurveOnSurface>` -> `Option<Arc<CurveOnSurface>>`
//!   (the kernel `base::proj_lib::adaptor::CurveOnSurface`, the sanctioned
//!   Adaptor3d_CurveOnSurface reuse),
//! - `occ::handle<Adaptor3d_Curve>` -> `Option<Arc<dyn Adaptor3dCurve>>`;
//!   the OCCT ctor discrimination `down_cast<Adaptor3d_CurveOnSurface>`
//!   (null / not-null) is carried by the [`CurveBoundary`] enum because the
//!   kernel adaptor trait exposes no downcast,
//! - `occ::handle<Geom2d_Curve>` -> `Option<Curve2d>`,
//! - `occ::handle<Adaptor2d_Curve2d>` -> `Option<Curve2dHandle>`,
//! - `occ::handle<Law_Function>` -> `Option<LawFunctionHandle>`,
//! - `GeomLProp_SLProps myLProp` -> GAP leaf (the package is untranslated);
//!   the `SetSurface` step keeps the OCCT GeomAdaptor_Surface check through
//!   the `kernel_surface()` bridge and the stored basis surface, and
//!   [`CurveConstraint::lprop_surf`] preserves the OCCT failure path.
//! - `Approx_Curve2d` (the myHCurve2d branch of Curve2dOnSurf) -> GAP leaf
//!   (the Approx/AdvApprox package is untranslated).
//!
//! gp_Pnt/gp_Pnt2d/gp_Vec/gp_Vec2d -> glam DVec3/DVec2 (architecture).

use std::sync::Arc;

use glam::DVec3;

use rcad_kernel::base::proj_lib::adaptor::{
    Adaptor3dCurve, Curve2dHandle, CurveHandle, CurveOnSurface, SurfaceHandle,
};
use rcad_kernel::geom::Curve2d;

use crate::geomalgo::gcpnts_abscissa_point::gcpnts_length_3d;
use crate::geomalgo::gcpnts_curve::{GCPntsCurve, GCPntsCurve3dHandle};
use crate::geomalgo::law::law_function::LawFunctionHandle;

/// OCCT ctor-1 parameter `const handle<Adaptor3d_Curve>& Boundary` — the
/// handle may hold a curve-on-surface (the down_cast succeeds) or another
/// 3D adaptor curve (the down_cast is null).
#[derive(Clone)]
pub enum CurveBoundary {
    /// OCCT: the handle holds an Adaptor3d_CurveOnSurface.
    OnSurface(Arc<CurveOnSurface>),
    /// OCCT: the down_cast<Adaptor3d_CurveOnSurface> is null.
    Curve(CurveHandle),
}

/// OCCT GeomPlate_CurveConstraint.
#[derive(Clone)]
pub struct CurveConstraint {
    /// hxx L146: handle(Adaptor3d_CurveOnSurface) myFrontiere.
    my_frontiere: Option<Arc<CurveOnSurface>>,
    /// hxx L147: Standard_Integer myNbPoints.
    my_nb_points: i32,
    /// hxx L148: Standard_Integer myOrder.
    my_order: i32,
    /// hxx L149: handle(Adaptor3d_Curve) my3dCurve.
    my3d_curve: Option<CurveHandle>,
    /// hxx L150: Standard_Integer myTang.
    my_tang: i32,
    /// hxx L151: handle(Geom2d_Curve) my2dCurve.
    my2d_curve: Option<Curve2d>,
    /// hxx L152: handle(Adaptor2d_Curve2d) myHCurve2d.
    my_hcurve2d: Option<Curve2dHandle>,
    /// hxx L153-155: handle(Law_Function) myG0Crit / myG1Crit / myG2Crit.
    my_g0_crit: Option<LawFunctionHandle>,
    my_g1_crit: Option<LawFunctionHandle>,
    my_g2_crit: Option<LawFunctionHandle>,
    /// hxx L156-158: bool myConstG0 / myConstG1 / myConstG2.
    my_const_g0: bool,
    my_const_g1: bool,
    my_const_g2: bool,
    /// hxx L159: GeomLProp_SLProps myLProp — GAP leaf; the surface set by
    /// the OCCT `myLProp.SetSurface(Surf)` step is kept for the LPropSurf
    /// consumers (LocalAnalysis chains, themselves staged).
    my_lprop_surface: Option<rcad_kernel::geom::Surface3>,
    /// hxx L160-162: double myTolDist / myTolAng / myTolCurv.
    my_tol_dist: f64,
    my_tol_ang: f64,
    my_tol_curv: f64,
    /// hxx L163-164: double myTolU / myTolV.
    my_tolu: f64,
    my_tolv: f64,
}

impl CurveConstraint {
    /// OCCT GeomPlate_CurveConstraint() — the empty constructor
    /// (GeomPlate_CurveConstraint.cxx L40-54).
    pub fn new_empty() -> Self {
        CurveConstraint {
            my_frontiere: None,
            my_nb_points: 0,
            my_order: 0,
            my3d_curve: None,
            my_tang: 0,
            my2d_curve: None,
            my_hcurve2d: None,
            my_g0_crit: None,
            my_g1_crit: None,
            my_g2_crit: None,
            my_const_g0: false,
            my_const_g1: false,
            my_const_g2: false,
            my_lprop_surface: None,
            my_tol_dist: 0.0,
            my_tol_ang: 0.0,
            my_tol_curv: 0.0,
            my_tolu: 0.0,
            my_tolv: 0.0,
        }
    }

    /// OCCT GeomPlate_CurveConstraint(Boundary, Tang, NPt = 10,
    /// TolDist = 0.0001, TolAng = 0.01, TolCurv = 0.1)
    /// (GeomPlate_CurveConstraint.cxx L59-115) — the constructor with a
    /// curve on surface.
    pub fn new(
        boundary: CurveBoundary,
        tang: i32,
        npt: i32,
        toldist: f64,
        tolang: f64,
        tolcurv: f64,
    ) -> Self {
        // myLProp(2, TolDist), myTolDist, myTolAng, myTolCurv (init list).
        let mut r = CurveConstraint::new_empty();
        r.my_tol_dist = toldist;
        r.my_tol_ang = tolang;
        r.my_tol_curv = tolcurv;

        r.my_order = tang;
        if (tang < -1) || (tang > 2) {
            panic!("GeomPlate : The continuity is not G0 G1 or G2");
        }
        r.my_nb_points = npt;
        r.my_const_g0 = true;
        r.my_const_g1 = true;
        r.my_const_g2 = true;

        // myFrontiere = down_cast<Adaptor3d_CurveOnSurface>(Boundary).
        match boundary {
            CurveBoundary::OnSurface(frontiere) => {
                r.my_frontiere = Some(frontiere.clone());

                // if (myFrontiere.IsNull()) ... else
                {
                    // handle(GeomAdaptor_Surface) GS1 =
                    //   down_cast<GeomAdaptor_Surface>(myFrontiere->GetSurface());
                    // if (!GS1.IsNull()) Surf = GS1->Surface();
                    // else throw Standard_Failure("... must be GeomAdaptor_Surface");
                    let gs1: &SurfaceHandle = frontiere.get_surface();
                    let surf = match gs1.kernel_surface() {
                        Some(s) => s.clone(),
                        None => panic!(
                            "GeomPlate_CurveConstraint : Surface must be GeomAdaptor_Surface"
                        ),
                    };

                    // myLProp.SetSurface(Surf).
                    r.my_lprop_surface = Some(surf);
                }
            }
            CurveBoundary::Curve(curve) => {
                // myFrontiere.IsNull() -> my3dCurve = Boundary.
                r.my3d_curve = Some(curve);
            }
        }

        // my2dCurve.Nullify(); myHCurve2d.Nullify();
        // myTolU = 0.; myTolV = 0.;
        // myG0Crit/myG1Crit/myG2Crit.Nullify().
        r.my2d_curve = None;
        r.my_hcurve2d = None;
        r.my_tolu = 0.0;
        r.my_tolv = 0.0;
        r.my_g0_crit = None;
        r.my_g1_crit = None;
        r.my_g2_crit = None;
        r
    }

    /// OCCT FirstParameter (GeomPlate_CurveConstraint.cxx L120-134).
    pub fn first_parameter(&self) -> f64 {
        if self.my_hcurve2d.is_some() {
            self.my_hcurve2d.as_ref().unwrap().first_parameter()
        } else if self.my3d_curve.is_none() {
            <CurveOnSurface as Adaptor3dCurve>::first_parameter(self.my_frontiere.as_ref().unwrap())
        } else {
            self.my3d_curve.as_ref().unwrap().first_parameter()
        }
    }

    /// OCCT LastParameter (GeomPlate_CurveConstraint.cxx L139-153).
    pub fn last_parameter(&self) -> f64 {
        if self.my_hcurve2d.is_some() {
            self.my_hcurve2d.as_ref().unwrap().last_parameter()
        } else if self.my3d_curve.is_none() {
            <CurveOnSurface as Adaptor3dCurve>::last_parameter(self.my_frontiere.as_ref().unwrap())
        } else {
            self.my3d_curve.as_ref().unwrap().last_parameter()
        }
    }

    /// OCCT Length (GeomPlate_CurveConstraint.cxx L158-176).
    pub fn length(&self) -> f64 {
        if self.my3d_curve.is_none() {
            // GCPnts_AbscissaPoint::Length(*myFrontiere).
            gcpnts_length_3d(self.my_frontiere.as_ref().unwrap().as_ref())
        } else {
            // GCPnts_AbscissaPoint::Length(*my3dCurve).
            let c = GCPntsCurve3dHandle(self.my3d_curve.as_ref().unwrap().clone());
            gcpnts_length_3d(&c)
        }
    }

    /// OCCT D0 (GeomPlate_CurveConstraint.cxx L181-194).
    pub fn d0(&self, u: f64) -> DVec3 {
        if self.my3d_curve.is_none() {
            let p2d = self.my_frontiere.as_ref().unwrap().get_curve().value(u);
            self.my_frontiere
                .as_ref()
                .unwrap()
                .get_surface()
                .value(p2d.x, p2d.y)
        } else {
            self.my3d_curve.as_ref().unwrap().value(u)
        }
    }

    /// OCCT D1 (GeomPlate_CurveConstraint.cxx L199-209) — returns
    /// (P, V1, V2).
    pub fn d1(&self, u: f64) -> (DVec3, DVec3, DVec3) {
        if self.my3d_curve.is_some() {
            panic!("GeomPlate_CurveConstraint.cxx : Curve must be on a Surface");
        }

        let p2d = self.my_frontiere.as_ref().unwrap().get_curve().value(u);
        let (p, v1, v2) = self
            .my_frontiere
            .as_ref()
            .unwrap()
            .get_surface()
            .d1(p2d.x, p2d.y);
        (p, v1, v2)
    }

    /// OCCT D2 (GeomPlate_CurveConstraint.cxx L214-230) — returns
    /// (P, V1, V2, V3, V4, V5).
    #[allow(clippy::type_complexity)]
    pub fn d2(&self, u: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        if self.my3d_curve.is_some() {
            panic!("GeomPlate_CurveConstraint.cxx : Curve must be on a Surface");
        }

        let p2d = self.my_frontiere.as_ref().unwrap().get_curve().value(u);
        self.my_frontiere
            .as_ref()
            .unwrap()
            .get_surface()
            .d2(p2d.x, p2d.y)
    }

    /// OCCT SetG0Criterion (GeomPlate_CurveConstraint.cxx L235-239).
    pub fn set_g0_criterion(&mut self, g0crit: LawFunctionHandle) {
        self.my_g0_crit = Some(g0crit);
        self.my_const_g0 = false;
    }

    /// OCCT SetG1Criterion (GeomPlate_CurveConstraint.cxx L244-252).
    pub fn set_g1_criterion(&mut self, g1crit: LawFunctionHandle) {
        if self.my3d_curve.is_some() {
            panic!("GeomPlate_CurveConstraint.cxx : Curve must be on a Surface");
        }
        self.my_g1_crit = Some(g1crit);
        self.my_const_g1 = false;
    }

    /// OCCT SetG2Criterion (GeomPlate_CurveConstraint.cxx L257-265).
    pub fn set_g2_criterion(&mut self, g2crit: LawFunctionHandle) {
        if self.my3d_curve.is_some() {
            panic!("GeomPlate_CurveConstraint.cxx : Curve must be on a Surface");
        }
        self.my_g2_crit = Some(g2crit);
        self.my_const_g2 = false;
    }

    /// OCCT G0Criterion (GeomPlate_CurveConstraint.cxx L270-280).
    pub fn g0_criterion(&self, u: f64) -> f64 {
        if self.my_const_g0 {
            self.my_tol_dist
        } else {
            self.my_g0_crit.as_ref().unwrap().borrow_mut().value(u)
        }
    }

    /// OCCT G1Criterion (GeomPlate_CurveConstraint.cxx L285-299).
    pub fn g1_criterion(&self, u: f64) -> f64 {
        if self.my3d_curve.is_some() {
            panic!("GeomPlate_CurveConstraint.cxx : Curve must be on a Surface");
        }
        if self.my_const_g1 {
            self.my_tol_ang
        } else {
            self.my_g1_crit.as_ref().unwrap().borrow_mut().value(u)
        }
    }

    /// OCCT G2Criterion (GeomPlate_CurveConstraint.cxx L304-318).
    pub fn g2_criterion(&self, u: f64) -> f64 {
        if self.my3d_curve.is_some() {
            panic!("GeomPlate_CurveConstraint.cxx : Curve must be on a Surface");
        }
        if self.my_const_g2 {
            self.my_tol_curv
        } else {
            self.my_g2_crit.as_ref().unwrap().borrow_mut().value(u)
        }
    }

    /// OCCT Curve2dOnSurf (GeomPlate_CurveConstraint.cxx L323-346).
    pub fn curve2d_on_surf(&self) -> Option<Curve2d> {
        if self.my2d_curve.is_none() && self.my_hcurve2d.is_some() {
            // GeomAbs_Shape Continuity = GeomAbs_C1;
            // int MaxDegree = 10;
            // int MaxSeg = 20 + myHCurve2d->NbIntervals(GeomAbs_C3);
            // Approx_Curve2d appr(myHCurve2d, first, last, myTolU, myTolV,
            //                     Continuity, MaxDegree, MaxSeg);
            // C2d = appr.Curve();
            //
            // GAP leaf: Approx_Curve2d rides on AdvApprox_ApproxAFunction
            // (untranslated Approx/AdvApprox package).  The OCCT call is
            // preserved at this anchor; until it lands the projected-curve
            // 2d approximation path fails here exactly where OCCT would
            // have run the approximation.
            unimplemented!(
                "Approx_Curve2d (GeomPlate_CurveConstraint::Curve2dOnSurf myHCurve2d branch) is not translated"
            );
        } else {
            self.my2d_curve.clone()
        }
    }

    /// OCCT SetCurve2dOnSurf (GeomPlate_CurveConstraint.cxx L351-354).
    pub fn set_curve2d_on_surf(&mut self, curve: Option<Curve2d>) {
        self.my2d_curve = curve;
    }

    /// OCCT ProjectedCurve (GeomPlate_CurveConstraint.cxx L359-362).
    pub fn projected_curve(&self) -> Option<Curve2dHandle> {
        self.my_hcurve2d.clone()
    }

    /// OCCT SetProjectedCurve (GeomPlate_CurveConstraint.cxx L367-374).
    pub fn set_projected_curve(&mut self, curve: Option<Curve2dHandle>, tolu: f64, tolv: f64) {
        self.my_hcurve2d = curve;
        self.my_tolu = tolu;
        self.my_tolv = tolv;
    }

    /// OCCT Curve3d (GeomPlate_CurveConstraint.cxx L379-389) — the boundary
    /// as a plain 3D adaptor handle (the frontiere coerces to it).
    pub fn curve3d(&self) -> Option<CurveHandle> {
        if self.my3d_curve.is_none() {
            // return occ::handle<Adaptor3d_Curve>(myFrontiere).
            match &self.my_frontiere {
                Some(f) => Some(f.clone() as Arc<dyn Adaptor3dCurve>),
                None => None,
            }
        } else {
            self.my3d_curve.clone()
        }
    }

    /// OCCT NbPoints (GeomPlate_CurveConstraint.cxx L394-397).
    pub fn nb_points(&self) -> i32 {
        self.my_nb_points
    }

    /// OCCT Order (GeomPlate_CurveConstraint.cxx L402-405).
    pub fn order(&self) -> i32 {
        self.my_order
    }

    /// OCCT SetNbPoints (GeomPlate_CurveConstraint.cxx L410-413).
    pub fn set_nb_points(&mut self, new_nb: i32) {
        self.my_nb_points = new_nb;
    }

    /// OCCT SetOrder (GeomPlate_CurveConstraint.cxx L418-421).
    pub fn set_order(&mut self, order: i32) {
        self.my_order = order;
    }

    /// OCCT LPropSurf (GeomPlate_CurveConstraint.cxx L426-435).
    ///
    /// GAP leaf: the returned GeomLProp_SLProps is an untranslated package
    /// (TKGeomBase/GeomLProp).  The OCCT flow is preserved up to the
    /// SetParameters step; the consumers (EcartContraintesMil case 2 /
    /// VerifPoints case 2 through LocalAnalysis_SurfaceContinuity) are
    /// themselves staged, so this anchor keeps the OCCT failure path.
    pub fn lprop_surf(&mut self, u: f64) -> GeomLPropSlProps {
        if self.my_frontiere.is_none() {
            panic!("GeomPlate_CurveConstraint.cxx : Curve must be on a Surface");
        }
        let p2d = self.my_frontiere.as_ref().unwrap().get_curve().value(u);
        GeomLPropSlProps {
            u: p2d.x,
            v: p2d.y,
        }
    }
}

/// GAP carrier for the OCCT GeomLProp_SLProps (TKGeomBase/GeomLProp) — the
/// surface-props object only carries the parameters set by LPropSurf; every
/// derivative evaluation preserves the untranslated-dependency failure path.
#[derive(Debug, Clone, Copy)]
pub struct GeomLPropSlProps {
    pub u: f64,
    pub v: f64,
}

impl GeomLPropSlProps {
    /// OCCT SetParameters — the step LPropSurf performs.
    pub fn set_parameters(&mut self, u: f64, v: f64) {
        self.u = u;
        self.v = v;
    }
}
