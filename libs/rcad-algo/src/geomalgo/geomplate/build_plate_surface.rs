//! OCCT GeomPlate_BuildPlateSurface (TKGeomAlgo/GeomPlate) — 1:1 port.
//!
//! The full class is ported across this file and its continuation module
//! `b` (build_plate_surface_b.rs, the #[path]-declared child module holding
//! the curve-path member bodies; the split keeps both files under the
//! 2000-line limit).  Anchor coverage:
//! - ctors 1/2/3 (L80-217), Init (L388), LoadInitSurface (L401),
//!   Add (L413/L429), SetNbBounds (L418), ProjectPoint (L352),
//! - Perform (L438-746) including the curve branch (L636-698) and the
//!   MakeApprox fallback chain (L535-581),
//! - ComputeSurfInit (L1454-1904) including the half-space selection
//!   (L1474-1648), the metrics comparison (L1746-1802) and the non-linear
//!   fallback (L1804-1903),
//! - LoadPoint (L2522-2583) G0/G1/G2, VerifPoints (L2737-2782) cases 0/1,
//!   ComputeAnisotropie (L2784), IsOrderG1 (L2797), TrierTab (L232),
//! - in `b`: ProjectCurve (L254-303), ProjectedCurve (L307-349),
//!   CourbeJointive (L1349-1445), Intersect (L1913-2145), Discretise
//!   (L2158-2358), CalculNbPtsInit (L2366-2400), LoadCurve (L2407-2516),
//!   VerifSurface (L2588-2732), EcartContraintesMil (L751-866),
//!   Disc2dContour (L871-1019), Disc3dContour (L1024-1163), Curves2d
//!   (L1206-1218), G0/G1/G2Error(Index) (L1260-1322).
//!
//! GAP leaves (untranslated dependencies, OCCT failure path preserved):
//! - EcartContraintesMil case 2 / VerifPoints case 2 need
//!   `LocalAnalysis_SurfaceContinuity` + `GeomLProp_SLProps`
//!   (TKGeomAlgo/LocalAnalysis, TKGeomBase/GeomLProp),
//! - a null-plate MakeApprox construction (Perform L537 with the
//!   null myGeomPlateSurface) panics exactly where OCCT dereferences null.
//!
//! Architecture adaptations:
//! - OCCT `Extrema_ExtPS myProj` -> rcad `closest_point_on_surface`
//!   (stateless; the Initialize sites only compute myTolU/myTolV).
//! - OCCT `Message_ProgressScope` is omitted (rcad Plate::solve_ti takes no
//!   progress range).
//! - `Geom_RectangularTrimmedSurface` -> rcad `Surface3::Trimmed`.
//! - OCCT leaves myTolU/myTolV uninitialized in the ctors; Rust requires an
//!   initializer, so they are zeroed until Perform assigns them.
//! - `NCollection_HArray1<NCollection_Sequence<double>>` (PntInter /
//!   PntG1G1 / myParCont / myPlateCont) -> `Vec<Vec<f64>>`.

#![allow(dead_code)]

use std::sync::Arc;

use glam::{DVec2, DVec3};

use rcad_kernel::base::geom_api::project::closest_point_on_surface;
use rcad_kernel::base::proj_lib::adaptor::{
    Adaptor3dCurve, CurveHandle, CurveOnSurface, SurfaceHandle,
};
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomSurfaceAdaptor;
use rcad_kernel::geom::{Surface3, SurfaceEval, TrimmedSurface};
use rcad_kernel::math::GeomAbsShape;

use super::build_average_plane::{Aij, BuildAveragePlane};
use super::curve_constraint::{CurveBoundary, CurveConstraint};
use super::make_approx::MakeApprox;
use super::point_constraint::PointConstraint;
use super::surface::GeomPlateSurface;
use crate::geomalgo::plate::{
    FreeGtoCConstraint, GtoCConstraint, PinpointConstraint, Plate, PlateD1, PlateD2,
};
use crate::geomalgo::proj_lib_h_comp_projected_curve::CompProjectedCurve;

#[path = "build_plate_surface_b.rs"]
pub(crate) mod b;

pub(crate) use b::vec3_angle;

/// OCCT GeomPlate_BuildPlateSurface.
#[derive(Clone)]
pub struct BuildPlateSurface {
    // OCCT member order (hxx L247-285).
    my_lin_cont: Vec<CurveConstraint>,
    my_par_cont: Option<Vec<Vec<f64>>>,
    my_plate_cont: Option<Vec<Vec<f64>>>,
    my_pnt_cont: Vec<PointConstraint>,
    my_surf_init: Option<Surface3>,
    my_planar_surf_init: Option<Surface3>,
    my_geom_plate_surface: Option<GeomPlateSurface>,
    my_plate: Plate,
    my_prev_plate: Plate,
    my_anisotropie: bool,
    my_sense: Option<Vec<i32>>,
    my_degree: i32,
    my_init_order: Option<Vec<i32>>,
    my_g0_error: f64,
    my_g1_error: f64,
    my_g2_error: f64,
    my_nb_pts_on_cur: i32,
    my_surf_init_is_give: bool,
    my_nb_iter: i32,
    // OCCT Extrema_ExtPS myProj — replaced by closest_point_on_surface.
    my_tol2d: f64,
    my_tol3d: f64,
    my_tolang: f64,
    my_tolu: f64,
    my_tolv: f64,
    my_nb_bounds: i32,
    my_is_linear: bool,
    my_free: bool,
}

/// OCCT static TrierTab (GeomPlate_BuildPlateSurface.cxx L232-250) — reorders
/// the table of transformations to preserve the initial order.
fn trier_tab(tab: &mut Vec<i32>) {
    let nb = tab.len();
    let mut tab_tri = vec![0i32; nb];
    // NCollection_Array1::SetValue(theItem, theIndex): TabTri(Tab(i)) = i.
    for i in 1..=nb {
        tab_tri[tab[i - 1] as usize - 1] = i as i32;
    }
    *tab = tab_tri;
}

/// OCCT GeomAdaptor_Surface::UResolution/VResolution — for a plane the
/// resolution equals the 3d tolerance.  OCCT load() recurses through
/// Geom_RectangularTrimmedSurface to the basis surface (GeomAdaptor_Surface.cxx
/// L423-431); only planes are in the point-path anchor scope.
fn surface_resolution(surf: &Surface3, tol3d: f64) -> f64 {
    match surf {
        Surface3::Plane(_) => tol3d,
        Surface3::Trimmed(trimmed) => surface_resolution(&trimmed.basis, tol3d),
        _ => unimplemented!(
            "GeomAdaptor_Surface UResolution/VResolution is only ported for planes"
        ),
    }
}

/// OCCT gp_Vec2d::Angle (gp_Vec2d.cxx L47-88) — the signed angle in
/// [-PI, PI] with the acos/asin precision switch.
pub(crate) fn vec2d_angle(v: DVec2, the_other: DVec2) -> f64 {
    let a_norm = v.length();
    let an_other_norm = the_other.length();
    if a_norm <= f64::MIN_POSITIVE || an_other_norm <= f64::MIN_POSITIVE {
        panic!("gp_VectorWithNullMagnitude");
    }
    let a_d = a_norm * an_other_norm;
    let a_cosinus = v.dot(the_other) / a_d;
    let a_sinus = (v.x * the_other.y - v.y * the_other.x) / a_d;
    const A_COS_45_DEG: f64 = std::f64::consts::FRAC_1_SQRT_2;
    if a_cosinus > -A_COS_45_DEG && a_cosinus < A_COS_45_DEG {
        // For angles near +/-90 degrees, use acos for better precision.
        if a_sinus > 0.0 {
            a_cosinus.acos()
        } else {
            -a_cosinus.acos()
        }
    } else {
        // For angles near 0 degrees or +/-180 degrees, use asin.
        if a_cosinus > 0.0 {
            a_sinus.asin()
        } else {
            -a_sinus.asin()
        }
    }
}

/// OCCT gp_Vec::IsOpposite (gp_Vec.hxx L130-134):
/// PI - Angle(theOther) <= theAngularTolerance.
pub(crate) fn vec_is_opposite(v: DVec3, the_other: DVec3, ang_tol: f64) -> bool {
    let an_ang = std::f64::consts::PI - vec3_angle(v, the_other);
    an_ang <= ang_tol
}

/// OCCT gp_Vec::IsEqual (gp_Vec.cxx L32-52).
pub(crate) fn vec_is_equal(v: DVec3, the_other: DVec3, lin_tol: f64, ang_tol: f64) -> bool {
    let a_magnitude = v.length();
    let an_other_magnitude = the_other.length();
    let a_val = (a_magnitude - an_other_magnitude).abs();
    if a_magnitude <= lin_tol || an_other_magnitude <= lin_tol {
        a_val <= lin_tol
    } else {
        a_val <= lin_tol && vec3_angle(v, the_other) <= ang_tol
    }
}

impl BuildPlateSurface {
    /// OCCT ctor 3 (GeomPlate_BuildPlateSurface.cxx L186-217) — the 8-parameter
    /// "Constructor with degree".  OCCT default arguments:
    /// Degree=3, NbPtsOnCur=10, NbIter=3, Tol2d=1e-5, Tol3d=1e-4, TolAng=0.01,
    /// TolCurv=0.1 (unused), Anisotropie=false.
    pub fn new(
        degree: i32,
        nb_pts_on_cur: i32,
        nb_iter: i32,
        tol2d: f64,
        tol3d: f64,
        tolang: f64,
        _tolcurv: f64,
        anisotropie: bool,
    ) -> Self {
        if nb_iter < 1 {
            panic!("GeomPlate :  Number of iteration must be >= 1");
        }
        if degree < 2 {
            panic!("GeomPlate : the degree resolution must be upper of 2");
        }
        BuildPlateSurface {
            my_lin_cont: Vec::new(),
            my_par_cont: None,
            my_plate_cont: None,
            my_pnt_cont: Vec::new(),
            my_surf_init: None,
            my_planar_surf_init: None,
            my_geom_plate_surface: None,
            my_plate: Plate::new(),
            my_prev_plate: Plate::new(),
            my_anisotropie: anisotropie,
            my_sense: None,
            my_degree: degree,
            my_init_order: None,
            my_g0_error: 0.0,
            my_g1_error: 0.0,
            my_g2_error: 0.0,
            my_nb_pts_on_cur: nb_pts_on_cur,
            my_surf_init_is_give: false,
            my_nb_iter: nb_iter,
            my_tol2d: tol2d,
            my_tol3d: tol3d,
            my_tolang: tolang,
            my_tolu: 0.0,
            my_tolv: 0.0,
            my_nb_bounds: 0,
            my_is_linear: true,
            my_free: false,
        }
    }

    /// OCCT ctor 2 (L145-183) — "Constructor with initial surface and degree".
    pub fn new_with_surface(
        surf: Surface3,
        degree: i32,
        nb_pts_on_cur: i32,
        nb_iter: i32,
        tol2d: f64,
        tol3d: f64,
        tolang: f64,
        _tolcurv: f64,
        anisotropie: bool,
    ) -> Self {
        if nb_iter < 1 {
            panic!("GeomPlate :  Number of iteration must be >= 1");
        }
        if degree < 2 {
            panic!("GeomPlate : the degree must be above 2");
        }
        BuildPlateSurface {
            my_lin_cont: Vec::new(),
            my_par_cont: None,
            my_plate_cont: None,
            my_pnt_cont: Vec::new(),
            my_surf_init: Some(surf),
            my_planar_surf_init: None,
            my_geom_plate_surface: None,
            my_plate: Plate::new(),
            my_prev_plate: Plate::new(),
            my_anisotropie: anisotropie,
            my_sense: None,
            my_degree: degree,
            my_init_order: None,
            my_g0_error: 0.0,
            my_g1_error: 0.0,
            my_g2_error: 0.0,
            my_nb_pts_on_cur: nb_pts_on_cur,
            my_surf_init_is_give: true,
            my_nb_iter: nb_iter,
            my_tol2d: tol2d,
            my_tol3d: tol3d,
            my_tolang: tolang,
            my_tolu: 0.0,
            my_tolv: 0.0,
            my_nb_bounds: 0,
            my_is_linear: true,
            my_free: false,
        }
    }

    /// OCCT ctor 1 (L80-143) — the curve-constraint constructor ("compatible
    /// with the old version").  `tab_curve` holds the OCCT
    /// `handle<NCollection_HArray1<handle(Adaptor3d_Curve)>>` bounds; each
    /// element carries the OCCT down_cast discrimination
    /// (CurveBoundary::OnSurface / ::Curve).
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_curves(
        npoints: &[i32],
        tab_curve: &[CurveBoundary],
        tang: &[i32],
        degree: i32,
        nb_iter: i32,
        tol2d: f64,
        tol3d: f64,
        tolang: f64,
        _tolcurv: f64,
        anisotropie: bool,
    ) -> Self {
        let mut r = BuildPlateSurface {
            my_lin_cont: Vec::new(),
            my_par_cont: None,
            my_plate_cont: None,
            my_pnt_cont: Vec::new(),
            my_surf_init: None,
            my_planar_surf_init: None,
            my_geom_plate_surface: None,
            my_plate: Plate::new(),
            my_prev_plate: Plate::new(),
            my_anisotropie: anisotropie,
            my_sense: None,
            my_degree: degree,
            my_init_order: None,
            my_g0_error: 0.0,
            my_g1_error: 0.0,
            my_g2_error: 0.0,
            my_nb_pts_on_cur: 0,
            my_surf_init_is_give: false,
            my_nb_iter: nb_iter,
            my_tol2d: tol2d,
            my_tol3d: tol3d,
            my_tolang: tolang,
            my_tolu: 0.0,
            my_tolv: 0.0,
            my_nb_bounds: 0,
            my_is_linear: true,
            my_free: false,
        };

        let ntcurve = tab_curve.len(); // Number of linear constraints
        // myNbPtsOnCur = 0 — different calculation of the number of points
        // depending on the length (L101).
        if r.my_nb_iter < 1 {
            panic!("GeomPlate :  Number of iteration must be >= 1");
        }
        if ntcurve == 0 {
            panic!("GeomPlate : the bounds Array is null");
        }
        if tang.is_empty() {
            panic!("GeomPlate : the constraints Array is null");
        }
        let mut nbp = 0i32;
        for np in npoints.iter().take(ntcurve) {
            nbp += *np;
        }
        if nbp == 0 {
            panic!("GeomPlate : the resolution is impossible if the number of constraints points is 0");
        }
        if r.my_degree < 2 {
            panic!("GeomPlate ; the degree resolution must be upper of 2");
        }
        // Filling fields passing from the old constructor to the new one.
        for i in 1..=ntcurve {
            // new GeomPlate_CurveConstraint(TabCurve(i), Tang(i), NPoints(i))
            // — the ctor-1 defaults TolDist = 0.0001, TolAng = 0.01,
            // TolCurv = 0.1 apply.
            let cont = CurveConstraint::new(
                tab_curve[i - 1].clone(),
                tang[i - 1],
                npoints[i - 1],
                0.0001,
                0.01,
                0.1,
            );
            r.my_lin_cont.push(cont);
        }
        r.my_surf_init_is_give = false;
        r.my_is_linear = true;
        r.my_free = false;
        r
    }

    /// OCCT Init (L388-395) — resets all constraints.
    pub fn init(&mut self) {
        self.my_lin_cont.clear();
        self.my_pnt_cont.clear();
        self.my_pnt_cont = Vec::new();
        self.my_lin_cont = Vec::new();
    }

    /// OCCT LoadInitSurface (L401-406).
    pub fn load_init_surface(&mut self, surf: Surface3) {
        self.my_surf_init = Some(surf);
        self.my_surf_init_is_give = true;
    }

    /// OCCT Add (L413-416) — adds a linear constraint.
    pub fn add_curve_constraint(&mut self, cont: CurveConstraint) {
        self.my_lin_cont.push(cont);
    }

    /// OCCT SetNbBounds (L418-421).
    pub fn set_nb_bounds(&mut self, nb_bounds: i32) {
        self.my_nb_bounds = nb_bounds;
    }

    /// OCCT Add (L429-432) — adds a point constraint.
    pub fn add_point_constraint(&mut self, cont: PointConstraint) {
        self.my_pnt_cont.push(cont);
    }

    /// OCCT Perform (L438-746) — calculates the surface filled with the
    /// loaded constraints.
    pub fn perform(&mut self) {
        // myGeomPlateSurface.Nullify().
        self.my_geom_plate_surface = None;

        if self.my_nb_bounds == 0 {
            self.my_nb_bounds = self.my_lin_cont.len() as i32;
        }

        self.my_plate.init();

        let ntlincont = self.my_lin_cont.len();
        let ntpntcont = self.my_pnt_cont.len();
        let mut nb_boucle: i32 = 0;
        // OCCT keeps `bool Fini = true;` and overwrites it inside the loop.
        #[allow(unused_assignments)]
        let mut fini = true;
        if (ntlincont + ntpntcont) == 0 {
            // OCCT prints a debug warning; then returns with myGeomPlateSurface
            // null (finding #30 semantics).
            return;
        }

        // Initial Surface (L473-504).
        if !self.my_surf_init_is_give {
            self.compute_surf_init();
        } else {
            if ntlincont >= 2 {
                // Table of transformations to preserve the initial order,
                // see TrierTab (L485-498).
                self.my_init_order = Some(vec![0i32; ntlincont]);
                for l in 1..=ntlincont {
                    self.my_init_order.as_mut().unwrap()[l - 1] = l as i32;
                }
                if !self.courbe_jointive(self.my_tol3d) {
                    // throw Standard_Failure("Curves are not joined") — the
                    // OCCT throw is commented out; the warning-only path runs.
                }
                // TrierTab(myInitOrder) — reorder the table of transformations.
                let order = self.my_init_order.as_mut().unwrap();
                trier_tab(order);
            } else if ntlincont > 0 {
                // Patch (L499-503).
                self.my_sense = Some(vec![0i32; ntlincont]);
                self.my_init_order = Some(vec![1i32; ntlincont]);
            }
        }

        if self.my_surf_init.is_none() {
            return;
        }

        // Bounds + GeomAdaptor resolution + myProj.Initialize (L511-516).
        {
            let surf = self.my_surf_init.as_ref().unwrap();
            let d = surf.default_domain();
            let (_u1, _v1, _u2, _v2) = (d[0], d[2], d[1], d[3]);
            self.my_tolu = surface_resolution(surf, self.my_tol3d);
            self.my_tolv = surface_resolution(surf, self.my_tol3d);
            // OCCT myProj.Initialize(aSurfInit, u1, v1, u2, v2, myTolU, myTolV):
            // the rcad projector (closest_point_on_surface) is stateless; the
            // (u1, v1, u2, v2) bounds are kept by the trimmed surface itself.
        }

        // Projection of curves (L518-534).
        let mut ok = true;
        for i in 1..=ntlincont {
            if self.my_lin_cont[i - 1].curve2d_on_surf().is_none() {
                let curve3d = self.my_lin_cont[i - 1].curve3d();
                let curve2d = match &curve3d {
                    Some(c) => self.project_curve(c),
                    None => None,
                };
                if curve2d.is_none() {
                    ok = false;
                    break;
                }
                self.my_lin_cont[i - 1].set_curve2d_on_surf(curve2d);
            }
        }
        if !ok {
            // GeomPlate_MakeApprox fallback chain (L535-573).
            if self.my_geom_plate_surface.is_none() {
                // OCCT constructs GeomPlate_MakeApprox from the null
                // myGeomPlateSurface and dereferences it (RealBounds) — a null
                // handle crash; the anchor preserves that failure path.
                panic!("GeomPlate BuildPlateSurface::Perform — MakeApprox fallback on a null myGeomPlateSurface (OCCT null-handle dereference)");
            }
            // OCCT L537-538: GeomPlate_MakeApprox App(myGeomPlateSurface, ...);
            // mySurfInit = App.Surface(); (the Geom_BSplineSurface upcast to
            // the Geom_Surface handle).
            let app = MakeApprox::new(
                self.my_geom_plate_surface.as_ref().unwrap(),
                self.my_tol3d,
                1,
                3,
                15.0 * self.my_tol3d,
                -1,
                GeomAbsShape::C0,
                1.3,
            );
            self.my_surf_init = app
                .surface()
                .map(|s| Surface3::BSpline(s.clone()));

            {
                let surf = self.my_surf_init.as_ref().unwrap();
                let d = surf.default_domain();
                let (_u1, _v1, _u2, _v2) = (d[0], d[2], d[1], d[3]);
                self.my_tolu = surface_resolution(surf, self.my_tol3d);
                self.my_tolv = surface_resolution(surf, self.my_tol3d);
            }

            ok = true;
            for i in 1..=ntlincont {
                let curve3d = self.my_lin_cont[i - 1].curve3d();
                let curve2d = match &curve3d {
                    Some(c) => self.project_curve(c),
                    None => None,
                };
                if curve2d.is_none() {
                    ok = false;
                    break;
                }
                self.my_lin_cont[i - 1].set_curve2d_on_surf(curve2d);
            }
            if !ok {
                // mySurfInit = myPlanarSurfInit (L559).
                self.my_surf_init = self.my_planar_surf_init.clone();
                {
                    let surf = self.my_surf_init.as_ref().unwrap();
                    self.my_tolu = surface_resolution(surf, self.my_tol3d);
                    self.my_tolv = surface_resolution(surf, self.my_tol3d);
                }
                for i in 1..=ntlincont {
                    let curve3d = self.my_lin_cont[i - 1].curve3d();
                    let curve2d = match &curve3d {
                        Some(c) => self.project_curve(c),
                        None => None,
                    };
                    self.my_lin_cont[i - 1].set_curve2d_on_surf(curve2d);
                }
            } else {
                // Project the points (L572-580).
                for i in 1..=ntpntcont {
                    let p = self.my_pnt_cont[i - 1].d0();
                    let p2d = self.project_point(p);
                    self.my_pnt_cont[i - 1].set_pnt2d_on_surf(p2d);
                }
            }
        }

        // Projection of points (L583-594).
        for i in 1..=ntpntcont {
            if !self.my_pnt_cont[i - 1].has_pnt2d_on_surf() {
                let p = self.my_pnt_cont[i - 1].d0();
                let p2d = self.project_point(p);
                self.my_pnt_cont[i - 1].set_pnt2d_on_surf(p2d);
            }
        }

        // Number of points by curve (L596-602).
        if (ntlincont != 0) && (self.my_nb_pts_on_cur != 0) {
            self.calcul_nb_pts_init();
        }

        // Management of incompatibilities between curves (L604-614).
        let mut pnt_inter: Vec<Vec<f64>> = vec![Vec::new(); ntlincont];
        let mut pnt_g1g1: Vec<Vec<f64>> = vec![Vec::new(); ntlincont];
        if ntlincont != 0 {
            self.intersect(&mut pnt_inter, &mut pnt_g1g1);
        }

        // Loop to obtain a better surface (L616-728).
        self.my_free = !self.my_is_linear;

        loop {
            nb_boucle += 1;
            if ntlincont != 0 {
                // Calculate the total number of points and the maximum of
                // points by curve (L637-647).
                let mut n_point_max = 0i32;
                for i in 1..=ntlincont {
                    if self.my_lin_cont[i - 1].nb_points() > n_point_max {
                        n_point_max = self.my_lin_cont[i - 1].nb_points();
                    }
                }
                let _ = n_point_max;

                // Discretization of curves (L651).
                self.discretise(&pnt_inter, &pnt_g1g1);

                // Preparation of constraint points for plate (L655).
                self.load_curve(nb_boucle, 2);
                if !self.my_pnt_cont.is_empty() {
                    self.load_point(nb_boucle, 2);
                }

                // Construction of the surface (L662-682).
                let anisotropie = self.compute_anisotropie();
                self.my_plate.solve_ti(self.my_degree, anisotropie);

                if !self.my_plate.is_done() {
                    return;
                }

                let plate_surface = GeomPlateSurface::new(
                    self.my_surf_init.clone().unwrap(),
                    self.my_plate.clone(),
                );
                self.my_geom_plate_surface = Some(plate_surface);
                let mut umin = 0.0;
                let mut umax = 0.0;
                let mut vmin = 0.0;
                let mut vmax = 0.0;
                self.my_plate.uv_box(&mut umin, &mut umax, &mut vmin, &mut vmax);
                self.my_geom_plate_surface
                    .as_mut()
                    .unwrap()
                    .set_bounds(umin, umax, vmin, vmax);

                fini = self.verif_surface(nb_boucle);
                if (nb_boucle >= self.my_nb_iter) && (!fini) {
                    fini = true;
                }

                if (ntpntcont != 0) && (fini) {
                    let mut di = 0.0;
                    let mut an = 0.0;
                    let mut cu = 0.0;
                    self.verif_points(&mut di, &mut an, &mut cu);
                    let _ = (di, an, cu);
                }
            } else {
                self.load_point(nb_boucle, 2);
                // Construction of the surface.
                let anisotropie = self.compute_anisotropie();
                self.my_plate.solve_ti(self.my_degree, anisotropie);

                if !self.my_plate.is_done() {
                    return;
                }

                let plate_surface = GeomPlateSurface::new(
                    self.my_surf_init.clone().unwrap(),
                    self.my_plate.clone(),
                );
                self.my_geom_plate_surface = Some(plate_surface);
                let mut umin = 0.0;
                let mut umax = 0.0;
                let mut vmin = 0.0;
                let mut vmax = 0.0;
                self.my_plate.uv_box(&mut umin, &mut umax, &mut vmin, &mut vmax);
                self.my_geom_plate_surface
                    .as_mut()
                    .unwrap()
                    .set_bounds(umin, umax, vmin, vmax);
                fini = true;
                let mut di = 0.0;
                let mut an = 0.0;
                let mut cu = 0.0;
                self.verif_points(&mut di, &mut an, &mut cu);
                let _ = (di, an, cu);
            }
            if fini {
                break;
            }
        } // End loop for better surface
    }

    /// OCCT ProjectPoint (L352-377) — projects a point on the initial surface.
    ///
    /// Architecture adaptation: OCCT iterates `myProj` extrema and keeps the
    /// nearest; rcad `closest_point_on_surface` returns the nearest directly.
    fn project_point(&self, p3d: DVec3) -> DVec2 {
        let proj = closest_point_on_surface(self.my_surf_init.as_ref().unwrap(), p3d, 8);
        DVec2::new(proj.params.0, proj.params.1)
    }

    /// OCCT ComputeSurfInit (L1454-1904) — computes the initial surface when
    /// none was given.
    fn compute_surf_init(&mut self) {
        let mut nopt = 2;
        let popt = 2;
        let mut np = 1usize;
        let mut is_half_space = true;
        let lin_tol = 0.001;
        let ang_tol = 0.001; // AngTol = 0.0001; //LinTol = 0.0001

        let ntlincont = self.my_lin_cont.len();
        let ntpntcont = self.my_pnt_cont.len();

        // Table of transformation to preserve the initial order (L1462-1470).
        if ntlincont != 0 {
            let mut init_order = vec![0i32; ntlincont];
            for i in 1..=ntlincont {
                init_order[i - 1] = i as i32;
            }
            self.my_init_order = Some(init_order);
        }

        let courbe_joint = (ntlincont != 0) && self.courbe_jointive(self.my_tol3d);
        if courbe_joint && self.is_order_g1() {
            nopt = 3;
            // Table contains the cloud of points for calculation of the plane
            // (L1477-1648).
            let nb_point = 20usize;
            let discr = nb_point / 4;
            let mut pnum = 0usize;
            let mut pts: Vec<DVec3> = vec![DVec3::ZERO; (nb_point + 1) * ntlincont + ntpntcont];
            let mut vecs: Vec<DVec3> = Vec::new();
            let mut new_vecs: Vec<DVec3> = Vec::new();
            let mut aset: Vec<Aij> = Vec::new();
            let mut last_vec = DVec3::ZERO;
            for i in 1..=ntlincont {
                let order = self.my_lin_cont[i - 1].order();

                new_vecs.clear();

                let mut uinit = self.my_lin_cont[i - 1].first_parameter();
                let ufinal = self.my_lin_cont[i - 1].last_parameter();
                let mut uif = ufinal - uinit;
                if self.my_sense.as_ref().unwrap()[i - 1] == 1 {
                    uinit = ufinal;
                    uif = -uif;
                }

                let mut to_reverse = false;
                // OCCT declares Vec1/Vec2/Normal outside the j-loop; Normal
                // keeps the last computed value for LastVec.
                let mut normal_last = DVec3::ZERO;
                if i > 1 && order >= 1 {
                    // GeomAbs_G1.
                    let (p, vec1, vec2) = self.my_lin_cont[i - 1].d1(uinit);
                    let _ = p;
                    let normal = vec1.cross(vec2);
                    if vec_is_opposite(last_vec, normal, ang_tol) {
                        to_reverse = true;
                    }
                }

                for j in 0..=nb_point {
                    // Number of points per curve = 20, linear distribution.
                    let inter = j as f64 * uif / (nb_point as f64);
                    if order < 1 || j % discr != 0 {
                        let p = self.my_lin_cont[i - 1].d0(uinit + inter);
                        pts[pnum] = p;
                        pnum += 1;
                    } else {
                        let (_p, vec1, vec2) = self.my_lin_cont[i - 1].d1(uinit + inter);
                        let mut normal = vec1.cross(vec2);
                        normal = normal.normalize();
                        if to_reverse {
                            normal = -normal;
                        }
                        normal_last = normal;
                        let mut is_new = true;
                        let mut k = 1usize;
                        while k <= vecs.len() {
                            if vec_is_equal(vecs[k - 1], normal, lin_tol, ang_tol) {
                                is_new = false;
                                break;
                            }
                            k += 1;
                        }
                        if is_new {
                            let mut k = 1usize;
                            while k <= new_vecs.len() {
                                if vec_is_equal(new_vecs[k - 1], normal, lin_tol, ang_tol) {
                                    is_new = false;
                                    break;
                                }
                                k += 1;
                            }
                        }
                        if is_new {
                            new_vecs.push(normal);
                        }
                    }
                }
                if order >= 1 {
                    is_half_space = BuildAveragePlane::half_space(
                        &new_vecs,
                        &mut vecs,
                        &mut aset,
                        lin_tol,
                        ang_tol,
                    );
                    if !is_half_space {
                        break;
                    }
                    last_vec = normal_last;
                }
            } // for (i = 1; i <= NTLinCont; i++)

            if is_half_space {
                for i in 1..=ntpntcont {
                    let order = self.my_pnt_cont[i - 1].order();

                    new_vecs.clear();
                    if order < 1 {
                        let p = self.my_pnt_cont[i - 1].d0();
                        pts[pnum] = p;
                        pnum += 1;
                    } else {
                        let (_p, vec1, vec2) = self.my_pnt_cont[i - 1].d1();
                        let mut normal = vec1.cross(vec2);
                        normal = normal.normalize();
                        let mut is_new = true;
                        for k in 1..=vecs.len() {
                            if vec_is_equal(vecs[k - 1], normal, lin_tol, ang_tol) {
                                is_new = false;
                                break;
                            }
                        }
                        if is_new {
                            new_vecs.push(normal);
                            is_half_space = BuildAveragePlane::half_space(
                                &new_vecs,
                                &mut vecs,
                                &mut aset,
                                lin_tol,
                                ang_tol,
                            );
                            if !is_half_space {
                                new_vecs[0] = -new_vecs[0];
                                is_half_space = BuildAveragePlane::half_space(
                                    &new_vecs,
                                    &mut vecs,
                                    &mut aset,
                                    lin_tol,
                                    ang_tol,
                                );
                            }
                            if !is_half_space {
                                break;
                            }
                        }
                    }
                } // for (i = 1; i <= NTPntCont; i++)

                if is_half_space {
                    loop {
                        let mut null_exist = false;
                        for i in 1..=vecs.len() {
                            if vecs[i - 1].length_squared() == 0.0 {
                                null_exist = true;
                                vecs.remove(i - 1);
                                break;
                            }
                        }
                        if !null_exist {
                            break;
                        }
                    }
                    // GeomPlate_BuildAveragePlane BAP(Vecs, Pts).
                    let bap = BuildAveragePlane::from_normals(&vecs, pts);
                    let mut u1 = 0.0;
                    let mut u2 = 0.0;
                    let mut v1 = 0.0;
                    let mut v2 = 0.0;
                    bap.min_max_box(&mut u1, &mut u2, &mut v1, &mut v2);
                    // The space is greater for projections.
                    let du = u2 - u1;
                    let dv = v2 - v1;
                    let (u1, u2, v1, v2) = (u1 - du, u2 + du, v1 - dv, v2 + dv);
                    // mySurfInit = new Geom_RectangularTrimmedSurface(BAP.Plane(),
                    //                                                u1, u2, v1, v2);
                    let plane = *bap.plane().unwrap();
                    self.my_surf_init = Some(Surface3::Trimmed(TrimmedSurface::new(
                        Surface3::Plane(plane),
                        u1,
                        u2,
                        v1,
                        v2,
                    )));
                }
            } // if (isHalfSpace)
            if !is_half_space {
                self.my_is_linear = false;
                nopt = 2;
            }
        } // if (NTLinCont != 0 && (CourbeJoint = CourbeJointive( myTol3d )) && IsOrderG1())

        if ntlincont != 0 {
            // TrierTab(myInitOrder) — reorder the table of transformations.
            self.trier_tab_my_init_order();
        }

        if nopt != 3 {
            let mut nopt = nopt;
            if ntpntcont != 0 {
                nopt = 1; // Calculate by the method of plane of inertia
            } else if !courbe_joint || ntlincont != self.my_nb_bounds as usize {
                // throw Standard_Failure("Curves are not joined") — the OCCT
                // throw is commented out; the inertia path runs.
                nopt = 1;
            }

            // Curve-length bookkeeping (L1671-1689).
            let mut len_t = 0.0;
            let mut npt = 0;
            let nt_point = 20 * ntlincont;
            for i in 1..=ntlincont {
                len_t += self.my_lin_cont[i - 1].length();
            }
            for i in 1..=ntlincont {
                let nb_point =
                    (nt_point as f64 * (self.my_lin_cont[i - 1].length()) / len_t) as i32;
                let nb_point = if nb_point < 10 { 10 } else { nb_point };
                npt += nb_point;
            }
            let _ = npt;

            // Table containing a cloud of points for the plane (L1691-1719).
            let mut pts: Vec<DVec3> = vec![DVec3::ZERO; 20 * ntlincont + ntpntcont];
            let nb_point = 20usize;
            for i in 1..=ntlincont {
                let mut uinit = self.my_lin_cont[i - 1].first_parameter();
                let ufinal = self.my_lin_cont[i - 1].last_parameter();
                let mut uif = ufinal - uinit;
                if self.my_sense.as_ref().unwrap()[i - 1] == 1 {
                    uinit = ufinal;
                    uif = -uif;
                }
                for j in 0..nb_point {
                    // Number of points per curve = 20, linear distribution.
                    let inter = j as f64 * uif / (nb_point as f64);
                    let p = self.my_lin_cont[i - 1].d0(uinit + inter);
                    pts[np - 1] = p;
                    np += 1;
                }
            }
            for i in 1..=ntpntcont {
                let p = self.my_pnt_cont[i - 1].d0();
                pts[np - 1] = p;
                np += 1;
            }
            if !courbe_joint {
                self.my_nb_bounds = 0;
            }
            let bap = BuildAveragePlane::new(
                pts,
                nb_point * self.my_nb_bounds as usize,
                self.my_tol3d / 1000.0,
                popt,
                nopt,
            );
            if !bap.is_plane() {
                return;
            }
            let mut u1 = 0.0;
            let mut u2 = 0.0;
            let mut v1 = 0.0;
            let mut v2 = 0.0;
            bap.min_max_box(&mut u1, &mut u2, &mut v1, &mut v2);
            // The space is greater for projections.
            let du = u2 - u1;
            let dv = v2 - v1;
            let (u1, u2, v1, v2) = (u1 - du, u2 + du, v1 - dv, v2 + dv);
            // mySurfInit = new Geom_RectangularTrimmedSurface(BAP.Plane(),
            //                                                u1, u2, v1, v2);
            let plane = *bap.plane().unwrap();
            self.my_surf_init = Some(Surface3::Trimmed(TrimmedSurface::new(
                Surface3::Plane(plane),
                u1,
                u2,
                v1,
                v2,
            )));
        } // if (nopt != 3)

        // Comparing metrics of curves and projected curves (L1746-1802).
        if ntlincont != 0 && self.my_is_linear {
            // InitPlane = down_cast<Geom_RectangularTrimmedSurface>
            //                (mySurfInit)->BasisSurface().
            let init_plane = match self.my_surf_init.as_ref().unwrap() {
                Surface3::Trimmed(trimmed) => (*trimmed.basis).clone(),
                _ => panic!("GeomPlate ComputeSurfInit: down_cast<Geom_RectangularTrimmedSurface>(mySurfInit)"),
            };
            let mut ratio = 0.0f64;
            let r1 = 2.0f64; // R1 = 3, R2 = 0.5; R1 = 1.4, R2 = 0.8; R1 = 5., R2 = 0.2;
            let r2 = 0.6f64;
            let hsur: SurfaceHandle = Arc::new(GeomSurfaceAdaptor::new(init_plane));
            let nb_point = 20usize;

            for i in 1..=ntlincont {
                if !self.my_is_linear {
                    break;
                }
                let first_par = self.my_lin_cont[i - 1].first_parameter();
                let last_par = self.my_lin_cont[i - 1].last_parameter();
                let uif = (last_par - first_par) / (nb_point as f64);

                // handle(Adaptor3d_Curve) Curve = myLinCont->Value(i)->Curve3d().
                let curve = self.my_lin_cont[i - 1].curve3d().expect("Curve3d");
                // occ::handle<ProjLib_HCompProjectedCurve> ProjCurve =
                //     new ProjLib_HCompProjectedCurve(hsur, Curve, myTol3d, myTol3d);
                let proj_curve = CompProjectedCurve::new(
                    hsur.clone(),
                    curve.clone(),
                    self.my_tol3d,
                    self.my_tol3d,
                );
                // Adaptor3d_CurveOnSurface AProj(ProjCurve, hsur).
                let a_proj = CurveOnSurface::new(Arc::new(proj_curve), hsur.clone());

                for j in 1..nb_point {
                    // OCCT L1778: the inner loop guard is
                    // `j < NbPoint && myIsLinear`.
                    if !self.my_is_linear {
                        break;
                    }
                    let inter = first_par + (j as f64) * uif;
                    let (_p, der_c) = curve.d1(inter);
                    let (_p, der_cproj) = a_proj.d1(inter);

                    let a1 = der_c.length();
                    let a2 = der_cproj.length();
                    if a2 <= 1.0e-20 {
                        ratio = 1.0e20;
                    } else {
                        ratio = a1 / a2;
                    }
                    if ratio > r1 || ratio < r2 {
                        self.my_is_linear = false;
                        break;
                    }
                } // for (int j = 1; j < NbPoint && myIsLinear; j++)
            }
            let _ = ratio;
        } // comparing metrics of curves and projected curves

        if !self.my_is_linear {
            // Free-form fallback (L1804-1903).
            self.my_planar_surf_init = self.my_surf_init.clone();
            {
                let surf = self.my_surf_init.as_ref().unwrap();
                self.my_tolu = surface_resolution(surf, self.my_tol3d);
                self.my_tolv = surface_resolution(surf, self.my_tol3d);
            }

            // Projection of curves (L1818-1824).
            for i in 1..=ntlincont {
                if self.my_lin_cont[i - 1].curve2d_on_surf().is_none() {
                    let curve3d = self.my_lin_cont[i - 1].curve3d();
                    let curve2d = match &curve3d {
                        Some(c) => self.project_curve(c),
                        None => None,
                    };
                    self.my_lin_cont[i - 1].set_curve2d_on_surf(curve2d);
                }
            }

            // Projection of points (L1829-1837).
            for i in 1..=ntpntcont {
                let p = self.my_pnt_cont[i - 1].d0();
                if !self.my_pnt_cont[i - 1].has_pnt2d_on_surf() {
                    let p2d = self.project_point(p);
                    self.my_pnt_cont[i - 1].set_pnt2d_on_surf(p2d);
                }
            }

            // Number of points by curve (L1842-1845).
            if (ntlincont != 0) && (self.my_nb_pts_on_cur != 0) {
                self.calcul_nb_pts_init();
            }

            // Management of incompatibilities between curves (L1847-1857).
            let mut pnt_inter: Vec<Vec<f64>> = vec![Vec::new(); ntlincont];
            let mut pnt_g1g1: Vec<Vec<f64>> = vec![Vec::new(); ntlincont];
            if ntlincont != 0 {
                self.intersect(&mut pnt_inter, &mut pnt_g1g1);
            }

            // Discretization of curves (L1862).
            self.discretise(&pnt_inter, &pnt_g1g1);

            // Preparation of points of constraint for plate (L1866-1870).
            self.load_curve(0, 0);
            if !self.my_pnt_cont.is_empty() {
                self.load_point(0, 0);
            }

            // Construction of the surface (L1874-1887).
            let anisotropie = self.compute_anisotropie();
            self.my_plate.solve_ti(2, anisotropie);

            if !self.my_plate.is_done() {
                return;
            }

            self.my_geom_plate_surface = Some(GeomPlateSurface::new(
                self.my_surf_init.clone().unwrap(),
                self.my_plate.clone(),
            ));

            // GeomPlate_MakeApprox App(myGeomPlateSurface, myTol3d, 1, 3,
            //     15 * myTol3d, -1, GeomAbs_C0);  mySurfInit = App.Surface().
            let app = MakeApprox::new(
                self.my_geom_plate_surface.as_ref().unwrap(),
                self.my_tol3d,
                1,
                3,
                15.0 * self.my_tol3d,
                -1,
                GeomAbsShape::C0,
                1.3,
            );
            // OCCT: mySurfInit = App.Surface(); (the Geom_BSplineSurface
            // upcast to the Geom_Surface handle).
            self.my_surf_init = app.surface().map(|s| Surface3::BSpline(s.clone()));

            self.my_surf_init_is_give = true;
            self.my_plate.init(); // Reset

            for i in 1..=ntlincont {
                self.my_lin_cont[i - 1].set_curve2d_on_surf(None);
            }
        }
    }

    /// TrierTab applied to myInitOrder (OCCT L1653-1656).
    fn trier_tab_my_init_order(&mut self) {
        if let Some(order) = self.my_init_order.as_mut() {
            trier_tab(order);
        }
    }

    /// OCCT IsOrderG1 (L2797-2809).
    pub fn is_order_g1(&self) -> bool {
        let mut result = true;
        for i in 1..=self.my_lin_cont.len() {
            if self.my_lin_cont[i - 1].order() < 1 {
                result = false;
                break;
            }
        }
        result
    }

    /// OCCT LoadPoint (L2522-2583) — loading of the point constraints.
    fn load_point(&mut self, _nb_boucle: i32, order_max: i32) {
        let ntpntcont = self.my_pnt_cont.len();
        // Loading of points of point constraints.
        for i in 1..=ntpntcont {
            let p3d = self.my_pnt_cont[i - 1].d0();
            let p2d = self.my_pnt_cont[i - 1].pnt2d_on_surf();
            let pp = self.my_surf_init.as_ref().unwrap().point_at(p2d.x, p2d.y);
            let pdif = DVec3::new(-pp.x + p3d.x, -pp.y + p3d.y, -pp.z + p3d.z);
            let pc = PinpointConstraint::new(p2d, pdif, 0, 0);
            self.my_plate.load_pinpoint(pc);
            let tang = self.my_pnt_cont[i - 1].order().min(order_max);
            if tang == 1 {
                // ==1 (L2542-2559).
                let (pi, v1, v2) = self.my_pnt_cont[i - 1].d1();
                let _ = pi;
                let (pp, v3, v4) = self
                    .my_surf_init
                    .as_ref()
                    .unwrap()
                    .derivatives(p2d.x, p2d.y);
                let d1final = PlateD1::new(v1, v2);
                let d1init = PlateD1::new(v3, v4);
                if !self.my_free {
                    let gcc = GtoCConstraint::new(p2d, &d1init, &d1final);
                    self.my_plate.load_gto_c(&gcc);
                } else {
                    let free_gcc = FreeGtoCConstraint::new(p2d, &d1init, &d1final, 1.0, 0);
                    self.my_plate.load_free_gto_c(&free_gcc);
                }
            }
            // Loading of points G2 (L2560-2581).
            if tang == 2 {
                // ==2
                let (_pi, v1, v2, v5, v6, v7) = self.my_pnt_cont[i - 1].d2();
                let (_pp, v3, v4, v8, v9, v10) = self
                    .my_surf_init
                    .as_ref()
                    .unwrap()
                    .derivatives2(p2d.x, p2d.y);
                let d1final = PlateD1::new(v1, v2);
                let d1init = PlateD1::new(v3, v4);
                let d2final = PlateD2::new(v5, v6, v7);
                let d2init = PlateD2::new(v8, v9, v10);
                let gcc = GtoCConstraint::new_g2(p2d, &d1init, &d1final, &d2init, &d2final);
                self.my_plate.load_gto_c(&gcc);
            }
        }
    }

    /// OCCT VerifPoints (L2737-2782) — the results are given through the
    /// output parameters (Dist, Ang, Curv).  Cases 0/1 are complete; case 2
    /// is a GAP leaf (LocalAnalysis_SurfaceContinuity).
    pub fn verif_points(&self, dist: &mut f64, ang: &mut f64, curv: &mut f64) {
        let ntpntcont = self.my_pnt_cont.len();
        *ang = 0.0;
        *dist = 0.0;
        *curv = 0.0;
        for i in 1..=ntpntcont {
            let pnt_cont = &self.my_pnt_cont[i - 1];
            match pnt_cont.order() {
                0 => {
                    // case 0 (L2751-2756).
                    let p2d = pnt_cont.pnt2d_on_surf();
                    let pi = pnt_cont.d0();
                    let pf = self
                        .my_geom_plate_surface
                        .as_ref()
                        .unwrap()
                        .eval_d0(p2d.x, p2d.y);
                    *dist = (pf - pi).length();
                }
                1 => {
                    // case 1 (L2757-2769).
                    let (pi, v1i, v2i) = pnt_cont.d1();
                    let p2d = pnt_cont.pnt2d_on_surf();
                    let (pf, v1f, v2f) = self
                        .my_geom_plate_surface
                        .as_ref()
                        .unwrap()
                        .eval_d1(p2d.x, p2d.y);
                    *dist = (pf - pi).length();
                    let v3i = v1i.cross(v2i);
                    let v3f = v1f.cross(v2f);
                    *ang = vec3_angle(v3f, v3i);
                    if *ang > (std::f64::consts::PI / 2.0) {
                        *ang = std::f64::consts::PI - *ang;
                    }
                }
                2 => {
                    // case 2 (L2770-2779): GeomLProp_SLProps +
                    // LocalAnalysis_SurfaceContinuity — untranslated
                    // dependencies; the anchor preserves the failure path.
                    unimplemented!(
                        "VerifPoints case 2 needs GeomLProp_SLProps + LocalAnalysis_SurfaceContinuity (untranslated)"
                    );
                }
                _ => {}
            }
        }
    }

    /// OCCT ComputeAnisotropie (L2784-2795).
    pub fn compute_anisotropie(&self) -> f64 {
        if self.my_anisotropie {
            // Temporary
            1.0
        } else {
            1.0
        }
    }

    /// OCCT IsDone (L1168-1171).
    pub fn is_done(&self) -> bool {
        self.my_plate.is_done()
    }

    /// OCCT Surface (L1176-1179) — the computation result (null when no
    /// result is available).
    pub fn surface(&self) -> Option<&GeomPlateSurface> {
        self.my_geom_plate_surface.as_ref()
    }

    /// OCCT SurfInit (L1184-1187).
    pub fn surf_init(&self) -> Option<&Surface3> {
        self.my_surf_init.as_ref()
    }

    /// OCCT Sense (L1192-1201).
    pub fn sense(&self) -> Vec<i32> {
        let ntcurve = self.my_lin_cont.len();
        let mut sens = vec![0i32; ntcurve];
        let sense_ref = self.my_sense.as_ref().unwrap();
        let order_ref = self.my_init_order.as_ref().unwrap();
        for i in 1..=ntcurve {
            sens[i - 1] = sense_ref[(order_ref[i - 1] - 1) as usize];
        }
        sens
    }

    /// OCCT Curves2d (L1206-1218).
    pub fn curves2d(&self) -> Vec<Option<rcad_kernel::geom::Curve2d>> {
        let ntcurve = self.my_lin_cont.len();
        let mut c2dfin: Vec<Option<rcad_kernel::geom::Curve2d>> = Vec::new();
        for i in 1..=ntcurve {
            c2dfin.push(self.my_lin_cont[i - 1].curve2d_on_surf());
        }
        c2dfin
    }

    /// OCCT Order (L1223-1231).
    pub fn order(&self) -> Vec<i32> {
        let mut result = vec![0i32; self.my_lin_cont.len()];
        let order_ref = self.my_init_order.as_ref().unwrap();
        for i in 1..=self.my_lin_cont.len() {
            result[(order_ref[i - 1] - 1) as usize] = i as i32;
        }
        result
    }

    /// OCCT G0Error() (L1237-1240).
    pub fn g0_error(&self) -> f64 {
        self.my_g0_error
    }

    /// OCCT G1Error() (L1245-1248).
    pub fn g1_error(&self) -> f64 {
        self.my_g1_error
    }

    /// OCCT G2Error() (L1253-1256).
    pub fn g2_error(&self) -> f64 {
        self.my_g2_error
    }

    /// OCCT CurveConstraint(order) (L1324-1327).
    pub fn curve_constraint(&self, order: usize) -> &CurveConstraint {
        &self.my_lin_cont[order - 1]
    }

    /// OCCT PointConstraint(order) (L1330-1333).
    pub fn point_constraint(&self, order: usize) -> &PointConstraint {
        &self.my_pnt_cont[order - 1]
    }
}

impl Default for BuildPlateSurface {
    /// OCCT ctor 3 default arguments (Degree=3, NbPtsOnCur=10, NbIter=3,
    /// Tol2d=1e-5, Tol3d=1e-4, TolAng=0.01, TolCurv=0.1, Anisotropie=false).
    fn default() -> Self {
        Self::new(3, 10, 3, 1.0e-5, 1.0e-4, 0.01, 0.1, false)
    }
}
