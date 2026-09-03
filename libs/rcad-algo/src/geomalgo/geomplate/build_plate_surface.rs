//! OCCT GeomPlate_BuildPlateSurface (TKGeomAlgo/GeomPlate) — 1:1 port of the
//! point-constraint path (anchor: GeomPlate_BuildPlateSurface_Test.cxx).
//!
//! Fully ported: ctor 2/3 (L145-217), Init (L388), LoadInitSurface (L401),
//! Add (L413/L429), SetNbBounds (L418), ProjectPoint (L352), Perform
//! (L438-746) point branch, ComputeSurfInit (L1454-1911) point branch,
//! LoadPoint (L2522-2583) G0 branch, VerifPoints (L2737-2782) G0/G1-less
//! case 0, ComputeAnisotropie (L2784), IsOrderG1 (L2797), TrierTab (L232),
//! and the scalar accessors.
//!
//! Curve-constraint machinery (ProjectCurve/ProjectedCurve/CourbeJointive/
//! Intersect/Discretise/LoadCurve/CalculNbPtsInit/VerifSurface/
//! EcartContraintesMil/Disc2dContour/Disc3dContour and the G2 verification)
//! is anchor-out-of-scope and follows the ThruSections precedent:
//! `unimplemented!()` skeletons, backfilled later.
//!
//! Architecture adaptations:
//! - OCCT `Extrema_ExtPS myProj` -> rcad `closest_point_on_surface`
//!   (stateless; the Initialize sites only compute myTolU/myTolV).
//! - OCCT `Message_ProgressScope` is omitted (rcad Plate::solve_ti takes no
//!   progress range).
//! - `Geom_RectangularTrimmedSurface` -> rcad `Surface3::Trimmed`.
//! - OCCT leaves myTolU/myTolV uninitialized in the ctors; Rust requires an
//!   initializer, so they are zeroed until Perform assigns them.

// Curve-path members and unported skeletons are declared but never called on
// the point-path anchor; see the ThruSections precedent.
#![allow(dead_code)]

use glam::{DVec2, DVec3};

use rcad_kernel::base::geom_api::project::closest_point_on_surface;
use rcad_kernel::geom::{Surface3, SurfaceEval, TrimmedSurface};

use super::build_average_plane::BuildAveragePlane;
use super::curve_constraint::CurveConstraint;
use super::point_constraint::PointConstraint;
use super::surface::GeomPlateSurface;
use crate::geomalgo::plate::{PinpointConstraint, Plate};

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
/// resolution equals the 3d tolerance.  Only the plane is in the point-path
/// anchor scope.
fn surface_resolution(surf: &Surface3, tol3d: f64) -> f64 {
    match surf {
        Surface3::Plane(_) => tol3d,
        _ => unimplemented!(
            "GeomAdaptor_Surface UResolution/VResolution is only ported for planes"
        ),
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

    /// OCCT ctor 1 (L69-143) — curve-constraint constructor, not ported
    /// (curve path, anchor-out-of-scope).
    pub fn new_with_curves(
        _npoints: &[i32],
        _tab_curve: &[CurveConstraint],
        _tang: &[i32],
        _degree: i32,
        _nb_iter: i32,
        _tol2d: f64,
        _tol3d: f64,
        _tolang: f64,
        _tolcurv: f64,
        _anisotropie: bool,
    ) -> Self {
        unimplemented!(
            "GeomPlate_BuildPlateSurface curve constructor is not ported (curve path)"
        );
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

        // Initial Surface (L477-501).
        if !self.my_surf_init_is_give {
            self.compute_surf_init();
        } else if ntlincont >= 2 {
            // Table of transformations to preserve the initial order
            // (L483-500) — curve path, anchor-out-of-scope.
            unimplemented!("Perform CourbeJointive/TrierTab branch is not ported (curve path)");
        } else if ntlincont > 0 {
            // Patch (L501-505).
            self.my_sense = Some(vec![0i32; ntlincont]);
            self.my_init_order = Some(vec![1i32; ntlincont]);
        }

        if self.my_surf_init.is_none() {
            return;
        }

        // Bounds + GeomAdaptor resolution + myProj.Initialize (L507-514).
        let (u1, v1, u2, v2) = {
            let surf = self.my_surf_init.as_ref().unwrap();
            let d = surf.default_domain();
            (d[0], d[2], d[1], d[3])
        };
        {
            let surf = self.my_surf_init.as_ref().unwrap();
            self.my_tolu = surface_resolution(surf, self.my_tol3d);
            self.my_tolv = surface_resolution(surf, self.my_tol3d);
        }
        // OCCT myProj.Initialize(aSurfInit, u1, v1, u2, v2, myTolU, myTolV):
        // the rcad projector (closest_point_on_surface) is stateless; the
        // (u1, v1, u2, v2) bounds are kept by the trimmed surface itself.
        let _ = (u1, v1, u2, v2);

        // Projection of curves (L518-535) — the loop body never runs without
        // curve constraints; ProjectCurve is an unported skeleton.
        #[allow(unused_mut)]
        let mut ok = true;
        for _i in 1..=ntlincont {
            unimplemented!("Perform curve projection loop is not ported (curve path)");
        }
        if !ok {
            // GeomPlate_MakeApprox fallback chain (L537-573) — curve path.
            unimplemented!("Perform MakeApprox fallback is not ported (curve path)");
        }

        // Projection of points (L598-606).
        for i in 1..=ntpntcont {
            if !self.my_pnt_cont[i - 1].has_pnt2d_on_surf() {
                let p = self.my_pnt_cont[i - 1].d0();
                let p2d = self.project_point(p);
                self.my_pnt_cont[i - 1].set_pnt2d_on_surf(p2d);
            }
        }

        // Number of points by curve (L611-614).
        if (ntlincont != 0) && (self.my_nb_pts_on_cur != 0) {
            self.calcul_nb_pts_init();
        }

        // Management of incompatibilities between curves (L618-628).
        if ntlincont != 0 {
            self.intersect();
        }

        // Loop to obtain a better surface (L631-744).
        self.my_free = !self.my_is_linear;

        loop {
            nb_boucle += 1;
            if ntlincont != 0 {
                // Curve branch (L655-722): NPointMax / Discretise / LoadCurve
                // / LoadPoint / SolveTI / VerifSurface — anchor-out-of-scope.
                unimplemented!("Perform curve branch is not ported (curve path)");
            } else {
                self.load_point(nb_boucle, 2);
                // Construction of the surface.
                let anisotropie = self.compute_anisotropie();
                self.my_plate.solve_ti(self.my_degree, anisotropie);

                if !self.my_plate.is_done() {
                    return;
                }

                let plate_surface =
                    GeomPlateSurface::new(self.my_surf_init.clone().unwrap(), self.my_plate.clone());
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

    /// OCCT ComputeSurfInit (L1454-1911) — computes the initial surface when
    /// none was given.
    fn compute_surf_init(&mut self) {
        let mut nopt = 2;
        let popt = 2;
        let mut np = 1usize;
        // OCCT mutates isHalfSpace inside the (curve-path) nopt==3 branch
        // only; kept for form.
        let mut is_half_space = true;
        let lin_tol = 0.001;
        let ang_tol = 0.001;
        let _ = (lin_tol, ang_tol);
        let _ = &mut is_half_space;

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
            // nopt = 3 — discretized curve normals + HalfSpace selection
            // (L1474-1648) — curve path, anchor-out-of-scope.
            unimplemented!("ComputeSurfInit half-space branch is not ported (curve path)");
        } // if (NTLinCont != 0 && CourbeJointive && IsOrderG1())

        if ntlincont != 0 {
            self.trier_tab_my_init_order();
        }

        if nopt != 3 {
            if ntpntcont != 0 {
                nopt = 1; // Calculate by the method of plane of inertia
            } else if !courbe_joint || ntlincont != self.my_nb_bounds as usize {
                nopt = 1;
            }

            // Curve-length bookkeeping (L1668-1682) — empty without curves.
            let mut len_t = 0.0;
            let mut npt = 0;
            let nt_point = 20 * ntlincont;
            for i in 1..=ntlincont {
                len_t += self.my_lin_cont[i - 1].length();
            }
            for i in 1..=ntlincont {
                let nb_point = (nt_point as f64 * (self.my_lin_cont[i - 1].length()) / len_t) as i32;
                let nb_point = if nb_point < 10 { 10 } else { nb_point };
                npt += nb_point;
            }
            let _ = npt;

            // Table containing a cloud of points for the plane (L1718-1741).
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
            u1 -= du;
            u2 += du;
            v1 -= dv;
            v2 += dv;
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

        // Comparing metrics of curves and projected curves (L1786-1849) —
        // curve path, anchor-out-of-scope.
        if ntlincont != 0 && self.my_is_linear {
            unimplemented!("ComputeSurfInit metrics comparison is not ported (curve path)");
        }

        if !self.my_is_linear {
            // Free-form fallback (L1851-1910): projections, discretisation,
            // Plate solve and MakeApprox — curve path, anchor-out-of-scope.
            unimplemented!("ComputeSurfInit non-linear fallback is not ported (curve path)");
        }
    }

    /// TrierTab applied to myInitOrder (OCCT L1653-1656).
    fn trier_tab_my_init_order(&mut self) {
        if let Some(order) = self.my_init_order.as_mut() {
            trier_tab(order);
        }
    }

    /// OCCT IsOrderG1 (L2797-2807).
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

    /// OCCT CourbeJointive (L1349-1451) — curve path, anchor-out-of-scope.
    fn courbe_jointive(&self, _tolerance: f64) -> bool {
        unimplemented!("CourbeJointive is not ported (curve path)");
    }

    /// OCCT CalculNbPtsInit (L2366-2403) — curve path, anchor-out-of-scope.
    fn calcul_nb_pts_init(&mut self) {
        unimplemented!("CalculNbPtsInit is not ported (curve path)");
    }

    /// OCCT Intersect (L1913-2154) — curve path, anchor-out-of-scope.
    fn intersect(&mut self) {
        unimplemented!("Intersect is not ported (curve path)");
    }

    /// OCCT LoadPoint (L2522-2583) — loading of the point constraints.
    fn load_point(&mut self, _nb_boucle: i32, order_max: i32) {
        let ntpntcont = self.my_pnt_cont.len();
        // Loading of points of point constraints.
        for i in 1..=ntpntcont {
            let p3d = self.my_pnt_cont[i - 1].d0();
            let p2d = self.my_pnt_cont[i - 1].pnt2d_on_surf();
            let pp = self
                .my_surf_init
                .as_ref()
                .unwrap()
                .point_at(p2d.x, p2d.y);
            let pdif = p3d - pp;
            let pc = PinpointConstraint::new(p2d, pdif, 0, 0);
            self.my_plate.load_pinpoint(pc);
            let tang = self.my_pnt_cont[i - 1].order().min(order_max);
            if tang == 1 {
                // OCCT builds Plate_GtoCConstraint / FreeGtoCConstraint from
                // PointConstraint::D1 — not ported (needs GeomLProp_SLProps).
                unimplemented!("LoadPoint G1 branch is not ported (PointConstraint::D1)");
            }
            if tang == 2 {
                unimplemented!("LoadPoint G2 branch is not ported (PointConstraint::D2)");
            }
        }
    }

    /// OCCT VerifPoints (L2737-2782) — the results are given through the
    /// output parameters (Dist, Ang, Curv).  Only the Order=0 case is fully
    /// ported; Order=1/2 need PointConstraint::D1/D2 (anchor-out-of-scope).
    pub fn verif_points(&self, dist: &mut f64, ang: &mut f64, curv: &mut f64) {
        let ntpntcont = self.my_pnt_cont.len();
        *ang = 0.0;
        *dist = 0.0;
        *curv = 0.0;
        for i in 1..=ntpntcont {
            let pnt_cont = &self.my_pnt_cont[i - 1];
            match pnt_cont.order() {
                0 => {
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
                    unimplemented!("VerifPoints case 1 needs PointConstraint::D1 (not ported)");
                }
                2 => {
                    unimplemented!(
                        "VerifPoints case 2 needs LocalAnalysis_SurfaceContinuity (not ported)"
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

    /// OCCT Curves2d (L1206-1218) — curve path, anchor-out-of-scope.
    pub fn curves2d(&self) {
        unimplemented!("Curves2d is not ported (curve path)");
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

    /// OCCT G0Error(Index) (L1260-1277) — curve path, anchor-out-of-scope.
    pub fn g0_error_index(&mut self, _index: i32) -> f64 {
        unimplemented!("G0Error(Index) is not ported");
    }

    /// OCCT G1Error(Index) (L1282-1299) — curve path, anchor-out-of-scope.
    pub fn g1_error_index(&mut self, _index: i32) -> f64 {
        unimplemented!("G1Error(Index) is not ported");
    }

    /// OCCT G2Error(Index) (L1304-1321) — curve path, anchor-out-of-scope.
    pub fn g2_error_index(&mut self, _index: i32) -> f64 {
        unimplemented!("G2Error(Index) is not ported");
    }

    /// OCCT CurveConstraint(order) (L1324-1327).
    pub fn curve_constraint(&self, order: usize) -> &CurveConstraint {
        &self.my_lin_cont[order - 1]
    }

    /// OCCT PointConstraint(order) (L1330-1333).
    pub fn point_constraint(&self, order: usize) -> &PointConstraint {
        &self.my_pnt_cont[order - 1]
    }

    /// OCCT EcartContraintesMil (L751-...) — curve path, anchor-out-of-scope.
    fn ecart_contraintes_mil(&mut self) {
        unimplemented!("EcartContraintesMil is not ported (curve path)");
    }

    /// OCCT Disc2dContour (L871-...) — curve path, anchor-out-of-scope.
    pub fn disc2d_contour(&self) {
        unimplemented!("Disc2dContour is not ported (curve path)");
    }

    /// OCCT Disc3dContour (L1024-...) — curve path, anchor-out-of-scope.
    pub fn disc3d_contour(&self) {
        unimplemented!("Disc3dContour is not ported (curve path)");
    }

    /// OCCT VerifSurface (L2588-2733) — curve path, anchor-out-of-scope.
    fn verif_surface(&mut self, _nb_boucle: i32) -> bool {
        unimplemented!("VerifSurface is not ported (curve path)");
    }

    /// OCCT ProjectCurve (L254-303) — curve path, anchor-out-of-scope.
    fn project_curve(&self) {
        unimplemented!("ProjectCurve is not ported (curve path)");
    }

    /// OCCT ProjectedCurve (L307-349) — curve path, anchor-out-of-scope.
    fn projected_curve(&self) {
        unimplemented!("ProjectedCurve is not ported (curve path)");
    }

    /// OCCT LoadCurve (L2407-2518) — curve path, anchor-out-of-scope.
    fn load_curve(&mut self, _nb_boucle: i32, _order_max: i32) {
        unimplemented!("LoadCurve is not ported (curve path)");
    }
}

impl Default for BuildPlateSurface {
    /// OCCT ctor 3 default arguments (Degree=3, NbPtsOnCur=10, NbIter=3,
    /// Tol2d=1e-5, Tol3d=1e-4, TolAng=0.01, TolCurv=0.1, Anisotropie=false).
    fn default() -> Self {
        Self::new(3, 10, 3, 1.0e-5, 1.0e-4, 0.01, 0.1, false)
    }
}
