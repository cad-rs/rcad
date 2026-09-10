//! OCCT GeomAPI_ProjectPointOnSurf (GeomAPI_ProjectPointOnSurf.hxx/.cxx) —
//! orthogonal projection of a point on a surface within a UV domain, with the
//! nearest-solution queries.
//!
//! 1:1 translation of GeomAPI_ProjectPointOnSurf.cxx (whole file, L25-314).
//! The point-surface extrema run on the rcad `Extrema_ExtPS` port
//! (`crate::base::extrema::ExtPS` — the Grad-style grid + Newton search;
//! OCCT's Extrema_ExtAlgo_Tree variant is not ported, the default
//! Extrema_ExtAlgo_Grad is the one used by the GTests).

use glam::DVec3;

use crate::base::extrema::ExtPS;
use crate::core::precision::{CONFUSION, PCONFUSION};
use crate::geom::{Surface3, SurfaceEval};

/// OCCT GeomAPI_ProjectPointOnSurf — projection of a point on a surface.
pub struct ProjectPointOnSurf {
    is_done: bool,
    index: usize,
    surface: Option<Surface3>,
    umin: f64,
    usup: f64,
    vmin: f64,
    vsup: f64,
    tol_u: f64,
    tol_v: f64,
    ext: Option<ExtPS>,
}

impl ProjectPointOnSurf {
    /// OCCT GeomAPI_ProjectPointOnSurf() (L25-29) — empty constructor.
    pub fn new() -> Self {
        ProjectPointOnSurf {
            is_done: false,
            index: 0,
            surface: None,
            umin: 0.0,
            usup: 0.0,
            vmin: 0.0,
            vsup: 0.0,
            tol_u: 0.0,
            tol_v: 0.0,
            ext: None,
        }
    }

    /// OCCT GeomAPI_ProjectPointOnSurf(P, Surface, Algo) (L33-38).
    pub fn new_point(point: DVec3, surface: &Surface3) -> Self {
        let mut p = ProjectPointOnSurf::new();
        p.init_point(point, surface, CONFUSION);
        p
    }

    /// OCCT GeomAPI_ProjectPointOnSurf(P, Surface, Tolerance, Algo) (L42-50).
    pub fn new_point_tol(point: DVec3, surface: &Surface3, tolerance: f64) -> Self {
        let mut p = ProjectPointOnSurf::new();
        p.init_point(point, surface, tolerance);
        p
    }

    /// OCCT GeomAPI_ProjectPointOnSurf(P, Surface, Umin, Usup, Vmin, Vsup, Algo)
    /// (L52-62) — the tolerance is Precision::PConfusion().
    pub fn new_point_domain(
        point: DVec3,
        surface: &Surface3,
        umin: f64,
        usup: f64,
        vmin: f64,
        vsup: f64,
    ) -> Self {
        let mut p = ProjectPointOnSurf::new();
        p.init_point_domain(point, surface, umin, usup, vmin, vsup, PCONFUSION);
        p
    }

    /// OCCT GeomAPI_ProjectPointOnSurf(P, Surface, Umin, Usup, Vmin, Vsup,
    /// Tolerance, Algo) (L66-77).
    pub fn new_point_domain_tol(
        point: DVec3,
        surface: &Surface3,
        umin: f64,
        usup: f64,
        vmin: f64,
        vsup: f64,
        tolerance: f64,
    ) -> Self {
        let mut p = ProjectPointOnSurf::new();
        p.init_point_domain(point, surface, umin, usup, vmin, vsup, tolerance);
        p
    }

    /// OCCT GeomAPI_ProjectPointOnSurf::Init() (L81-101) — refresh
    /// `is_done` from the extrema and locate the minimal-distance index.
    fn refresh(&mut self) {
        self.is_done = match &self.ext {
            Some(ext) => ext.is_done() && ext.nb_ext() > 0,
            None => false,
        };

        if self.is_done {
            let ext = self.ext.as_ref().unwrap();
            let mut dist2_min = ext.square_distance(1);
            self.index = 1;
            for i in 2..=ext.nb_ext() {
                let dist2 = ext.square_distance(i);
                if dist2 < dist2_min {
                    dist2_min = dist2;
                    self.index = i;
                }
            }
        }
    }

    /// OCCT GeomAPI_ProjectPointOnSurf::Init(P, Surface, Algo) (L105-111) —
    /// the domain is the surface's natural bounds, the tolerance is
    /// Precision::Confusion().
    pub fn init_point(&mut self, point: DVec3, surface: &Surface3, tolerance: f64) {
        let dom = surface.default_domain();
        self.init_point_domain(point, surface, dom[0], dom[1], dom[2], dom[3], tolerance);
    }

    /// OCCT GeomAPI_ProjectPointOnSurf::Init(P, Surface, Umin, Usup, Vmin,
    /// Vsup, Algo) (L136-153) — the tolerance is Precision::PConfusion().
    pub fn init_point_domain_pc(
        &mut self,
        point: DVec3,
        surface: &Surface3,
        umin: f64,
        usup: f64,
        vmin: f64,
        vsup: f64,
    ) {
        self.init_point_domain(point, surface, umin, usup, vmin, vsup, PCONFUSION);
    }

    /// OCCT GeomAPI_ProjectPointOnSurf::Init(P, Surface, Umin, Usup, Vmin,
    /// Vsup, Tolerance, Algo) (L157-173).
    pub fn init_point_domain(
        &mut self,
        point: DVec3,
        surface: &Surface3,
        umin: f64,
        usup: f64,
        vmin: f64,
        vsup: f64,
        tolerance: f64,
    ) {
        self.load_surface(surface, umin, usup, vmin, vsup);
        self.run_extrema(point);
        self.refresh();
    }

    /// OCCT GeomAPI_ProjectPointOnSurf::Init(Surface, Umin, Usup, Vmin, Vsup,
    /// Algo) (L177-192) — the tolerance is Precision::PConfusion(); no point
    /// is projected yet (use Perform(P)).
    pub fn init_surface(&mut self, surface: &Surface3, umin: f64, usup: f64, vmin: f64, vsup: f64) {
        self.tol_u = PCONFUSION;
        self.tol_v = PCONFUSION;
        self.load_surface(surface, umin, usup, vmin, vsup);
        self.is_done = false;
    }

    /// OCCT GeomAPI_ProjectPointOnSurf::Init(Surface, Umin, Usup, Vmin, Vsup,
    /// Tolerance, Algo) (L196-210).
    pub fn init_surface_tol(
        &mut self,
        surface: &Surface3,
        umin: f64,
        usup: f64,
        vmin: f64,
        vsup: f64,
        tolerance: f64,
    ) {
        self.tol_u = tolerance;
        self.tol_v = tolerance;
        self.load_surface(surface, umin, usup, vmin, vsup);
        self.is_done = false;
    }

    /// OCCT myGeomAdaptor.Load(Surface, Umin, Usup, Vmin, Vsup) — store the
    /// (possibly trimmed) surface with the given UV bounds.
    fn load_surface(&mut self, surface: &Surface3, umin: f64, usup: f64, vmin: f64, vsup: f64) {
        self.surface = Some(surface.clone());
        self.umin = umin;
        self.usup = usup;
        self.vmin = vmin;
        self.vsup = vsup;
        self.ext = None;
    }

    /// OCCT myExtPS.Initialize(...) + Perform(P) — run the extrema.
    fn run_extrema(&mut self, point: DVec3) {
        let surface = self.surface.clone().expect("surface loaded");
        self.ext = Some(ExtPS::with_domain(
            point,
            &surface,
            self.umin,
            self.usup,
            self.vmin,
            self.vsup,
            self.tol_u,
            self.tol_v,
        ));
    }

    /// OCCT GeomAPI_ProjectPointOnSurf::Perform(P) (L214-218).
    pub fn perform(&mut self, point: DVec3) {
        self.run_extrema(point);
        self.refresh();
    }

    /// OCCT GeomAPI_ProjectPointOnSurf::IsDone (L222-225).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// OCCT GeomAPI_ProjectPointOnSurf::NbPoints (L229-239).
    pub fn nb_points(&self) -> usize {
        if self.is_done {
            self.ext.as_ref().map_or(0, |e| e.nb_ext())
        } else {
            0
        }
    }

    /// OCCT GeomAPI_ProjectPointOnSurf::Point(Index) (L243-248) — 1-based,
    /// panics (Standard_OutOfRange) out of range.
    pub fn point(&self, index: usize) -> DVec3 {
        assert!(
            index >= 1 && index <= self.nb_points(),
            "GeomAPI_ProjectPointOnSurf::Point"
        );
        self.ext.as_ref().unwrap().point(index).point
    }

    /// OCCT GeomAPI_ProjectPointOnSurf::Parameters(Index, U, V) (L252-257).
    pub fn parameters(&self, index: usize) -> (f64, f64) {
        assert!(
            index >= 1 && index <= self.nb_points(),
            "GeomAPI_ProjectPointOnSurf::Parameters"
        );
        let p = self.ext.as_ref().unwrap().point(index);
        (p.u, p.v)
    }

    /// OCCT GeomAPI_ProjectPointOnSurf::Distance(Index) (L261-266).
    pub fn distance(&self, index: usize) -> f64 {
        assert!(
            index >= 1 && index <= self.nb_points(),
            "GeomAPI_ProjectPointOnSurf::Distance"
        );
        self.ext.as_ref().unwrap().square_distance(index).sqrt()
    }

    /// OCCT GeomAPI_ProjectPointOnSurf::NearestPoint (L270-275) — panics
    /// (StdFail_NotDone) when the projection failed.
    pub fn nearest_point(&self) -> DVec3 {
        assert!(self.is_done, "GeomAPI_ProjectPointOnSurf::NearestPoint");
        self.ext.as_ref().unwrap().point(self.index).point
    }

    /// OCCT GeomAPI_ProjectPointOnSurf::LowerDistanceParameters (L293-298).
    pub fn lower_distance_parameters(&self) -> (f64, f64) {
        assert!(
            self.is_done,
            "GeomAPI_ProjectPointOnSurf::LowerDistanceParameters"
        );
        let p = self.ext.as_ref().unwrap().point(self.index);
        (p.u, p.v)
    }

    /// OCCT GeomAPI_ProjectPointOnSurf::LowerDistance (L302-307).
    pub fn lower_distance(&self) -> f64 {
        assert!(self.is_done, "GeomAPI_ProjectPointOnSurf::LowerDistance");
        self.ext.as_ref().unwrap().square_distance(self.index).sqrt()
    }
}

impl Default for ProjectPointOnSurf {
    fn default() -> Self {
        Self::new()
    }
}
