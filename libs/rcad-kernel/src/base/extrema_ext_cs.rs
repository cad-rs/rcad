//! OCCT Extrema_ExtCS (TKGeomBase/Extrema/Extrema_ExtCS.hxx L38-128 and
//! Extrema_ExtCS.cxx L37-630) — all the extremum distances between a curve and
//! a surface.
//!
//! Companion classes of the same OCCT package consumed here:
//! - `Extrema_ExtElCS` (analytic curve-surface extrema) — translated in
//!   [`super::extrema_ext_el_cs`],
//! - `Extrema_ExtPS` (point-surface extrema) — the pre-existing kernel engine
//!   [`crate::base::extrema::ExtPS`],
//! - `Extrema_GenExtCS` / `Extrema_FuncExtCS` — GAP carriers, see
//!   [`ExtremaGenExtCS`] / [`ExtremaFuncExtCS`].
//!
//! GAP carriers (dependency + OCCT anchor + preserved failure path):
//! - `BndLib_AddSurface::Add(S, u1, u2, v1, v2, tol, box)` + `Bnd_Box::Get`
//!   (TKGeomBase/BndLib) have no adaptor-range form in rcad, so the line-range
//!   clipping of cxx L159-194 is skipped.  Preserved OCCT path: the
//!   `if (!(IsInfinite(ufirst) || ...))` guard of cxx L159-161 evaluating
//!   false — `cfirst`/`clast` keep the curve parameters, exactly as the carrier
//!   leaves them.
//! - `Extrema_GenExtCS` (cxx L228-256, L418-433) — its PSO seed generators
//!   (`math_PSO`, `math_PSOParticlesPool`, `Extrema_GlobOptFuncCS` /
//!   `GlobOptFuncConicS` / `GlobOptFuncCQuadric`, `Extrema_GenLocateExtPS`,
//!   `GeomGridEval_Surface`) are not translated, so the carrier answers
//!   `IsDone() == false`.  Preserved OCCT path: `myDone = Ext.IsDone()`
//!   (cxx L435-436) — Extrema_ExtCS then reports not-done exactly as it does
//!   when the generic search finds nothing.

use glam::DVec3;

use crate::base::extrema::POnCurve;
use crate::base::extrema::POnSurface;
use crate::base::extrema_curve_tool::ExtremaCurveTool;
use crate::base::extrema_ext_el_cs::{ExtElCSSurface, ExtremaExtElCS};
use crate::base::extrema_ext_elc::{elclib_in_period, REAL_LAST};
use crate::base::proj_lib::adaptor::{Adaptor3dSurface, GeomAbsSurfaceType};
use crate::base::proj_lib::CurveType;
use crate::core::precision::{is_infinite_value, CONFUSION, PCONFUSION};
use crate::geom::Surface3;
use crate::math::el::{
    elslib_cone_parameters, elslib_cone_value, elslib_cylinder_value, elslib_cylinder_parameters, elslib_plane_parameters,
    elslib_plane_value, elslib_sphere_parameters, elslib_sphere_value, elslib_torus_parameters,
    elslib_torus_value,
};

/// OCCT Extrema_ExtCS (hxx L38-128).
pub struct ExtremaExtCS<'a> {
    /// hxx L112: const Adaptor3d_Surface* myS.
    my_s: Option<&'a dyn Adaptor3dSurface>,
    /// hxx L113: bool myDone.
    my_done: bool,
    /// hxx L114: bool myIsPar.
    my_is_par: bool,
    /// hxx L115: Extrema_ExtElCS myExtElCS.
    my_ext_el_cs: ExtremaExtElCS,
    /// hxx L116: NCollection_Sequence<Extrema_POnSurf> myPOnS.
    my_p_on_s: Vec<POnSurface>,
    /// hxx L117: NCollection_Sequence<Extrema_POnCurv> myPOnC.
    my_p_on_c: Vec<POnCurve>,
    /// hxx L118-121: myuinf / myusup / myvinf / myvsup.
    my_uinf: f64,
    my_usup: f64,
    my_vinf: f64,
    my_vsup: f64,
    /// hxx L122: double mytolC.
    my_tol_c: f64,
    /// hxx L123: double mytolS.
    my_tol_s: f64,
    /// hxx L124: double myucinf.
    my_ucinf: f64,
    /// hxx L125: double myucsup.
    my_ucsup: f64,
    /// hxx L126: NCollection_Sequence<double> mySqDist.
    my_sq_dist: Vec<f64>,
    /// hxx L127: GeomAbs_SurfaceType myStype.
    my_stype: GeomAbsSurfaceType,
}

impl<'a> ExtremaExtCS<'a> {
    /// OCCT Extrema_ExtCS() (cxx L37-51).
    pub fn new() -> Self {
        ExtremaExtCS {
            my_s: None,
            my_done: false,
            my_is_par: false,
            my_ext_el_cs: ExtremaExtElCS::new(),
            my_p_on_s: Vec::new(),
            my_p_on_c: Vec::new(),
            my_uinf: 0.0,
            my_usup: 0.0,
            my_vinf: 0.0,
            my_vsup: 0.0,
            my_tol_c: 0.0,
            my_tol_s: 0.0,
            my_ucinf: 0.0,
            my_ucsup: 0.0,
            my_sq_dist: Vec::new(),
            my_stype: GeomAbsSurfaceType::OtherSurface,
        }
    }

    /// OCCT Extrema_ExtCS(C, S, TolC, TolS) (cxx L53-61).
    pub fn new_curve_surface(
        c: &dyn ExtremaCurveTool,
        s: &'a dyn Adaptor3dSurface,
        tol_c: f64,
        tol_s: f64,
    ) -> Self {
        let mut this = ExtremaExtCS::new();
        this.initialize(s, tol_c, tol_s);
        let (first, last) = (c.first_parameter(), c.last_parameter());
        this.perform(c, first, last);
        this
    }

    /// OCCT Extrema_ExtCS::Initialize(S, TolC, TolS) (cxx L79-88).
    pub fn initialize(&mut self, s: &'a dyn Adaptor3dSurface, tol_c: f64, tol_s: f64) {
        let (u1, u2, v1, v2) = (
            s.first_u_parameter(),
            s.last_u_parameter(),
            s.first_v_parameter(),
            s.last_v_parameter(),
        );
        self.initialize_range(s, u1, u2, v1, v2, tol_c, tol_s);
    }

    /// OCCT Extrema_ExtCS::Initialize(S, Uinf, Usup, Vinf, Vsup, TolC, TolS)
    /// (cxx L90-107).
    #[allow(clippy::too_many_arguments)]
    pub fn initialize_range(
        &mut self,
        s: &'a dyn Adaptor3dSurface,
        uinf: f64,
        usup: f64,
        vinf: f64,
        vsup: f64,
        tol_c: f64,
        tol_s: f64,
    ) {
        self.my_s = Some(s);
        self.my_is_par = false;
        self.my_uinf = uinf;
        self.my_usup = usup;
        self.my_vinf = vinf;
        self.my_vsup = vsup;
        self.my_tol_c = tol_c;
        self.my_tol_s = tol_s;
        self.my_stype = s.get_type();
    }

    /// The kernel surface the adaptor wraps (the bridge the Extrema engine
    /// needs for the point-surface calls; the analytic branches take the
    /// elementary payloads instead).
    fn kernel_sur(&self) -> &'a Surface3 {
        self.my_s
            .expect("Extrema_ExtCS: myS")
            .kernel_surface()
            .expect("Extrema_ExtCS: adaptor without a kernel surface")
    }

    /// OCCT Extrema_ExtCS::Perform(C, Uinf, Usup) (cxx L109-528).
    pub fn perform(&mut self, c: &dyn ExtremaCurveTool, uinf: f64, usup: f64) {
        self.my_ucinf = uinf;
        self.my_ucsup = usup;
        self.my_p_on_s.clear();
        self.my_p_on_c.clear();
        self.my_sq_dist.clear();

        // cxx L118-119.
        let nb_t = 12;
        let mut nb_u = 10;
        let mut nb_v = 10;
        let my_ctype = c.get_type();

        self.my_done = false;
        // cxx L124.
        let mut is_compute_analytic = true;

        // cxx L126-295: the switch on the curve type and the surface type.
        match my_ctype {
            CurveType::Line => match self.my_stype {
                GeomAbsSurfaceType::Sphere => {
                    let Surface3::Sphere(sph) = self.kernel_sur() else {
                        panic!("Extrema_ExtCS: surface kind mismatch")
                    };
                    self.my_ext_el_cs =
                        ExtremaExtElCS::with_line(&c.line(), ExtElCSSurface::Sphere(sph));
                }
                GeomAbsSurfaceType::Cylinder => {
                    let Surface3::Cylinder(cyl) = self.kernel_sur() else {
                        panic!("Extrema_ExtCS: surface kind mismatch")
                    };
                    self.my_ext_el_cs =
                        ExtremaExtElCS::with_line(&c.line(), ExtElCSSurface::Cylinder(cyl));
                }
                GeomAbsSurfaceType::Plane => {
                    let Surface3::Plane(plane) = self.kernel_sur() else {
                        panic!("Extrema_ExtCS: surface kind mismatch")
                    };
                    self.my_ext_el_cs =
                        ExtremaExtElCS::with_line(&c.line(), ExtElCSSurface::Plane(plane));
                    // cxx L141-144.
                    if !self.my_ext_el_cs.is_parallel() {
                        self.line_generic_branch(c, nb_t, &mut nb_u, &mut nb_v);
                        return;
                    }
                }
                _ => {
                    // cxx L147-258: Torus / Cone / Bezier / BSpline /
                    // SurfaceOfRevolution / SurfaceOfExtrusion / OffsetSurface /
                    // OtherSurface share the generic line branch.
                    self.line_generic_branch(c, nb_t, &mut nb_u, &mut nb_v);
                    return;
                }
            },
            CurveType::Circle => match self.my_stype {
                // cxx L264-281.
                GeomAbsSurfaceType::Cylinder => {
                    let Surface3::Cylinder(cyl) = self.kernel_sur() else {
                        panic!("Extrema_ExtCS: surface kind mismatch")
                    };
                    self.my_ext_el_cs =
                        ExtremaExtElCS::with_circle(&c.circle(), ExtElCSSurface::Cylinder(cyl));
                }
                GeomAbsSurfaceType::Plane => {
                    let Surface3::Plane(plane) = self.kernel_sur() else {
                        panic!("Extrema_ExtCS: surface kind mismatch")
                    };
                    self.my_ext_el_cs =
                        ExtremaExtElCS::with_circle(&c.circle(), ExtElCSSurface::Plane(plane));
                }
                GeomAbsSurfaceType::Sphere => {
                    let Surface3::Sphere(sph) = self.kernel_sur() else {
                        panic!("Extrema_ExtCS: surface kind mismatch")
                    };
                    self.my_ext_el_cs =
                        ExtremaExtElCS::with_circle(&c.circle(), ExtElCSSurface::Sphere(sph));
                }
                _ => is_compute_analytic = false,
            },
            CurveType::Hyperbola => {
                // cxx L282-289.
                if self.my_stype == GeomAbsSurfaceType::Plane {
                    let Surface3::Plane(plane) = self.kernel_sur() else {
                        panic!("Extrema_ExtCS: surface kind mismatch")
                    };
                    self.my_ext_el_cs =
                        ExtremaExtElCS::with_hyperbola(&c.hyperbola(), ExtElCSSurface::Plane(plane));
                } else {
                    is_compute_analytic = false;
                }
            }
            _ => is_compute_analytic = false,
        }

        if is_compute_analytic {
            // cxx L297-415.
            if self.my_ext_el_cs.is_done() {
                self.my_done = true;
                self.my_is_par = self.my_ext_el_cs.is_parallel();
                if self.my_is_par {
                    // cxx L303-307.
                    let d = self.my_ext_el_cs.square_distance(1);
                    self.my_sq_dist.push(d);
                } else {
                    // cxx L309-319.
                    let nb_ext = self.my_ext_el_cs.nb_ext();
                    for i in 1..=nb_ext {
                        let mut pc = POnCurve {
                            param: 0.0,
                            point: DVec3::ZERO,
                        };
                        let mut ps = POnSurface {
                            u: 0.0,
                            v: 0.0,
                            point: DVec3::ZERO,
                        };
                        self.my_ext_el_cs.points(i, &mut pc, &mut ps);
                        let u_curve = pc.param;
                        let (u, v) = (ps.u, ps.v);
                        let sd = self.my_ext_el_cs.square_distance(i);
                        self.add_solution(c, u_curve, u, v, pc.point, ps.point, sd);
                    }

                    // cxx L321-411.
                    if self.my_sq_dist.is_empty() && nb_ext > 0 {
                        self.add_curve_end_solutions(c);
                    }
                }
                return;
            }
        }

        // cxx L417-433: the elementary extrema are not done, try the generic
        // solution.  GAP carrier: Extrema_GenExtCS is not translated (see the
        // module header), so the tool reports not-done and `myDone` follows
        // `Ext.IsDone()` (cxx L435) — the OCCT state whenever the generic
        // search finds nothing.
        let s = self.my_s.expect("Extrema_ExtCS::Perform: myS");
        let ext = ExtremaGenExtCS::initialize_and_perform(s, nb_u, nb_v, self.my_tol_s, c, nb_t, self.my_tol_c);
        self.my_done = ext.is_done();
    }

    /// The OCCT line × generic-surface branch (cxx L155-258).
    fn line_generic_branch(
        &mut self,
        c: &dyn ExtremaCurveTool,
        nb_t: i32,
        nb_u: &mut i32,
        nb_v: &mut i32,
    ) {
        let s = self.my_s.expect("Extrema_ExtCS::Perform: myS");
        // cxx L155-157.
        let cfirst = self.my_ucinf;
        let clast = self.my_ucsup;
        let ufirst = s.first_u_parameter();
        let ulast = s.last_u_parameter();
        let vfirst = s.first_v_parameter();
        let vlast = s.last_v_parameter();

        // cxx L159-194: the surface box clipping.  GAP carrier: the
        // adaptor-form `BndLib_AddSurface::Add` + `Bnd_Box::Get` are not
        // translated, so the guard is treated as false and cfirst/clast keep
        // the curve parameters (see the module header).
        if !(is_infinite_value(ufirst)
            || is_infinite_value(ulast)
            || is_infinite_value(vfirst)
            || is_infinite_value(vlast))
        {
            // (clipping skipped — see the module header)
        }

        // cxx L196-203.
        if s.is_u_periodic() {
            *nb_u = 13;
        }
        if s.is_v_periodic() {
            *nb_v = 13;
        }

        // cxx L205-226.
        if clast - cfirst <= CONFUSION {
            let a_c_par = (cfirst + clast) / 2.0;
            let a_pm = c.value(a_c_par);
            let an_ext_ps = crate::base::extrema::ExtPS::with_domain(
                a_pm,
                self.kernel_sur(),
                ufirst,
                ulast,
                vfirst,
                vlast,
                self.my_tol_s,
                self.my_tol_s,
            );
            self.my_done = an_ext_ps.is_done();
            if self.my_done {
                // cxx L214-224.
                let nb_ext = an_ext_ps.nb_ext();
                let t = a_c_par;
                for i in 1..=nb_ext {
                    let (u, v, p, sd) = {
                        let ps = an_ext_ps.point(i);
                        (ps.u, ps.v, ps.point, an_ext_ps.square_distance(i))
                    };
                    self.add_solution(c, t, u, v, DVec3::ZERO, p, sd);
                }
            }
            return;
        }

        // cxx L228-256: Extrema_GenExtCS Ext(C, *myS, NbT, NbU, NbV, cfirst,
        // clast, ufirst, ulast, vfirst, vlast, mytolC, mytolS).  GAP carrier
        // (see the module header): not-done, so the `if (myDone)` block of
        // cxx L243-257 is skipped and `myDone` is false.
        let _ = (nb_t, nb_u, nb_v, cfirst, clast);
        let _ = (ufirst, ulast, vfirst, vlast);
        self.my_done = false;
    }

    /// OCCT Extrema_ExtCS::Perform cxx L321-411 — the curve-extremity fallback
    /// when the analytical extrema fall outside the boundaries.
    fn add_curve_end_solutions(&mut self, c: &dyn ExtremaCurveTool) {
        let a_t = [self.my_ucinf, self.my_ucsup];
        let mut a_p_on_c = [DVec3::ZERO; 2];
        let mut a_p_on_s = [DVec3::ZERO; 2];
        let mut u = [0.0f64; 2];
        let mut v = [0.0f64; 2];
        let mut a_dist = [-1.0f64, -1.0];

        for i in 0..2 {
            if is_infinite_value(a_t[i]) {
                continue;
            }

            a_p_on_c[i] = c.value(a_t[i]);
            match self.my_stype {
                GeomAbsSurfaceType::Plane => {
                    let Surface3::Plane(plane) = self.kernel_sur() else {
                        continue;
                    };
                    let (uu, vv) = elslib_plane_parameters(
                        a_p_on_c[i],
                        plane.origin,
                        plane.u_dir,
                        plane.v_dir,
                    );
                    u[i] = uu;
                    v[i] = vv;
                    a_p_on_s[i] =
                        elslib_plane_value(uu, vv, plane.origin, plane.u_dir, plane.v_dir);
                }
                GeomAbsSurfaceType::Sphere => {
                    let Surface3::Sphere(sphere) = self.kernel_sur() else {
                        continue;
                    };
                    let x = sphere.ref_dir;
                    let y = sphere.axis.cross(sphere.ref_dir);
                    let (uu, vv) =
                        elslib_sphere_parameters(a_p_on_c[i], sphere.center, x, y, sphere.axis);
                    u[i] = uu;
                    v[i] = vv;
                    a_p_on_s[i] = elslib_sphere_value(
                        uu,
                        vv,
                        sphere.center,
                        sphere.axis,
                        x,
                        sphere.radius,
                    );
                }
                GeomAbsSurfaceType::Cylinder => {
                    let Surface3::Cylinder(cyl) = self.kernel_sur() else {
                        continue;
                    };
                    let x = cyl.ref_dir;
                    let y = cyl.axis.cross(cyl.ref_dir);
                    let (uu, vv) = elslib_cylinder_parameters(
                        a_p_on_c[i],
                        cyl.origin,
                        x,
                        y,
                        cyl.axis,
                        cyl.radius,
                    );
                    u[i] = uu;
                    v[i] = vv;
                    a_p_on_s[i] =
                        elslib_cylinder_value(uu, vv, cyl.origin, cyl.axis, x, cyl.radius);
                }
                GeomAbsSurfaceType::Torus => {
                    let Surface3::Torus(tor) = self.kernel_sur() else {
                        continue;
                    };
                    let x = tor.ref_dir;
                    let y = tor.axis.cross(tor.ref_dir);
                    let (uu, vv) = elslib_torus_parameters(
                        a_p_on_c[i],
                        tor.center,
                        x,
                        y,
                        tor.axis,
                        tor.major_radius,
                        tor.minor_radius,
                    );
                    u[i] = uu;
                    v[i] = vv;
                    a_p_on_s[i] = elslib_torus_value(
                        uu,
                        vv,
                        tor.center,
                        tor.axis,
                        tor.major_radius,
                        tor.minor_radius,
                    );
                }
                GeomAbsSurfaceType::Cone => {
                    let Surface3::Cone(cone) = self.kernel_sur() else {
                        continue;
                    };
                    let x = cone.ref_dir;
                    let y = cone.axis.cross(cone.ref_dir);
                    let (uu, vv) = elslib_cone_parameters(
                        a_p_on_c[i],
                        cone.apex,
                        x,
                        y,
                        cone.axis,
                        cone.radius,
                        cone.half_angle_rad,
                    );
                    u[i] = uu;
                    v[i] = vv;
                    // OCCT: ElSLib::Value(U, V, myS->Cone()).
                    a_p_on_s[i] = elslib_cone_value(
                        uu,
                        vv,
                        cone.apex,
                        cone.axis,
                        cone.half_angle_rad,
                        cone.radius,
                    );
                }
                _ => continue,
            }

            a_dist[i] = (a_p_on_c[i] - a_p_on_s[i]).length_squared();
        }

        // cxx L370-410: choose the solution(s) to add.
        let mut b_add = [false, false];
        if a_dist[0] >= 0.0 && a_dist[1] >= 0.0 {
            let a_diff = a_dist[0] - a_dist[1];
            if a_diff.abs() < CONFUSION {
                b_add[0] = true;
                b_add[1] = true;
            } else if a_diff < 0.0 {
                b_add[0] = true;
            } else {
                b_add[1] = true;
            }
        } else if a_dist[0] >= 0.0 {
            b_add[0] = true;
        } else if a_dist[1] >= 0.0 {
            b_add[1] = true;
        }

        for i in 0..2 {
            if b_add[i] {
                self.add_solution(c, a_t[i], u[i], v[i], a_p_on_c[i], a_p_on_s[i], a_dist[i]);
            }
        }
    }

    /// OCCT Extrema_ExtCS::AddSolution (cxx L576-630).
    #[allow(clippy::too_many_arguments)]
    fn add_solution(
        &mut self,
        the_curve: &dyn ExtremaCurveTool,
        a_t: f64,
        a_u: f64,
        a_v: f64,
        point_on_curve: DVec3,
        point_on_surf: DVec3,
        square_dist: f64,
    ) -> bool {
        let s = self.my_s.expect("Extrema_ExtCS::AddSolution: myS");
        let mut added = false;

        let mut t = a_t;
        let mut u = a_u;
        let mut v = a_v;

        // cxx L588-599.
        if the_curve.is_periodic() {
            t = elclib_in_period(t, self.my_ucinf, self.my_ucinf + the_curve.period());
        }
        if s.is_u_periodic() {
            u = elclib_in_period(u, self.my_uinf, self.my_uinf + s.u_period());
        }
        if s.is_v_periodic() {
            v = elclib_in_period(v, self.my_vinf, self.my_vinf + s.v_period());
        }

        // cxx L603-604.
        if (self.my_ucinf - t) <= self.my_tol_c
            && (t - self.my_ucsup) <= self.my_tol_c
            && (self.my_uinf - u) <= self.my_tol_s
            && (u - self.my_usup) <= self.my_tol_s
            && (self.my_vinf - v) <= self.my_tol_s
            && (v - self.my_vsup) <= self.my_tol_s
        {
            let mut is_new_solution = true;
            // cxx L607-619.
            for j in 1..=self.my_sq_dist.len() {
                let tj = self.my_p_on_c[j - 1].param;
                let uj = self.my_p_on_s[j - 1].u;
                let vj = self.my_p_on_s[j - 1].v;
                if (t - tj).abs() <= self.my_tol_c
                    && (u - uj).abs() <= self.my_tol_s
                    && (v - vj).abs() <= self.my_tol_s
                {
                    is_new_solution = false;
                    break;
                }
            }
            if is_new_solution {
                self.my_sq_dist.push(square_dist);
                self.my_p_on_c.push(POnCurve {
                    param: t,
                    point: point_on_curve,
                });
                self.my_p_on_s.push(POnSurface {
                    u,
                    v,
                    point: point_on_surf,
                });
                added = true;
            }
        }
        added
    }

    /// OCCT IsDone() (cxx L530-533).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT IsParallel() (cxx L535-543).
    pub fn is_parallel(&self) -> bool {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_is_par
    }

    /// OCCT NbExt() (cxx L555-563).
    pub fn nb_ext(&self) -> usize {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_sq_dist.len()
    }

    /// OCCT SquareDistance(N) (cxx L545-553) — 1-based.
    pub fn square_distance(&self, n: usize) -> f64 {
        if n < 1 || n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.my_sq_dist[n - 1]
    }

    /// OCCT Points(N, P1, P2) (cxx L565-574) — 1-based.
    pub fn points(&self, n: usize, p1: &mut POnCurve, p2: &mut POnSurface) {
        if n < 1 || n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        *p1 = self.my_p_on_c[n - 1].clone();
        *p2 = self.my_p_on_s[n - 1].clone();
    }
}

impl Default for ExtremaExtCS<'_> {
    fn default() -> Self {
        ExtremaExtCS::new()
    }
}

/// OCCT Extrema_GenExtCS (hxx L33-163) — GAP carrier.
///
/// Dependency: the class needs `math_PSO` and `math_PSOParticlesPool`
/// (TKMath/math), `Extrema_GlobOptFuncCS` / `Extrema_GlobOptFuncConicS` /
/// `Extrema_GlobOptFuncCQuadric` and `Extrema_GenLocateExtPS`
/// (TKGeomBase/Extrema) plus `GeomGridEval_Surface` — none translated.
///
/// Preserved OCCT failure path: the carrier constructs, initializes and
/// performs as a not-done tool, so `Extrema_ExtCS` reads `Ext.IsDone() == false`
/// (Extrema_ExtCS.cxx L435) and reports not-done — the OCCT state whenever the
/// generic search produces no solution.
pub struct ExtremaGenExtCS;

impl ExtremaGenExtCS {
    /// OCCT Extrema_GenExtCS(C, S, NbT, NbU, NbV, tmin, tsup, Umin, Usup,
    /// Vmin, Vsup, Tol1, Tol2) (cxx L174-190) reduced to the not-done carrier.
    #[allow(clippy::too_many_arguments)]
    pub fn initialize_and_perform(
        _s: &dyn Adaptor3dSurface,
        _nb_u: i32,
        _nb_v: i32,
        _tol2: f64,
        _c: &dyn ExtremaCurveTool,
        _nb_t: i32,
        _tol1: f64,
    ) -> Self {
        ExtremaGenExtCS
    }

    /// OCCT IsDone() (cxx L827-830).
    pub fn is_done(&self) -> bool {
        false
    }
}

/// OCCT Extrema_FuncExtCS (hxx L37-99, cxx L48-238) — the function set
/// { F1(t,u,v) = (C(t)-S(u,v)).Dtc(t) }, { F2 = (C(t)-S(u,v)).Dus(u,v) },
/// { F3 = (C(t)-S(u,v)).Dvs(u,v) } whose roots are the curve-surface
/// extremalities.
pub struct ExtremaFuncExtCS;

impl ExtremaFuncExtCS {
    /// OCCT NbVariables() (cxx L81-84).
    pub fn nb_variables() -> i32 {
        3
    }

    /// OCCT NbEquations() (cxx L88-91).
    pub fn nb_equations() -> i32 {
        3
    }

    /// OCCT Values(UV, F, Df) (cxx L132-167) — the residual and the Jacobian of
    /// the curve-surface extremality system.
    pub fn values(
        c_d2: (DVec3, DVec3, DVec3),
        s_d2: (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3),
        f: &mut [f64; 3],
        df: &mut [[f64; 3]; 3],
    ) {
        let (p1, dtc, dttc) = c_d2;
        let (p2, dus, dvs, duus, dvvs, duvs) = s_d2;

        let p1p2 = p1 - p2;

        // cxx L150-152.
        f[0] = p1p2.dot(dtc);
        f[1] = p1p2.dot(dus);
        f[2] = p1p2.dot(dvs);

        // cxx L154-164.
        df[0][0] = dtc.length_squared() + p1p2.dot(dttc);
        df[0][1] = -dus.dot(dtc);
        df[0][2] = -dvs.dot(dtc);
        df[1][0] = -df[0][1];
        df[1][1] = -dus.length_squared() + p1p2.dot(duus);
        df[1][2] = -dvs.dot(dus) + p1p2.dot(duvs);
        df[2][0] = -df[0][2];
        df[2][1] = df[1][2];
        df[2][2] = -dvs.length_squared() + p1p2.dot(dvvs);
    }

    /// OCCT GetStateNumber() (cxx L171-198) — deduplicates by the curve
    /// parameter with the OCCT SquarePConfusion tolerance and records the
    /// solution.
    pub fn record_solution(
        t: f64,
        u: f64,
        v: f64,
        p1: DVec3,
        p2: DVec3,
        sq_dist: &mut Vec<f64>,
        points_on_curve: &mut Vec<POnCurve>,
        points_on_surface: &mut Vec<POnSurface>,
    ) {
        let tol2d = PCONFUSION * PCONFUSION;
        let nb_sol = sq_dist.len();
        let mut i = 1usize;
        while i <= nb_sol {
            let mut a_t = points_on_curve[i - 1].param;
            a_t -= t;
            a_t *= a_t;
            if a_t <= tol2d {
                break;
            }
            i += 1;
        }
        if i <= nb_sol {
            return;
        }
        sq_dist.push((p1 - p2).length_squared());
        points_on_curve.push(POnCurve { param: t, point: p1 });
        points_on_surface.push(POnSurface { u, v, point: p2 });
    }
}

/// OCCT RealLast() reuse marker for the table of squared distances.
#[allow(dead_code)]
fn _real_last_marker() -> f64 {
    REAL_LAST
}
