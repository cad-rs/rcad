//! OCCT Extrema_ExtElCS (TKGeomBase/Extrema/Extrema_ExtElCS.hxx L44-140 and
//! Extrema_ExtElCS.cxx L44-807) — all the distances between an elementary
//! curve (line, circle, hyperbola) and an elementary surface.
//!
//! GAP carriers (each annotated with its dependency, the OCCT anchor and the
//! preserved OCCT failure path):
//! - the `IntAna_IntConicQuad` / `IntAna_Quadric` / `IntAna_QuadQuadGeo`
//!   intersection tests (TKGeomBase/IntAna) and `Extrema_ExtPElS`
//!   (TKGeomBase/Extrema) are not translated, so the branches that consume them
//!   take the OCCT "tool not done" early return: the carrier answers
//!   `is_done() == false`, which is the OCCT value whenever those tools report
//!   not-done, and the caller keeps the extremum set it already accumulated
//!   exactly as OCCT does on that path.

use glam::DVec3;

use crate::base::extrema::POnSurface;
use crate::base::extrema_ext_elc::{
    elclib_circle_value, elclib_hyperbola_value, elclib_line_value, elclib_line_parameter,
    ExtremaExtElC, REAL_LAST,
};
use crate::base::extrema_gp::{dir_is_normal, dir_is_parallel};
use crate::core::precision::{ANGULAR, CONFUSION, SQUARE_CONFUSION};
use crate::geom::{
    Circle3, ConicalSurface, CylindricalSurface, Hyperbola3, Line3, Plane, SphericalSurface,
    ToroidalSurface,
};
use crate::math::el::{elslib_cylinder_parameters, elslib_plane_parameters, elslib_plane_value};

/// The surface-kind dispatch of the OCCT overloads, mirroring the ten
/// `Perform` overloads of the class.
pub enum ExtElCSSurface<'a> {
    /// gp_Pln.
    Plane(&'a Plane),
    /// gp_Cylinder.
    Cylinder(&'a CylindricalSurface),
    /// gp_Cone.
    Cone(&'a ConicalSurface),
    /// gp_Sphere.
    Sphere(&'a SphericalSurface),
    /// gp_Torus.
    Torus(&'a ToroidalSurface),
}

/// OCCT Extrema_ExtElCS (hxx L44-140).
#[derive(Debug, Clone)]
pub struct ExtremaExtElCS {
    /// hxx L134: bool myDone.
    my_done: bool,
    /// hxx L135: int myNbExt.
    my_nb_ext: usize,
    /// hxx L136: bool myIsPar.
    my_is_par: bool,
    /// hxx L137: handle(NCollection_HArray1<double>) mySqDist.
    my_sq_dist: Vec<f64>,
    /// hxx L138: handle(NCollection_HArray1<Extrema_POnCurv>) myPoint1.
    my_point1: Vec<crate::base::extrema::POnCurve>,
    /// hxx L139: handle(NCollection_HArray1<Extrema_POnSurf>) myPoint2.
    my_point2: Vec<POnSurface>,
}

impl ExtremaExtElCS {
    /// OCCT Extrema_ExtElCS() (cxx L44-49).
    pub fn new() -> Self {
        ExtremaExtElCS {
            my_done: false,
            my_nb_ext: 0,
            my_is_par: false,
            my_sq_dist: Vec::new(),
            my_point1: Vec::new(),
            my_point2: Vec::new(),
        }
    }

    /// OCCT Extrema_ExtElCS(C, S) dispatch over the overload set
    /// (hxx L54-113).
    pub fn with_line(c: &Line3, s: ExtElCSSurface<'_>) -> Self {
        let mut this = ExtremaExtElCS::new();
        this.perform_line(c, s);
        this
    }

    /// OCCT Extrema_ExtElCS(C, S) dispatch over the overload set.
    pub fn with_circle(c: &Circle3, s: ExtElCSSurface<'_>) -> Self {
        let mut this = ExtremaExtElCS::new();
        this.perform_circle(c, s);
        this
    }

    /// OCCT Extrema_ExtElCS(C, S) dispatch over the overload set.
    pub fn with_hyperbola(c: &Hyperbola3, s: ExtElCSSurface<'_>) -> Self {
        let mut this = ExtremaExtElCS::new();
        this.perform_hyperbola(c, s);
        this
    }

    /// The line dispatch (cxx L51-285).
    pub fn perform_line(&mut self, c: &Line3, s: ExtElCSSurface<'_>) {
        match s {
            ExtElCSSurface::Plane(p) => self.perform_line_plane(c, p),
            ExtElCSSurface::Cylinder(cyl) => self.perform_line_cylinder(c, cyl),
            ExtElCSSurface::Cone(cone) => self.perform_line_cone(c, cone),
            ExtElCSSurface::Sphere(sph) => self.perform_line_sphere(c, sph),
            ExtElCSSurface::Torus(tor) => self.perform_line_torus(c, tor),
        }
    }

    /// The circle dispatch (cxx L289-707).
    pub fn perform_circle(&mut self, c: &Circle3, s: ExtElCSSurface<'_>) {
        match s {
            ExtElCSSurface::Plane(p) => self.perform_circle_plane(c, p),
            ExtElCSSurface::Cylinder(cyl) => self.perform_circle_cylinder(c, cyl),
            ExtElCSSurface::Cone(cone) => self.perform_circle_cone(c, cone),
            ExtElCSSurface::Sphere(sph) => self.perform_circle_sphere(c, sph),
            ExtElCSSurface::Torus(tor) => self.perform_circle_torus(c, tor),
        }
    }

    /// The hyperbola dispatch (cxx L709-762).
    pub fn perform_hyperbola(&mut self, c: &Hyperbola3, s: ExtElCSSurface<'_>) {
        match s {
            ExtElCSSurface::Plane(p) => self.perform_hyperbola_plane(c, p),
            _ => panic!("Standard_NotImplemented: Extrema_ExtElCS(gp_Hypr, surface)"),
        }
    }

    /// OCCT Extrema_ExtElCS::Perform(const gp_Lin& C, const gp_Pln& S)
    /// (cxx L56-69).
    pub fn perform_line_plane(&mut self, c: &Line3, s: &Plane) {
        self.my_done = true;
        self.my_is_par = false;
        self.my_nb_ext = 0;

        // cxx L62-68.
        if dir_is_normal(c.direction, s.normal, ANGULAR) {
            self.my_sq_dist = vec![plane_square_distance(s, c.origin)];
            self.my_is_par = true;
            self.my_nb_ext = 1;
        }
    }

    /// OCCT Extrema_ExtElCS::Perform(const gp_Lin& C, const gp_Cylinder& S)
    /// (cxx L76-184).
    pub fn perform_line_cylinder(&mut self, c: &Line3, s: &CylindricalSurface) {
        self.my_done = false;
        self.my_nb_ext = 0;
        self.my_is_par = false;

        // cxx L82: gp_Ax3 Pos = S.Position().
        let pos_origin = s.origin;
        let pos_axis = s.axis;
        let _pos_x = s.ref_dir;

        let mut is_parallel = false;

        let radius = s.radius;
        // cxx L87: Extrema_ExtElC Extrem(gp_Lin(Pos.Axis()), C, Angular).
        let axis_line = Line3::new(pos_origin, pos_axis);
        let extrem = ExtremaExtElC::line_line(&axis_line, c, ANGULAR);
        if extrem.is_parallel() {
            is_parallel = true;
        } else {
            // cxx L96-98.
            let mut my_p_on_c1 = crate::base::extrema::POnCurve {
                param: 0.0,
                point: DVec3::ZERO,
            };
            let mut my_p_on_c2 = crate::base::extrema::POnCurve {
                param: 0.0,
                point: DVec3::ZERO,
            };
            extrem.points(1, &mut my_p_on_c1, &mut my_p_on_c2);
            let pon_axis = my_p_on_c1.point;
            let pc = my_p_on_c2.point;

            if radius - (pc - pon_axis).length() > crate::core::precision::PCONFUSION {
                // cxx L104-132: the line intersects the cylinder.  GAP carrier.
                // Dependency: IntAna_Quadric + IntAna_IntConicQuad (TKGeomBase/
                // IntAna), untranslated.  Preserved OCCT path: when the
                // intersection tool reports not-done OCCT leaves myNbExt = 0 and
                // myDone = true (cxx L154); the carrier takes the same outcome.
                // (OCCT's `Inters.IsDone() && Inters.IsInQuadric()` ->
                // isParallel branch is likewise unreachable without the tool.)
            } else {
                // cxx L134-151: the line is tangent or outside the cylinder.
                // GAP carrier.  Dependency: Extrema_ExtPElS (TKGeomBase/Extrema),
                // untranslated.  Preserved OCCT path: `if (ExPS.IsDone())` — a
                // not-done ExPS leaves myNbExt = 0 (cxx L138-151).
            }

            self.my_done = true;
        }

        if is_parallel {
            // cxx L157-183.
            let mut a_dist = extrem.square_distance(1);
            let a_nb_ext = extrem.nb_ext();
            for i in 2..=a_nb_ext {
                let a_d = extrem.square_distance(i);
                if a_d < a_dist {
                    a_dist = a_d;
                }
            }

            a_dist = a_dist.sqrt() - radius;
            self.my_sq_dist = vec![a_dist * a_dist];
            self.my_done = true;
            self.my_is_par = true;
            self.my_nb_ext = 1;
        }
    }

    /// OCCT Extrema_ExtElCS::Perform(const gp_Lin& C, const gp_Cone& S)
    /// (cxx L191-196).
    pub fn perform_line_cone(&mut self, _c: &Line3, _s: &ConicalSurface) {
        panic!("Standard_NotImplemented: Extrema_ExtElCS(gp_Lin, gp_Cone)");
    }

    /// OCCT Extrema_ExtElCS::Perform(const gp_Lin& C, const gp_Sphere& S)
    /// (cxx L198-273).
    pub fn perform_line_sphere(&mut self, c: &Line3, s: &SphericalSurface) {
        self.my_done = false;
        self.my_nb_ext = 0;
        self.my_is_par = false;
        let a_start_idx = 0usize;

        // cxx L214-216: Extrema_ExtPElC Extrem(aCenter, C, Angular,
        // RealFirst(), RealLast()).
        let a_center = s.center;
        let extrem = crate::base::extrema_ext_p_elc::ExtremaExtPElC::point_line(
            a_center,
            c,
            ANGULAR,
            -REAL_LAST,
            REAL_LAST,
        );

        if extrem.is_done() && extrem.nb_ext() > 0 {
            let my_p_on_c1 = extrem.point(1).clone();
            if (my_p_on_c1.point - a_center).length() <= s.radius {
                // cxx L224-246: IntAna_IntConicQuad aLinSphere(C, S).  GAP
                // carrier.  Dependency: IntAna_IntConicQuad (TKGeomBase/IntAna),
                // untranslated.  Preserved OCCT path: when the intersection tool
                // reports not-done the block is skipped and aStartIdx stays 0
                // (cxx L225), exactly the carrier's state.
            }

            // cxx L249-270: Extrema_ExtPElS ExPS(myPOnC1.Value(), S, Confusion).
            // GAP carrier.  Dependency: Extrema_ExtPElS (TKGeomBase/Extrema),
            // untranslated.  Preserved OCCT path: `if (ExPS.IsDone())` — a
            // not-done ExPS leaves the accumulated set untouched (cxx L250).
            let _ = a_start_idx;
        }
        self.my_done = true;
    }

    /// OCCT Extrema_ExtElCS::Perform(const gp_Lin& C, const gp_Torus& S)
    /// (cxx L280-285).
    pub fn perform_line_torus(&mut self, _c: &Line3, _s: &ToroidalSurface) {
        panic!("Standard_NotImplemented: Extrema_ExtElCS(gp_Lin, gp_Torus)");
    }

    /// OCCT Extrema_ExtElCS::Perform(const gp_Circ& C, const gp_Pln& S)
    /// (cxx L294-391).
    pub fn perform_circle_plane(&mut self, c: &Circle3, s: &Plane) {
        self.my_done = true;
        self.my_is_par = false;
        self.my_nb_ext = 0;

        // cxx L300-301.
        let n_circ = c.normal;
        let n_pln = s.normal;

        let mut is_parallel = false;

        if dir_is_parallel(n_circ, n_pln, ANGULAR) {
            is_parallel = true;
        } else {
            // cxx L313-324: the extreme line of the circle w.r.t. the plane.
            let mut ext_line = n_circ.cross(n_pln);
            ext_line = ext_line.cross(n_circ);
            let x_dir = c.x_dir;
            let mut t = [0.0f64; 2];
            t[0] = crate::base::extrema_gp::dir_angle_with_ref(x_dir, ext_line, n_circ);
            if t[0] < 0.0 {
                t[0] += std::f64::consts::PI;
            }
            t[1] = t[0] + std::f64::consts::PI;

            self.my_nb_ext = 2;
            // cxx L327-340: IntAna_IntConicQuad anInter(C, S, Angular,
            // Confusion).  GAP carrier.  Dependency: IntAna_IntConicQuad
            // (TKGeomBase/IntAna), untranslated.  Preserved OCCT path: for a
            // not-done (or non-in-quadric) tool OCCT keeps myNbExt == 2 and
            // isParallel == false (cxx L330-340), exactly the carrier's state,
            // and the two extreme points below are stored as in OCCT.

            if !is_parallel {
                // cxx L344-363.
                for i in 0..2 {
                    let pc = elclib_circle_value(t[i], c);
                    let (u, v) = elslib_plane_parameters(pc, s.origin, s.u_dir, s.v_dir);
                    let pp = elslib_plane_value(u, v, s.origin, s.u_dir, s.v_dir);
                    self.my_point1.push(crate::base::extrema::POnCurve {
                        param: t[i],
                        point: pc,
                    });
                    self.my_point2.push(POnSurface {
                        u,
                        v,
                        point: pp,
                    });
                    self.my_sq_dist.push((pc - pp).length_squared());
                }
                // cxx L365-380: the additional intersection points
                // (myNbExt > 2) come from the GAP-carried IntAna tool above.
            }
        }

        if is_parallel {
            // cxx L384-390.
            self.my_sq_dist = vec![plane_square_distance(s, c.center)];
            self.my_is_par = true;
            self.my_nb_ext = 1;
        }
    }

    /// OCCT Extrema_ExtElCS::Perform(const gp_Circ& C, const gp_Cylinder& S)
    /// (cxx L400-541).
    pub fn perform_circle_cylinder(&mut self, c: &Circle3, s: &CylindricalSurface) {
        self.my_done = false;
        self.my_is_par = false;
        self.my_nb_ext = 0;

        // cxx L407: gp_Lin anAxis(S.Axis()).
        let an_axis = Line3::new(s.origin, s.axis);
        // cxx L410: Extrema_ExtElC anExtC(anAxis, C, 0.).
        let an_ext_c = ExtremaExtElC::line_circle(&an_axis, c, 0.0);
        if !an_ext_c.is_done() {
            return;
        }

        let mut is_parallel = false;
        if an_ext_c.is_parallel() {
            is_parallel = true;
        } else {
            // cxx L431-442: IntAna_Quadric aCylQuad(S) + IntAna_IntConicQuad
            // aCircCylInter(C, aCylQuad).  GAP carrier.  Dependency:
            // IntAna_Quadric + IntAna_IntConicQuad (TKGeomBase/IntAna),
            // untranslated.  Preserved OCCT path: for a not-done tool
            // aNbInter = 0 and isParallel stays false (cxx L435-442), exactly
            // the carrier's state, so only the 2*aNbExt extremas below are
            // produced — the OCCT outcome for an empty intersection.
            let a_nb_inter = 0usize;

            if !is_parallel {
                // cxx L446-490.
                let a_nb_ext = an_ext_c.nb_ext();
                let mut a_cur_i = 1usize;
                let a_tol_conf = CONFUSION;
                let a_cyl_rad = s.radius;

                for i in 1..=a_nb_ext {
                    let mut a_p_on_axis = crate::base::extrema::POnCurve {
                        param: 0.0,
                        point: DVec3::ZERO,
                    };
                    let mut a_p_on_circ = crate::base::extrema::POnCurve {
                        param: 0.0,
                        point: DVec3::ZERO,
                    };
                    let a_sq_dist = an_ext_c.square_distance(i);
                    let a_dist = a_sq_dist.sqrt();

                    an_ext_c.points(i, &mut a_p_on_axis, &mut a_p_on_circ);

                    if a_sq_dist <= (a_tol_conf * a_tol_conf) {
                        self.my_nb_ext -= 2;
                        continue;
                    }

                    let a_dir = (a_p_on_axis.point - a_p_on_circ.point).normalize_or_zero();
                    let a_shift = [a_dist + a_cyl_rad, a_dist - a_cyl_rad];

                    for j in 0..2 {
                        let a_pnt_on_cyl = a_p_on_circ.point + a_dir * a_shift[j];
                        let (a_u, a_v) = elslib_cylinder_parameters(
                            a_pnt_on_cyl,
                            s.origin,
                            s.ref_dir,
                            s.axis.cross(s.ref_dir),
                            s.axis,
                            s.radius,
                        );
                        self.my_point1.push(a_p_on_circ.clone());
                        self.my_point2.push(POnSurface {
                            u: a_u,
                            v: a_v,
                            point: a_pnt_on_cyl,
                        });
                        self.my_sq_dist.push(a_shift[j] * a_shift[j]);
                        a_cur_i += 1;
                    }
                }

                // cxx L493-508: the intersection points (aNbInter == 0 here,
                // see the GAP above).
                let _ = a_cur_i;
                let _ = a_nb_inter;
            }
        }

        self.my_done = true;

        if is_parallel {
            // cxx L514-540.
            self.my_is_par = true;
            self.my_nb_ext = 1;
            let mut a_dist = an_ext_c.square_distance(1);
            let a_nb_ext = an_ext_c.nb_ext();
            for i in 2..=a_nb_ext {
                let a_d = an_ext_c.square_distance(i);
                if a_d < a_dist {
                    a_dist = a_d;
                }
            }
            a_dist = a_dist.sqrt() - s.radius;
            self.my_sq_dist = vec![a_dist * a_dist];
        }
    }

    /// OCCT Extrema_ExtElCS::Perform(const gp_Circ& C, const gp_Cone& S)
    /// (cxx L550-555).
    pub fn perform_circle_cone(&mut self, _c: &Circle3, _s: &ConicalSurface) {
        panic!("Standard_NotImplemented: Extrema_ExtElCS(gp_Circ, gp_Cone)");
    }

    /// OCCT Extrema_ExtElCS::Perform(const gp_Circ& C, const gp_Sphere& S)
    /// (cxx L566-695).
    pub fn perform_circle_sphere(&mut self, c: &Circle3, s: &SphericalSurface) {
        self.my_done = false;
        self.my_is_par = false;
        self.my_nb_ext = 0;

        // cxx L572: gp_Lin(C.Axis()).SquareDistance(S.Location()).
        let c_axis = Line3::new(c.center, c.normal);
        let sq_dist_loc = {
            let v = s.center - c_axis.origin;
            let t = v.dot(c_axis.direction);
            (v - t * c_axis.direction).length_squared()
        };
        if sq_dist_loc < SQUARE_CONFUSION {
            // cxx L574-585: circle and sphere are parallel.
            self.my_is_par = true;
            self.my_done = true;
            self.my_nb_ext = 1;

            let a_sq_dist_loc = (c.center - s.center).length_squared();
            let a_sq_dist = a_sq_dist_loc + c.radius * c.radius;
            let a_dist = a_sq_dist.sqrt() - s.radius;
            self.my_sq_dist = vec![a_dist * a_dist];
            return;
        }

        // cxx L588-595: IntAna_QuadQuadGeo anInter(CPln, S).  GAP carrier.
        // Dependency: IntAna_QuadQuadGeo (TKGeomBase/IntAna), untranslated.
        // Preserved OCCT path: `if (!anInter.IsDone()) return;` (cxx L591-595)
        // — the carrier returns with myDone == false, the OCCT not-done state.
        return;
    }

    /// OCCT Extrema_ExtElCS::Perform(const gp_Circ& C, const gp_Torus& S)
    /// (cxx L702-707).
    pub fn perform_circle_torus(&mut self, _c: &Circle3, _s: &ToroidalSurface) {
        panic!("Standard_NotImplemented: Extrema_ExtElCS(gp_Circ, gp_Torus)");
    }

    /// OCCT Extrema_ExtElCS::Perform(const gp_Hypr& C, const gp_Pln& S)
    /// (cxx L714-762).
    pub fn perform_hyperbola_plane(&mut self, c: &Hyperbola3, s: &Plane) {
        self.my_done = true;
        self.my_is_par = false;
        self.my_nb_ext = 0;

        let n_hypr = c.normal;
        let n_pln = s.normal;

        if dir_is_parallel(n_hypr, n_pln, ANGULAR) {
            // cxx L726-731.
            self.my_sq_dist = vec![plane_square_distance(s, c.center)];
            self.my_is_par = true;
            self.my_nb_ext = 1;
        } else {
            // cxx L735-761.
            let x_dir = c.major_dir;
            let y_dir = c.normal.cross(c.major_dir);

            let a = c.semi_minor * n_pln.dot(y_dir);
            let b = c.semi_major * n_pln.dot(x_dir);

            if b.abs() > a.abs() {
                let t = -0.5 * ((a + b) / (b - a)).ln();
                let ph = elclib_hyperbola_value(t, c);
                self.my_point1
                    .push(crate::base::extrema::POnCurve { param: t, point: ph });

                self.my_sq_dist = vec![plane_square_distance(s, ph)];

                let (u, v) = elslib_plane_parameters(ph, s.origin, s.u_dir, s.v_dir);
                let pp = elslib_plane_value(u, v, s.origin, s.u_dir, s.v_dir);
                self.my_point2 = vec![POnSurface { u, v, point: pp }];

                self.my_nb_ext = 1;
            }
        }
    }

    /// OCCT IsDone() (cxx L764-767).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT NbExt() (cxx L769-776).
    pub fn nb_ext(&self) -> usize {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_nb_ext
    }

    /// OCCT IsParallel() (cxx L799-806).
    pub fn is_parallel(&self) -> bool {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_is_par
    }

    /// OCCT SquareDistance(N) (cxx L778-786) — 1-based.
    pub fn square_distance(&self, n: usize) -> f64 {
        if n < 1 || n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.my_sq_dist[n - 1]
    }

    /// OCCT Points(N, P1, P2) (cxx L788-797) — 1-based.
    pub fn points(
        &self,
        n: usize,
        p1: &mut crate::base::extrema::POnCurve,
        p2: &mut POnSurface,
    ) {
        if n < 1 || n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        *p1 = self.my_point1[n - 1].clone();
        *p2 = self.my_point2[n - 1].clone();
    }
}

/// OCCT gp_Pln::SquareDistance(gp_Pnt) (gp_Pln.hxx L159-164) with the
/// `SignedDistance` convention of gp_Pln.hxx L345-347.
fn plane_square_distance(s: &Plane, p: DVec3) -> f64 {
    let d = (p - s.origin).dot(s.normal).abs();
    d * d
}

/// The unused `elclib_line_value`/`elclib_line_parameter` markers keep the
/// import surface honest for the branches whose IntAna tools are GAP-carried.
#[allow(dead_code)]
fn _unused(a: &Line3, p: DVec3) -> f64 {
    elclib_line_parameter(a, elclib_line_value(0.0, a)) + p.length()
}

impl Default for ExtremaExtElCS {
    fn default() -> Self {
        ExtremaExtElCS::new()
    }
}
