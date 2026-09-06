// OCCT Contap_SurfFunction (TKHLR) — F(u,v) on a parametric surface for
// the contour problem (math_FunctionSetWithDerivatives, 2 variables,
// 1 equation).
//
// Contap_SurfFunction.hxx L35-118 + .cxx L29-287 + .lxx L19-106.
// The OCCT `occ::handle<Adaptor3d_Surface> mySurf` maps to
// `Option<SurfaceHandle>`; the sampling of Set() goes through
// Contap_HContTool and Contap_SurfProps exactly as in OCCT.

use rcad_kernel::geom::Point3;

use crate::hlr::contap::h_cont_tool as hcont;
use crate::hlr::contap::surf_props;
use crate::hlr::contap::surface_adaptor::SurfaceHandle;
use crate::hlr::contap::t_function::TFunction;

/// OCCT gp::Resolution().
const GP_RESOLUTION: f64 = 1e-15;

/// OCCT Contap_SurfFunction.
#[derive(Clone)]
pub struct SurfFunction {
    my_surf: Option<SurfaceHandle>,
    my_mean: f64,
    my_type: TFunction,
    my_dir: Vec3Shape,
    my_eye: Point3,
    my_ang: f64,
    my_cos_ang: f64,
    tol: f64,
    pub(crate) solpt: Point3,
    valf: f64,
    pub(crate) usol: f64,
    pub(crate) vsol: f64,
    fpu: f64,
    fpv: f64,
    d2d: (f64, f64),
    d3d: Vec3Shape,
    tangent: bool,
    computed: bool,
    derived: bool,
}

// rcad keeps OCCT gp_Dir/gp_Vec semantics with plain DVec3 (unit invariant
// is preserved by construction at the Set sites).
type Vec3Shape = glam::DVec3;

impl SurfFunction {
    /// OCCT Contap_SurfFunction() (cxx L29-45).
    pub fn new() -> Self {
        SurfFunction {
            my_surf: None,
            my_mean: 1.0,
            my_type: TFunction::ContourStd,
            my_dir: glam::DVec3::new(0.0, 0.0, 1.0), // gp_Dir::D::Z
            my_eye: Point3::ZERO,
            my_ang: 0.0,
            my_cos_ang: 0.0, // PI/2 - Angle de depouille
            tol: 1.0e-6,
            solpt: Point3::ZERO,
            valf: 0.0,
            usol: 0.0,
            vsol: 0.0,
            fpu: 0.0,
            fpv: 0.0,
            d2d: (0.0, 0.0),
            d3d: Vec3Shape::ZERO,
            tangent: false,
            computed: false,
            derived: false,
        }
    }

    /// OCCT Set(const handle<Adaptor3d_Surface>& S) (cxx L47-69).
    pub fn set(&mut self, s: SurfaceHandle) {
        self.my_surf = Some(s.clone());
        let nbs = hcont::nb_sample_points(s.as_ref());
        if nbs > 0 {
            self.my_mean = 0.0;
            for i in 1..=nbs {
                let (u, v) = hcont::sample_point(s.as_ref(), i);
                let mut norm = Vec3Shape::ZERO;
                let mut solpt = Point3::ZERO;
                surf_props::normale(s.as_ref(), u, v, &mut solpt, &mut norm);
                self.solpt = solpt;
                self.my_mean += norm.length();
            }
            self.my_mean /= nbs as f64;
        }
        self.computed = false;
        self.derived = false;
    }

    /// OCCT Set(const gp_Pnt& Eye) (lxx L19-24).
    pub fn set_eye(&mut self, eye: Point3) {
        self.my_type = TFunction::ContourPrs; // pers
        self.my_eye = eye;
        self.my_ang = 0.0;
    }

    /// OCCT Set(const gp_Dir& Direction) (lxx L26-31).
    pub fn set_dir(&mut self, direction: Vec3Shape) {
        self.my_type = TFunction::ContourStd; // Contour app
        self.my_dir = direction;
        self.my_ang = 0.0;
    }

    /// OCCT Set(const gp_Dir& Direction, const double Angle) (lxx L33-39).
    pub fn set_dir_angle(&mut self, direction: Vec3Shape, angle: f64) {
        self.my_type = TFunction::DraftStd; // Contour vu
        self.my_dir = direction;
        self.my_ang = angle;
        self.my_cos_ang = (std::f64::consts::PI / 2.0 + angle).cos();
    }

    /// OCCT Set(const gp_Pnt& Eye, const double Angle) (lxx L41-47).
    pub fn set_eye_angle(&mut self, eye: Point3, angle: f64) {
        self.my_type = TFunction::DraftPrs; // Contour vu "conique"...
        self.my_eye = eye;
        self.my_ang = angle;
        self.my_cos_ang = (std::f64::consts::PI / 2.0 + angle).cos();
    }

    /// OCCT Set(const double Tolerance) (lxx L49-52).
    pub fn set_tolerance(&mut self, tolerance: f64) {
        self.tol = tolerance.max(1.0e-12);
    }

    /// OCCT NbVariables (cxx L71-74).
    pub fn nb_variables(&self) -> i32 {
        2
    }

    /// OCCT NbEquations (cxx L76-79).
    pub fn nb_equations(&self) -> i32 {
        1
    }

    /// OCCT Value (cxx L81-109) — X is the [u, v] pair, F the [f] value.
    pub fn value(&mut self, x: &[f64; 2]) -> Option<f64> {
        self.usol = x[0];
        self.vsol = x[1];
        let s = self.my_surf.as_deref().expect("Contap_SurfFunction::Set");
        let mut norm = Vec3Shape::ZERO;
        let mut solpt = Point3::ZERO;
        surf_props::normale(s, self.usol, self.vsol, &mut solpt, &mut norm);
        self.solpt = solpt;
        match self.my_type {
            TFunction::ContourStd => {
                self.valf = norm.dot(self.my_dir) / self.my_mean;
            }
            TFunction::ContourPrs => {
                self.valf = norm.dot(self.solpt - self.my_eye) / self.my_mean;
            }
            TFunction::DraftStd => {
                self.valf = (norm.dot(self.my_dir) - self.my_cos_ang * norm.length()) / self.my_mean;
            }
            _ => {}
        }
        self.computed = false;
        self.derived = false;
        Some(self.valf)
    }

    /// OCCT Derivatives (cxx L111-156) — Grad is the [1 x 2] row.
    pub fn derivatives(&mut self, x: &[f64; 2]) -> Option<[f64; 2]> {
        self.usol = x[0];
        self.vsol = x[1];
        let s = self.my_surf.as_deref().expect("Contap_SurfFunction::Set");
        let mut norm = Vec3Shape::ZERO;
        let mut dnu = Vec3Shape::ZERO;
        let mut dnv = Vec3Shape::ZERO;
        let mut solpt = Point3::ZERO;
        surf_props::norm_and_dn(s, self.usol, self.vsol, &mut solpt, &mut norm, &mut dnu, &mut dnv);
        self.solpt = solpt;

        match self.my_type {
            TFunction::ContourStd => {
                self.fpu = dnu.dot(self.my_dir) / self.my_mean;
                self.fpv = dnv.dot(self.my_dir) / self.my_mean;
            }
            TFunction::ContourPrs => {
                let ep = self.solpt - self.my_eye;
                self.fpu = dnu.dot(ep) / self.my_mean;
                self.fpv = dnv.dot(ep) / self.my_mean;
            }
            TFunction::DraftStd => {
                let normunit = norm.normalize();
                self.fpu = (dnu.dot(self.my_dir) - self.my_cos_ang * dnu.dot(normunit)) / self.my_mean;
                self.fpv = (dnv.dot(self.my_dir) - self.my_cos_ang * dnv.dot(normunit)) / self.my_mean;
            }
            _ => {}
        }
        let grad = [self.fpu, self.fpv];
        self.computed = false;
        self.derived = true;
        Some(grad)
    }

    /// OCCT Values (cxx L158-212) — returns (F(1), Grad(1,1), Grad(1,2)).
    pub fn values(&mut self, x: &[f64; 2]) -> Option<(f64, [f64; 2])> {
        self.usol = x[0];
        self.vsol = x[1];
        let s = self.my_surf.as_deref().expect("Contap_SurfFunction::Set");
        let mut norm = Vec3Shape::ZERO;
        let mut dnu = Vec3Shape::ZERO;
        let mut dnv = Vec3Shape::ZERO;
        let mut solpt = Point3::ZERO;
        surf_props::norm_and_dn(s, self.usol, self.vsol, &mut solpt, &mut norm, &mut dnu, &mut dnv);
        self.solpt = solpt;

        let f;
        match self.my_type {
            TFunction::ContourStd => {
                f = norm.dot(self.my_dir) / self.my_mean;
                self.fpu = dnu.dot(self.my_dir) / self.my_mean;
                self.fpv = dnv.dot(self.my_dir) / self.my_mean;
            }
            TFunction::ContourPrs => {
                let ep = self.solpt - self.my_eye;
                f = norm.dot(ep) / self.my_mean;
                self.fpu = dnu.dot(ep) / self.my_mean;
                self.fpv = dnv.dot(ep) / self.my_mean;
            }
            TFunction::DraftStd => {
                f = (norm.dot(self.my_dir) - self.my_cos_ang * norm.length()) / self.my_mean;
                let normunit = norm.normalize();
                self.fpu = (dnu.dot(self.my_dir) - self.my_cos_ang * dnu.dot(normunit)) / self.my_mean;
                self.fpv = (dnv.dot(self.my_dir) - self.my_cos_ang * dnv.dot(normunit)) / self.my_mean;
            }
            _ => {
                f = self.valf;
            }
        }
        self.valf = f;
        let grad = [self.fpu, self.fpv];
        self.computed = false;
        self.derived = true;
        Some((f, grad))
    }

    /// OCCT IsTangent (cxx L214-287) — also refreshes d2d/d3d when not
    /// tangent.
    pub fn is_tangent(&mut self) -> bool {
        if !self.computed {
            self.computed = true;
            if !self.derived {
                let s = self.my_surf.as_deref().expect("Contap_SurfFunction::Set");
                let mut norm = Vec3Shape::ZERO;
                let mut dnu = Vec3Shape::ZERO;
                let mut dnv = Vec3Shape::ZERO;
                let mut solpt = Point3::ZERO;
                surf_props::norm_and_dn(
                    s,
                    self.usol,
                    self.vsol,
                    &mut solpt,
                    &mut norm,
                    &mut dnu,
                    &mut dnv,
                );
                self.solpt = solpt;

                match self.my_type {
                    TFunction::ContourStd => {
                        self.fpu = dnu.dot(self.my_dir) / self.my_mean;
                        self.fpv = dnv.dot(self.my_dir) / self.my_mean;
                    }
                    TFunction::ContourPrs => {
                        let ep = self.solpt - self.my_eye;
                        self.fpu = dnu.dot(ep) / self.my_mean;
                        self.fpv = dnv.dot(ep) / self.my_mean;
                    }
                    TFunction::DraftStd => {
                        let normunit = norm.normalize();
                        self.fpu =
                            (dnu.dot(self.my_dir) - self.my_cos_ang * dnu.dot(normunit)) / self.my_mean;
                        self.fpv =
                            (dnv.dot(self.my_dir) - self.my_cos_ang * dnv.dot(normunit)) / self.my_mean;
                    }
                    _ => {}
                }
                self.derived = true;
            }
            self.tangent = false;
            let d = (self.fpu * self.fpu + self.fpv * self.fpv).sqrt();

            if d <= GP_RESOLUTION {
                self.tangent = true;
            } else {
                // OCCT cxx L271: `d2d = gp_Dir2d(-Fpv, Fpu)` — the d2d member
                // is a gp_Dir2d (hxx L113), normalized by construction.
                self.d2d = (-self.fpv / d, self.fpu / d);
                let s = self.my_surf.as_deref().expect("Contap_SurfFunction::Set");
                let (pt, d1u, d1v) = s.d1(self.usol, self.vsol); // ajout jag 02.95
                self.solpt = pt;

                // gp_XYZ d3dxyz(-Fpv * d1u.XYZ()); d3dxyz.Add(Fpu * d1v.XYZ());
                self.d3d = -self.fpv * d1u + self.fpu * d1v;

                // jag 940616    if (d3d.Magnitude() <= Tolpetit) {
                if self.d3d.length() <= self.tol {
                    self.tangent = true;
                }
            }
        }
        self.tangent
    }

    /// OCCT Root (lxx L59-62) — the value of the function at the solution.
    pub fn root(&self) -> f64 {
        self.valf
    }

    /// OCCT Tolerance (lxx L64-67).
    pub fn tolerance(&self) -> f64 {
        self.tol
    }

    /// OCCT Point (lxx L54-57) — the solution point on the surface.
    pub fn point(&self) -> Point3 {
        self.solpt
    }

    /// OCCT Direction3d (lxx L69-74).
    pub fn direction_3d(&mut self) -> Vec3Shape {
        if self.is_tangent() {
            panic!("StdFail_UndefinedDerivative");
        }
        self.d3d
    }

    /// OCCT Direction2d (lxx L76-81).
    pub fn direction_2d(&mut self) -> glam::DVec2 {
        if self.is_tangent() {
            panic!("StdFail_UndefinedDerivative");
        }
        self.d2d.into()
    }

    /// OCCT Surface (lxx L83-86).
    pub fn surface(&self) -> &SurfaceHandle {
        self.my_surf.as_ref().expect("Contap_SurfFunction::Surface")
    }

    /// OCCT PSurface — method entered for compatibility with
    /// IntPatch_TheSurfFunction; returns Surface() (lxx L93-96).
    pub fn p_surface(&self) -> &SurfaceHandle {
        self.surface()
    }

    /// OCCT Eye (lxx L88-91).
    pub fn eye(&self) -> Point3 {
        self.my_eye
    }

    /// OCCT Direction (lxx L93-96).
    pub fn direction(&self) -> Vec3Shape {
        self.my_dir
    }

    /// OCCT Angle (lxx L98-101).
    pub fn angle(&self) -> f64 {
        self.my_ang
    }

    /// OCCT FunctionType (lxx L103-106).
    pub fn function_type(&self) -> TFunction {
        self.my_type
    }
}

impl Default for SurfFunction {
    fn default() -> Self {
        Self::new()
    }
}

impl crate::geomalgo::int_patch::imp_prm::function_set_root::FunctionSetWithDerivatives2
    for SurfFunction
{
    fn value(&mut self, x: &[f64; 2]) -> Option<f64> {
        SurfFunction::value(self, x)
    }
    fn values(&mut self, x: &[f64; 2]) -> Option<(f64, [f64; 2])> {
        SurfFunction::values(self, x)
    }
    fn tolerance(&self) -> f64 {
        SurfFunction::tolerance(self)
    }
}

impl crate::geomalgo::int_patch::imp_prm::i_walking::IWFunction for SurfFunction {
    fn derivatives(&mut self, x: &[f64; 2]) -> Option<[f64; 2]> {
        SurfFunction::derivatives(self, x)
    }
    fn root(&self) -> f64 {
        SurfFunction::root(self)
    }
    fn is_tangent(&mut self) -> bool {
        SurfFunction::is_tangent(self)
    }
    fn direction_3d(&mut self) -> glam::DVec3 {
        SurfFunction::direction_3d(self)
    }
    fn direction_2d(&mut self) -> glam::DVec2 {
        SurfFunction::direction_2d(self)
    }
    fn point(&self) -> glam::DVec3 {
        SurfFunction::point(self)
    }
    /// OCCT Func.Set(Caro) — the Contap function is bound to the adaptor
    /// (Contap_Contour.cxx L152), and the walking re-Set (gxx L247) hits
    /// the same surface: the sampling state (myMean) is recomputed from
    /// the stored adapter, exactly the OCCT result.
    fn set_surface(&mut self, _s: &rcad_kernel::geom::Surface3) {
        let surf = self.my_surf.clone().expect("Contap_SurfFunction::Set");
        self.set(surf);
    }
    fn surface_value(&self, u: f64, v: f64) -> glam::DVec3 {
        self.surface().value(u, v)
    }
    /// OCCT ThePSurfaceTool::UResolution(Func.PSurface(), R3d) — the
    /// analytic resolution of the bound adaptor.
    fn u_resolution(&self, r3d: f64) -> f64 {
        self.surface().u_resolution(r3d)
    }
    /// OCCT ThePSurfaceTool::VResolution(Func.PSurface(), R3d).
    fn v_resolution(&self, r3d: f64) -> f64 {
        self.surface().v_resolution(r3d)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::{CylindricalSurface, Surface3, ToroidalSurface};

    use crate::hlr::contap::surface_adaptor::GeomSurfaceAdapter;

    /// TEMP-DEBUG: probe the solver from the seam point.
    #[test]
    fn temp_probe_solver_seam() {
        use crate::geomalgo::int_patch::imp_prm::function_set_root::FunctionSetRoot;
        let tor = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::new(0.0, 0.0, 1.0),
            major_radius: 30.0,
            minor_radius: 10.0,
            ref_dir: DVec3::new(1.0, 0.0, 0.0),
        };
        let mut f = SurfFunction::new();
        f.set_dir(DVec3::new(1.0, -1.0, 1.0).normalize());
        f.set(std::sync::Arc::new(GeomSurfaceAdapter::new(
            Surface3::Torus(tor),
        )));
        f.set_tolerance(1e-7);
        let mut solver = FunctionSetRoot::new(&mut f, [1.5915494309189534e-8, 1.5915494309189534e-8]);
        let start = [0.0, 5.497515982079354];
        solver.perform(&mut f, start, [0.0, 0.0], [6.283185307179586, 6.283185307179586]);
        println!(
            "SOLVER done={} root={:?} F={:.3e}",
            solver.is_done(),
            solver.root(),
            f.root()
        );
        let (val, grad) = f.values(&[0.0, 5.497515982079354]).unwrap();
        println!("SOLVER F(start)={val:.3e} grad={grad:?}");
    }

    /// TEMP-DEBUG: probe the torus gradient at the departure point.
    #[test]
    fn temp_probe_torus_gradient() {
        let tor = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::new(0.0, 0.0, 1.0),
            major_radius: 30.0,
            minor_radius: 10.0,
            ref_dir: DVec3::new(1.0, 0.0, 0.0),
        };
        let adapter = GeomSurfaceAdapter::new(Surface3::Torus(tor));
        let mut p = Point3::ZERO;
        let mut n = Vec3Shape::ZERO;
        let mut dnu = Vec3Shape::ZERO;
        let mut dnv = Vec3Shape::ZERO;
        super::surf_props::norm_and_dn(
            &adapter,
            0.7853981633974486,
            6.283185307179586,
            &mut p,
            &mut n,
            &mut dnu,
            &mut dnv,
        );
        println!("PROBE p={p:?} n={n:?} dnu={dnu:?} dnv={dnv:?}");
        // finite differences for reference
        let h = 1e-4;
        let u0 = 0.7853981633974486;
        let v0 = 6.283185307179586;
        let (p0, du_fd, dv_fd) = rcad_kernel::geom::SurfaceEval::derivatives(
            &rcad_kernel::geom::Surface3::Torus(ToroidalSurface {
                center: DVec3::ZERO,
                axis: DVec3::new(0.0, 0.0, 1.0),
                major_radius: 30.0,
                minor_radius: 10.0,
                ref_dir: DVec3::new(1.0, 0.0, 0.0),
            }),
            u0,
            v0,
        );
        let _ = p0;
        println!("PROBE d1 du={du_fd:?} dv={dv_fd:?}");
        let (_, _, _, puu, puv, pvv) = rcad_kernel::geom::SurfaceEval::derivatives2(
            &rcad_kernel::geom::Surface3::Torus(ToroidalSurface {
                center: DVec3::ZERO,
                axis: DVec3::new(0.0, 0.0, 1.0),
                major_radius: 30.0,
                minor_radius: 10.0,
                ref_dir: DVec3::new(1.0, 0.0, 0.0),
            }),
            u0,
            v0,
        );
        println!("PROBE d2 puu={puu:?} puv={puv:?} pvv={pvv:?}");
        let mut f = SurfFunction::new();
        f.set_dir(DVec3::new(1.0, -1.0, 1.0).normalize());
        f.set(std::sync::Arc::new(GeomSurfaceAdapter::new(
            Surface3::Torus(tor),
        )));
        let (val, grad) = f.values(&[0.7853981633974486, 6.283185307179586]).unwrap();
        println!("PROBE F={val} grad={grad:?} mean={}", f.my_mean);
        let d = (grad[0] * grad[0] + grad[1] * grad[1]).sqrt();
        println!("PROBE d2d=({},{})", -grad[1] / d, grad[0] / d);
        assert!(!f.is_tangent());
        let d2d = f.direction_2d();
        let d3d = f.direction_3d();
        println!("PROBE d2d={d2d:?} d3d={d3d:?}");
    }

    /// OCCT anchor: on a unit cylinder along Z with the contour direction
    /// +X, F(u,v) = cos(u); the silhouette is at u = pi/2.  myMean from
    /// the 5 HContTool samples on a finite window equals 1 (normals are
    /// unit), so F is exactly cos(u) (cxx L47-69, L81-109).
    #[test]
    fn surf_function_cylinder_values() {
        let surf = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
            ref_dir: DVec3::new(1.0, 0.0, 0.0),
            y_dir: None,
        });
        let mut f = SurfFunction::new();
        f.set_dir(DVec3::new(1.0, 0.0, 0.0));
        f.set(std::sync::Arc::new(GeomSurfaceAdapter::with_domain(
            surf,
            [0.0, std::f64::consts::TAU, -1.0, 1.0],
        )));
        assert!((f.my_mean - 1.0).abs() < 1e-12);

        // F(pi/2, 0) = cos(pi/2) = 0 — the silhouette.
        let v = f.value(&[std::f64::consts::FRAC_PI_2, 0.0]).unwrap();
        assert!(v.abs() < 1e-12, "F={}", v);

        // F(0, 0) = 1.
        let v = f.value(&[0.0, 0.0]).unwrap();
        assert!((v - 1.0).abs() < 1e-12);

        // Values at (pi/2, 0): gradient = (-sin u, 0)/myMean * ... =>
        // Grad(1,1) = dnu.Dir = -sin(u) = -1.
        let (_, grad) = f.values(&[std::f64::consts::FRAC_PI_2, 0.0]).unwrap();
        assert!((grad[0] + 1.0).abs() < 1e-12, "g={}", grad[0]);
        assert!(grad[1].abs() < 1e-12);

        // Not tangent on the silhouette: the 3d tangent -Fpv*d1u + Fpu*d1v
        // = d1v = Z (magnitude 1 > tol).
        assert!(!f.is_tangent());
        let dir2d = f.direction_2d();
        assert!(dir2d.x.abs() < 1e-12 && (dir2d.y + 1.0).abs() < 1e-12);
    }
}
