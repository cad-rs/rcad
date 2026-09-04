// OCCT Contap_ArcFunction (TKHLR) — F(t) = f(P(t)) on a restriction arc
// of the surface (math_FunctionWithDerivative).
//
// Contap_ArcFunction.hxx L29-84 + .cxx L24-209 + .lxx L17-62.

use rcad_kernel::geom::Point3;
use rcad_kernel::math::root::{FunctionValue, FunctionWithDerivative};

use crate::geomalgo::int_patch::transitions::Transition;
use crate::geomalgo::int_surf::quadric::Quadric;

use crate::hlr::contap::h_cont_tool as hcont;
use crate::hlr::contap::h_curve2d_tool as hcurve2d;
use crate::hlr::contap::point::Arc;
use crate::hlr::contap::surf_props;
use crate::hlr::contap::surface_adaptor::SurfaceHandle;
use crate::hlr::contap::t_function::TFunction;

/// OCCT Contap_ArcFunction.
#[derive(Clone)]
pub struct ArcFunction {
    my_arc: Option<Arc>,
    my_surf: Option<SurfaceHandle>,
    my_mean: f64,
    my_type: TFunction,
    my_dir: glam::DVec3,
    my_cos_ang: f64,
    my_eye: Point3,
    pub(crate) solpt: Point3,
    seqpt: Vec<Point3>,
    my_quad: Quadric,
}

impl ArcFunction {
    /// OCCT Contap_ArcFunction() (cxx L24-30).
    pub fn new() -> Self {
        ArcFunction {
            my_arc: None,
            my_surf: None,
            my_mean: 1.0,
            my_type: TFunction::ContourStd,
            my_dir: glam::DVec3::new(0.0, 0.0, 1.0), // gp_Dir::D::Z
            my_cos_ang: 0.0,
            my_eye: Point3::ZERO,
            solpt: Point3::ZERO,
            seqpt: Vec::new(),
            my_quad: Quadric::new(),
        }
    }

    /// OCCT Set(const handle<Adaptor3d_Surface>& S) (cxx L32-53).
    pub fn set(&mut self, s: SurfaceHandle) {
        self.my_surf = Some(s.clone());
        let nbs = hcont::nb_sample_points(s.as_ref());
        if nbs > 0 {
            self.my_mean = 0.0;
            for i in 1..=nbs {
                let (u, v) = hcont::sample_point(s.as_ref(), i);
                let mut norm = glam::DVec3::ZERO;
                let mut solpt = Point3::ZERO;
                surf_props::normale(s.as_ref(), u, v, &mut solpt, &mut norm);
                self.solpt = solpt;
                self.my_mean += norm.length();
            }
            self.my_mean /= nbs as f64;
        }
    }

    /// OCCT Set(const gp_Dir& Direction, const double Angle) (lxx L17-22).
    pub fn set_dir_angle(&mut self, direction: glam::DVec3, angle: f64) {
        self.my_type = TFunction::DraftStd;
        self.my_dir = direction;
        self.my_cos_ang = (std::f64::consts::PI / 2.0 + angle).cos();
    }

    /// OCCT Set(const gp_Pnt& Eye, const double Angle) (lxx L24-29).
    pub fn set_eye_angle(&mut self, eye: Point3, angle: f64) {
        self.my_type = TFunction::DraftPrs;
        self.my_eye = eye;
        self.my_cos_ang = (std::f64::consts::PI / 2.0 + angle).cos();
    }

    /// OCCT Set(const gp_Dir& Direction) (lxx L31-35).
    pub fn set_dir(&mut self, direction: glam::DVec3) {
        self.my_type = TFunction::ContourStd;
        self.my_dir = direction;
    }

    /// OCCT Set(const gp_Pnt& Eye) (lxx L37-41).
    pub fn set_eye(&mut self, eye: Point3) {
        self.my_type = TFunction::ContourPrs;
        self.my_eye = eye;
    }

    /// OCCT Set(const handle<Adaptor2d_Curve2d>& A) (lxx L43-47).
    pub fn set_arc(&mut self, a: Arc) {
        self.my_arc = Some(a);
        self.seqpt.clear();
    }

    /// OCCT Value(U, F) (cxx L55-83).
    fn value_f(&mut self, u: f64, f: &mut f64) -> bool {
        let arc = self.my_arc.clone().expect("Contap_ArcFunction::Set(A)");
        let pt2d = hcurve2d::value(arc.as_ref(), u);
        let s = self.my_surf.as_deref().expect("Contap_ArcFunction::Set(S)");
        let mut norm = glam::DVec3::ZERO;
        let mut solpt = Point3::ZERO;
        surf_props::normale(s, pt2d.x, pt2d.y, &mut solpt, &mut norm);
        self.solpt = solpt;

        match self.my_type {
            TFunction::ContourStd => {
                *f = norm.dot(self.my_dir) / self.my_mean;
            }
            TFunction::ContourPrs => {
                *f = norm.dot(self.solpt - self.my_eye) / self.my_mean;
            }
            TFunction::DraftStd => {
                *f = (norm.dot(self.my_dir) - self.my_cos_ang * norm.length()) / self.my_mean;
            }
            _ => {}
        }
        true
    }

    /// OCCT Derivative(U, D) (cxx L85-132).
    fn derivative_f(&mut self, u: f64, d: &mut f64) -> bool {
        let arc = self.my_arc.clone().expect("Contap_ArcFunction::Set(A)");
        let (pt2d, d2d) = hcurve2d::d1(arc.as_ref(), u);
        let s = self.my_surf.as_deref().expect("Contap_ArcFunction::Set(S)");
        let mut norm = glam::DVec3::ZERO;
        let mut dnu = glam::DVec3::ZERO;
        let mut dnv = glam::DVec3::ZERO;
        let mut solpt = Point3::ZERO;
        surf_props::norm_and_dn(s, pt2d.x, pt2d.y, &mut solpt, &mut norm, &mut dnu, &mut dnv);
        self.solpt = solpt;

        let (mut dfu, mut dfv) = (0.0, 0.0);
        match self.my_type {
            TFunction::ContourStd => {
                dfu = dnu.dot(self.my_dir) / self.my_mean;
                dfv = dnv.dot(self.my_dir) / self.my_mean;
            }
            TFunction::ContourPrs => {
                let ep = self.solpt - self.my_eye;
                dfu = dnu.dot(ep) / self.my_mean;
                dfv = dnv.dot(ep) / self.my_mean;
            }
            TFunction::DraftStd => {
                let normunit = norm.normalize();
                dfu = (dnu.dot(self.my_dir) - self.my_cos_ang * dnu.dot(normunit)) / self.my_mean;
                dfv = (dnv.dot(self.my_dir) - self.my_cos_ang * dnv.dot(normunit)) / self.my_mean;
            }
            _ => {}
        }
        *d = d2d.x * dfu + d2d.y * dfv;
        true
    }

    /// OCCT Values(U, F, D) (cxx L134-185).
    fn values_f(&mut self, u: f64, f: &mut f64, d: &mut f64) -> bool {
        let arc = self.my_arc.clone().expect("Contap_ArcFunction::Set(A)");
        let (pt2d, d2d) = hcurve2d::d1(arc.as_ref(), u);
        let s = self.my_surf.as_deref().expect("Contap_ArcFunction::Set(S)");
        let mut norm = glam::DVec3::ZERO;
        let mut dnu = glam::DVec3::ZERO;
        let mut dnv = glam::DVec3::ZERO;
        let mut solpt = Point3::ZERO;
        surf_props::norm_and_dn(s, pt2d.x, pt2d.y, &mut solpt, &mut norm, &mut dnu, &mut dnv);
        self.solpt = solpt;

        let (mut dfu, mut dfv) = (0.0, 0.0);
        match self.my_type {
            TFunction::ContourStd => {
                *f = norm.dot(self.my_dir) / self.my_mean;
                dfu = dnu.dot(self.my_dir) / self.my_mean;
                dfv = dnv.dot(self.my_dir) / self.my_mean;
            }
            TFunction::ContourPrs => {
                let ep = self.solpt - self.my_eye;
                *f = norm.dot(ep) / self.my_mean;
                dfu = dnu.dot(ep) / self.my_mean;
                dfv = dnv.dot(ep) / self.my_mean;
            }
            TFunction::DraftStd => {
                *f = (norm.dot(self.my_dir) - self.my_cos_ang * norm.length()) / self.my_mean;
                let normunit = norm.normalize();
                dfu = (dnu.dot(self.my_dir) - self.my_cos_ang * dnu.dot(normunit)) / self.my_mean;
                dfv = (dnv.dot(self.my_dir) - self.my_cos_ang * dnv.dot(normunit)) / self.my_mean;
            }
            _ => {}
        }

        *d = d2d.x * dfu + d2d.y * dfv;
        true
    }

    /// OCCT NbSamples (cxx L193-198).
    pub fn nb_samples(&self) -> i32 {
        let arc = self.my_arc.as_deref().expect("Contap_ArcFunction::Set(A)");
        let surf = self.my_surf.as_deref().expect("Contap_ArcFunction::Set(S)");
        hcont::nb_samples_u(surf, 0.0, 0.0)
            .max(hcont::nb_samples_v(surf, 0.0, 0.0))
            .max(hcont::nb_samples_on_arc(arc))
    }

    /// OCCT GetStateNumber (cxx L187-191) — appends the computed point and
    /// returns the sequence length.
    pub fn get_state_number(&mut self) -> i32 {
        self.seqpt.push(self.solpt);
        self.seqpt.len() as i32
    }

    /// OCCT Valpoint(Index) (lxx L49-52) — 1-based.
    pub fn valpoint(&self, index: i32) -> Point3 {
        self.seqpt[(index - 1) as usize]
    }

    /// OCCT Quadric (cxx L204-207).
    pub fn quadric(&self) -> &Quadric {
        &self.my_quad
    }

    /// OCCT Surface (lxx L54-57).
    pub fn surface(&self) -> &SurfaceHandle {
        self.my_surf.as_ref().expect("Contap_ArcFunction::Surface")
    }

    /// OCCT LastComputedPoint (lxx L59-62).
    pub fn last_computed_point(&self) -> Point3 {
        self.solpt
    }
}

impl Default for ArcFunction {
    fn default() -> Self {
        Self::new()
    }
}

impl FunctionValue for ArcFunction {
    fn value(&mut self, x: f64) -> Option<f64> {
        let mut f = 0.0;
        if self.value_f(x, &mut f) {
            Some(f)
        } else {
            None
        }
    }
}

impl FunctionWithDerivative for ArcFunction {
    fn derivative(&mut self, x: f64) -> Option<f64> {
        let mut d = 0.0;
        if self.derivative_f(x, &mut d) {
            Some(d)
        } else {
            None
        }
    }
    fn values(&mut self, x: f64) -> Option<(f64, f64)> {
        let mut f = 0.0;
        let mut d = 0.0;
        if self.values_f(x, &mut f, &mut d) {
            Some((f, d))
        } else {
            None
        }
    }
    fn get_state_number(&mut self) -> i32 {
        ArcFunction::get_state_number(self)
    }
}

// The Transition import documents the SetArc partner type used by callers
// (IntSurf_Transition flows through Contap_Point::SetArc).
#[allow(unused)]
fn _transition_shape() -> Transition {
    Transition::new()
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::{DVec2, DVec3};
    use rcad_kernel::geom::{Circle2d, Curve2d, CylindricalSurface, Surface3};

    use crate::hlr::contap::surface_adaptor::GeomSurfaceAdapter;

    /// OCCT anchor: on the unit cylinder with direction +X, the
    /// restriction line u=t (the pcurve t -> (t, v0)) has
    /// F(t) = cos(t); the zero is the silhouette t = pi/2 (cxx L55-83).
    #[test]
    fn arc_function_on_cylinder_silhouette_line() {
        let surf = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
            ref_dir: DVec3::new(1.0, 0.0, 0.0),
            y_dir: None,
        });
        let mut f = ArcFunction::new();
        f.set_dir(DVec3::new(1.0, 0.0, 0.0));
        f.set(std::sync::Arc::new(GeomSurfaceAdapter::with_domain(
            surf,
            [0.0, std::f64::consts::TAU, -1.0, 1.0],
        )));
        // Restriction arc: u-iso line at v = 0.2, parameter = angle t.
        let arc = Curve2d::Line(rcad_kernel::geom::Line2d {
            origin: DVec2::new(0.0, 0.2),
            direction: DVec2::new(1.0, 0.0),
        });
        // NB: the pcurve coordinates are (u, v) in the (origin, direction)
        // frame of the Line2d — a horizontal line at v = 0.2 evaluated as
        // P(t) = (t, 0.2).
        f.set_arc(std::sync::Arc::new(arc));

        let v0 = f.value(0.0).unwrap();
        assert!((v0 - 1.0).abs() < 1e-12, "F(0)={}", v0);
        let vs = f.value(std::f64::consts::FRAC_PI_2).unwrap();
        assert!(vs.abs() < 1e-12, "F(pi/2)={}", vs);

        // Derivative F'(t) = -sin(t): at pi/2 it is -1.
        let d = f.derivative(std::f64::consts::FRAC_PI_2).unwrap();
        assert!((d + 1.0).abs() < 1e-9, "D={}", d);

        // GetStateNumber stores the last computed point.
        let _ = f.value(0.0);
        let n = f.get_state_number();
        assert_eq!(n, 1);
        assert_eq!(f.valpoint(1), f.last_computed_point());
    }
}
