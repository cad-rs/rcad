// OCCT Contap_Contour (TKHLR) — the contour engine for one surface.
//
// Contap_Contour.hxx L32-116 + .lxx L19-54 + .cxx L46-239 (ctors, Init,
// Perform dispatch) — the body lives in the submodule files:
//   - [`functions`]     the file-static helpers (Recadre, LineConstructor,
//     ComputeTangency, ProcessSegments, ComputeInternalPoints...).
//   - [`perform`]       Perform(Domain) — the walking path (cxx L1541-1971).
//   - [`perform_ana`]   PerformAna(Domain) — the analytic path
//     (cxx L2156-2389).
//
// File constants (cxx L46-48): Tolpetit = 1.0e-10 (versus 1.e-8 in
// Contap_ContAna — the conflicting constants of the plan §5 are kept
// verbatim) and tole = 5.0e-6.

mod functions;
mod perform;
mod perform_ana;

use perform::perform_domain;
use perform_ana::perform_ana;

use rcad_kernel::geom::Point3;

use crate::hlr::contap::arc_function::ArcFunction;
use crate::hlr::contap::domain::ContapDomain;
use crate::hlr::contap::line::Line;
use crate::hlr::contap::surf_function::SurfFunction;
use crate::hlr::contap::the_search::TheSearch;
use crate::hlr::contap::the_search_inside::TheSearchInside;

/// OCCT Contap_Contour.cxx L46: `static const double Tolpetit = 1.0e-10;`
pub(crate) const TOLPETIT: f64 = 1.0e-10;
/// OCCT Contap_Contour.cxx L48: `static const double tole = 5.0e-6;`
pub(crate) const TOLE: f64 = 5.0e-6;

/// OCCT Contap_Contour.
#[derive(Clone)]
pub struct Contour {
    done: bool,
    slin: Vec<Line>,
    solrst: TheSearch,
    solins: TheSearchInside,
    my_sfunc: SurfFunction,
    my_afunc: ArcFunction,
    modeset: bool,
}

impl Contour {
    /// OCCT Contap_Contour() (cxx L50-54).
    pub fn new() -> Self {
        Contour {
            done: false,
            slin: Vec::new(),
            solrst: TheSearch::new(),
            solins: TheSearchInside::new(),
            my_sfunc: SurfFunction::new(),
            my_afunc: ArcFunction::new(),
            modeset: false,
        }
    }

    /// OCCT Contap_Contour(const gp_Vec& Direction) (cxx L56-64).
    pub fn with_direction(direction: glam::DVec3) -> Self {
        let mut c = Contour::new();
        c.done = false;
        c.modeset = true;
        c.my_sfunc.set_dir(direction);
        c.my_afunc.set_dir(direction);
        c
    }

    /// OCCT Contap_Contour(const gp_Vec& Direction, const double Angle)
    /// (cxx L66-74).
    pub fn with_direction_angle(direction: glam::DVec3, angle: f64) -> Self {
        let mut c = Contour::new();
        c.done = false;
        c.modeset = true;
        c.my_sfunc.set_dir_angle(direction, angle);
        c.my_afunc.set_dir_angle(direction, angle);
        c
    }

    /// OCCT Contap_Contour(const gp_Pnt& Eye) (cxx L76-84).
    pub fn with_eye(eye: Point3) -> Self {
        let mut c = Contour::new();
        c.done = false;
        c.modeset = true;
        c.my_sfunc.set_eye(eye);
        c.my_afunc.set_eye(eye);
        c
    }

    /// OCCT Init(const gp_Vec& Direction) (cxx L120-127).
    pub fn init(&mut self, direction: glam::DVec3) {
        self.done = false;
        self.modeset = true;
        self.my_sfunc.set_dir(direction);
        self.my_afunc.set_dir(direction);
    }

    /// OCCT Init(const gp_Vec& Direction, const double Angle) (cxx L129-135).
    pub fn init_angle(&mut self, direction: glam::DVec3, angle: f64) {
        self.done = false;
        self.modeset = true;
        self.my_sfunc.set_dir_angle(direction, angle);
        self.my_afunc.set_dir_angle(direction, angle);
    }

    /// OCCT Init(const gp_Pnt& Eye) (cxx L137-143).
    pub fn init_eye(&mut self, eye: Point3) {
        self.done = false;
        self.modeset = true;
        self.my_sfunc.set_eye(eye);
        self.my_afunc.set_eye(eye);
    }

    /// OCCT Perform(Surf, Domain) (cxx L145-171) — dispatch on the surface
    /// type.
    pub fn perform(
        &mut self,
        surf: &crate::hlr::contap::surface_adaptor::SurfaceHandle,
        domain: &mut dyn ContapDomain,
    ) {
        self.my_sfunc.set(surf.clone());
        self.my_afunc.set(surf.clone());

        let typ_s = surf.get_type();
        match typ_s {
            crate::geomalgo::int_patch::GeomAbsSurfaceType::Plane
            | crate::geomalgo::int_patch::GeomAbsSurfaceType::Sphere
            | crate::geomalgo::int_patch::GeomAbsSurfaceType::Cylinder
            | crate::geomalgo::int_patch::GeomAbsSurfaceType::Cone => {
                perform_ana(self, domain); // Surf,Domain,Direction,0.,gp_Pnt(0.,0.,0.),1);
            }
            _ => {
                perform_domain(self, domain); // Surf,Domain,Direction,0.,gp_Pnt(0.,0.,0.),1);
            }
        }
    }

    /// OCCT Perform(Surf, Domain, Direction) (cxx L173-180).
    pub fn perform_direction(
        &mut self,
        surf: &crate::hlr::contap::surface_adaptor::SurfaceHandle,
        domain: &mut dyn ContapDomain,
        direction: glam::DVec3,
    ) {
        self.init(direction);
        self.perform(surf, domain);
    }

    /// OCCT Perform(Surf, Domain, Direction, Angle) (cxx L182-190).
    pub fn perform_direction_angle(
        &mut self,
        surf: &crate::hlr::contap::surface_adaptor::SurfaceHandle,
        domain: &mut dyn ContapDomain,
        direction: glam::DVec3,
        angle: f64,
    ) {
        self.init_angle(direction, angle);
        self.perform(surf, domain);
    }

    /// OCCT Perform(Surf, Domain, Eye) (cxx L192-199).
    pub fn perform_eye(
        &mut self,
        surf: &crate::hlr::contap::surface_adaptor::SurfaceHandle,
        domain: &mut dyn ContapDomain,
        eye: Point3,
    ) {
        self.init_eye(eye);
        self.perform(surf, domain);
    }

    /// OCCT IsDone (.lxx L19-22).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT IsEmpty (.lxx L24-27).
    pub fn is_empty(&self) -> bool {
        self.nb_lines() == 0
    }

    /// OCCT NbLines (.lxx L29-35).
    pub fn nb_lines(&self) -> usize {
        if !self.done {
            panic!("StdFail_NotDone: Contap_Contour::NbLines");
        }
        self.slin.len()
    }

    /// OCCT Line(Index) (.lxx L38-45) — 1-based.
    pub fn line(&self, index: usize) -> &Line {
        if !self.done {
            panic!("StdFail_NotDone: Contap_Contour::Line");
        }
        &self.slin[index - 1]
    }

    /// OCCT SurfaceFunction (.lxx L47-54) — a reference on the internal
    /// SurfaceFunction, used to compute tangents on the lines.
    pub fn surface_function(&mut self) -> &mut SurfFunction {
        if !self.done {
            panic!("StdFail_NotDone: Contap_Contour::SurfaceFunction");
        }
        &mut self.my_sfunc
    }

    // ---- internal (private-member) accessors for the submodules ----

    pub(crate) fn slin(&self) -> &Vec<Line> {
        &self.slin
    }
    pub(crate) fn slin_mut(&mut self) -> &mut Vec<Line> {
        &mut self.slin
    }
    pub(crate) fn solrst(&self) -> &TheSearch {
        &self.solrst
    }
    pub(crate) fn solrst_mut(&mut self) -> &mut TheSearch {
        &mut self.solrst
    }
    pub(crate) fn solins(&self) -> &TheSearchInside {
        &self.solins
    }
    pub(crate) fn solins_mut(&mut self) -> &mut TheSearchInside {
        &mut self.solins
    }
    pub(crate) fn my_sfunc(&mut self) -> &mut SurfFunction {
        &mut self.my_sfunc
    }
    pub(crate) fn sfunc_ref(&self) -> &SurfFunction {
        &self.my_sfunc
    }
    pub(crate) fn my_afunc(&mut self) -> &mut ArcFunction {
        &mut self.my_afunc
    }
    pub(crate) fn set_done(&mut self, done: bool) {
        self.done = done;
    }

    /// Split borrow of the internal members (the OCCT members are all
    /// mutable through `this`).
    #[allow(clippy::type_complexity)]
    pub(crate) fn parts_mut(
        &mut self,
    ) -> (
        &mut SurfFunction,
        &mut ArcFunction,
        &mut TheSearch,
        &mut TheSearchInside,
        &mut Vec<Line>,
    ) {
        (
            &mut self.my_sfunc,
            &mut self.my_afunc,
            &mut self.solrst,
            &mut self.solins,
            &mut self.slin,
        )
    }
}

impl Default for Contour {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::{CylindricalSurface, Surface3};

    use crate::hlr::contap::geom_tool::GeomTool;
    use crate::hlr::contap::surface_adaptor::GeomSurfaceAdapter;
    use crate::topalgo::adaptor3d::topol_tool::TopolTool;

    /// OCCT anchor: a unit cylinder along Z viewed along +X — the two
    /// silhouette generatrices at y = +/-1 (ContAna cylinder branch),
    /// oriented transitions, plus the restriction-solution treatment of
    /// PerformSolRst (cxx L2156-2389 through LineConstructor).
    #[test]
    fn contour_cylinder_silhouettes() {
        let surf = GeomSurfaceAdapter::with_domain(
            Surface3::Cylinder(CylindricalSurface {
                origin: DVec3::ZERO,
                axis: DVec3::new(0.0, 0.0, 1.0),
                radius: 1.0,
                ref_dir: DVec3::new(1.0, 0.0, 0.0),
                y_dir: None,
            }),
            [0.0, 2.0 * std::f64::consts::PI, -2.0, 2.0],
        );

        let mut contour = Contour::with_direction(DVec3::new(1.0, 0.0, 0.0));
        let mut domain = TopolTool::<'_, _, GeomTool>::new(&surf);
        let handle: crate::hlr::contap::surface_adaptor::SurfaceHandle = std::sync::Arc::new(surf.clone());
        contour.perform(&handle, &mut domain);

        assert!(contour.is_done());
        assert!(contour.nb_lines() >= 2, "nb_lines={}", contour.nb_lines());

        // The two silhouette lines: origin (0, +/-1, 0), direction Z.
        let mut silhouette = 0;
        let mut transitions = Vec::new();
        for i in 1..=contour.nb_lines() {
            let l = contour.line(i);
            if l.type_contour() == crate::hlr::contap::i_type::IType::Lin {
                silhouette += 1;
                let lin = l.line();
                assert!(lin.direction.z.abs() - 1.0 < 1e-9);
                assert!(lin.origin.x.abs() < 1e-9, "ox={}", lin.origin.x);
                assert!((lin.origin.y.abs() - 1.0).abs() < 1e-9, "oy={}", lin.origin.y);
                assert!(l.nb_vertex() >= 2, "nb_vertex={}", l.nb_vertex());
                transitions.push(l.transition_on_s());
            }
        }
        assert_eq!(silhouette, 2, "silhouette lines={}", silhouette);
        // OCCT-literal transitions: the +Y line has det=+1 (Out); the -Y
        // line falls in the `det < RealEpsilon()` arm of
        // ComputeTransitionOnLine (Contap_Contour.cxx L945-950, the
        // "revoir le test jag 940620" quirk) and stays Undecided.
        assert_eq!(
            transitions,
            vec![
                crate::geomalgo::int_patch::transitions::TypeTrans::Out,
                crate::geomalgo::int_patch::transitions::TypeTrans::Undecided
            ]
        );
    }
}
