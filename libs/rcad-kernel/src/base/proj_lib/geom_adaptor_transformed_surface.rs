//! OCCT GeomAdaptor_TransformedSurface (TKG3d/GeomAdaptor) — the 1:1
//! translation of `GeomAdaptor_TransformedSurface.hxx` (L36-261) +
//! `GeomAdaptor_TransformedSurface.cxx` (L34-421).
//!
//! OCCT inheritance `Adaptor3d_Surface <- GeomAdaptor_TransformedSurface`
//! becomes composition: the struct owns the [`GeomSurfaceAdaptor`] base and
//! every inherited member is an explicit delegation (the OCCT `mySurf.X()`
//! one-liners of the hxx).  The lazy `mutable std::optional` transformed
//! cache becomes a `RefCell<Option<...>>` (the OCCT mutable cache is
//! non-thread-safe by design, hxx L246-247).

use std::cell::RefCell;
use std::sync::Arc;

use glam::DVec3;

use super::adaptor::{Adaptor3dCurve, Adaptor3dSurface, GeomAbsSurfaceType};
use super::geom_adaptor_curve::GeomCurveAdaptor;
use super::geom_adaptor_surface::GeomSurfaceAdaptor;
use super::proj_lib_projected_curve::Adaptor3dSurfaceGeom;
use crate::geom::{
    ConicalSurface, CylindricalSurface, Plane, SphericalSurface, Surface3, ToroidalSurface,
};
use crate::math::gp::{Trsf, TrsfForm};
use crate::math::GeomAbsShape;

// =========================================================================
// OCCT GeomAdaptor_TransformedSurface (hxx L36-261)
// =========================================================================

/// OCCT GeomAdaptor_TransformedSurface — an adaptor for surfaces with an
/// applied transformation.
#[derive(Debug)]
pub struct TransformedSurfaceAdaptor {
    /// OCCT: GeomAdaptor_Surface mySurf (hxx L258).
    pub my_surf: GeomSurfaceAdaptor,
    /// OCCT: gp_Trsf myTrsf (hxx L259).
    pub my_trsf: Trsf,
    /// OCCT: mutable std::optional<GeomAdaptor_Surface> myTransformedAdaptor
    /// (hxx L260) — the lazy transformed cache.
    pub my_transformed_adaptor: RefCell<Option<GeomSurfaceAdaptor>>,
}

impl Clone for TransformedSurfaceAdaptor {
    /// The OCCT copy semantics (the ShallowCopy member copy of
    /// myTransformedAdaptor, GeomAdaptor_TransformedSurface.cxx L74).
    fn clone(&self) -> Self {
        TransformedSurfaceAdaptor {
            my_surf: self.my_surf.clone(),
            my_trsf: self.my_trsf,
            my_transformed_adaptor: RefCell::new(self.my_transformed_adaptor.borrow().clone()),
        }
    }
}

impl TransformedSurfaceAdaptor {
    /// OCCT GeomAdaptor_TransformedSurface() (cxx L34) — the undefined
    /// surface with identity transformation.
    pub fn new() -> Self {
        TransformedSurfaceAdaptor {
            my_surf: GeomSurfaceAdaptor::empty(),
            my_trsf: Trsf::identity(),
            my_transformed_adaptor: RefCell::new(None),
        }
    }

    /// OCCT GeomAdaptor_TransformedSurface(theSurface, theTrsf) (cxx
    /// L38-45): mySurf(theSurface), myTrsf(theTrsf),
    /// invalidateTransformedCache().
    pub fn with_surface_trsf(the_surface: Surface3, the_trsf: Trsf) -> Self {
        TransformedSurfaceAdaptor {
            my_surf: GeomSurfaceAdaptor::new(the_surface),
            my_trsf: the_trsf,
            my_transformed_adaptor: RefCell::new(None),
        }
    }

    /// OCCT GeomAdaptor_TransformedSurface(S, U1, U2, V1, V2, Trsf, TolU,
    /// TolV) (cxx L49-62).
    pub fn with_window_trsf(
        the_surface: Surface3,
        the_u_first: f64,
        the_u_last: f64,
        the_v_first: f64,
        the_v_last: f64,
        the_trsf: Trsf,
        the_tol_u: f64,
        the_tol_v: f64,
    ) -> Self {
        TransformedSurfaceAdaptor {
            my_surf: GeomSurfaceAdaptor::with_window(
                the_surface,
                the_u_first,
                the_u_last,
                the_v_first,
                the_v_last,
                the_tol_u,
                the_tol_v,
            ),
            my_trsf: the_trsf,
            my_transformed_adaptor: RefCell::new(None),
        }
    }

    /// OCCT ShallowCopy (cxx L66-77).
    pub fn shallow_copy_of(&self) -> Self {
        self.clone()
    }

    /// OCCT Load(theSurface, theTrsf) (cxx L81-87).
    pub fn load(&mut self, the_surface: Surface3, the_trsf: Trsf) {
        self.my_surf.load(the_surface);
        self.my_trsf = the_trsf;
        self.invalidate_transformed_cache();
    }

    /// OCCT Load(theSurface, U1, U2, V1, V2, theTrsf, TolU, TolV) (cxx
    /// L91-103).
    #[allow(clippy::too_many_arguments)]
    pub fn load_with_window(
        &mut self,
        the_surface: Surface3,
        the_u_first: f64,
        the_u_last: f64,
        the_v_first: f64,
        the_v_last: f64,
        the_trsf: Trsf,
        the_tol_u: f64,
        the_tol_v: f64,
    ) {
        self.my_surf.load_with_tols(
            the_surface,
            the_u_first,
            the_u_last,
            the_v_first,
            the_v_last,
            the_tol_u,
            the_tol_v,
        );
        self.my_trsf = the_trsf;
        self.invalidate_transformed_cache();
    }

    /// OCCT SetTrsf(theTrsf) (cxx L107-111).
    pub fn set_trsf(&mut self, the_trsf: Trsf) {
        self.my_trsf = the_trsf;
        self.invalidate_transformed_cache();
    }

    /// OCCT HasTrsf() (hxx L98): myTrsf.Form() != gp_Identity.
    pub fn has_trsf(&self) -> bool {
        self.my_trsf.form != TrsfForm::Identity
    }

    /// OCCT ToleranceU() (hxx L169): mySurf.ToleranceU().
    pub fn tolerance_u(&self) -> f64 {
        self.my_surf.tolerance_u()
    }

    /// OCCT ToleranceV() (hxx L172): mySurf.ToleranceV().
    pub fn tolerance_v(&self) -> f64 {
        self.my_surf.tolerance_v()
    }

    /// OCCT Trsf() (hxx L101).
    pub fn trsf(&self) -> &Trsf {
        &self.my_trsf
    }

    /// OCCT Surface() / AdaptorSurfaceOriginal() (hxx L106 / L109).
    pub fn adaptor_surface_original(&self) -> &GeomSurfaceAdaptor {
        &self.my_surf
    }

    /// OCCT AdaptorSurfaceTransformed() (cxx L404-407).  The OCCT const&
    /// return becomes an owned copy — the RefCell cache cannot lend out a
    /// reference past its guard.
    pub fn adaptor_surface_transformed(&self) -> GeomSurfaceAdaptor {
        self.transformed_adaptor()
    }

    /// OCCT GeomSurfaceOriginal() (hxx L116): mySurf.Surface().
    pub fn geom_surface_original(&self) -> &Surface3 {
        self.my_surf.surface()
    }

    /// OCCT GeomSurfaceTransformed() (cxx L115-124).
    pub fn geom_surface_transformed(&self) -> Surface3 {
        if self.my_trsf.form == TrsfForm::Identity {
            return self.my_surf.surface().clone();
        }
        self.ensure_transformed_cache();
        self.my_transformed_adaptor
            .borrow()
            .as_ref()
            .expect("transformed cache")
            .surface()
            .clone()
    }

    /// OCCT invalidateTransformedCache (cxx L128-131).
    fn invalidate_transformed_cache(&mut self) {
        *self.my_transformed_adaptor.borrow_mut() = None;
    }

    /// OCCT ensureTransformedCache (cxx L135-142).
    fn ensure_transformed_cache(&self) {
        if self.my_trsf.form == TrsfForm::Identity
            || self.my_transformed_adaptor.borrow().is_some()
        {
            return;
        }
        self.init_transformed_cache();
    }

    /// OCCT initTransformedCache (cxx L146-199): the per-type transformed
    /// surface construction.
    fn init_transformed_cache(&self) {
        // OCCT L148-152: the null-surface early out cannot be triggered with
        // the rcad value encoding.
        let is_identity = self.my_trsf.form == TrsfForm::Identity;
        let a_trsf = &self.my_trsf;
        // OCCT L157-190: the type switch.  Each analytic arm transforms the
        // gp payload (Pln/Cylinder/.../Torus .Transformed(aTrsf)); the rcad
        // encoding transforms the kernel payload frame (location apply +
        // direction transform, the gp_Dir::Transform normalization
        // semantics).  The default arm transforms the whole surface.
        let a_surface = match self.my_surf.get_type() {
            GeomAbsSurfaceType::Plane => {
                if is_identity {
                    self.my_surf.surface().clone()
                } else {
                    let Surface3::Plane(pl) = self.my_surf.surface() else {
                        unreachable!()
                    };
                    Surface3::Plane(transform_plane(pl, a_trsf))
                }
            }
            GeomAbsSurfaceType::Cylinder => {
                if is_identity {
                    self.my_surf.surface().clone()
                } else {
                    let Surface3::Cylinder(cy) = self.my_surf.surface() else {
                        unreachable!()
                    };
                    Surface3::Cylinder(CylindricalSurface {
                        origin: a_trsf.apply(cy.origin),
                        axis: a_trsf.transform_dir(cy.axis),
                        radius: cy.radius,
                        ref_dir: a_trsf.transform_dir(cy.ref_dir),
                        y_dir: cy.y_dir.map(|y| a_trsf.transform_dir(y)),
                    })
                }
            }
            GeomAbsSurfaceType::Cone => {
                if is_identity {
                    self.my_surf.surface().clone()
                } else {
                    let Surface3::Cone(co) = self.my_surf.surface() else {
                        unreachable!()
                    };
                    Surface3::Cone(ConicalSurface {
                        apex: a_trsf.apply(co.apex),
                        axis: a_trsf.transform_dir(co.axis),
                        radius: co.radius,
                        half_angle_rad: co.half_angle_rad,
                        ref_dir: a_trsf.transform_dir(co.ref_dir),
                    })
                }
            }
            GeomAbsSurfaceType::Sphere => {
                if is_identity {
                    self.my_surf.surface().clone()
                } else {
                    let Surface3::Sphere(sp) = self.my_surf.surface() else {
                        unreachable!()
                    };
                    Surface3::Sphere(SphericalSurface {
                        center: a_trsf.apply(sp.center),
                        axis: a_trsf.transform_dir(sp.axis),
                        radius: sp.radius,
                        ref_dir: a_trsf.transform_dir(sp.ref_dir),
                    })
                }
            }
            GeomAbsSurfaceType::Torus => {
                if is_identity {
                    self.my_surf.surface().clone()
                } else {
                    let Surface3::Torus(to) = self.my_surf.surface() else {
                        unreachable!()
                    };
                    Surface3::Torus(ToroidalSurface {
                        center: a_trsf.apply(to.center),
                        axis: a_trsf.transform_dir(to.axis),
                        ref_dir: a_trsf.transform_dir(to.ref_dir),
                        major_radius: to.major_radius,
                        minor_radius: to.minor_radius,
                    })
                }
            }
            // OCCT L186-189: default -> mySurf.Surface()->Transformed(aTrsf).
            _ => {
                if is_identity {
                    self.my_surf.surface().clone()
                } else {
                    crate::geom::transform_surface(self.my_surf.surface(), &a_trsf.to_daffine3())
                }
            }
        };

        // OCCT L192-198: the windowed adaptor over the transformed surface.
        *self.my_transformed_adaptor.borrow_mut() = Some(GeomSurfaceAdaptor::with_window(
            a_surface,
            self.my_surf.first_u_parameter(),
            self.my_surf.last_u_parameter(),
            self.my_surf.first_v_parameter(),
            self.my_surf.last_v_parameter(),
            self.my_surf.tolerance_u(),
            self.my_surf.tolerance_v(),
        ));
    }

    /// OCCT transformedAdaptor (cxx L411-420): the original adaptor for the
    /// identity transformation, the transformed cache otherwise.  The OCCT
    /// const& return becomes an owned copy — the RefCell cache cannot lend
    /// out a reference past its guard.
    fn transformed_adaptor(&self) -> GeomSurfaceAdaptor {
        if self.my_trsf.form == TrsfForm::Identity {
            return self.my_surf.clone();
        }
        self.ensure_transformed_cache();
        self.my_transformed_adaptor
            .borrow()
            .as_ref()
            .expect("transformed cache")
            .clone()
    }

    /// OCCT UIntervals (cxx L203-207): mySurf.UIntervals(theT, theS).
    pub fn u_intervals_of(&self, s: crate::math::GeomAbsShape) -> Vec<f64> {
        self.my_surf.u_intervals_of(s)
    }

    /// OCCT VIntervals (cxx L211-215).
    pub fn v_intervals_of(&self, s: crate::math::GeomAbsShape) -> Vec<f64> {
        self.my_surf.v_intervals_of(s)
    }

    /// OCCT UTrim (cxx L219-224):
    /// transformedAdaptor().UTrim(theFirst, theLast, theTol).
    pub fn u_trim_of(&self, first: f64, last: f64, tol: f64) -> GeomSurfaceAdaptor {
        self.transformed_adaptor().u_trim_of(first, last, tol)
    }

    /// OCCT VTrim (cxx L228-233).
    pub fn v_trim_of(&self, first: f64, last: f64, tol: f64) -> GeomSurfaceAdaptor {
        self.transformed_adaptor().v_trim_of(first, last, tol)
    }

    /// OCCT EvalD0 (cxx L237-244): the point transform.
    pub fn value_at(&self, u: f64, v: f64) -> DVec3 {
        if self.my_trsf.form == TrsfForm::Identity {
            return self.my_surf.value_at(u, v);
        }
        self.my_trsf.apply(self.my_surf.value_at(u, v))
    }

    /// OCCT EvalD1 (cxx L248-260): the point and the derivative vectors
    /// transform (gp_Pnt::Transformed / gp_Vec::Transform).
    pub fn d1_at(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        if self.my_trsf.form == TrsfForm::Identity {
            return self.my_surf.d1_at(u, v);
        }
        let (p, d1u, d1v) = self.my_surf.d1_at(u, v);
        (
            self.my_trsf.apply(p),
            self.my_trsf.transform_vec(d1u),
            self.my_trsf.transform_vec(d1v),
        )
    }

    /// OCCT EvalD2 (cxx L264-279).
    pub fn d2_at(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        if self.my_trsf.form == TrsfForm::Identity {
            return self.my_surf.d2_at(u, v);
        }
        let (p, d1u, d1v, d2u, d2v, d2uv) = self.my_surf.d2_at(u, v);
        let t = &self.my_trsf;
        (
            t.apply(p),
            t.transform_vec(d1u),
            t.transform_vec(d1v),
            t.transform_vec(d2u),
            t.transform_vec(d2v),
            t.transform_vec(d2uv),
        )
    }

    /// OCCT EvalD3 (cxx L283-302).
    #[allow(clippy::type_complexity)]
    pub fn d3_at(
        &self,
        u: f64,
        v: f64,
    ) -> (
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
    ) {
        if self.my_trsf.form == TrsfForm::Identity {
            return self.my_surf.d3_at(u, v);
        }
        let (p, d1u, d1v, d2u, d2v, d2uv, d3u, d3v, d3uuv, d3uvv) = self.my_surf.d3_at(u, v);
        let t = &self.my_trsf;
        (
            t.apply(p),
            t.transform_vec(d1u),
            t.transform_vec(d1v),
            t.transform_vec(d2u),
            t.transform_vec(d2v),
            t.transform_vec(d2uv),
            t.transform_vec(d3u),
            t.transform_vec(d3v),
            t.transform_vec(d3uuv),
            t.transform_vec(d3uvv),
        )
    }

    /// OCCT EvalDN (cxx L306-316).
    pub fn dn_at(&self, u: f64, v: f64, nu: i32, nv: i32) -> DVec3 {
        if self.my_trsf.form == TrsfForm::Identity {
            return self.my_surf.dn_at(u, v, nu, nv);
        }
        self.my_trsf
            .transform_vec(self.my_surf.dn_at(u, v, nu, nv))
    }

    /// OCCT Plane (cxx L320-323).
    pub fn plane(&self) -> Plane {
        self.transformed_adaptor().plane()
    }

    /// OCCT Cylinder (cxx L327-330).
    pub fn cylinder(&self) -> CylindricalSurface {
        self.transformed_adaptor().cylinder()
    }

    /// OCCT Cone (cxx L334-337).
    pub fn cone(&self) -> ConicalSurface {
        self.transformed_adaptor().cone()
    }

    /// OCCT Sphere (cxx L341-344).
    pub fn sphere(&self) -> SphericalSurface {
        self.transformed_adaptor().sphere()
    }

    /// OCCT Torus (cxx L348-351).
    pub fn torus(&self) -> ToroidalSurface {
        self.transformed_adaptor().torus()
    }

    /// OCCT Bezier (cxx L355-358).
    pub fn bezier(&self) -> crate::geom::BezierSurface {
        self.transformed_adaptor().bezier()
    }

    /// OCCT BSpline (cxx L362-365).
    pub fn bspline(&self) -> crate::geom::BSplineSurface {
        self.transformed_adaptor().bspline()
    }

    /// OCCT AxeOfRevolution (cxx L369-372).
    pub fn axe_of_revolution(&self) -> (DVec3, DVec3) {
        self.transformed_adaptor().axe_of_revolution()
    }

    /// OCCT Direction (cxx L376-379).
    pub fn direction(&self) -> DVec3 {
        self.transformed_adaptor().direction()
    }

    /// OCCT BasisCurve (cxx L383-386).
    pub fn basis_curve_of(&self) -> GeomCurveAdaptor {
        self.transformed_adaptor().basis_curve_of()
    }

    /// OCCT BasisSurface (cxx L390-393).
    pub fn basis_surface_of(&self) -> GeomSurfaceAdaptor {
        self.transformed_adaptor().basis_surface_of()
    }

    /// OCCT OffsetValue (cxx L397-400).
    pub fn offset_value_of(&self) -> f64 {
        self.transformed_adaptor().offset_value_of()
    }
}

impl Default for TransformedSurfaceAdaptor {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT gp_Pln::Transformed(theT) — the plane frame transformation (location
/// apply + the gp_Dir::Transform of the axes).
fn transform_plane(pl: &Plane, t: &Trsf) -> Plane {
    Plane {
        origin: t.apply(pl.origin),
        normal: t.transform_dir(pl.normal),
        u_dir: t.transform_dir(pl.u_dir),
        v_dir: t.transform_dir(pl.v_dir),
    }
}

// =========================================================================
// The trait implementations — the inherited-member delegation
// =========================================================================

impl Adaptor3dSurface for TransformedSurfaceAdaptor {
    /// hxx L126.
    fn first_u_parameter(&self) -> f64 {
        self.my_surf.first_u_parameter()
    }

    /// hxx L128.
    fn last_u_parameter(&self) -> f64 {
        self.my_surf.last_u_parameter()
    }

    /// hxx L130.
    fn first_v_parameter(&self) -> f64 {
        self.my_surf.first_v_parameter()
    }

    /// hxx L132.
    fn last_v_parameter(&self) -> f64 {
        self.my_surf.last_v_parameter()
    }

    /// hxx L134.
    fn u_continuity(&self) -> GeomAbsShape {
        self.my_surf.u_continuity()
    }

    /// hxx L136.
    fn v_continuity(&self) -> GeomAbsShape {
        self.my_surf.v_continuity()
    }

    /// hxx L138.
    fn nb_u_intervals(&self, s: crate::math::GeomAbsShape) -> usize {
        self.my_surf.nb_u_intervals(s)
    }

    /// hxx L140.
    fn nb_v_intervals(&self, s: crate::math::GeomAbsShape) -> usize {
        self.my_surf.nb_v_intervals(s)
    }

    /// cxx L203-207.
    fn u_intervals(&self, s: crate::math::GeomAbsShape) -> Vec<f64> {
        self.u_intervals_of(s)
    }

    /// cxx L211-215.
    fn v_intervals(&self, s: crate::math::GeomAbsShape) -> Vec<f64> {
        self.v_intervals_of(s)
    }

    /// cxx L219-224.
    fn u_trim(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor3dSurface> {
        Arc::new(self.u_trim_of(first, last, tol))
    }

    /// cxx L228-233.
    fn v_trim(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor3dSurface> {
        Arc::new(self.v_trim_of(first, last, tol))
    }

    /// hxx L156.
    fn is_u_closed(&self) -> bool {
        self.my_surf.is_u_closed()
    }

    /// hxx L158.
    fn is_v_closed(&self) -> bool {
        self.my_surf.is_v_closed()
    }

    /// hxx L160.
    fn is_u_periodic(&self) -> bool {
        self.my_surf.is_u_periodic()
    }

    /// hxx L162.
    fn u_period(&self) -> f64 {
        self.my_surf.u_period()
    }

    /// hxx L164.
    fn is_v_periodic(&self) -> bool {
        self.my_surf.is_v_periodic()
    }

    /// hxx L166.
    fn v_period(&self) -> f64 {
        self.my_surf.v_period()
    }

    /// hxx L195 (UResolution).
    fn u_resolution(&self, r3d: f64) -> f64 {
        self.my_surf.u_resolution(r3d)
    }

    /// hxx L197 (VResolution).
    fn v_resolution(&self, r3d: f64) -> f64 {
        self.my_surf.v_resolution(r3d)
    }

    /// hxx L199.
    fn get_type(&self) -> GeomAbsSurfaceType {
        self.my_surf.get_type()
    }

    /// cxx L237-244.
    fn value(&self, u: f64, v: f64) -> DVec3 {
        self.value_at(u, v)
    }

    /// cxx L248-260.
    fn d1(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        self.d1_at(u, v)
    }

    /// cxx L264-279.
    fn d2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        self.d2_at(u, v)
    }

    /// cxx L283-302.
    #[allow(clippy::type_complexity)]
    fn d3(
        &self,
        u: f64,
        v: f64,
    ) -> (
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
    ) {
        self.d3_at(u, v)
    }

    /// cxx L306-316.
    fn dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> DVec3 {
        self.dn_at(u, v, nu, nv)
    }

    /// hxx L211.
    fn u_degree(&self) -> usize {
        self.my_surf.u_degree()
    }

    /// hxx L213.
    fn nb_u_poles(&self) -> usize {
        self.my_surf.nb_u_poles()
    }

    /// hxx L215.
    fn v_degree(&self) -> usize {
        self.my_surf.v_degree()
    }

    /// hxx L217.
    fn nb_v_poles(&self) -> usize {
        self.my_surf.nb_v_poles()
    }

    /// hxx L219.
    fn nb_u_knots(&self) -> usize {
        self.my_surf.nb_u_knots()
    }

    /// hxx L221.
    fn nb_v_knots(&self) -> usize {
        self.my_surf.nb_v_knots()
    }

    /// hxx L223.
    fn is_u_rational(&self) -> bool {
        self.my_surf.is_u_rational()
    }

    /// hxx L225.
    fn is_v_rational(&self) -> bool {
        self.my_surf.is_v_rational()
    }

    /// cxx L362-365.
    fn bspline(&self) -> crate::geom::BSplineSurface {
        self.bspline()
    }

    /// cxx L376-379.
    fn direction(&self) -> DVec3 {
        self.direction()
    }

    /// cxx L383-386.
    fn basis_curve(&self) -> Arc<dyn Adaptor3dCurve> {
        Arc::new(self.basis_curve_of())
    }

    /// cxx L390-393.
    fn basis_surface(&self) -> Arc<dyn Adaptor3dSurface> {
        Arc::new(self.basis_surface_of())
    }

    /// cxx L397-400.
    fn offset_value(&self) -> f64 {
        self.offset_value_of()
    }

    /// cxx L66-77.
    fn shallow_copy(&self) -> Arc<dyn Adaptor3dSurface> {
        Arc::new(self.shallow_copy_of())
    }

    /// The rcad routing bridge — the untransformed kernel surface.
    fn kernel_surface(&self) -> Option<&Surface3> {
        self.my_surf.kernel_surface()
    }
}

impl Adaptor3dSurfaceGeom for TransformedSurfaceAdaptor {
    /// cxx L320-323.
    fn plane(&self) -> Plane {
        self.plane()
    }

    /// cxx L327-330.
    fn cylinder(&self) -> CylindricalSurface {
        self.cylinder()
    }

    /// cxx L334-337.
    fn cone(&self) -> ConicalSurface {
        self.cone()
    }

    /// cxx L341-344.
    fn sphere(&self) -> SphericalSurface {
        self.sphere()
    }

    /// cxx L348-351.
    fn torus(&self) -> ToroidalSurface {
        self.torus()
    }

    /// cxx L369-372.
    fn axe_of_revolution(&self) -> (DVec3, DVec3) {
        self.axe_of_revolution()
    }
}
