//! OCCT GeomAdaptor_TransformedCurve (TKG3d/GeomAdaptor) — the 1:1
//! translation of `GeomAdaptor_TransformedCurve.hxx` (L31-210) +
//! `GeomAdaptor_TransformedCurve.cxx` (L25-325).
//!
//! OCCT inheritance `Adaptor3d_Curve <- GeomAdaptor_TransformedCurve`
//! becomes composition: the struct owns the [`GeomCurveAdaptor`] base
//! (`myCurve`), the optional curve-on-surface adaptor (`myConSurf`) and the
//! transformation; every inherited member is an explicit delegation of the
//! OCCT `myConSurf.IsNull() ? myCurve.X() : myConSurf->X()` one-liners
//! (hxx L103-198).  The `GetType()` virtual lands on the
//! [`Adaptor3dCurveGeom`] companion trait (see the geom_adaptor_curve module
//! header for the placement rationale).

use std::sync::Arc;

use glam::DVec3;

use super::adaptor::Adaptor3dCurve;
use super::adaptor::CurveOnSurface;
use super::geom_adaptor_curve::{Adaptor3dCurveGeom, GeomCurveAdaptor};
use super::CurveType;
use crate::geom::{
    Circle3, Curve3, Ellipse3, Hyperbola3, Line3, Parabola3,
};
use crate::math::gp::{Trsf, TrsfForm};
use crate::math::{GeomAbsCurveType, GeomAbsShape};

// =========================================================================
// OCCT GeomAdaptor_TransformedCurve (hxx L31-210)
// =========================================================================

/// OCCT GeomAdaptor_TransformedCurve — an adaptor for curves with an applied
/// transformation.
#[derive(Clone)]
pub struct TransformedCurveAdaptor {
    /// OCCT: GeomAdaptor_Curve myCurve (hxx L207).
    pub my_curve: GeomCurveAdaptor,
    /// OCCT: occ::handle<Adaptor3d_CurveOnSurface> myConSurf (hxx L208).
    pub my_con_surf: Option<Arc<CurveOnSurface>>,
    /// OCCT: gp_Trsf myTrsf (hxx L209).
    pub my_trsf: Trsf,
}

impl TransformedCurveAdaptor {
    /// OCCT GeomAdaptor_TransformedCurve() (cxx L25) — the undefined curve
    /// with identity transformation.
    pub fn new() -> Self {
        TransformedCurveAdaptor {
            my_curve: empty_curve_adaptor(),
            my_con_surf: None,
            my_trsf: Trsf::identity(),
        }
    }

    /// OCCT GeomAdaptor_TransformedCurve(theCurve, theTrsf) (cxx L29-34).
    pub fn with_curve_trsf(the_curve: Curve3, the_trsf: Trsf) -> Self {
        TransformedCurveAdaptor {
            my_curve: GeomCurveAdaptor::new(the_curve),
            my_con_surf: None,
            my_trsf: the_trsf,
        }
    }

    /// OCCT GeomAdaptor_TransformedCurve(theCurve, theFirst, theLast,
    /// theTrsf) (cxx L38-45).
    pub fn with_range_trsf(the_curve: Curve3, the_first: f64, the_last: f64, the_trsf: Trsf) -> Self {
        TransformedCurveAdaptor {
            my_curve: GeomCurveAdaptor::with_range(the_curve, the_first, the_last),
            my_con_surf: None,
            my_trsf: the_trsf,
        }
    }

    /// OCCT ShallowCopy (cxx L49-64).
    pub fn shallow_copy_of(&self) -> Self {
        self.clone()
    }

    /// OCCT Load(theCurve) (hxx L59): myCurve.Load(theCurve).
    pub fn load(&mut self, the_curve: Curve3) {
        self.my_curve.load(the_curve);
    }

    /// OCCT Load(theCurve, theFirst, theLast) (hxx L65-68).
    pub fn load_with_range(&mut self, the_curve: Curve3, the_first: f64, the_last: f64) {
        self.my_curve.load_with_range(the_curve, the_first, the_last);
    }

    /// OCCT LoadCurveOnSurface(theConSurf) (hxx L72-75).
    pub fn load_curve_on_surface(&mut self, the_con_surf: Arc<CurveOnSurface>) {
        self.my_con_surf = Some(the_con_surf);
    }

    /// OCCT SetTrsf(theTrsf) (hxx L79).
    pub fn set_trsf(&mut self, the_trsf: Trsf) {
        self.my_trsf = the_trsf;
    }

    /// OCCT Trsf() (hxx L82).
    pub fn trsf(&self) -> &Trsf {
        &self.my_trsf
    }

    /// OCCT Is3DCurve() (hxx L85): myConSurf.IsNull().
    pub fn is_3d_curve(&self) -> bool {
        self.my_con_surf.is_none()
    }

    /// OCCT IsCurveOnSurface() (hxx L88): !myConSurf.IsNull().
    pub fn is_curve_on_surface(&self) -> bool {
        self.my_con_surf.is_some()
    }

    /// OCCT Curve() (hxx L91) — the underlying GeomAdaptor_Curve.
    pub fn curve(&self) -> &GeomCurveAdaptor {
        &self.my_curve
    }

    /// OCCT ChangeCurve() (hxx L94) — the underlying GeomAdaptor_Curve for
    /// modification.
    pub fn change_curve(&mut self) -> &mut GeomCurveAdaptor {
        &mut self.my_curve
    }

    /// OCCT CurveOnSurface() (hxx L97) — the CurveOnSurface adaptor.
    pub fn curve_on_surface(&self) -> &CurveOnSurface {
        self.my_con_surf
            .as_ref()
            .expect("Standard_NullObject: GeomAdaptor_TransformedCurve::CurveOnSurface")
    }

    /// OCCT GeomCurve() (hxx L100): myCurve.Curve() — the underlying Geom
    /// curve.
    pub fn geom_curve(&self) -> &Curve3 {
        &self.my_curve.curve
    }

    // ---------------------------------------------------------------------
    // The member bodies (cxx).
    // ---------------------------------------------------------------------

    /// OCCT Intervals(theT, theS) (cxx L68-79).
    pub fn intervals_of(&self, s: GeomAbsShape) -> Vec<f64> {
        if self.my_con_surf.is_none() {
            self.my_curve.intervals_of(s)
        } else {
            self.my_con_surf
                .as_ref()
                .expect("con surf")
                .intervals(s)
        }
    }

    /// OCCT Trim(theFirst, theLast, theTol) (cxx L83-99).
    pub fn trim_of(&self, the_first: f64, the_last: f64, the_tol: f64) -> Self {
        let mut a_copy = TransformedCurveAdaptor::new();
        if self.my_con_surf.is_none() {
            // OCCT L90: aCopy->myCurve.Load(myCurve.Curve(), theFirst,
            // theLast).
            a_copy
                .my_curve
                .load_with_range(self.my_curve.curve.clone(), the_first, the_last);
        } else {
            // OCCT L94-95: the downcast of myConSurf->Trim(theFirst, theLast,
            // theTol) to Adaptor3d_CurveOnSurface — the rcad encoding builds
            // the CurveOnSurface Trim body directly (the OCCT Trim of
            // Adaptor3d_CurveOnSurface.cxx L1133-1141: Load(mySurface) +
            // Load(myCurve->Trim(...))).
            let con = self.my_con_surf.as_ref().expect("con surf");
            a_copy.my_con_surf = Some(Arc::new(CurveOnSurface::new(
                con.my2d_curve.trim(the_first, the_last, the_tol),
                con.my_surface.clone(),
            )));
        }
        // OCCT L97: aCopy->myTrsf = myTrsf.
        a_copy.my_trsf = self.my_trsf;
        a_copy
    }

    /// OCCT Line (cxx L103-116).
    pub fn line(&self) -> Line3 {
        let mut a_l = if self.my_con_surf.is_none() {
            self.my_curve.line()
        } else {
            // OCCT L112: myConSurf->Line() — raises for the non-line
            // curve-on-surface composites.
            panic!("Standard_NoSuchObject: Adaptor3d_CurveOnSurface::Line(): curve is not a line")
        };
        // OCCT L114: aL.Transform(myTrsf) — the gp_Lin transformation
        // (location apply + direction transform).
        a_l.origin = self.my_trsf.apply(a_l.origin);
        a_l.direction = self.my_trsf.transform_dir(a_l.direction);
        a_l
    }

    /// OCCT Circle (cxx L120-133).
    pub fn circle(&self) -> Circle3 {
        let mut a_c = if self.my_con_surf.is_none() {
            self.my_curve.circle()
        } else {
            panic!("Standard_NoSuchObject: Adaptor3d_CurveOnSurface::Line(): curve is not a circle")
        };
        // OCCT L130: aC.Transform(myTrsf).
        a_c.center = self.my_trsf.apply(a_c.center);
        a_c.x_dir = self.my_trsf.transform_dir(a_c.x_dir);
        a_c.y_dir = self.my_trsf.transform_dir(a_c.y_dir);
        a_c.normal = self.my_trsf.transform_dir(a_c.normal);
        a_c
    }

    /// OCCT Ellipse (cxx L137-150).
    pub fn ellipse(&self) -> Ellipse3 {
        let mut a_e = if self.my_con_surf.is_none() {
            self.my_curve.ellipse()
        } else {
            panic!("Standard_NoSuchObject: Adaptor3d_CurveOnSurface::Ellipse")
        };
        // OCCT L147: aE.Transform(myTrsf).
        a_e.center = self.my_trsf.apply(a_e.center);
        a_e.normal = self.my_trsf.transform_dir(a_e.normal);
        a_e.major_dir = self.my_trsf.transform_dir(a_e.major_dir);
        a_e
    }

    /// OCCT Hyperbola (cxx L154-167).
    pub fn hyperbola(&self) -> Hyperbola3 {
        let mut a_h = if self.my_con_surf.is_none() {
            self.my_curve.hyperbola()
        } else {
            panic!("Standard_NoSuchObject: Adaptor3d_CurveOnSurface::Hyperbola")
        };
        // OCCT L164: aH.Transform(myTrsf).
        a_h.center = self.my_trsf.apply(a_h.center);
        a_h.normal = self.my_trsf.transform_dir(a_h.normal);
        a_h.major_dir = self.my_trsf.transform_dir(a_h.major_dir);
        a_h
    }

    /// OCCT Parabola (cxx L171-184).
    pub fn parabola(&self) -> Parabola3 {
        let mut a_p = if self.my_con_surf.is_none() {
            self.my_curve.parabola()
        } else {
            panic!("Standard_NoSuchObject: Adaptor3d_CurveOnSurface::Parabola")
        };
        // OCCT L181: aP.Transform(myTrsf).
        a_p.vertex = self.my_trsf.apply(a_p.vertex);
        a_p.normal = self.my_trsf.transform_dir(a_p.normal);
        a_p.axis_dir = self.my_trsf.transform_dir(a_p.axis_dir);
        a_p
    }

    /// OCCT Bezier (cxx L188-201): the payload, transformed when the
    /// transformation is not the identity (the poles transform; the OCCT
    /// downcast to Geom_BezierCurve is the payload type itself).
    pub fn bezier(&self) -> crate::geom::BezierCurve3 {
        let a_bc = if self.my_con_surf.is_none() {
            self.my_curve.bezier()
        } else {
            // OCCT L196: myConSurf->Bezier() — raises for the non-plane
            // surfaces (the Adaptor3d_CurveOnSurface requirement).
            panic!("Standard_NoSuchObject: Adaptor3d_CurveOnSurface : Bezier")
        };
        if self.my_trsf.form == TrsfForm::Identity {
            a_bc
        } else {
            let mut b = a_bc;
            b.control_points = b
                .control_points
                .iter()
                .map(|p| self.my_trsf.apply(*p))
                .collect();
            b
        }
    }

    /// OCCT BSpline (cxx L205-218).
    pub fn bspline(&self) -> crate::geom::BSplineCurve3 {
        let a_bs = if self.my_con_surf.is_none() {
            self.my_curve.bspline()
        } else {
            // OCCT L214: myConSurf->BSpline() — raises for the non-plane
            // surfaces.
            panic!("Standard_NoSuchObject: Adaptor3d_CurveOnSurface : BSpline")
        };
        if self.my_trsf.form == TrsfForm::Identity {
            a_bs
        } else {
            let mut b = a_bs;
            b.control_points = b
                .control_points
                .iter()
                .map(|p| self.my_trsf.apply(*p))
                .collect();
            b
        }
    }

    /// OCCT OffsetCurve (cxx L222-233).
    pub fn offset_curve(&self) -> crate::geom::OffsetCurve3 {
        // OCCT L224: raise when !Is3DCurve() or GetType() != OffsetCurve.
        if !self.is_3d_curve() || self.my_curve.my_type_curve != GeomAbsCurveType::OffsetCurve {
            panic!("Standard_NoSuchObject: GeomAdaptor_TransformedCurve::OffsetCurve");
        }
        let an_off_c = self.my_curve.offset_curve();
        if self.my_trsf.form == TrsfForm::Identity {
            return an_off_c;
        }
        // OCCT L232: anOffC->Transformed(myTrsf) — the basis curve and the
        // offset direction transform (Geom_OffsetCurve::Transform).
        let basis_t = an_off_c.basis.as_ref().transformed_payload(&self.my_trsf);
        crate::geom::OffsetCurve3 {
            basis: Box::new(basis_t),
            offset_distance: an_off_c.offset_distance,
            offset_dir: self.my_trsf.transform_dir(an_off_c.offset_dir),
        }
    }

    /// OCCT EvalD0 (cxx L237-250).
    pub fn value_at(&self, u: f64) -> DVec3 {
        let a_p = if self.my_con_surf.is_none() {
            self.my_curve.value_at(u)
        } else {
            // OCCT L246: myConSurf->D0(theU, aP).
            self.my_con_surf.as_ref().expect("con surf").value(u)
        };
        self.my_trsf.apply(a_p)
    }

    /// OCCT EvalD1 (cxx L254-268).
    pub fn d1_at(&self, u: f64) -> (DVec3, DVec3) {
        let (p, d1) = if self.my_con_surf.is_none() {
            self.my_curve.d1_at(u)
        } else {
            // OCCT L263: myConSurf->D1(theU, aRes.Point, aRes.D1).
            let con = self.my_con_surf.as_ref().expect("con surf");
            let (p, d1) = con.d1(u);
            (p, d1)
        };
        (self.my_trsf.apply(p), self.my_trsf.transform_vec(d1))
    }

    /// OCCT EvalD2 (cxx L272-287).
    pub fn d2_at(&self, u: f64) -> (DVec3, DVec3, DVec3) {
        let (p, d1, d2) = if self.my_con_surf.is_none() {
            self.my_curve.d2_at(u)
        } else {
            let con = self.my_con_surf.as_ref().expect("con surf");
            let (p, d1, d2) = con.d2(u);
            (p, d1, d2)
        };
        (
            self.my_trsf.apply(p),
            self.my_trsf.transform_vec(d1),
            self.my_trsf.transform_vec(d2),
        )
    }

    /// OCCT EvalD3 (cxx L291-307).
    pub fn d3_at(&self, u: f64) -> (DVec3, DVec3, DVec3, DVec3) {
        let (p, d1, d2, d3) = if self.my_con_surf.is_none() {
            self.my_curve.d3_at(u)
        } else {
            let con = self.my_con_surf.as_ref().expect("con surf");
            // OCCT L300: myConSurf->D3(theU, Point, D1, D2, D3) — the
            // curve-on-surface D3.
            let (p, d1, d2) = con.d2(u);
            let d3 = con.d3(u).3;
            (p, d1, d2, d3)
        };
        (
            self.my_trsf.apply(p),
            self.my_trsf.transform_vec(d1),
            self.my_trsf.transform_vec(d2),
            self.my_trsf.transform_vec(d3),
        )
    }

    /// OCCT EvalDN (cxx L311-324).
    pub fn dn_at(&self, u: f64, n: i32) -> DVec3 {
        let a_v = if self.my_con_surf.is_none() {
            self.my_curve.dn_at(u, n)
        } else {
            self.my_con_surf.as_ref().expect("con surf").dn(u, n)
        };
        self.my_trsf.transform_vec(a_v)
    }

    // ---------------------------------------------------------------------
    // The delegated members (hxx L103-198 one-liners).
    // ---------------------------------------------------------------------

    /// hxx L103-106.
    pub fn first_parameter(&self) -> f64 {
        match &self.my_con_surf {
            None => self.my_curve.first_parameter(),
            Some(c) => c.first_parameter(),
        }
    }

    /// hxx L108-111.
    pub fn last_parameter(&self) -> f64 {
        match &self.my_con_surf {
            None => self.my_curve.last_parameter(),
            Some(c) => c.last_parameter(),
        }
    }

    /// hxx L113-116.
    pub fn continuity(&self) -> GeomAbsShape {
        match &self.my_con_surf {
            None => self.my_curve.continuity(),
            Some(c) => c.continuity(),
        }
    }

    /// hxx L118-121.
    pub fn nb_intervals_of(&self, s: GeomAbsShape) -> usize {
        match &self.my_con_surf {
            None => self.my_curve.nb_intervals_of(s),
            Some(c) => c.nb_intervals(s),
        }
    }

    /// hxx L130-133.
    pub fn is_closed(&self) -> bool {
        match &self.my_con_surf {
            None => self.my_curve.is_closed(),
            Some(c) => c.is_closed(),
        }
    }

    /// hxx L135-138.
    pub fn is_periodic(&self) -> bool {
        match &self.my_con_surf {
            None => self.my_curve.is_periodic(),
            Some(c) => c.is_periodic(),
        }
    }

    /// hxx L140-143.
    pub fn period(&self) -> f64 {
        match &self.my_con_surf {
            None => self.my_curve.period(),
            Some(c) => c.period(),
        }
    }

    /// hxx L160-163.
    pub fn resolution(&self, r3d: f64) -> f64 {
        match &self.my_con_surf {
            None => self.my_curve.resolution(r3d),
            Some(c) => c.resolution(r3d),
        }
    }

    /// hxx L165-168 — the curve kind (the Adaptor3dCurveGeom encoding).
    pub fn get_type_of(&self) -> CurveType {
        match &self.my_con_surf {
            None => self.my_curve.get_type(),
            Some(c) => c.get_type(),
        }
    }

    /// hxx L180-183.
    pub fn degree(&self) -> usize {
        match &self.my_con_surf {
            None => self.my_curve.degree(),
            Some(c) => c.degree(),
        }
    }

    /// hxx L185-188.
    pub fn is_rational(&self) -> bool {
        match &self.my_con_surf {
            None => self.my_curve.is_rational(),
            Some(c) => c.is_rational(),
        }
    }

    /// hxx L190-193.
    pub fn nb_poles(&self) -> usize {
        match &self.my_con_surf {
            None => self.my_curve.nb_poles(),
            Some(c) => c.nb_poles(),
        }
    }

    /// hxx L195-198.
    pub fn nb_knots(&self) -> usize {
        match &self.my_con_surf {
            None => self.my_curve.nb_knots(),
            Some(c) => c.nb_knots(),
        }
    }
}

impl Default for TransformedCurveAdaptor {
    fn default() -> Self {
        Self::new()
    }
}

/// The undefined GeomAdaptor_Curve (the OCCT default constructor
/// GeomAdaptor_Curve.hxx L94-99: myTypeCurve = GeomAbs_OtherCurve, myFirst =
/// myLast = 0).
fn empty_curve_adaptor() -> GeomCurveAdaptor {
    GeomCurveAdaptor {
        curve: Curve3::Line(Line3::new(DVec3::ZERO, glam::DVec3::X)),
        my_type_curve: GeomAbsCurveType::OtherCurve,
        first: 0.0,
        last: 0.0,
    }
}

/// The kernel Curve3 transform used by the OffsetCurve transformation (the
/// Geom_OffsetCurve::Transform encoding: the basis curve transforms through
/// the point/direction rules).
trait CurveTransformPayload {
    fn transformed_payload(&self, t: &Trsf) -> Curve3;
}

impl CurveTransformPayload for Curve3 {
    fn transformed_payload(&self, t: &Trsf) -> Curve3 {
        match self {
            Curve3::Line(l) => Curve3::Line(Line3::new(t.apply(l.origin), t.transform_dir(l.direction))),
            Curve3::Circle(c) => Curve3::Circle(Circle3 {
                center: t.apply(c.center),
                normal: t.transform_dir(c.normal),
                x_dir: t.transform_dir(c.x_dir),
                y_dir: t.transform_dir(c.y_dir),
                radius: c.radius,
            }),
            _ => crate::geom::transform_curve(self, &t.to_daffine3()),
        }
    }
}

// =========================================================================
// The trait implementations
// =========================================================================

impl Adaptor3dCurve for TransformedCurveAdaptor {
    /// hxx L103-106.
    fn first_parameter(&self) -> f64 {
        self.first_parameter()
    }

    /// hxx L108-111.
    fn last_parameter(&self) -> f64 {
        self.last_parameter()
    }

    /// cxx L237-250.
    fn value(&self, u: f64) -> DVec3 {
        self.value_at(u)
    }

    /// cxx L254-268.
    fn d1(&self, u: f64) -> (DVec3, DVec3) {
        self.d1_at(u)
    }

    /// cxx L272-287.
    fn d2(&self, u: f64) -> (DVec3, DVec3, DVec3) {
        self.d2_at(u)
    }

    /// cxx L291-307.
    fn d3(&self, u: f64) -> (DVec3, DVec3, DVec3, DVec3) {
        self.d3_at(u)
    }

    /// cxx L311-324.
    fn dn(&self, u: f64, n: i32) -> DVec3 {
        self.dn_at(u, n)
    }

    /// hxx L113-116.
    fn continuity(&self) -> GeomAbsShape {
        self.continuity()
    }

    /// hxx L118-121.
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        self.nb_intervals_of(s)
    }

    /// cxx L68-79.
    fn intervals(&self, s: GeomAbsShape) -> Vec<f64> {
        self.intervals_of(s)
    }

    /// cxx L83-99.
    fn trim(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor3dCurve> {
        Arc::new(self.trim_of(first, last, tol))
    }

    /// hxx L130-133.
    fn is_closed(&self) -> bool {
        self.is_closed()
    }

    /// hxx L135-138.
    fn is_periodic(&self) -> bool {
        self.is_periodic()
    }

    /// hxx L140-143.
    fn period(&self) -> f64 {
        self.period()
    }

    /// hxx L160-163.
    fn resolution(&self, r3d: f64) -> f64 {
        self.resolution(r3d)
    }

    /// hxx L180-183.
    fn degree(&self) -> usize {
        self.degree()
    }

    /// hxx L185-188.
    fn is_rational(&self) -> bool {
        self.is_rational()
    }

    /// hxx L190-193.
    fn nb_poles(&self) -> usize {
        self.nb_poles()
    }

    /// hxx L195-198.
    fn nb_knots(&self) -> usize {
        self.nb_knots()
    }

    /// cxx L188-201.
    fn bezier(&self) -> crate::geom::BezierCurve3 {
        self.bezier()
    }

    /// cxx L205-218.
    fn bspline(&self) -> crate::geom::BSplineCurve3 {
        self.bspline()
    }

    /// cxx L222-233.
    fn offset_curve(&self) -> crate::geom::OffsetCurve3 {
        self.offset_curve()
    }

    /// cxx L49-64.
    fn shallow_copy(&self) -> Arc<dyn Adaptor3dCurve> {
        Arc::new(self.shallow_copy_of())
    }
}

impl Adaptor3dCurveGeom for TransformedCurveAdaptor {
    /// hxx L165-168.
    fn get_type(&self) -> CurveType {
        self.get_type_of()
    }

    /// cxx L103-116.
    fn line(&self) -> Line3 {
        self.line()
    }

    /// cxx L120-133.
    fn circle(&self) -> Circle3 {
        self.circle()
    }

    /// cxx L137-150.
    fn ellipse(&self) -> Ellipse3 {
        self.ellipse()
    }

    /// cxx L154-167.
    fn hyperbola(&self) -> Hyperbola3 {
        self.hyperbola()
    }

    /// cxx L171-184.
    fn parabola(&self) -> Parabola3 {
        self.parabola()
    }

    /// cxx L83-99.
    fn trim_geom(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor3dCurveGeom> {
        Arc::new(self.trim_of(first, last, tol))
    }
}
