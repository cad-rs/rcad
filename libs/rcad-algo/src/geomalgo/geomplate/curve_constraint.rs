//! OCCT GeomPlate_CurveConstraint (TKGeomAlgo/GeomPlate) — API skeleton.
//!
//! The curve-constraint machinery (ProjectCurve / Approx_CurveOnSurface /
//! Geom2dInt_GInter / Discretise / LoadCurve / Intersect) is outside the
//! point-constraint anchor scope (GeomPlate_BuildPlateSurface_Test) and
//! follows the ThruSections precedent: declared API with `unimplemented!()`,
//! backfilled by a later unit.

use rcad_kernel::geom::Curve3;

/// OCCT GeomPlate_CurveConstraint.
#[derive(Debug, Clone)]
pub struct CurveConstraint {
    // OCCT members (hxx) — not ported: myAdpCurv, myD2d, myNbPoints,
    // myCurv2dOnSurf, myTypan, myTol, myOrder.
}

impl CurveConstraint {
    /// OCCT ctor 1 (GeomPlate_CurveConstraint.cxx) — curve path,
    /// anchor-out-of-scope.
    pub fn new(_curve: &Curve3, _order: i32, _nb_points: i32) -> Self {
        unimplemented!(
            "GeomPlate_CurveConstraint is not ported (curve path, anchor-out-of-scope)"
        );
    }

    /// OCCT Order().
    pub fn order(&self) -> i32 {
        unimplemented!("GeomPlate_CurveConstraint::Order is not ported");
    }

    /// OCCT NbPoints().
    pub fn nb_points(&self) -> i32 {
        unimplemented!("GeomPlate_CurveConstraint::NbPoints is not ported");
    }

    /// OCCT FirstParameter().
    pub fn first_parameter(&self) -> f64 {
        unimplemented!("GeomPlate_CurveConstraint::FirstParameter is not ported");
    }

    /// OCCT LastParameter().
    pub fn last_parameter(&self) -> f64 {
        unimplemented!("GeomPlate_CurveConstraint::LastParameter is not ported");
    }

    /// OCCT Length().
    pub fn length(&self) -> f64 {
        unimplemented!("GeomPlate_CurveConstraint::Length is not ported");
    }

    /// OCCT Curve3d().
    pub fn curve3d(&self) -> &Curve3 {
        unimplemented!("GeomPlate_CurveConstraint::Curve3d is not ported");
    }

    /// OCCT Curve2dOnSurf().
    pub fn curve2d_on_surf(&self) {
        unimplemented!("GeomPlate_CurveConstraint::Curve2dOnSurf is not ported");
    }

    /// OCCT SetCurve2dOnSurf().
    pub fn set_curve2d_on_surf(&mut self) {
        unimplemented!("GeomPlate_CurveConstraint::SetCurve2dOnSurf is not ported");
    }

    /// OCCT D0(U, P).
    pub fn d0(&self, _u: f64) -> glam::DVec3 {
        unimplemented!("GeomPlate_CurveConstraint::D0 is not ported");
    }

    /// OCCT D1(U, P, V1, V2).
    pub fn d1(&self, _u: f64) -> (glam::DVec3, glam::DVec3) {
        unimplemented!("GeomPlate_CurveConstraint::D1 is not ported");
    }
}
