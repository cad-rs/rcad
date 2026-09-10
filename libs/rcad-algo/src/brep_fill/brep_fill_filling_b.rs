//! OCCT BRepFill_Filling — GAP placeholder types and split helpers
//! (TKTopAlgo / TKGeomAlgo dependencies not translated yet, plan §0.6).
//!
//! The main file ([`super::brep_fill_filling`]) keeps the OCCT control flow
//! 1:1; the classes below carry the OCCT constructor/method surface so the
//! call sites stay line-by-line, and panic until the owning translation units
//! land:
//! - BRepFill_CurveConstraint (TKBool/BRepFill — own file, separate unit);
//! - Adaptor3d_CurveOnSurface + GeomAdaptor_Surface + Geom2dAdaptor_Curve
//!   (TKGeomBase/adaptor3d hierarchy);
//! - GeomPlate_PlateG0Criterion / GeomPlate_MakeApprox (TKGeomAlgo — the
//!   fillet module carries the same two gaps);
//! - GeomPlate_PointConstraint(U,V,surf,...) — the UV ctor of the point
//!   constraint (the rcad GeomPlate port covers the point ctor only);
//! - GeomPlate curve-constraint consumption (Add / Curves2d()->Value) — the
//!   rcad GeomPlate port is point-constraint-path only (see the geomplate
//!   module header; the curve machinery follows the ThruSections precedent).

use glam::{DVec2, DVec3};

use rcad_kernel::geom::{Curve2d, Curve3, Surface3};
use rcad_kernel::topo::topods::{BRep, GeomAbsShape, Shape};

use crate::geomalgo::geomplate::build_plate_surface::BuildPlateSurface;
use crate::geomalgo::geomplate::point_constraint::PointConstraint;

// ---------------------------------------------------------------------------
// Adaptor3d hierarchy (TKGeomBase / TKG3d)
// ---------------------------------------------------------------------------

/// OCCT BRepAdaptor_Curve — mapped to the stored 3d curve of the edge and its
/// range (the direct `TEdgeData` mapping used across brep_fill).
pub struct BRepAdaptorCurve {
    /// The 3d curve (`TEdgeData::curve`).
    pub curve: Curve3,
    /// OCCT FirstParameter().
    pub first: f64,
    /// OCCT LastParameter().
    pub last: f64,
}

impl BRepAdaptorCurve {
    /// OCCT BRepAdaptor_Curve(); Initialize(E).
    pub fn initialize(brep: &BRep, e: &Shape) -> Self {
        let ed = brep.edge(e.clone());
        BRepAdaptorCurve {
            curve: ed.curve.clone().expect("BRepAdaptor_Curve: no 3d curve"),
            first: ed.range[0],
            last: ed.range[1],
        }
    }
}

/// OCCT GeomAdaptor_Surface — GAP: the adaptor hierarchy is not translated
/// (plan §0.6).  The wrapped `Surface3` is kept so the data path stays typed.
pub struct GeomAdaptorSurface {
    #[allow(dead_code)]
    pub surface: Surface3,
}

impl GeomAdaptorSurface {
    /// OCCT GeomAdaptor_Surface(S).
    pub fn new(surface: Surface3) -> Self {
        GeomAdaptorSurface { surface }
    }
}

/// OCCT Geom2dAdaptor_Curve — GAP (plan §0.6).
pub struct Geom2dAdaptorCurve {
    #[allow(dead_code)]
    pub curve: Curve2d,
    #[allow(dead_code)]
    pub first: f64,
    #[allow(dead_code)]
    pub last: f64,
}

impl Geom2dAdaptorCurve {
    /// OCCT Geom2dAdaptor_Curve(C2d).
    pub fn new(curve: Curve2d, first: f64, last: f64) -> Self {
        Geom2dAdaptorCurve {
            curve,
            first,
            last,
        }
    }
}

/// OCCT Adaptor3d_CurveOnSurface — GAP: not translated (plan §0.6).  The
/// wrapped (pcurve, surface) pair is kept so the data path stays typed.
pub struct Adaptor3dCurveOnSurface {
    #[allow(dead_code)]
    pub curve_2d: Geom2dAdaptorCurve,
    #[allow(dead_code)]
    pub surface: GeomAdaptorSurface,
}

impl Adaptor3dCurveOnSurface {
    /// OCCT Adaptor3d_CurveOnSurface(Curve2d, Surf).
    pub fn new(curve_2d: Geom2dAdaptorCurve, surface: GeomAdaptorSurface) -> Self {
        Adaptor3dCurveOnSurface {
            curve_2d,
            surface,
        }
    }
}

/// OCCT handle<Adaptor3d_Curve> (the polymorphic curve handle).
pub enum Adaptor3dCurve {
    /// BRepAdaptor_Curve over an edge.
    BRep(BRepAdaptorCurve),
    /// Adaptor3d_CurveOnSurface.
    OnSurface(Adaptor3dCurveOnSurface),
}

// ---------------------------------------------------------------------------
// BRepFill_CurveConstraint (TKBool/BRepFill — separate translation unit)
// ---------------------------------------------------------------------------

/// OCCT BRepFill_CurveConstraint — GAP: own file in TKBool/BRepFill, not part
/// of this batch (plan §0.6).
pub struct BRepFillCurveConstraint;

impl BRepFillCurveConstraint {
    /// OCCT BRepFill_CurveConstraint(HCurve, Order, NbPts, Tol3d).
    pub fn new(
        _hcurve: &Adaptor3dCurve,
        _order: GeomAbsShape,
        _nb_pts: i32,
        _tol3d: f64,
    ) -> Self {
        panic!(
            "GAP: BRepFill_CurveConstraint (TKBool/BRepFill) is not translated — see file header"
        )
    }

    /// OCCT BRepFill_CurveConstraint(HCurvOnSurf, Order, NbPts, Tol3d, TolAng, TolCurv).
    pub fn new_on_surface(
        _hcurv_on_surf: &Adaptor3dCurveOnSurface,
        _order: GeomAbsShape,
        _nb_pts: i32,
        _tol3d: f64,
        _tolang: f64,
        _tolcurv: f64,
    ) -> Self {
        panic!(
            "GAP: BRepFill_CurveConstraint (TKBool/BRepFill) is not translated — see file header"
        )
    }
}

/// OCCT handle<GeomPlate_CurveConstraint> Constr — the polymorphic handle
/// assigned in AddConstraints (BRepFill_Filling.cxx L303-389).
pub enum CurveConstraintHandle {
    /// OCCT: new BRepFill_CurveConstraint(...) — GAP.
    BRepFill(BRepFillCurveConstraint),
    /// OCCT: new GeomPlate_CurveConstraint(HCurvOnSurf, Order, NbPts, Tol3d,
    /// TolAng, TolCurv) — GAP: the rcad GeomPlate port covers the
    /// point-constraint path only.
    GeomPlate {
        #[allow(dead_code)]
        curv_on_surf: Adaptor3dCurveOnSurface,
        #[allow(dead_code)]
        order: GeomAbsShape,
        #[allow(dead_code)]
        nb_pts: i32,
        #[allow(dead_code)]
        tol3d: f64,
        #[allow(dead_code)]
        tolang: f64,
        #[allow(dead_code)]
        tolcurv: f64,
    },
}

impl CurveConstraintHandle {
    /// OCCT Constr->SetCurve2dOnSurf(Curve2d) (L384) — the rcad GeomPlate
    /// CurveConstraint carries only the anchor-out-of-scope form.
    pub fn set_curve2d_on_surf(&mut self, _curve2d: Curve2d) {
        panic!(
            "GAP: GeomPlate_CurveConstraint::SetCurve2dOnSurf (curve path) is not ported — \
             see the geomplate module header"
        )
    }
}

/// OCCT myBuilder->Add(Constr) (L387) — the rcad GeomPlate port covers the
/// point-constraint path only; the curve-constraint machinery (Add included)
/// follows the ThruSections anchor-out-of-scope precedent.
pub fn builder_add_curve_constraint(
    _builder: &mut BuildPlateSurface,
    _cont: CurveConstraintHandle,
) {
    panic!(
        "GAP: GeomPlate_BuildPlateSurface::Add(CurveConstraint) (curve path) is not ported — \
         see the geomplate module header"
    )
}

/// OCCT myBuilder->Curves2d()->Value(i) (L720/732) — the rcad GeomPlate port
/// keeps the point-constraint path only (Curves2d is curve path).
pub fn curves_on_plate_value(_i: usize) -> Curve2d {
    panic!(
        "GAP: GeomPlate_BuildPlateSurface::Curves2d (curve path) is not ported — \
         see the geomplate module header"
    )
}

// ---------------------------------------------------------------------------
// GeomPlate UV point constraint + approximation stack (TKGeomAlgo)
// ---------------------------------------------------------------------------

/// OCCT GeomPlate_PointConstraint(U, V, Surface, Order, Tol3d, TolAng, TolCurv)
/// — the UV ctor; the rcad GeomPlate port covers the point ctor only.
#[allow(clippy::too_many_arguments)]
pub fn point_constraint_uv(
    _u: f64,
    _v: f64,
    _surface: Surface3,
    _order: GeomAbsShape,
    _tol3d: f64,
    _tolang: f64,
    _tolcurv: f64,
) -> PointConstraint {
    panic!(
        "GAP: GeomPlate_PointConstraint(U,V,Surf,...) (UV ctor) is not ported — \
         see the geomplate module header"
    )
}

/// OCCT GeomPlate_PlateG0Criterion — GAP: not translated (plan §0.6; the
/// fillet module carries the same gap).
pub struct GeomPlatePlateG0Criterion;

impl GeomPlatePlateG0Criterion {
    /// OCCT GeomPlate_PlateG0Criterion(S2d, S3d, Seuil).
    pub fn new(_s2d: &[DVec2], _s3d: &[DVec3], _seuil: f64) -> Self {
        panic!(
            "GAP: GeomPlate_PlateG0Criterion (TKGeomAlgo) is not translated — see file header"
        )
    }
}

/// OCCT GeomPlate_MakeApprox — GAP: not translated (plan §0.6).
pub struct GeomPlateMakeApprox;

impl GeomPlateMakeApprox {
    /// OCCT GeomPlate_MakeApprox(GPlate, Criterion, Tol3d, NbMax, MaxDeg)
    /// (L705) — the two-ctor overload set maps to `new_criterion` / `new_dmax`.
    pub fn new_criterion(
        _gplate: &Surface3,
        _criterion: &GeomPlatePlateG0Criterion,
        _tol3d: f64,
        _nb_max: i32,
        _max_deg: i32,
    ) -> Self {
        panic!(
            "GAP: GeomPlate_MakeApprox (TKGeomAlgo) is not translated — see file header"
        )
    }

    /// OCCT GeomPlate_MakeApprox(GPlate, Tol3d, NbMax, MaxDeg, dmax, 0) (L711).
    pub fn new_dmax(
        _gplate: &Surface3,
        _tol3d: f64,
        _nb_max: i32,
        _max_deg: i32,
        _dmax: f64,
        _ordering: i32,
    ) -> Self {
        panic!(
            "GAP: GeomPlate_MakeApprox (TKGeomAlgo) is not translated — see file header"
        )
    }

    /// OCCT Surface().
    pub fn surface(&self) -> Surface3 {
        panic!(
            "GAP: GeomPlate_MakeApprox (TKGeomAlgo) is not translated — see file header"
        )
    }
}
