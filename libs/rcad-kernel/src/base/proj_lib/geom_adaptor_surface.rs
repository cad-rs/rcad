//! OCCT GeomAdaptor_Surface (TKG3d/GeomAdaptor) — the 1:1 translation of
//! `GeomAdaptor_Surface.hxx` (L55-389) + `GeomAdaptor_Surface.cxx`
//! (L72-2351).
//!
//! The rcad type keeps the OCCT member set (`mySurface`, `myUFirst`,
//! `myULast`, `myVFirst`, `myVLast`, `myTolU`, `myTolV`, `mySurfaceType`).
//! Architecture differences (rcad encodings):
//!   - `occ::handle<Geom_Surface> mySurface` maps to the kernel surface
//!     value [`crate::geom::Surface3`];
//!   - the `SurfaceDataVariant` (gp_Pln / gp_Cylinder / ... / ExtrusionData
//!     / RevolutionData / OffsetData / BezierData / BSplineData) is carried
//!     by the `Surface3` variant payloads — the OCCT variant exists to avoid
//!     downcasts in C++;
//!   - `Geom_RectangularTrimmedSurface` is unwrapped in `load` exactly as
//!     OCCT does (L423-430 recurses on the basis surface, keeping the
//!     window);
//!   - the BSplSLib_Cache / GeomEval_RepSurfaceDesc evaluation caches and
//!     the `hasEvalRep` branches are C++-internal performance state: the
//!     kernel engine evaluates the same poles directly, so the
//!     `EvalD0/D1/D2/D3` arms share the kernel evaluation (results equal);
//!     the `IfUVBound`/`Span` span locators (L2244-2350) only choose between
//!     the OCCT LocalD1 and the cached global D1 — both evaluate identical
//!     values — and have no observable counterpart; the OCCT boundary
//!     snapping on `myTolU`/`myTolV` IS observable and is translated
//!     verbatim in the EvalD1/D2/D3/DN composition.
//!   - the `GeomAdaptor_Curve` temporaries of the interval members ride the
//!     canonical [`GeomCurveAdaptor`].

use glam::DVec3;
use std::sync::Arc;

use super::adaptor::{Adaptor3dCurve, Adaptor3dSurface, GeomAbsSurfaceType};
use super::geom_adaptor_curve::{Adaptor3dCurveGeom, GeomCurveAdaptor};
use super::proj_lib_projected_curve::TWO_PI;
use super::CurveType;
use crate::core::precision;
use crate::geom::{
    ConicalSurface, CylindricalSurface, Plane, SphericalSurface, Surface3, SurfaceEval,
    ToroidalSurface,
};
use crate::math::bspl::surface_rational_flags;
use crate::math::bspl_lib::{intervals as bspl_intervals, locate_parameter_knots_mults};
use crate::math::GeomAbsShape;

/// OCCT `static const double PosTol = Precision::PConfusion() * 0.5`
/// (GeomAdaptor_Surface.cxx L72).
const POS_TOL: f64 = precision::PCONFUSION * 0.5;

/// 1-based array read (the OCCT `NCollection_Array1` access convention).
fn at(v: &[f64], i: i32) -> f64 {
    v[(i - 1) as usize]
}

/// 1-based array read of the i32 multiplicity array.
fn ati(v: &[i32], i: i32) -> i32 {
    v[(i - 1) as usize]
}

// =========================================================================
// OCCT file-static LocalContinuity (GeomAdaptor_Surface.cxx L110-165)
// =========================================================================

/// OCCT `static GeomAbs_Shape LocalContinuity(Degree, Nb, TK, TM, PFirst,
/// PLast, IsPeriodic)` (GeomAdaptor_Surface.cxx L110-165).
fn local_continuity(
    degree: usize,
    nb: i32,
    tk: &[f64],
    tm: &[i32],
    p_first: f64,
    p_last: f64,
    is_periodic: bool,
) -> GeomAbsShape {
    let mut index1 = 0i32;
    let mut index2 = 0i32;
    let mut new_first = 0.0f64;
    let mut new_last = 0.0f64;
    // OCCT L122-123: BSplCLib::LocateParameter for both bounds.
    locate_parameter_knots_mults(
        degree,
        tk,
        tm,
        p_first,
        is_periodic,
        1,
        nb,
        &mut index1,
        &mut new_first,
    );
    locate_parameter_knots_mults(
        degree,
        tk,
        tm,
        p_last,
        is_periodic,
        1,
        nb,
        &mut index2,
        &mut new_last,
    );
    // OCCT L124-132: the knot-proximity index adjustments.
    let eps_knot = precision::PCONFUSION;
    if (new_first - at(tk, index1 + 1)).abs() < eps_knot {
        index1 += 1;
    }
    if (new_last - at(tk, index2)).abs() < eps_knot {
        index2 -= 1;
    }
    // OCCT L134-137: the periodic wrap.
    if is_periodic && index1 == nb {
        index1 = 1;
    }

    if index2 != index1 {
        // OCCT L141-149: the max multiplicity between the bounds.
        let mut multmax = ati(tm, index1 + 1);
        for i in index1 + 1..=index2 {
            if ati(tm, i) > multmax {
                multmax = ati(tm, i);
            }
        }
        let multmax = degree as i32 - multmax;
        // OCCT L150-153.
        if multmax <= 0 {
            return GeomAbsShape::C0;
        }
        // OCCT L154-162.
        match multmax {
            1 => return GeomAbsShape::C1,
            2 => return GeomAbsShape::C2,
            3 => return GeomAbsShape::C3,
            _ => {}
        }
    }
    // OCCT L164.
    GeomAbsShape::CN
}

// =========================================================================
// OCCT GeomAdaptor_Surface (GeomAdaptor_Surface.hxx L55-389)
// =========================================================================

/// OCCT GeomAdaptor_Surface — the Adaptor3d_Surface instance over a
/// Geom_Surface.  Constructor forms: [`GeomSurfaceAdaptor::empty`]
/// (hxx L117-126), [`GeomSurfaceAdaptor::new`] = `GeomAdaptor_Surface(S)`
/// (L128-133 + Load L150-160), [`GeomSurfaceAdaptor::with_window`] =
/// `GeomAdaptor_Surface(S, U1, U2, V1, V2, TolU, TolV)` (L136-145).
#[derive(Clone, Debug)]
pub struct GeomSurfaceAdaptor {
    /// OCCT: handle(Geom_Surface) mySurface.
    pub surface: Surface3,
    /// OCCT: Standard_Real myUFirst.
    pub u_first: f64,
    /// OCCT: Standard_Real myULast.
    pub u_last: f64,
    /// OCCT: Standard_Real myVFirst.
    pub v_first: f64,
    /// OCCT: Standard_Real myVLast.
    pub v_last: f64,
    /// OCCT: Standard_Real myTolU.
    pub my_tol_u: f64,
    /// OCCT: Standard_Real myTolV.
    pub my_tol_v: f64,
    /// OCCT: GeomAbs_SurfaceType mySurfaceType.
    pub my_surface_type: GeomAbsSurfaceType,
}

impl GeomSurfaceAdaptor {
    /// OCCT GeomAdaptor_Surface() (hxx L117-126) — the undefined adaptor.
    pub fn empty() -> Self {
        GeomSurfaceAdaptor {
            surface: Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z)),
            u_first: 0.0,
            u_last: 0.0,
            v_first: 0.0,
            v_last: 0.0,
            my_tol_u: 0.0,
            my_tol_v: 0.0,
            my_surface_type: GeomAbsSurfaceType::OtherSurface,
        }
    }

    /// OCCT GeomAdaptor_Surface(theSurf) -> Load(theSurf)
    /// (hxx L128-133 + L150-160): the surface natural bounds domain.  (The
    /// OCCT Standard_NullObject raise of Load for a null handle cannot be
    /// triggered — the rcad surface value is never null.)
    pub fn new(surface: Surface3) -> Self {
        let mut a = GeomSurfaceAdaptor::empty();
        a.load(surface);
        a
    }

    /// OCCT GeomAdaptor_Surface(S, U1, U2, V1, V2, TolU, TolV) (hxx
    /// L136-145): the restricted window; raises Standard_ConstructionError
    /// when U1 > U2 or V1 > V2 (hxx L162-181).
    pub fn with_window(
        surface: Surface3,
        u_first: f64,
        u_last: f64,
        v_first: f64,
        v_last: f64,
        tol_u: f64,
        tol_v: f64,
    ) -> Self {
        assert!(
            !(u_first > u_last || v_first > v_last),
            "Standard_ConstructionError: GeomAdaptor_Surface::Load"
        );
        let mut a = GeomSurfaceAdaptor::empty();
        a.load_with_tols(surface, u_first, u_last, v_first, v_last, tol_u, tol_v);
        a
    }

    /// OCCT Load(theSurf) (hxx L150-160): Bounds + load.
    pub fn load(&mut self, surface: Surface3) {
        let [u1, u2, v1, v2] = SurfaceEval::default_domain(&surface);
        self.load_with_tols(surface, u1, u2, v1, v2, 0.0, 0.0);
    }

    /// OCCT Load(theSurf, theUFirst, theULast, theVFirst, theVLast, TolU,
    /// TolV) (hxx L163-181).
    pub fn load_with_tols(
        &mut self,
        surface: Surface3,
        u_first: f64,
        u_last: f64,
        v_first: f64,
        v_last: f64,
        tol_u: f64,
        tol_v: f64,
    ) {
        assert!(
            !(u_first > u_last || v_first > v_last),
            "Standard_ConstructionError: GeomAdaptor_Surface::Load"
        );
        self.load_inner(surface, u_first, u_last, v_first, v_last, tol_u, tol_v);
    }

    /// OCCT Load(theSurf, U1, U2, V1, V2) — the restricted window form
    /// without the tolerance arguments (the consumer entry kept from the
    /// consumed subset: Load(S, U1, U2, V1, V2)).
    pub fn load_with_window(&mut self, surface: Surface3, u1: f64, u2: f64, v1: f64, v2: f64) {
        self.load_with_tols(surface, u1, u2, v1, v2, 0.0, 0.0);
    }

    /// OCCT `void GeomAdaptor_Surface::load(S, UFirst, ULast, VFirst, VLast,
    /// TolU, TolV)` (GeomAdaptor_Surface.cxx L402-531).
    fn load_inner(
        &mut self,
        s: Surface3,
        u_first: f64,
        u_last: f64,
        v_first: f64,
        v_last: f64,
        tol_u: f64,
        tol_v: f64,
    ) {
        self.my_tol_u = tol_u;
        self.my_tol_v = tol_v;
        self.u_first = u_first;
        self.u_last = u_last;
        self.v_first = v_first;
        self.v_last = v_last;

        // OCCT L417: `if (mySurface != S)` — the rcad value encoding always
        // (re)loads the dispatched type.
        // OCCT L419: mySurface = S.
        // OCCT L423-430: Geom_RectangularTrimmedSurface -> Load(basis, ...).
        let s = match &s {
            Surface3::Trimmed(t) => t.basis.as_ref().clone(),
            _ => s,
        };
        self.surface = s;

        // OCCT L431-529: the DynamicType dispatch.
        self.my_surface_type = match &self.surface {
            Surface3::Plane(_) => GeomAbsSurfaceType::Plane,
            Surface3::Cylinder(_) => GeomAbsSurfaceType::Cylinder,
            Surface3::Cone(_) => GeomAbsSurfaceType::Cone,
            Surface3::Sphere(_) => GeomAbsSurfaceType::Sphere,
            Surface3::Torus(_) => GeomAbsSurfaceType::Torus,
            Surface3::Revolution(_) => GeomAbsSurfaceType::SurfaceOfRevolution,
            Surface3::LinearExtrusion(_) => GeomAbsSurfaceType::SurfaceOfExtrusion,
            Surface3::Bezier(_) => GeomAbsSurfaceType::BezierSurface,
            Surface3::BSpline(_) => GeomAbsSurfaceType::BSplineSurface,
            Surface3::Offset(_) => GeomAbsSurfaceType::OffsetSurface,
            _ => GeomAbsSurfaceType::OtherSurface,
        };
    }

    /// OCCT Surface() (hxx L183) — the underlying surface.
    pub fn surface(&self) -> &Surface3 {
        &self.surface
    }

    /// OCCT Bounds(U1, U2, V1, V2) (hxx L198-204).
    pub fn bounds(&self) -> [f64; 4] {
        [
            self.first_u_parameter(),
            self.last_u_parameter(),
            self.first_v_parameter(),
            self.last_v_parameter(),
        ]
    }

    /// OCCT ToleranceU() (hxx L207).
    pub fn tolerance_u(&self) -> f64 {
        self.my_tol_u
    }

    /// OCCT ToleranceV() (hxx L210).
    pub fn tolerance_v(&self) -> f64 {
        self.my_tol_v
    }

    // ---------------------------------------------------------------------
    // The evaluation dispatch (GeomAdaptor_Surface.cxx EvalD0/D1/D2/D3/DN).
    // The per-type arms (ElSLib for the quadrics, BSplSLib_Cache for the
    // polynomial kinds, Geom_*Utils for extrusion/revolution/offset, and
    // the mySurface->Eval* fall-through) share the kernel engine encoding;
    // the myTolU/myTolV boundary snap of EvalD1/D2/D3/DN is translated
    // verbatim.
    // ---------------------------------------------------------------------

    /// The OCCT EvalD1/D2/D3/DN boundary snap (L1127-1148, L1271-1292,
    /// L1478-1499, L1706-1727): USide/VSide and the snapped (u, v).
    fn snapped_uv(&self, u: f64, v: f64) -> (f64, f64) {
        let mut us = u;
        let mut vs = v;
        if (u - self.u_first).abs() <= self.my_tol_u {
            us = self.u_first;
        } else if (u - self.u_last).abs() <= self.my_tol_u {
            us = self.u_last;
        }
        if (v - self.v_first).abs() <= self.my_tol_v {
            vs = self.v_first;
        } else if (v - self.v_last).abs() <= self.my_tol_v {
            vs = self.v_last;
        }
        (us, vs)
    }

    /// OCCT GeomAdaptor_Surface::EvalD0 (L1026-1118) — the point at (U, V).
    pub fn value_at(&self, u: f64, v: f64) -> DVec3 {
        // OCCT L1115-1117 default: mySurface->EvalD0 — with the quadric /
        // polynomial / composite arms answering the same values through the
        // kernel engine encoding.
        SurfaceEval::point_at(&self.surface, u, v)
    }

    /// OCCT GeomAdaptor_Surface::EvalD1 (L1122-1262).
    pub fn d1_at(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        let (u, v) = self.snapped_uv(u, v);
        SurfaceEval::derivatives(&self.surface, u, v)
    }

    /// OCCT GeomAdaptor_Surface::EvalD2 (L1266-1469).  The Plane arm
    /// (L1296-1301) answers the D1 pair with zero second partials exactly as
    /// OCCT.
    pub fn d2_at(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        let (u, v) = self.snapped_uv(u, v);
        if self.my_surface_type == GeomAbsSurfaceType::Plane {
            // OCCT L1297-1300: ElSLib::D1 + the zero D2U/D2V/D2UV.
            let (p, d1u, d1v) = SurfaceEval::derivatives(&self.surface, u, v);
            return (p, d1u, d1v, DVec3::ZERO, DVec3::ZERO, DVec3::ZERO);
        }
        SurfaceEval::derivatives2(&self.surface, u, v)
    }

    /// OCCT GeomAdaptor_Surface::EvalD3 (L1473-1693).  The Plane arm
    /// (L1503-1512) answers the D1 pair with zero second and third partials
    /// exactly as OCCT; the remaining kinds ride the kernel finite-difference
    /// scaffold for the third partials (the C++ arms evaluate the exact
    /// ElSLib/BSplSLib third orders — recorded leaf gap for the quadric
    /// third partials).
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
        let (u, v) = self.snapped_uv(u, v);
        if self.my_surface_type == GeomAbsSurfaceType::Plane {
            // OCCT L1504-1511: ElSLib::D1 + the zero partials.
            let (p, d1u, d1v) = SurfaceEval::derivatives(&self.surface, u, v);
            return (
                p,
                d1u,
                d1v,
                DVec3::ZERO,
                DVec3::ZERO,
                DVec3::ZERO,
                DVec3::ZERO,
                DVec3::ZERO,
                DVec3::ZERO,
                DVec3::ZERO,
            );
        }
        // The kernel finite-difference scaffold for D3U/D3V/D3UUV/D3UVV
        // (established encoding of the consumed subset, kept).
        let h = 1e-4;
        let base = SurfaceEval::derivatives2(&self.surface, u, v);
        let d2u = SurfaceEval::derivatives2(&self.surface, u + h, v);
        let d2v = SurfaceEval::derivatives2(&self.surface, u, v + h);
        let puu = base.3;
        let puv = base.4;
        let pvv = base.5;
        (
            base.0,
            base.1,
            base.2,
            puu,
            puv,
            pvv,
            (d2u.3 - puu) / h,
            (d2u.4 - puv) / h,
            (d2v.5 - pvv) / h,
            (d2v.4 - puv) / h,
        )
    }

    /// OCCT GeomAdaptor_Surface::EvalDN (L1697-1814).  Orders 1/2 answer the
    /// D1/D2 partials (the OCCT DN of those orders); the higher orders ride
    /// the untranslated ElSLib::DN / BSplCLib DN leaf (recorded leaf gap).
    pub fn dn_at(&self, u: f64, v: f64, nu: i32, nv: i32) -> DVec3 {
        let (u, v) = self.snapped_uv(u, v);
        if nu == 1 && nv == 0 {
            return self.d1_at(u, v).1;
        }
        if nu == 0 && nv == 1 {
            return self.d1_at(u, v).2;
        }
        panic!("Standard_NotImplemented: GeomAdaptor_Surface::EvalDN ({nu},{nv})")
    }

    // ---------------------------------------------------------------------
    // The global members (GeomAdaptor_Surface.cxx L539-972, L1818-2235).
    // ---------------------------------------------------------------------

    /// OCCT GeomAdaptor_Surface::UContinuity (L539-587).
    pub fn u_continuity_of(&self) -> GeomAbsShape {
        match self.my_surface_type {
            GeomAbsSurfaceType::BSplineSurface => {
                let Surface3::BSpline(a_bspl) = &self.surface else {
                    unreachable!()
                };
                let (tk, tm) = bspline_knots_mults_u(a_bspl);
                local_continuity(
                    a_bspl.degree_u,
                    tm.len() as i32,
                    &tk,
                    &tm,
                    self.u_first,
                    self.u_last,
                    self.is_u_periodic(),
                )
            }
            // OCCT L555-571: the offset basis downshift.  The G1/G2 arms of
            // the OCCT switch cannot be spelled with the rcad 5-level
            // GeomAbsShape and cannot be reached.
            GeomAbsSurfaceType::OffsetSurface => match self.basis_surface_of().u_continuity_of() {
                GeomAbsShape::CN | GeomAbsShape::C3 => GeomAbsShape::CN,
                GeomAbsShape::C2 => GeomAbsShape::C1,
                GeomAbsShape::C1 | GeomAbsShape::C0 => GeomAbsShape::C0,
            },
            GeomAbsSurfaceType::SurfaceOfExtrusion => {
                self.basis_curve_of().continuity()
            }
            // OCCT L575-576: throw Standard_NoSuchObject.
            GeomAbsSurfaceType::OtherSurface => {
                panic!("Standard_NoSuchObject: GeomAdaptor_Surface::UContinuity")
            }
            // OCCT L577-586: the CN kinds.
            _ => GeomAbsShape::CN,
        }
    }

    /// OCCT GeomAdaptor_Surface::VContinuity (L591-639).
    pub fn v_continuity_of(&self) -> GeomAbsShape {
        match self.my_surface_type {
            GeomAbsSurfaceType::BSplineSurface => {
                let Surface3::BSpline(a_bspl) = &self.surface else {
                    unreachable!()
                };
                let (tk, tm) = bspline_knots_mults_v(a_bspl);
                local_continuity(
                    a_bspl.degree_v,
                    tm.len() as i32,
                    &tk,
                    &tm,
                    self.v_first,
                    self.v_last,
                    self.is_v_periodic(),
                )
            }
            GeomAbsSurfaceType::OffsetSurface => match self.basis_surface_of().v_continuity_of() {
                GeomAbsShape::CN | GeomAbsShape::C3 => GeomAbsShape::CN,
                GeomAbsShape::C2 => GeomAbsShape::C1,
                GeomAbsShape::C1 | GeomAbsShape::C0 => GeomAbsShape::C0,
            },
            GeomAbsSurfaceType::SurfaceOfRevolution => {
                self.basis_curve_of().continuity()
            }
            GeomAbsSurfaceType::OtherSurface => {
                panic!("Standard_NoSuchObject: GeomAdaptor_Surface::VContinuity")
            }
            _ => GeomAbsShape::CN,
        }
    }

    /// OCCT GeomAdaptor_Surface::NbUIntervals (L643-697).  The BSpline arm
    /// builds the first-V-knot isocurve adaptor; the isocurve shares the
    /// surface U knot structure, so the interval count rides the same
    /// BSplCLib::Intervals evaluation on the U knots (the pole-magnitude
    /// eps term rides the surface resolution — the kernel encoding of the
    /// isocurve construction).
    pub fn nb_u_intervals_of(&self, s: GeomAbsShape) -> usize {
        match self.my_surface_type {
            GeomAbsSurfaceType::BSplineSurface => {
                let Surface3::BSpline(a_bspl) = &self.surface else {
                    unreachable!()
                };
                let (knots, mults) = bspline_knots_mults_u(a_bspl);
                let cont = self.u_continuity_of();
                // OCCT GeomAdaptor_Curve::NbIntervals L376 on the isocurve.
                if (!self.is_u_periodic() && (s as u8) <= (cont as u8)) || s == GeomAbsShape::C0 {
                    return 1;
                }
                let a_cont: i32 = match s {
                    GeomAbsShape::C1 => 1,
                    GeomAbsShape::C2 => 2,
                    GeomAbsShape::C3 => 3,
                    GeomAbsShape::CN => a_bspl.degree_u as i32,
                    _ => panic!("Standard_DomainError: GeomAdaptor_Surface::NbUIntervals"),
                };
                let an_eps = self
                    .u_resolution(precision::CONFUSION)
                    .min(precision::PCONFUSION);
                return bspl_intervals(
                    &knots,
                    &mults,
                    a_bspl.degree_u,
                    self.is_u_periodic(),
                    a_cont,
                    self.u_first,
                    self.u_last,
                    an_eps,
                )
                .len()
                    - 1;
            }
            GeomAbsSurfaceType::SurfaceOfExtrusion => {
                let my_basis_curve = self.basis_curve_of();
                if my_basis_curve.get_type() == CurveType::BSpline {
                    return my_basis_curve.nb_intervals_of(s);
                }
                1
            }
            GeomAbsSurfaceType::OffsetSurface => {
                // OCCT L664-685: the basis shift + recursion.
                let base_s = intervals_base_shape(s);
                return self.basis_surface_of().nb_u_intervals_of(base_s);
            }
            _ => 1,
        }
    }

    /// OCCT GeomAdaptor_Surface::NbVIntervals (L701-755).
    pub fn nb_v_intervals_of(&self, s: GeomAbsShape) -> usize {
        match self.my_surface_type {
            GeomAbsSurfaceType::BSplineSurface => {
                let Surface3::BSpline(a_bspl) = &self.surface else {
                    unreachable!()
                };
                let (knots, mults) = bspline_knots_mults_v(a_bspl);
                let cont = self.v_continuity_of();
                if (!self.is_v_periodic() && (s as u8) <= (cont as u8)) || s == GeomAbsShape::C0 {
                    return 1;
                }
                let a_cont: i32 = match s {
                    GeomAbsShape::C1 => 1,
                    GeomAbsShape::C2 => 2,
                    GeomAbsShape::C3 => 3,
                    GeomAbsShape::CN => a_bspl.degree_v as i32,
                    _ => panic!("Standard_DomainError: GeomAdaptor_Surface::NbVIntervals"),
                };
                let an_eps = self
                    .v_resolution(precision::CONFUSION)
                    .min(precision::PCONFUSION);
                return bspl_intervals(
                    &knots,
                    &mults,
                    a_bspl.degree_v,
                    self.is_v_periodic(),
                    a_cont,
                    self.v_first,
                    self.v_last,
                    an_eps,
                )
                .len()
                    - 1;
            }
            GeomAbsSurfaceType::SurfaceOfRevolution => {
                let my_basis_curve = self.basis_curve_of();
                if my_basis_curve.get_type() == CurveType::BSpline {
                    return my_basis_curve.nb_intervals_of(s);
                }
                1
            }
            GeomAbsSurfaceType::OffsetSurface => {
                let base_s = intervals_base_shape(s);
                return self.basis_surface_of().nb_v_intervals_of(base_s);
            }
            _ => 1,
        }
    }

    /// OCCT GeomAdaptor_Surface::UIntervals (L759-819) — the T values.
    pub fn u_intervals_of(&self, s: GeomAbsShape) -> Vec<f64> {
        match self.my_surface_type {
            GeomAbsSurfaceType::BSplineSurface => {
                let Surface3::BSpline(a_bspl) = &self.surface else {
                    unreachable!()
                };
                let (knots, mults) = bspline_knots_mults_u(a_bspl);
                let cont = self.u_continuity_of();
                if (!self.is_u_periodic() && (s as u8) <= (cont as u8)) || s == GeomAbsShape::C0 {
                    return vec![self.u_first, self.u_last];
                }
                let a_cont: i32 = match s {
                    GeomAbsShape::C1 => 1,
                    GeomAbsShape::C2 => 2,
                    GeomAbsShape::C3 => 3,
                    GeomAbsShape::CN => a_bspl.degree_u as i32,
                    _ => panic!("Standard_DomainError: GeomAdaptor_Surface::UIntervals"),
                };
                let an_eps = self
                    .u_resolution(precision::CONFUSION)
                    .min(precision::PCONFUSION);
                bspl_intervals(
                    &knots,
                    &mults,
                    a_bspl.degree_u,
                    self.is_u_periodic(),
                    a_cont,
                    self.u_first,
                    self.u_last,
                    an_eps,
                )
            }
            GeomAbsSurfaceType::SurfaceOfExtrusion => {
                let my_basis_curve = self.basis_curve_of();
                if my_basis_curve.get_type() == CurveType::BSpline {
                    return my_basis_curve.intervals_of(s);
                }
                vec![self.u_first, self.u_last]
            }
            GeomAbsSurfaceType::OffsetSurface => {
                let base_s = intervals_base_shape(s);
                return self.basis_surface_of().u_intervals_of(base_s);
            }
            // OCCT L817-819.
            _ => vec![self.u_first, self.u_last],
        }
    }

    /// OCCT GeomAdaptor_Surface::VIntervals (L823-882).
    pub fn v_intervals_of(&self, s: GeomAbsShape) -> Vec<f64> {
        match self.my_surface_type {
            GeomAbsSurfaceType::BSplineSurface => {
                let Surface3::BSpline(a_bspl) = &self.surface else {
                    unreachable!()
                };
                let (knots, mults) = bspline_knots_mults_v(a_bspl);
                let cont = self.v_continuity_of();
                if (!self.is_v_periodic() && (s as u8) <= (cont as u8)) || s == GeomAbsShape::C0 {
                    return vec![self.v_first, self.v_last];
                }
                let a_cont: i32 = match s {
                    GeomAbsShape::C1 => 1,
                    GeomAbsShape::C2 => 2,
                    GeomAbsShape::C3 => 3,
                    GeomAbsShape::CN => a_bspl.degree_v as i32,
                    _ => panic!("Standard_DomainError: GeomAdaptor_Surface::VIntervals"),
                };
                let an_eps = self
                    .v_resolution(precision::CONFUSION)
                    .min(precision::PCONFUSION);
                bspl_intervals(
                    &knots,
                    &mults,
                    a_bspl.degree_v,
                    self.is_v_periodic(),
                    a_cont,
                    self.v_first,
                    self.v_last,
                    an_eps,
                )
            }
            GeomAbsSurfaceType::SurfaceOfRevolution => {
                let my_basis_curve = self.basis_curve_of();
                if my_basis_curve.get_type() == CurveType::BSpline {
                    return my_basis_curve.intervals_of(s);
                }
                vec![self.v_first, self.v_last]
            }
            GeomAbsSurfaceType::OffsetSurface => {
                let base_s = intervals_base_shape(s);
                return self.basis_surface_of().v_intervals_of(base_s);
            }
            // OCCT L880-882.
            _ => vec![self.v_first, self.v_last],
        }
    }

    /// OCCT GeomAdaptor_Surface::UTrim (L886-892).
    pub fn u_trim_of(&self, first: f64, last: f64, tol: f64) -> GeomSurfaceAdaptor {
        GeomSurfaceAdaptor::with_window(
            self.surface.clone(),
            first,
            last,
            self.v_first,
            self.v_last,
            tol,
            self.my_tol_v,
        )
    }

    /// OCCT GeomAdaptor_Surface::VTrim (L896-902).
    pub fn v_trim_of(&self, first: f64, last: f64, tol: f64) -> GeomSurfaceAdaptor {
        GeomSurfaceAdaptor::with_window(
            self.surface.clone(),
            self.u_first,
            self.u_last,
            first,
            last,
            self.my_tol_u,
            tol,
        )
    }

    /// OCCT GeomAdaptor_Surface::IsUClosed (L906-922).
    pub fn is_u_closed_of(&self) -> bool {
        // OCCT L908-911: mySurface->IsUClosed().
        if !SurfaceEval::is_u_closed(&self.surface) {
            return false;
        }
        // OCCT L913-921: the natural bounds comparison.
        let [u1, u2, _v1, _v2] = SurfaceEval::default_domain(&self.surface);
        if SurfaceEval::is_u_periodic(&self.surface) {
            return ((u1 - u2).abs() - (self.u_first - self.u_last).abs()).abs()
                < precision::PCONFUSION;
        }
        (u1 - self.u_first).abs() < precision::PCONFUSION
            && (u2 - self.u_last).abs() < precision::PCONFUSION
    }

    /// OCCT GeomAdaptor_Surface::IsVClosed (L926-942).
    pub fn is_v_closed_of(&self) -> bool {
        if !SurfaceEval::is_v_closed(&self.surface) {
            return false;
        }
        let [_u1, _u2, v1, v2] = SurfaceEval::default_domain(&self.surface);
        if SurfaceEval::is_v_periodic(&self.surface) {
            return ((v1 - v2).abs() - (self.v_first - self.v_last).abs()).abs()
                < precision::PCONFUSION;
        }
        (v1 - self.v_first).abs() < precision::PCONFUSION
            && (v2 - self.v_last).abs() < precision::PCONFUSION
    }

    /// OCCT GeomAdaptor_Surface::IsUPeriodic (L946-949).
    pub fn is_u_periodic_of(&self) -> bool {
        SurfaceEval::is_u_periodic(&self.surface)
    }

    /// OCCT GeomAdaptor_Surface::UPeriod (L953-957): raise when not
    /// periodic, else mySurface->UPeriod().
    pub fn u_period_of(&self) -> f64 {
        if !self.is_u_periodic_of() {
            panic!("Standard_NoSuchObject: GeomAdaptor_Surface::UPeriod");
        }
        // The Geom_*::UPeriod() dispatch: the quadric and revolution /
        // extrusion U periods are 2*pi; the offset passes to the basis; the
        // kernel BSplineSurface carries no periodic flag (never periodic).
        TWO_PI
    }

    /// OCCT GeomAdaptor_Surface::IsVPeriodic (L961-964).
    pub fn is_v_periodic_of(&self) -> bool {
        SurfaceEval::is_v_periodic(&self.surface)
    }

    /// OCCT GeomAdaptor_Surface::VPeriod (L968-972).
    pub fn v_period_of(&self) -> f64 {
        if !self.is_v_periodic_of() {
            panic!("Standard_NoSuchObject: GeomAdaptor_Surface::VPeriod");
        }
        match self.my_surface_type {
            // OCCT Geom_ToroidalSurface::VPeriod = 2*pi.
            GeomAbsSurfaceType::Torus => TWO_PI,
            // The offset passes to the basis surface.
            GeomAbsSurfaceType::OffsetSurface => self.basis_surface_of().v_period_of(),
            _ => TWO_PI,
        }
    }

    /// OCCT GeomAdaptor_Surface::UResolution (L1818-1896).
    pub fn u_resolution_of(&self, r3d: f64) -> f64 {
        let mut res = 0.0f64;
        match self.my_surface_type {
            GeomAbsSurfaceType::SurfaceOfExtrusion => {
                // OCCT L1824-1827: BasisCurve->Resolution(R3d).
                return self.basis_curve_of().resolution(r3d);
            }
            GeomAbsSurfaceType::Torus => {
                // OCCT L1828-1836: R = Major + Minor; Res = R3d / (2*R).
                let Surface3::Torus(s) = &self.surface else {
                    unreachable!()
                };
                let r = s.major_radius + s.minor_radius;
                if r > precision::CONFUSION {
                    res = r3d / (2.0 * r);
                }
            }
            GeomAbsSurfaceType::Sphere => {
                let Surface3::Sphere(s) = &self.surface else {
                    unreachable!()
                };
                let r = s.radius;
                if r > precision::CONFUSION {
                    res = r3d / (2.0 * r);
                }
            }
            GeomAbsSurfaceType::Cylinder => {
                let Surface3::Cylinder(s) = &self.surface else {
                    unreachable!()
                };
                let r = s.radius;
                if r > precision::CONFUSION {
                    res = r3d / (2.0 * r);
                }
            }
            GeomAbsSurfaceType::Cone => {
                // OCCT L1855-1868.
                if self.v_last - self.v_first > 1.0e10 {
                    // Not truly bounded => unknown resolution.
                    return precision::parametric_default(r3d);
                }
                let Surface3::Cone(s) = &self.surface else {
                    unreachable!()
                };
                // OCCT L1861-1865: Rayon = VIso(V)->Radius()
                //                = RefRadius + V * tan(SemiAngle).
                let rayon1 = (s.radius + self.v_last * s.half_angle_rad.tan()).abs();
                let rayon2 = (s.radius + self.v_first * s.half_angle_rad.tan()).abs();
                let r = if rayon1 > rayon2 { rayon1 } else { rayon2 };
                return if r > precision::CONFUSION { r3d / r } else { 0.0 };
            }
            GeomAbsSurfaceType::Plane => {
                // OCCT L1869-1871.
                return r3d;
            }
            GeomAbsSurfaceType::BezierSurface => {
                let Surface3::Bezier(b) = &self.surface else {
                    unreachable!()
                };
                // OCCT L1872-1876: Bezier Resolution -> Ures.
                let (ures, _vres) = crate::math::bspl::bezier_surface_resolution(b, r3d);
                return ures;
            }
            GeomAbsSurfaceType::BSplineSurface => {
                let Surface3::BSpline(b) = &self.surface else {
                    unreachable!()
                };
                // OCCT L1877-1881: BSpline Resolution -> Ures.
                let (ures, _vres) = crate::math::bspl::bspline_surface_resolution(b, r3d);
                return ures;
            }
            GeomAbsSurfaceType::OffsetSurface => {
                // OCCT L1882-1885: BasisAdaptor->UResolution.
                return self.basis_surface_of().u_resolution_of(r3d);
            }
            // OCCT L1886-1888: default -> Precision::Parametric.
            _ => return precision::parametric_default(r3d),
        }

        // OCCT L1890-1895: the angular conversion of the break arms.
        if res <= 1.0 {
            return 2.0 * res.asin();
        }
        2.0 * std::f64::consts::PI
    }

    /// OCCT GeomAdaptor_Surface::VResolution (L1900-1958).
    pub fn v_resolution_of(&self, r3d: f64) -> f64 {
        let mut res = 0.0f64;
        match self.my_surface_type {
            GeomAbsSurfaceType::SurfaceOfRevolution => {
                // OCCT L1906-1909: BasisCurve->Resolution(R3d).
                return self.basis_curve_of().resolution(r3d);
            }
            GeomAbsSurfaceType::Torus => {
                // OCCT L1910-1918: R = MinorRadius.
                let Surface3::Torus(s) = &self.surface else {
                    unreachable!()
                };
                let r = s.minor_radius;
                if r > precision::CONFUSION {
                    res = r3d / (2.0 * r);
                }
            }
            GeomAbsSurfaceType::Sphere => {
                let Surface3::Sphere(s) = &self.surface else {
                    unreachable!()
                };
                let r = s.radius;
                if r > precision::CONFUSION {
                    res = r3d / (2.0 * r);
                }
            }
            // OCCT L1928-1933: Extrusion / Cylinder / Cone / Plane -> R3d.
            GeomAbsSurfaceType::SurfaceOfExtrusion
            | GeomAbsSurfaceType::Cylinder
            | GeomAbsSurfaceType::Cone
            | GeomAbsSurfaceType::Plane => {
                return r3d;
            }
            GeomAbsSurfaceType::BezierSurface => {
                let Surface3::Bezier(b) = &self.surface else {
                    unreachable!()
                };
                let (_ures, vres) = crate::math::bspl::bezier_surface_resolution(b, r3d);
                return vres;
            }
            GeomAbsSurfaceType::BSplineSurface => {
                let Surface3::BSpline(b) = &self.surface else {
                    unreachable!()
                };
                let (_ures, vres) = crate::math::bspl::bspline_surface_resolution(b, r3d);
                return vres;
            }
            GeomAbsSurfaceType::OffsetSurface => {
                return self.basis_surface_of().v_resolution_of(r3d);
            }
            _ => return precision::parametric_default(r3d),
        }

        if res <= 1.0 {
            return 2.0 * res.asin();
        }
        2.0 * std::f64::consts::PI
    }

    // ---------------------------------------------------------------------
    // The geometry downcasts (GeomAdaptor_Surface.cxx L1962-2235).
    // ---------------------------------------------------------------------

    /// OCCT GeomAdaptor_Surface::Plane (L1962-1969) — the gp_Pln payload.
    pub fn plane(&self) -> Plane {
        if self.my_surface_type != GeomAbsSurfaceType::Plane {
            panic!("Standard_NoSuchObject: GeomAdaptor_Surface::Plane");
        }
        let Surface3::Plane(g) = &self.surface else {
            unreachable!()
        };
        *g
    }

    /// OCCT GeomAdaptor_Surface::Cylinder (L1973-1980).
    pub fn cylinder(&self) -> CylindricalSurface {
        if self.my_surface_type != GeomAbsSurfaceType::Cylinder {
            panic!("Standard_NoSuchObject: GeomAdaptor_Surface::Cylinder");
        }
        let Surface3::Cylinder(g) = &self.surface else {
            unreachable!()
        };
        *g
    }

    /// OCCT GeomAdaptor_Surface::Cone (L1984-1991).
    pub fn cone(&self) -> ConicalSurface {
        if self.my_surface_type != GeomAbsSurfaceType::Cone {
            panic!("Standard_NoSuchObject: GeomAdaptor_Surface::Cone");
        }
        let Surface3::Cone(g) = &self.surface else {
            unreachable!()
        };
        *g
    }

    /// OCCT GeomAdaptor_Surface::Sphere (L1995-2002).
    pub fn sphere(&self) -> SphericalSurface {
        if self.my_surface_type != GeomAbsSurfaceType::Sphere {
            panic!("Standard_NoSuchObject: GeomAdaptor_Surface::Sphere");
        }
        let Surface3::Sphere(g) = &self.surface else {
            unreachable!()
        };
        *g
    }

    /// OCCT GeomAdaptor_Surface::Torus (L2006-2013).
    pub fn torus(&self) -> ToroidalSurface {
        if self.my_surface_type != GeomAbsSurfaceType::Torus {
            panic!("Standard_NoSuchObject: GeomAdaptor_Surface::Torus");
        }
        let Surface3::Torus(g) = &self.surface else {
            unreachable!()
        };
        *g
    }

    /// OCCT GeomAdaptor_Surface::UDegree (L2017-2032).
    pub fn u_degree_of(&self) -> usize {
        if self.my_surface_type == GeomAbsSurfaceType::BSplineSurface {
            let Surface3::BSpline(b) = &self.surface else {
                unreachable!()
            };
            return b.degree_u;
        }
        if self.my_surface_type == GeomAbsSurfaceType::BezierSurface {
            let Surface3::Bezier(b) = &self.surface else {
                unreachable!()
            };
            return b.control_points.len() - 1;
        }
        if self.my_surface_type == GeomAbsSurfaceType::SurfaceOfExtrusion {
            return self.basis_curve_of().degree();
        }
        panic!("Standard_NoSuchObject: GeomAdaptor_Surface::UDegree");
    }

    /// OCCT GeomAdaptor_Surface::NbUPoles (L2036-2051).
    pub fn nb_u_poles_of(&self) -> usize {
        if self.my_surface_type == GeomAbsSurfaceType::BSplineSurface {
            let Surface3::BSpline(b) = &self.surface else {
                unreachable!()
            };
            return b.control_points.len();
        }
        if self.my_surface_type == GeomAbsSurfaceType::BezierSurface {
            let Surface3::Bezier(b) = &self.surface else {
                unreachable!()
            };
            return b.control_points.len();
        }
        if self.my_surface_type == GeomAbsSurfaceType::SurfaceOfExtrusion {
            return self.basis_curve_of().nb_poles();
        }
        panic!("Standard_NoSuchObject: GeomAdaptor_Surface::NbUPoles");
    }

    /// OCCT GeomAdaptor_Surface::VDegree (L2055-2070).
    pub fn v_degree_of(&self) -> usize {
        if self.my_surface_type == GeomAbsSurfaceType::BSplineSurface {
            let Surface3::BSpline(b) = &self.surface else {
                unreachable!()
            };
            return b.degree_v;
        }
        if self.my_surface_type == GeomAbsSurfaceType::BezierSurface {
            let Surface3::Bezier(b) = &self.surface else {
                unreachable!()
            };
            return b.control_points[0].len() - 1;
        }
        if self.my_surface_type == GeomAbsSurfaceType::SurfaceOfRevolution {
            return self.basis_curve_of().degree();
        }
        panic!("Standard_NoSuchObject: GeomAdaptor_Surface::VDegree");
    }

    /// OCCT GeomAdaptor_Surface::NbVPoles (L2074-2089).
    pub fn nb_v_poles_of(&self) -> usize {
        if self.my_surface_type == GeomAbsSurfaceType::BSplineSurface {
            let Surface3::BSpline(b) = &self.surface else {
                unreachable!()
            };
            return b.control_points[0].len();
        }
        if self.my_surface_type == GeomAbsSurfaceType::BezierSurface {
            let Surface3::Bezier(b) = &self.surface else {
                unreachable!()
            };
            return b.control_points[0].len();
        }
        if self.my_surface_type == GeomAbsSurfaceType::SurfaceOfRevolution {
            return self.basis_curve_of().nb_poles();
        }
        panic!("Standard_NoSuchObject: GeomAdaptor_Surface::NbVPoles");
    }

    /// OCCT GeomAdaptor_Surface::NbUKnots (L2093-2104).
    pub fn nb_u_knots_of(&self) -> usize {
        if self.my_surface_type == GeomAbsSurfaceType::BSplineSurface {
            let Surface3::BSpline(b) = &self.surface else {
                unreachable!()
            };
            return bspline_knots_mults_u(b).0.len();
        }
        if self.my_surface_type == GeomAbsSurfaceType::SurfaceOfExtrusion {
            return self.basis_curve_of().nb_knots();
        }
        panic!("Standard_NoSuchObject: GeomAdaptor_Surface::NbUKnots");
    }

    /// OCCT GeomAdaptor_Surface::NbVKnots (L2108-2115).
    pub fn nb_v_knots_of(&self) -> usize {
        if self.my_surface_type == GeomAbsSurfaceType::BSplineSurface {
            let Surface3::BSpline(b) = &self.surface else {
                unreachable!()
            };
            return bspline_knots_mults_v(b).0.len();
        }
        panic!("Standard_NoSuchObject: GeomAdaptor_Surface::NbVKnots");
    }

    /// OCCT GeomAdaptor_Surface::IsURational (L2119-2130).
    pub fn is_u_rational_of(&self) -> bool {
        if self.my_surface_type == GeomAbsSurfaceType::BSplineSurface {
            let Surface3::BSpline(b) = &self.surface else {
                unreachable!()
            };
            return surface_rational_flags(&b.weights).0;
        }
        if self.my_surface_type == GeomAbsSurfaceType::BezierSurface {
            let Surface3::Bezier(b) = &self.surface else {
                unreachable!()
            };
            return surface_rational_flags(&b.weights).0;
        }
        false
    }

    /// OCCT GeomAdaptor_Surface::IsVRational (L2134-2145).
    pub fn is_v_rational_of(&self) -> bool {
        if self.my_surface_type == GeomAbsSurfaceType::BSplineSurface {
            let Surface3::BSpline(b) = &self.surface else {
                unreachable!()
            };
            return surface_rational_flags(&b.weights).1;
        }
        if self.my_surface_type == GeomAbsSurfaceType::BezierSurface {
            let Surface3::Bezier(b) = &self.surface else {
                unreachable!()
            };
            return surface_rational_flags(&b.weights).1;
        }
        false
    }

    /// OCCT GeomAdaptor_Surface::Bezier (L2149-2156) — the
    /// Geom_BezierSurface payload (no copy).
    pub fn bezier(&self) -> crate::geom::BezierSurface {
        if self.my_surface_type != GeomAbsSurfaceType::BezierSurface {
            panic!("Standard_NoSuchObject: GeomAdaptor_Surface::Bezier");
        }
        let Surface3::Bezier(b) = &self.surface else {
            unreachable!()
        };
        b.clone()
    }

    /// OCCT GeomAdaptor_Surface::BSpline (L2160-2167) — the
    /// Geom_BSplineSurface payload (no copy).
    pub fn bspline(&self) -> crate::geom::BSplineSurface {
        if self.my_surface_type != GeomAbsSurfaceType::BSplineSurface {
            panic!("Standard_NoSuchObject: GeomAdaptor_Surface::BSpline");
        }
        let Surface3::BSpline(b) = &self.surface else {
            unreachable!()
        };
        b.clone()
    }

    /// OCCT GeomAdaptor_Surface::AxeOfRevolution (L2171-2178) — the
    /// (location, direction) pair of the gp_Ax1.
    pub fn axe_of_revolution(&self) -> (DVec3, DVec3) {
        if self.my_surface_type != GeomAbsSurfaceType::SurfaceOfRevolution {
            panic!("Standard_NoSuchObject: GeomAdaptor_Surface::AxeOfRevolution");
        }
        let Surface3::Revolution(r) = &self.surface else {
            unreachable!()
        };
        (r.axis_origin, r.axis_dir)
    }

    /// OCCT GeomAdaptor_Surface::Direction (L2182-2189) — the extrusion
    /// direction.
    pub fn direction(&self) -> DVec3 {
        if self.my_surface_type != GeomAbsSurfaceType::SurfaceOfExtrusion {
            panic!("Standard_NoSuchObject: GeomAdaptor_Surface::Direction");
        }
        let Surface3::LinearExtrusion(e) = &self.surface else {
            unreachable!()
        };
        e.direction
    }

    /// OCCT GeomAdaptor_Surface::BasisCurve (L2193-2209) — the basis curve
    /// adaptor of an extrusion / revolution surface.
    pub fn basis_curve_of(&self) -> GeomCurveAdaptor {
        let c = match &self.surface {
            // OCCT L2196-2199: extrusion -> BasisCurve().
            Surface3::LinearExtrusion(e) => e.profile.as_ref().clone(),
            // OCCT L2200-2203: revolution -> BasisCurve().
            Surface3::Revolution(r) => r.profile.as_ref().clone(),
            // OCCT L2204-2207: else throw.
            _ => panic!("Standard_NoSuchObject: GeomAdaptor_Surface::BasisCurve"),
        };
        GeomCurveAdaptor::new(c)
    }

    /// OCCT GeomAdaptor_Surface::BasisSurface (L2213-2224) — the basis
    /// surface adaptor of an offset surface (same window).
    pub fn basis_surface_of(&self) -> GeomSurfaceAdaptor {
        if self.my_surface_type != GeomAbsSurfaceType::OffsetSurface {
            panic!("Standard_NoSuchObject: GeomAdaptor_Surface::BasisSurface");
        }
        let Surface3::Offset(o) = &self.surface else {
            unreachable!()
        };
        GeomSurfaceAdaptor::with_window(
            o.basis.as_ref().clone(),
            self.u_first,
            self.u_last,
            self.v_first,
            self.v_last,
            self.my_tol_u,
            self.my_tol_v,
        )
    }

    /// OCCT GeomAdaptor_Surface::OffsetValue (L2228-2235).
    pub fn offset_value_of(&self) -> f64 {
        if self.my_surface_type != GeomAbsSurfaceType::OffsetSurface {
            panic!("Standard_NoSuchObject: GeomAdaptor_Surface::BasisSurface");
        }
        let Surface3::Offset(o) = &self.surface else {
            unreachable!()
        };
        o.offset_distance
    }

    /// OCCT GeomAdaptor_Surface::ShallowCopy (L331-398) — the value copy.
    pub fn shallow_copy_of(&self) -> Self {
        self.clone()
    }
}

// =========================================================================
// The BSpline knot-structure readers
// =========================================================================

/// The (knots, mults) U structure of the kernel BSpline surface (the OCCT
/// `UKnots()` / `UMultiplicities()` pair).
fn bspline_knots_mults_u(s: &crate::geom::BSplineSurface) -> (Vec<f64>, Vec<i32>) {
    compress_knots(&s.knots_u)
}

/// The (knots, mults) V structure of the kernel BSpline surface.
fn bspline_knots_mults_v(s: &crate::geom::BSplineSurface) -> (Vec<f64>, Vec<i32>) {
    compress_knots(&s.knots_v)
}

/// Compress a flat knot vector into the OCCT (knots, mults) form.
fn compress_knots(flat: &[f64]) -> (Vec<f64>, Vec<i32>) {
    let mut knots: Vec<f64> = Vec::new();
    let mut mults: Vec<i32> = Vec::new();
    for (i, k) in flat.iter().enumerate() {
        if i > 0 && *k == knots.last().copied().unwrap_or(f64::NAN) {
            *mults.last_mut().unwrap() += 1;
        } else {
            knots.push(*k);
            mults.push(1);
        }
    }
    (knots, mults)
}

/// OCCT NbUIntervals / NbVIntervals / UIntervals / VIntervals offset arm
/// (L665-684, L723-742, L784-803, L847-866): the S -> BaseS downshift.
fn intervals_base_shape(s: GeomAbsShape) -> GeomAbsShape {
    match s {
        // OCCT: G1/G2 throw Standard_DomainError — unreachable with the rcad
        // 5-level GeomAbsShape.
        GeomAbsShape::C0 => GeomAbsShape::C1,
        GeomAbsShape::C1 => GeomAbsShape::C2,
        GeomAbsShape::C2 => GeomAbsShape::C3,
        GeomAbsShape::C3 | GeomAbsShape::CN => GeomAbsShape::CN,
    }
}

// =========================================================================
// The trait implementations
// =========================================================================

impl Adaptor3dSurface for GeomSurfaceAdaptor {
    fn first_u_parameter(&self) -> f64 {
        self.u_first
    }

    fn last_u_parameter(&self) -> f64 {
        self.u_last
    }

    fn first_v_parameter(&self) -> f64 {
        self.v_first
    }

    fn last_v_parameter(&self) -> f64 {
        self.v_last
    }

    fn value(&self, u: f64, v: f64) -> DVec3 {
        self.value_at(u, v)
    }

    fn d1(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        self.d1_at(u, v)
    }

    fn d2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
        self.d2_at(u, v)
    }

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

    fn dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> DVec3 {
        self.dn_at(u, v, nu, nv)
    }

    fn u_resolution(&self, r3d: f64) -> f64 {
        self.u_resolution_of(r3d)
    }

    fn v_resolution(&self, r3d: f64) -> f64 {
        self.v_resolution_of(r3d)
    }

    fn get_type(&self) -> GeomAbsSurfaceType {
        self.my_surface_type
    }

    fn is_u_periodic(&self) -> bool {
        self.is_u_periodic_of()
    }

    fn u_period(&self) -> f64 {
        self.u_period_of()
    }

    fn is_v_periodic(&self) -> bool {
        self.is_v_periodic_of()
    }

    fn v_period(&self) -> f64 {
        self.v_period_of()
    }

    fn u_continuity(&self) -> GeomAbsShape {
        self.u_continuity_of()
    }

    fn v_continuity(&self) -> GeomAbsShape {
        self.v_continuity_of()
    }

    fn nb_u_intervals(&self, s: GeomAbsShape) -> usize {
        self.nb_u_intervals_of(s)
    }

    fn nb_v_intervals(&self, s: GeomAbsShape) -> usize {
        self.nb_v_intervals_of(s)
    }

    fn u_intervals(&self, s: GeomAbsShape) -> Vec<f64> {
        self.u_intervals_of(s)
    }

    fn v_intervals(&self, s: GeomAbsShape) -> Vec<f64> {
        self.v_intervals_of(s)
    }

    fn u_trim(&self, u1: f64, u2: f64, eps: f64) -> Arc<dyn Adaptor3dSurface> {
        Arc::new(self.u_trim_of(u1, u2, eps))
    }

    fn v_trim(&self, v1: f64, v2: f64, eps: f64) -> Arc<dyn Adaptor3dSurface> {
        Arc::new(self.v_trim_of(v1, v2, eps))
    }

    fn is_u_closed(&self) -> bool {
        self.is_u_closed_of()
    }

    fn is_v_closed(&self) -> bool {
        self.is_v_closed_of()
    }

    fn u_degree(&self) -> usize {
        self.u_degree_of()
    }

    fn nb_u_poles(&self) -> usize {
        self.nb_u_poles_of()
    }

    fn v_degree(&self) -> usize {
        self.v_degree_of()
    }

    fn nb_v_poles(&self) -> usize {
        self.nb_v_poles_of()
    }

    fn nb_u_knots(&self) -> usize {
        self.nb_u_knots_of()
    }

    fn nb_v_knots(&self) -> usize {
        self.nb_v_knots_of()
    }

    fn is_u_rational(&self) -> bool {
        self.is_u_rational_of()
    }

    fn is_v_rational(&self) -> bool {
        self.is_v_rational_of()
    }

    fn bezier(&self) -> crate::geom::BezierSurface {
        self.bezier()
    }

    fn bspline(&self) -> crate::geom::BSplineSurface {
        self.bspline()
    }

    fn direction(&self) -> DVec3 {
        self.direction()
    }

    fn basis_curve(&self) -> Arc<dyn Adaptor3dCurve> {
        Arc::new(self.basis_curve_of())
    }

    fn basis_surface(&self) -> Arc<dyn Adaptor3dSurface> {
        Arc::new(self.basis_surface_of())
    }

    fn offset_value(&self) -> f64 {
        self.offset_value_of()
    }

    fn shallow_copy(&self) -> Arc<dyn Adaptor3dSurface> {
        Arc::new(self.shallow_copy_of())
    }

    fn kernel_surface(&self) -> Option<&Surface3> {
        Some(&self.surface)
    }
}

/// The OCCT Adaptor3d_Surface geometry accessors (Plane() / Cylinder() /
/// Cone() / Sphere() / Torus() / AxeOfRevolution()) — the companion trait
/// consumed by the ProjLib Perform dispatch.
pub trait Adaptor3dSurfaceGeom: Adaptor3dSurface {
    /// OCCT Plane() — valid when GetType() == GeomAbs_Plane.
    fn plane(&self) -> Plane;
    /// OCCT Cylinder() — valid when GetType() == GeomAbs_Cylinder.
    fn cylinder(&self) -> CylindricalSurface;
    /// OCCT Cone() — valid when GetType() == GeomAbs_Cone.
    fn cone(&self) -> ConicalSurface;
    /// OCCT Sphere() — valid when GetType() == GeomAbs_Sphere.
    fn sphere(&self) -> SphericalSurface;
    /// OCCT Torus() — valid when GetType() == GeomAbs_Torus.
    fn torus(&self) -> ToroidalSurface;
    /// OCCT AxeOfRevolution() — the (location, direction) of the revolution
    /// axis (gp_Ax1).
    fn axe_of_revolution(&self) -> (DVec3, DVec3) {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::AxeOfRevolution")
    }
}

impl Adaptor3dSurfaceGeom for GeomSurfaceAdaptor {
    fn plane(&self) -> Plane {
        self.plane()
    }

    fn cylinder(&self) -> CylindricalSurface {
        self.cylinder()
    }

    fn cone(&self) -> ConicalSurface {
        self.cone()
    }

    fn sphere(&self) -> SphericalSurface {
        self.sphere()
    }

    fn torus(&self) -> ToroidalSurface {
        self.torus()
    }

    fn axe_of_revolution(&self) -> (DVec3, DVec3) {
        self.axe_of_revolution()
    }
}
