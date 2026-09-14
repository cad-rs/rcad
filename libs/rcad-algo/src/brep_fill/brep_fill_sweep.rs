//! OCCT BRepFill_Sweep (TKBool/BRepFill) — 1:1 translation of
//! BRepFill_Sweep.hxx (L43-189) + BRepFill_Sweep.cxx, part A (this file):
//! the class shell, the simple setters/accessors, SetBounds, BuildWire,
//! MergeVertex / UpdateVertex / EvalExtrapol and the self-contained file
//! statics (Translate / Box / BuildVertex / NullEdge / HasPCurves /
//! ReverseEdgeInFirstOrLastWire / ReverseModifiedEdges).
//!
//! Part B (re-dispatch, listed in the delivery report) — the remaining
//! translation units of BRepFill_Sweep.cxx: NumberOfPoles, Couture,
//! CheckSameParameter / CheckSameParameterExact / SameParameter /
//! CorrectSameParameter, Oriente, UpdateEdgeOnPlane, BuildFace, BuildEdge
//! (both forms), Filling, Substitute, KeepEdge, UpdateEdge, IsDegen,
//! CorrectApproxParameters, BuildShell, Build, PerformCorner,
//! RebuildTopOrBottomEdge.
//!
//! Architecture differences:
//! - `handle(BRepFill_SectionLaw) mySec` /
//!   `handle(BRepFill_LocationLaw) myLoc` map to
//!   `Rc<RefCell<dyn ...Ops>>` (the trait slots of the brep_fill law files).
//! - `NCollection_HArray2<TopoDS_Shape>` maps to
//!   [`crate::brep_fill::brep_fill_pipe_shell_b::ShapeHArray2`]
//!   (`Vec<Vec<Shape>>`, row-major, 1-based OCCT indexing -> `[r-1][c-1]`).
//! - `NCollection_Map / NCollection_DataMap` keyed by
//!   TopTools_ShapeMapHasher map to HashSet/HashMap keyed by `ShapeKey`.
//! - `BRep_Builder` maps to `rcad_kernel::topo::topods::BRepBuilder`; the
//!   owning pool is the leading `brep` argument.
//! - `BRep_Tool::Tolerance / Pnt` map to the stored TVertexData fields.
//! - `GeomAdaptor/BRepAdaptor` surface bounds map to
//!   `surface_adaptor_basis_and_bounds`; the UIso/VIso dispatch is the local
//!   ElSLib / BSpline re-host ([`surface_uiso`] / [`surface_viso`]).

use std::cell::RefCell;
use std::collections::HashSet;
use std::rc::Rc;

use glam::DVec3;

use rcad_kernel::base::extrema_ext_elc::epsilon_of;
use rcad_kernel::base::proj_lib::adaptor::{Adaptor3dSurface, GeomAbsSurfaceType};
use rcad_kernel::base::proj_lib::geom_adaptor_surface::GeomSurfaceAdaptor;
use rcad_kernel::geom::extrusion_utils;
use rcad_kernel::geom::offset_surface_utils::{
    offset_basis_and_value, offset_equivalent_surface, offset_surface_eval_d0,
    offset_surface_eval_d1, offset_surface_osculating,
};
use rcad_kernel::geom::osculating_surface::OsculatingSurface;
use rcad_kernel::geom::{
    BezierCurve3, BezierSurface, BSplineCurve3, BSplineSurface, Circle3, Curve3, CurveEval, Ellipse3,
    Line3, Surface3, TrimmedCurve3,
};
use rcad_kernel::base::proj_lib::elslib_iso::{
    elslib_cone_u_iso, elslib_cone_v_iso, elslib_cylinder_u_iso, elslib_cylinder_v_iso,
    elslib_plane_u_iso, elslib_plane_v_iso, elslib_sphere_u_iso, elslib_sphere_v_iso,
    elslib_torus_u_iso, elslib_torus_v_iso, Ax3View,
};
use rcad_kernel::math::adv_approx::{ApproxAFunction, EvaluatorFunction};
use rcad_kernel::math::bspl_lib::bspl_slib_iso;
use rcad_kernel::math::gp::GP_RESOLUTION;
use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::topo::topods::{surface_adaptor_basis_and_bounds, BRep, BRepBuilder, Shape};

use crate::brep_fill::brep_fill_location_law::BRepFillLocationLawOps;
use crate::brep_fill::brep_fill_pipe_shell_b::{
    ShapeHArray2, ShapeSet,
};
use crate::brep_fill::generator::{shape_reversed, ShapeKey};
use crate::geomalgo::geomfill::section_law::SectionLaw;
use crate::geomalgo::geomfill::sweep::GeomFillApproxStyle;

/// OCCT BRepFill_TransitionStyle (BRepFill_TransitionStyle.hxx).
pub use crate::brep_fill::brep_fill_pipe_shell_b::BRepFillTransitionStyle;

// ---------------------------------------------------------------------------
// UIso / VIso re-hosts (OCCT Geom_Surface::UIso / VIso dispatch)
// ---------------------------------------------------------------------------

/// OCCT ElSLib + Geom_BSplineSurface::UIso + Geom_BezierSurface::UIso over the
/// rcad Surface3 — the u-varying iso curve at `u`.  The Offset variant keeps
/// the OCCT failure path (the missing OCCT dependency is named at the arm),
/// and the rcad-only variants (Ellipsoid / Helicoid / Pipe / Ruled / Coons /
/// TriBezier — no Geom_Surface override exists in TKG3d/Geom) land on the
/// catch-all.  `pub(crate)`: Geom_Surface::UIso is a public Geom-level
/// operation; SplitSurf (ChFi3d_FilBuilder.cxx L2291-2292) consumes it over
/// the stored blend surface.
pub(crate) fn surface_uiso(surf: &Surface3, u: f64) -> Curve3 {
    // OCCT Geom_*::UIso are thin wrappers over the ElSLib U-iso constructors
    // (Geom_Plane.cxx L261-265, Geom_CylindricalSurface.cxx L294-298,
    // Geom_SphericalSurface.cxx L292-297, Geom_ConicalSurface.cxx L337-341,
    // Geom_ToroidalSurface.cxx L305-309).  Delegate to the single kernel
    // ElSLib translation (`base/proj_lib/elslib_iso.rs`) instead of
    // re-deriving the frames here.
    let ax3 = |location: DVec3, direction: DVec3, x_direction: DVec3| {
        Ax3View::from_axes(location, direction, x_direction)
    };
    match surf {
        // Geom_Plane::UIso = ElSLib::PlaneUIso.
        Surface3::Plane(pl) => Curve3::Line(elslib_plane_u_iso(
            &ax3(pl.origin, pl.normal, pl.u_dir),
            u,
        )),
        // Geom_CylindricalSurface::UIso = ElSLib::CylinderUIso: the ruling
        // line at longitude U, anchored at P(U, 0) with direction
        // ConeD1/CylinderD1's DV (= the axis).
        Surface3::Cylinder(cy) => {
            let pos = match cy.y_dir {
                Some(y) => Ax3View::with_y_dir(cy.origin, cy.axis, cy.ref_dir, y),
                None => ax3(cy.origin, cy.axis, cy.ref_dir),
            };
            Curve3::Line(elslib_cylinder_u_iso(&pos, cy.radius, u))
        }
        // Geom_SphericalSurface::UIso = ElSLib::SphereUIso wrapped in
        // Geom_TrimmedCurve(GC, -M_PI/2, M_PI/2) so that the iso is
        // parameterised by latitude.
        Surface3::Sphere(sp) => {
            let circ = elslib_sphere_u_iso(&ax3(sp.center, sp.axis, sp.ref_dir), sp.radius, u);
            Curve3::Trimmed(TrimmedCurve3 {
                curve: Box::new(Curve3::Circle(circ)),
                first: -0.5 * std::f64::consts::PI,
                last: 0.5 * std::f64::consts::PI,
            })
        }
        // Geom_ConicalSurface::UIso = ElSLib::ConeUIso: the ruling line at
        // longitude U.  The rcad payload's `apex` is the reference-circle
        // centre (= OCCT `pos.Location()`); `apex_point()` is the true apex.
        Surface3::Cone(co) => Curve3::Line(elslib_cone_u_iso(
            &ax3(co.apex, co.axis, co.ref_dir),
            co.radius,
            co.half_angle_rad,
            u,
        )),
        // Geom_ToroidalSurface::UIso = ElSLib::TorusUIso.
        Surface3::Torus(to) => Curve3::Circle(elslib_torus_u_iso(
            &ax3(to.center, to.axis, to.ref_dir),
            to.major_radius,
            to.minor_radius,
            u,
        )),
        // Geom_BSplineSurface::UIso: the poles of the iso are the
        // homogeneous De Boor evaluation of the u basis per V column.
        Surface3::BSpline(bs) => Curve3::BSpline(bspline_surface_uiso(bs, u)),
        // OCCT Geom_BezierSurface::UIso (Geom_BezierSurface.cxx L1769-1810).
        Surface3::Bezier(bz) => bezier_surface_uiso(bz, u),
        // OCCT Geom_SurfaceOfLinearExtrusion::UIso
        // (Geom_SurfaceOfLinearExtrusion.cxx L275-281): the ruling line
        // (basisCurve->Value(U), direction).
        Surface3::LinearExtrusion(le) => Curve3::Line(Line3 {
            origin: le.profile.point_at(u),
            direction: le.direction,
        }),
        // OCCT Geom_SurfaceOfRevolution::UIso (Geom_SurfaceOfRevolution.cxx
        // L372-381): the basis curve rotated by U about the axis
        // (C->Rotate(Ax1(loc, direction), U)).
        Surface3::Revolution(rev) => {
            curve3_rotated_about_axis(&rev.profile, rev.axis_origin, rev.axis_dir, u)
        }
        // OCCT Geom_RectangularTrimmedSurface::UIso
        // (Geom_RectangularTrimmedSurface.cxx L444-458): the basis UIso,
        // restricted to the v-trim range when isvtrimmed (the rcad flag maps
        // to "the stored v-trim bounds differ from the basis natural domain").
        Surface3::Trimmed(ts) => {
            let c = surface_uiso(&ts.basis, u);
            if basis_v_bounds_of(&ts.basis) != (ts.trim[2], ts.trim[3]) {
                Curve3::Trimmed(TrimmedCurve3 {
                    curve: Box::new(c),
                    first: ts.trim[2],
                    last: ts.trim[3],
                })
            } else {
                c
            }
        }
        // OCCT Geom_OffsetSurface::UIso (Geom_OffsetSurface.cxx L601-655).
        Surface3::Offset(of) => {
            // L603-604: `const handle(Geom_Surface) anEquivSurface =
            // directRepSurface(*this)` — the eval representation that
            // SetBasisSurface sets to makeFullSurfaceRep(Surface())
            // (L255-256).  The rcad OffsetSurface payload carries no eval
            // representation member, so the equivalent surface is recomputed
            // on demand (offset_equivalent_surface = Geom_OffsetSurface::
            // Surface(), L867-993).
            let an_equiv_surface =
                offset_equivalent_surface(of.basis.as_ref(), of.offset_distance);
            match an_equiv_surface {
                // L652: `return anEquivSurface->UIso(UU);`.
                Some(an_equiv) => surface_uiso(&an_equiv, u),
                None => {
                    // OCCT SetBasisSurface L246-253 keeps `basisSurf` as the
                    // (possibly re-wrapped) trimmed basis; the rcad value
                    // model uses the unwrapped basis of the
                    // offset_basis_and_value unwrap, whose evaluation the
                    // kernel Geom_Surface leaf accepts.
                    let (basis, offset_value) =
                        offset_basis_and_value(of.basis.as_ref(), of.offset_distance);
                    // L606: `GeomAdaptor_Surface aGAsurf(basisSurf);`.
                    let a_ga_surf = GeomSurfaceAdaptor::new(basis.clone());
                    if a_ga_surf.get_type() == GeomAbsSurfaceType::SurfaceOfExtrusion {
                        // L609: `handle(Geom_Curve) aL = basisSurf->UIso(UU);`.
                        let a_l = surface_uiso(&basis, u);
                        // L612-613: `basisSurf->D1(UU, 0., aP, aD1U, aD1V);`.
                        let a_d1 = extrusion_utils::surface_eval_d1(&basis, u, 0.0);
                        // L614: `gp_Vec aDir = aD1U.Crossed(aD1V);`.
                        let a_dir = a_d1.d1u.cross(a_d1.d1v);
                        // L615-618: `if (aDir.SquareMagnitude() <
                        // gp::Resolution()) return aL;`.
                        if a_dir.length_squared() < GP_RESOLUTION {
                            return a_l;
                        }
                        // L619-620: `aDir.Normalize(); aDir *= offsetValue;`.
                        let a_dir = a_dir.normalize() * offset_value;
                        // L622: `aL->Translate(aDir); return aL;`.
                        return curve3_translated(&a_l, a_dir);
                    }

                    // L625-635: the general approximation arm.
                    let num1 = 0;
                    let num2 = 0;
                    let num3 = 1;
                    // T3 = HArray1<double>(1, Num3) with
                    // T3->Init(Precision::Approximation()).
                    let t3 = [rcad_kernel::core::precision::APPROXIMATION];
                    // Bounds(U1, U2, V1, V2) — Geom_OffsetSurface::Bounds
                    // (L313-316) delegates to basisSurf.
                    let bounds = rcad_kernel::geom::SurfaceEval::default_domain(&basis);
                    let (_u1, _u2, v1, v2) = (bounds[0], bounds[1], bounds[2], bounds[3]);
                    let cont = GeomAbsShape::C1;
                    let max_seg = 100;
                    let max_deg = 14;

                    // L633: `Geom_OffsetSurface_UIsoEvaluator ev(*this, UU);`.
                    let mut ev = GeomOffsetSurfaceUIsoEvaluator {
                        basis,
                        offset: offset_value,
                        osc: offset_surface_osculating(of.basis.as_ref()),
                        my_iso_par: u,
                    };
                    // L634-635: the AdvApprox_ApproxAFunction ctor + Perform
                    // over [V1, V2].
                    let approx = ApproxAFunction::new(
                        num1,
                        num2,
                        num3,
                        None,
                        None,
                        Some(&t3),
                        v1,
                        v2,
                        cont,
                        max_deg,
                        max_seg,
                        &mut ev,
                    );
                    // L637: Standard_ConstructionError_Raise_if(!Approx
                    // .IsDone(), " Geom_OffsetSurface : UIso").
                    assert!(
                        approx.is_done(),
                        "Standard_ConstructionError: Geom_OffsetSurface : UIso"
                    );
                    // L639-650: the Geom_BSplineCurve of the approximation.
                    approx_to_bspline_curve(&approx)
                }
            }
        }
        // Ellipsoid / Helicoid / Pipe / Ruled / Coons / TriBezier: rcad-only
        // Surface3 variants with no Geom_Surface override in TKG3d/Geom, so
        // there is no OCCT body to translate (the OCCT-faithful failure path).
        _ => panic!(
            "GAP: Geom_*Surface::UIso (TKMath/TKG3d kernel re-host) is not \
             translated for this surface type — BRepFill_Sweep::BuildWire"
        ),
    }
}

/// The natural v-bounds of a basis surface (the OCCT `isvtrimmed` check
/// compares the trim bounds against the basis domain).
fn basis_v_bounds_of(surf: &Surface3) -> (f64, f64) {
    use rcad_kernel::geom::SurfaceEval;
    let d = surf.default_domain();
    (d[2], d[3])
}

/// OCCT ElSLib + Geom_BSplineSurface::VIso + Geom_BezierSurface::VIso over the
/// rcad Surface3 — the v-varying iso curve at `v`.  Same GAP note as
/// [`surface_uiso`] (the Offset arm names its missing OCCT dependency; the
/// rcad-only variants hit the catch-all).  `pub(crate)`: Geom_Surface::VIso is
/// a public Geom-level operation (ChFi3d_ComputeArete, ChFi3d_Builder_0.cxx
/// L2044).
pub(crate) fn surface_viso(surf: &Surface3, v: f64) -> Curve3 {
    // The same delegation to the kernel ElSLib V-iso bodies as
    // [`surface_uiso`] (Geom_Plane.cxx L269-273, Geom_CylindricalSurface.cxx
    // L302-306, Geom_SphericalSurface.cxx L301-305, Geom_ConicalSurface.cxx
    // L345-349, Geom_ToroidalSurface.cxx L314-318).
    match surf {
        // Geom_Plane::VIso = ElSLib::PlaneVIso: the line through P(0, V)
        // along the X direction.
        Surface3::Plane(pl) => Curve3::Line(elslib_plane_v_iso(
            &Ax3View::from_axes(pl.origin, pl.normal, pl.u_dir),
            v,
        )),
        // Geom_CylindricalSurface::VIso = ElSLib::CylinderVIso: the circle
        // of radius R on the frame Pos.Ax2(), translated by V along the axis.
        Surface3::Cylinder(cy) => {
            let pos = match cy.y_dir {
                Some(y) => Ax3View::with_y_dir(cy.origin, cy.axis, cy.ref_dir, y),
                None => Ax3View::from_axes(cy.origin, cy.axis, cy.ref_dir),
            };
            Curve3::Circle(elslib_cylinder_v_iso(&pos, cy.radius, v))
        }
        // Geom_SphericalSurface::VIso = ElSLib::SphereVIso: the parallel at
        // latitude V — the centre is lifted by R*sin(V) and the radius is
        // R*cos(V) (with the OCCT #23170 negative-radius direction flip).
        Surface3::Sphere(sp) => Curve3::Circle(elslib_sphere_v_iso(
            &Ax3View::from_axes(sp.center, sp.axis, sp.ref_dir),
            sp.radius,
            v,
        )),
        // Geom_ConicalSurface::VIso = ElSLib::ConeVIso: the parallel circle
        // Radius + V*sin(SAngle) at height V*cos(SAngle) above the
        // reference circle.
        Surface3::Cone(co) => Curve3::Circle(elslib_cone_v_iso(
            &Ax3View::from_axes(co.apex, co.axis, co.ref_dir),
            co.radius,
            co.half_angle_rad,
            v,
        )),
        // Geom_ToroidalSurface::VIso = ElSLib::TorusVIso.
        Surface3::Torus(to) => Curve3::Circle(elslib_torus_v_iso(
            &Ax3View::from_axes(to.center, to.axis, to.ref_dir),
            to.major_radius,
            to.minor_radius,
            v,
        )),
        // Geom_BSplineSurface::VIso.
        Surface3::BSpline(bs) => Curve3::BSpline(bspline_surface_viso_full(bs, v)),
        // OCCT Geom_BezierSurface::VIso (Geom_BezierSurface.cxx L1821-1863).
        Surface3::Bezier(bz) => bezier_surface_viso(bz, v),
        // OCCT Geom_SurfaceOfLinearExtrusion::VIso
        // (Geom_SurfaceOfLinearExtrusion.cxx L285-291): the basis curve
        // translated by V*direction (Vdir.Multiply(V); basis->Translated).
        Surface3::LinearExtrusion(le) => {
            curve3_translated(&le.profile, le.direction * v)
        }
        // OCCT Geom_SurfaceOfRevolution::VIso (Geom_SurfaceOfRevolution.cxx
        // L383-411): the parallel circle (Loc, axis, Rad) built at the basis
        // point of parameter V.
        Surface3::Revolution(rev) => {
            // Pnt Pc = basisCurve->Value(V);
            let pc = rev.profile.point_at(v);
            // gp_Lin L1(loc, direction); Rad = L1.Distance(Pc).
            let axis = rev.axis_dir;
            let rad = (pc - rev.axis_origin - (pc - rev.axis_origin).dot(axis) * axis).length();
            // Ax2 Rep — the circle frame: center C on the axis, XDir D the
            // radial direction (P = Pc - C normalized when Rad > Resolution).
            let c = rev.axis_origin + (pc - rev.axis_origin).dot(axis) * axis;
            let x_dir = if rad > GP_RESOLUTION {
                let d = pc - c;
                if d.length() > GP_RESOLUTION {
                    d.normalize_or_zero()
                } else {
                    // gp_Ax2(C, direction) default X — any orthogonal
                    // direction (the OCCT default picks a solver-dependent
                    // one; the radial projection is degenerate here).
                    DVec3::Z.cross(axis).normalize_or_zero()
                }
            } else {
                DVec3::Z.cross(axis).normalize_or_zero()
            };
            Curve3::Circle(Circle3 {
                center: c,
                normal: axis,
                x_dir,
                y_dir: axis.cross(x_dir).normalize_or_zero(),
                radius: rad,
            })
        }
        // OCCT Geom_RectangularTrimmedSurface::VIso
        // (Geom_RectangularTrimmedSurface.cxx L463-477): the basis VIso,
        // restricted to the u-trim range when isutrimmed.
        Surface3::Trimmed(ts) => {
            let c = surface_viso(&ts.basis, v);
            if basis_v_bounds_of(&ts.basis) != (ts.trim[0], ts.trim[1]) {
                Curve3::Trimmed(TrimmedCurve3 {
                    curve: Box::new(c),
                    first: ts.trim[0],
                    last: ts.trim[1],
                })
            } else {
                c
            }
        }
        // OCCT Geom_OffsetSurface::VIso (Geom_OffsetSurface.cxx L657-688).
        Surface3::Offset(of) => {
            // L659-660: `directRepSurface(*this)` — see the UIso arm for the
            // rcad recomputation of the eval representation.
            let an_equiv_surface =
                offset_equivalent_surface(of.basis.as_ref(), of.offset_distance);
            match an_equiv_surface {
                // L687: `return anEquivSurface->VIso(VV);`.
                Some(an_equiv) => surface_viso(&an_equiv, v),
                None => {
                    // There is no extrusion special case on this arm.
                    let (basis, offset_value) =
                        offset_basis_and_value(of.basis.as_ref(), of.offset_distance);
                    // L662-672: the general approximation arm.
                    let num1 = 0;
                    let num2 = 0;
                    let num3 = 1;
                    let t3 = [rcad_kernel::core::precision::APPROXIMATION];
                    let bounds = rcad_kernel::geom::SurfaceEval::default_domain(&basis);
                    let (u1, u2, _v1, _v2) = (bounds[0], bounds[1], bounds[2], bounds[3]);
                    let cont = GeomAbsShape::C1;
                    let max_seg = 100;
                    let max_deg = 14;

                    // L670: `Geom_OffsetSurface_VIsoEvaluator ev(*this, VV);`.
                    let mut ev = GeomOffsetSurfaceVIsoEvaluator {
                        basis,
                        offset: offset_value,
                        osc: offset_surface_osculating(of.basis.as_ref()),
                        my_iso_par: v,
                    };
                    // L671-672: the AdvApprox_ApproxAFunction ctor + Perform
                    // over [U1, U2].
                    let approx = ApproxAFunction::new(
                        num1,
                        num2,
                        num3,
                        None,
                        None,
                        Some(&t3),
                        u1,
                        u2,
                        cont,
                        max_deg,
                        max_seg,
                        &mut ev,
                    );
                    // L674: Standard_ConstructionError_Raise_if(!Approx
                    // .IsDone(), " Geom_OffsetSurface : VIso").
                    assert!(
                        approx.is_done(),
                        "Standard_ConstructionError: Geom_OffsetSurface : VIso"
                    );
                    // L676-685: the Geom_BSplineCurve of the approximation.
                    approx_to_bspline_curve(&approx)
                }
            }
        }
        // Same rcad-only variants as the UIso catch-all (no OCCT override).
        _ => panic!(
            "GAP: Geom_*Surface::VIso (TKMath/TKG3d kernel re-host) is not \
             translated for this surface type — BRepFill_Sweep::BuildWire"
        ),
    }
}

/// OCCT Geom_OffsetSurface_UIsoEvaluator (Geom_OffsetSurface.cxx L505-550) —
/// `mySurface` is the Geom_OffsetSurface itself; the rcad re-host carries the
/// value payload of that surface (its basis, the accumulated offset value and
/// the `myOscSurf` member rebuilt by [`offset_surface_osculating`]), and the
/// `Value` / `D1` calls route to the kernel `Geom_OffsetSurface::EvalD0` /
/// `EvalD1` re-hosts.
struct GeomOffsetSurfaceUIsoEvaluator {
    basis: Surface3,
    offset: f64,
    osc: Option<OsculatingSurface>,
    my_iso_par: f64,
}

impl EvaluatorFunction for GeomOffsetSurfaceUIsoEvaluator {
    /// OCCT Geom_OffsetSurface_UIsoEvaluator::Evaluate (cxx L526-550).
    fn evaluate(
        &mut self,
        _start_end: &[f64; 2],
        parameter: f64,
        derivative_request: i32,
        result: &mut [f64],
    ) -> i32 {
        if derivative_request == 0 {
            // OCCT: `P = mySurface.Value(myIsoPar, *Parameter);`.
            let p = offset_surface_eval_d0(
                &self.basis,
                self.offset,
                self.osc.as_ref(),
                self.my_iso_par,
                parameter,
            )
            .expect("Geom_UndefinedValue: Geom_OffsetSurface::EvalD0");
            result[0] = p.x;
            result[1] = p.y;
            result[2] = p.z;
        } else {
            // OCCT: `mySurface.D1(myIsoPar, *Parameter, P, DU, DV);` — the
            // derivative request answers DV.
            let d1 = offset_surface_eval_d1(
                &self.basis,
                self.offset,
                self.osc.as_ref(),
                self.my_iso_par,
                parameter,
            )
            .expect("Geom_UndefinedDerivative: Geom_OffsetSurface::EvalD1");
            result[0] = d1.d1v.x;
            result[1] = d1.d1v.y;
            result[2] = d1.d1v.z;
        }
        0
    }
}

/// OCCT Geom_OffsetSurface_VIsoEvaluator (Geom_OffsetSurface.cxx L552-597).
struct GeomOffsetSurfaceVIsoEvaluator {
    basis: Surface3,
    offset: f64,
    osc: Option<OsculatingSurface>,
    my_iso_par: f64,
}

impl EvaluatorFunction for GeomOffsetSurfaceVIsoEvaluator {
    /// OCCT Geom_OffsetSurface_VIsoEvaluator::Evaluate (cxx L573-597).
    fn evaluate(
        &mut self,
        _start_end: &[f64; 2],
        parameter: f64,
        derivative_request: i32,
        result: &mut [f64],
    ) -> i32 {
        if derivative_request == 0 {
            // OCCT: `P = mySurface.Value(*Parameter, myIsoPar);`.
            let p = offset_surface_eval_d0(
                &self.basis,
                self.offset,
                self.osc.as_ref(),
                parameter,
                self.my_iso_par,
            )
            .expect("Geom_UndefinedValue: Geom_OffsetSurface::EvalD0");
            result[0] = p.x;
            result[1] = p.y;
            result[2] = p.z;
        } else {
            // OCCT: `mySurface.D1(*Parameter, myIsoPar, P, DU, DV);` — the
            // derivative request answers DU.
            let d1 = offset_surface_eval_d1(
                &self.basis,
                self.offset,
                self.osc.as_ref(),
                parameter,
                self.my_iso_par,
            )
            .expect("Geom_UndefinedDerivative: Geom_OffsetSurface::EvalD1");
            result[0] = d1.d1u.x;
            result[1] = d1.d1u.y;
            result[2] = d1.d1u.z;
        }
        0
    }
}

/// OCCT Geom_OffsetSurface::UIso / VIso tail (cxx L639-650 / L676-685) — the
/// `AdvApprox_ApproxAFunction` result converted to the `Geom_BSplineCurve`
/// `new Geom_BSplineCurve(Poles, Knots, Mults, Approx.Degree())`.
fn approx_to_bspline_curve(approx: &ApproxAFunction) -> Curve3 {
    let degree = approx.degree().max(0) as usize;
    let knots = approx.knots_vec().to_vec();
    let mults = approx.multiplicities_vec().to_vec();
    // OCCT: `Approx.Poles(1, Poles)` — the 3D subspace #1.
    let poles = approx.poles_flat(1);
    let control_points: Vec<DVec3> = (0..poles.len() / 3)
        .map(|i| DVec3::new(poles[3 * i], poles[3 * i + 1], poles[3 * i + 2]))
        .collect();
    Curve3::BSpline(BSplineCurve3::from_knots_mults(
        degree,
        knots,
        mults,
        control_points,
    ))
}

/// OCCT gp_GTrsf-based `Geom_Curve::Translated(T)` over the rcad Curve3 —
/// supported for the analytic / poles-based variants (the OCCT operation
/// translates the geometry representation in place); the remaining curve
/// kinds keep the OCCT failure path.
fn curve3_translated(c: &Curve3, t: DVec3) -> Curve3 {
    match c {
        Curve3::Line(l) => Curve3::Line(Line3 {
            origin: l.origin + t,
            direction: l.direction,
        }),
        Curve3::Circle(ci) => Curve3::Circle(Circle3 {
            center: ci.center + t,
            ..ci.clone()
        }),
        Curve3::Ellipse(el) => Curve3::Ellipse(Ellipse3 {
            center: el.center + t,
            ..el.clone()
        }),
        Curve3::BSpline(bs) => Curve3::BSpline(rcad_kernel::geom::BSplineCurve3 {
            control_points: bs.control_points.iter().map(|p| p + t).collect(),
            ..bs.clone()
        }),
        Curve3::Bezier(bz) => Curve3::Bezier(rcad_kernel::geom::BezierCurve3 {
            control_points: bz.control_points.iter().map(|p| p + t).collect(),
            weights: bz.weights.clone(),
        }),
        Curve3::Trimmed(tr) => Curve3::Trimmed(TrimmedCurve3 {
            curve: Box::new(curve3_translated(&tr.curve, t)),
            first: tr.first,
            last: tr.last,
        }),
        _ => panic!(
            "GAP: Geom_Curve::Translated (TKG3d/Geom kernel re-host) for this \
             curve type — Geom_SurfaceOfLinearExtrusion::VIso"
        ),
    }
}

/// OCCT `Geom_Curve::Rotate(Ax1, Angle)` over the rcad Curve3 — the
/// rotation about the (origin, direction) axis applied to the geometry
/// representation; supported for the analytic / poles-based variants.
fn curve3_rotated_about_axis(c: &Curve3, origin: DVec3, direction: DVec3, angle: f64) -> Curve3 {
    use glam::DAffine3;
    let axis = DAffine3::from_rotation_translation(
        glam::DQuat::from_axis_angle(direction, angle),
        DVec3::ZERO,
    );
    let rot = |p: DVec3| origin + axis.transform_point3(p - origin);    match c {
        Curve3::Line(l) => Curve3::Line(Line3 {
            origin: rot(l.origin),
            direction: axis.transform_vector3(l.direction),
        }),
        Curve3::Circle(ci) => Curve3::Circle(Circle3 {
            center: rot(ci.center),
            normal: axis.transform_vector3(ci.normal),
            x_dir: axis.transform_vector3(ci.x_dir),
            y_dir: axis.transform_vector3(ci.y_dir),
            radius: ci.radius,
        }),
        Curve3::Ellipse(el) => Curve3::Ellipse(Ellipse3 {
            center: rot(el.center),
            normal: axis.transform_vector3(el.normal),
            major_dir: axis.transform_vector3(el.major_dir),
            major_radius: el.major_radius,
            minor_radius: el.minor_radius,
        }),
        Curve3::BSpline(bs) => Curve3::BSpline(rcad_kernel::geom::BSplineCurve3 {
            control_points: bs.control_points.iter().map(|p| rot(*p)).collect(),
            ..bs.clone()
        }),
        Curve3::Bezier(bz) => Curve3::Bezier(rcad_kernel::geom::BezierCurve3 {
            control_points: bz.control_points.iter().map(|p| rot(*p)).collect(),
            weights: bz.weights.clone(),
        }),
        Curve3::Trimmed(tr) => Curve3::Trimmed(TrimmedCurve3 {
            curve: Box::new(curve3_rotated_about_axis(&tr.curve, origin, direction, angle)),
            first: tr.first,
            last: tr.last,
        }),
        _ => panic!(
            "GAP: Geom_Curve::Rotate (TKG3d/Geom kernel re-host) for this \
             curve type — Geom_SurfaceOfRevolution::UIso"
        ),
    }
}

/// OCCT Geom_BSplineSurface::UIso (Geom_BSplineSurface_1.cxx L598-630) — the
/// v-varying iso curve at `u`: `BSplSLib::Iso` over the U direction
/// (BSplSLib.cxx L1617-1740, the rcad [`bspl_slib_iso`] re-host) wrapped in a
/// `Geom_BSplineCurve` carrying the opposite (V) knot vector.
fn bspline_surface_uiso(surf: &BSplineSurface, u: f64) -> rcad_kernel::geom::BSplineCurve3 {
    // OCCT L602: `if (myURational || myVRational)` selects `Weights()`,
    // otherwise `BSplSLib::NoWeights()`.  The rcad surface derives the two
    // constructor flags on read (`BSplineSurface::is_rational_u/v`).
    let weights = if surf.is_rational_u() || surf.is_rational_v() {
        Some(surf.weights.as_slice())
    } else {
        None
    };
    // OCCT L605-609: BSplSLib::Iso(U, true, myPoles, Weights(), myUFlatKnots,
    // BSplCLib::NoMults(), myUDeg, myUPeriodic, cpoles, &cweights).  The rcad
    // surface stores the multiplicity-expanded knot vector, which is the
    // `Mults == NoMults()` form.
    let (cpoles, cweights) = bspl_slib_iso(
        u,
        true,
        surf.degree_u,
        &surf.knots_u,
        &surf.control_points,
        weights,
        surf.is_periodic_u,
    );
    // OCCT L610: `new Geom_BSplineCurve(cpoles, cweights, myVKnots, myVMults,
    // myVDeg, myVPeriodic)` — the iso curve runs along V, so it carries the V
    // knot vector and the V periodic flag.  The non-rational arm (L617)
    // omits the weights; the rcad BSplineCurve3 always carries the weights
    // vector, and bspl_slib_iso returns 1.0 for every pole in that arm.
    BSplineCurve3 {
        degree: surf.degree_v,
        knots: surf.knots_v.clone(),
        control_points: cpoles,
        weights: cweights,
        is_periodic: surf.is_periodic_v,
    }
}

/// OCCT Geom_BSplineSurface::VIso (Geom_BSplineSurface_1.cxx L775-807) — the
/// u-varying iso curve at `v`: `BSplSLib::Iso` over the V direction
/// (BSplSLib.cxx L1617-1740, the rcad [`bspl_slib_iso`] re-host) wrapped in a
/// `Geom_BSplineCurve` carrying the opposite (U) knot vector.  The OCCT
/// consumers call this through `Geom_Surface::VIso` (BRepFill_NSections.cxx
/// L679/L687/L780/L872, GeomFill_NSections.cxx L300).
fn bspline_surface_viso_full(surf: &BSplineSurface, v: f64) -> rcad_kernel::geom::BSplineCurve3 {
    // OCCT L779: `if (myURational || myVRational)` selects `Weights()`,
    // otherwise `BSplSLib::NoWeights()`.
    let weights = if surf.is_rational_u() || surf.is_rational_v() {
        Some(surf.weights.as_slice())
    } else {
        None
    };
    // OCCT L782-786: BSplSLib::Iso(V, false, myPoles, Weights(), myVFlatKnots,
    // BSplCLib::NoMults(), myVDeg, myVPeriodic, cpoles, &cweights).
    let (cpoles, cweights) = bspl_slib_iso(
        v,
        false,
        surf.degree_v,
        &surf.knots_v,
        &surf.control_points,
        weights,
        surf.is_periodic_v,
    );
    // OCCT L787: `new Geom_BSplineCurve(cpoles, cweights, myUKnots, myUMults,
    // myUDeg, myUPeriodic)` — the iso curve runs along U, so it carries the U
    // knot vector and the U periodic flag.
    BSplineCurve3 {
        degree: surf.degree_u,
        knots: surf.knots_u.clone(),
        control_points: cpoles,
        weights: cweights,
        is_periodic: surf.is_periodic_u,
    }
}

/// OCCT static Rational(Weights, Urational, Vrational)
/// (Geom_BezierSurface.cxx L66-96) — the weight-variation flags of a surface
/// weight grid.  OCCT derives them once in the constructor; the rcad
/// `BezierSurface` is an immutable value, so the same derivation evaluates on
/// read.  The names are OCCT's verbatim: `Vrational` is set when two weights
/// adjacent along the first (U) pole index differ, `Urational` when two
/// weights adjacent along the second (V) pole index differ.
fn bezier_surface_rational(weights: &[Vec<f64>]) -> (bool, bool) {
    let nb_cols = weights.first().map(|row| row.len()).unwrap_or(0);
    // OCCT L68-80: Vrational over the column loop J, comparing
    // Weights(I, J) with Weights(I + 1, J).
    let mut v_rational = false;
    'v_cols: for jj in 0..nb_cols {
        for ii in 0..weights.len().saturating_sub(1) {
            let w = weights[ii][jj];
            let w_next = weights[ii + 1][jj];
            if (w - w_next).abs() > epsilon_of(w.abs()) {
                v_rational = true;
                break 'v_cols;
            }
        }
    }
    // OCCT L82-96: Urational over the row loop I, comparing Weights(I, J)
    // with Weights(I, J + 1).
    let mut u_rational = false;
    'u_rows: for row in weights.iter() {
        for jj in 0..row.len().saturating_sub(1) {
            let w = row[jj];
            let w_next = row[jj + 1];
            if (w - w_next).abs() > epsilon_of(w.abs()) {
                u_rational = true;
                break 'u_rows;
            }
        }
    }
    (u_rational, v_rational)
}

/// OCCT `Geom_BezierSurface::UKnotSequence()` / `VKnotSequence()`
/// (Geom_BezierSurface.cxx L2167-2193) — the multiplicity-expanded form of
/// the implicit Bezier knot vector, `BSplCLib::FlatBezierKnots(degree)` =
/// `[0; degree + 1] ++ [1; degree + 1]`.
fn bezier_flat_knots(degree: usize) -> Vec<f64> {
    std::iter::repeat_n(0.0, degree + 1)
        .chain(std::iter::repeat_n(1.0, degree + 1))
        .collect()
}

/// OCCT Geom_BezierSurface::UIso (Geom_BezierSurface.cxx L1769-1810) — the
/// iso curve at u: BSplSLib::Iso over the U direction, wrapped in a
/// `Geom_BezierCurve` of `myPoles.RowLength()` poles.
///
/// Architecture difference: the OCCT call passes the implicit Bezier knot
/// vector as the pair `UKnots()` = [0, 1] and `&UMultiplicities()` =
/// [nbU, nbU]; the rcad [`bspl_slib_iso`] takes the flat form (`Mults ==
/// BSplCLib::NoMults()`), so the equivalent `UKnotSequence()` is passed with
/// the same degree `(myPoles.ColLength() - 1)` and `Periodic = false` (the
/// same note as the Geom_BSplineSurface::UIso arm, which passes
/// `myUFlatKnots` + `myUDeg`).
fn bezier_surface_uiso(surf: &BezierSurface, u: f64) -> Curve3 {
    // OCCT L1771-1773: VCurvePoles(1, myPoles.RowLength()) and the degree
    // (myPoles.ColLength() - 1).
    let degree = surf.control_points.len() - 1;
    let (my_u_rational, my_v_rational) = bezier_surface_rational(&surf.weights);
    if my_u_rational || my_v_rational {
        // OCCT L1775-1785: BSplSLib::Iso(U, true, myPoles, &myWeights,
        // UKnots(), &UMultiplicities(), degree, false, VCurvePoles,
        // &VCurveWeights).
        let (v_curve_poles, v_curve_weights) = bspl_slib_iso(
            u,
            true,
            degree,
            &bezier_flat_knots(degree),
            &surf.control_points,
            Some(&surf.weights),
            false,
        );
        if my_u_rational {
            // OCCT L1786: new Geom_BezierCurve(VCurvePoles, VCurveWeights).
            Curve3::Bezier(BezierCurve3 {
                control_points: v_curve_poles,
                weights: v_curve_weights,
            })
        } else {
            // OCCT L1790: new Geom_BezierCurve(VCurvePoles) — non-rational.
            let nb_poles = v_curve_poles.len();
            Curve3::Bezier(BezierCurve3 {
                control_points: v_curve_poles,
                weights: vec![1.0; nb_poles],
            })
        }
    } else {
        // OCCT L1795-1803: BSplSLib::Iso with BSplSLib::NoWeights() and
        // PLib::NoWeights(); new Geom_BezierCurve(VCurvePoles).
        let (v_curve_poles, _) = bspl_slib_iso(
            u,
            true,
            degree,
            &bezier_flat_knots(degree),
            &surf.control_points,
            None,
            false,
        );
        let nb_poles = v_curve_poles.len();
        Curve3::Bezier(BezierCurve3 {
            control_points: v_curve_poles,
            weights: vec![1.0; nb_poles],
        })
    }
}

/// OCCT Geom_BezierSurface::VIso (Geom_BezierSurface.cxx L1821-1863) — the
/// v-isoparametric counterpart of [`bezier_surface_uiso`]: BSplSLib::Iso over
/// the V direction with the degree `(myPoles.RowLength() - 1)` and the
/// `&VMultiplicities()` implicit knot vector, wrapped in a `Geom_BezierCurve`
/// of `myPoles.ColLength()` poles.
fn bezier_surface_viso(surf: &BezierSurface, v: f64) -> Curve3 {
    // OCCT L1823-1824: VCurvePoles(1, myPoles.ColLength()) and the degree
    // (myPoles.RowLength() - 1).
    let nb_u_poles = surf.control_points.first().map(|row| row.len()).unwrap_or(0);
    let degree = nb_u_poles - 1;
    let (my_u_rational, my_v_rational) = bezier_surface_rational(&surf.weights);
    if my_v_rational || my_u_rational {
        // OCCT L1825-1837: BSplSLib::Iso(V, false, myPoles, &myWeights,
        // UKnots(), &VMultiplicities(), degree, false, VCurvePoles,
        // &VCurveWeights).
        let (v_curve_poles, v_curve_weights) = bspl_slib_iso(
            v,
            false,
            degree,
            &bezier_flat_knots(degree),
            &surf.control_points,
            Some(&surf.weights),
            false,
        );
        if my_v_rational {
            // OCCT L1838: new Geom_BezierCurve(VCurvePoles, VCurveWeights).
            Curve3::Bezier(BezierCurve3 {
                control_points: v_curve_poles,
                weights: v_curve_weights,
            })
        } else {
            // OCCT L1842: new Geom_BezierCurve(VCurvePoles) — non-rational.
            let nb_poles = v_curve_poles.len();
            Curve3::Bezier(BezierCurve3 {
                control_points: v_curve_poles,
                weights: vec![1.0; nb_poles],
            })
        }
    } else {
        // OCCT L1847-1855: BSplSLib::Iso with BSplSLib::NoWeights() and
        // PLib::NoWeights(); new Geom_BezierCurve(VCurvePoles).
        let (v_curve_poles, _) = bspl_slib_iso(
            v,
            false,
            degree,
            &bezier_flat_knots(degree),
            &surf.control_points,
            None,
            false,
        );
        let nb_poles = v_curve_poles.len();
        Curve3::Bezier(BezierCurve3 {
            control_points: v_curve_poles,
            weights: vec![1.0; nb_poles],
        })
    }
}

// ---------------------------------------------------------------------------
// File statics (part A)
// ---------------------------------------------------------------------------

/// OCCT static Translate (cxx L178-189) — Copy a column from one table to
/// another (1-based OCCT indexing -> `[i-1]`).
pub(super) fn translate(
    array_in: &ShapeHArray2,
    in_col: i32,
    array_out: &mut ShapeHArray2,
    out_col: i32,
) {
    let nb = array_out[0].len();
    for ii in 1..=nb {
        array_out[ii - 1][out_col as usize - 1] = array_in[ii - 1][in_col as usize - 1].clone();
    }
}

/// OCCT static Box (cxx L193-206) — the section poles bounding box (the
/// Bnd_Box SetVoid / Add / Get triple reduced to the accumulated min/max,
/// the pure-math Bnd_Box re-host).
pub(super) fn box_of_section(sec: &Rc<RefCell<dyn SectionLaw>>, u: f64) -> ([f64; 3], [f64; 3]) {
    let mut nb_poles = 0usize;
    let mut bid = 0usize;
    {
        let mut bid2 = 0usize;
        sec.borrow().section_shape(&mut nb_poles, &mut bid, &mut bid2);
    }
    let mut poles = vec![DVec3::ZERO; nb_poles];
    let mut w = vec![0.0f64; nb_poles];
    sec.borrow().d0(u, &mut poles, &mut w);
    let mut min = [f64::MAX; 3];
    let mut max = [f64::MIN; 3];
    for p in &poles {
        for k in 0..3 {
            min[k] = min[k].min(p[k]);
            max[k] = max[k].max(p[k]);
        }
    }
    (min, max)
}

/// OCCT static BuildVertex (cxx L1440-1458).
pub(super) fn build_vertex(
    brep: &mut BRep,
    iso: &Curve3,
    isfirst: bool,
    first: f64,
    last: f64,
    vertex: &mut Shape,
) {
    let mut b = BRepBuilder::new();
    let val = if isfirst { first } else { last };
    // B.MakeVertex(TopoDS::Vertex(Vertex), Iso->Value(val), Confusion());
    *vertex = b.add_vertex(brep, iso.point_at(val), rcad_kernel::core::precision::CONFUSION);
}

/// OCCT static NullEdge (cxx L1462-1472).
pub(super) fn null_edge(brep: &mut BRep, vertex: &mut Shape) -> Shape {
    use rcad_kernel::topo::topods::Orientation;
    let mut b = BRepBuilder::new();
    let mut e = b.add_edge(brep, None, vertex.clone(), shape_reversed(vertex), [0.0, 0.0]);
    // Vertex.Orientation(TopAbs_FORWARD); B.Add(E, Vertex); B.Add(E,
    // Vertex.Reversed()); B.Degenerated(E, true);
    vertex.orientation = Orientation::Forward;
    b.set_edge_degenerated(brep, e.clone(), true);
    e
}

/// OCCT static HasPCurves (cxx L156-172) — the edge carries a
/// curve-on-surface representation.
pub(super) fn has_pcurves(brep: &BRep, e: &Shape) -> bool {
    !brep.edge(e.clone()).pcurves.is_empty()
}

/// OCCT static ReverseEdgeInFirstOrLastWire (cxx L1877-1900).
pub(super) fn reverse_edge_in_first_or_last_wire(brep: &mut BRep, the_wire: &mut Shape, the_edge: &Shape) {
    // TopoDS_Iterator itw(theWire) — the wire-ordered children.
    let mut edge_to_reverse: Option<Shape> = None;
    for an_edge in crate::brep_fill::compatible_wires::wire_edges(brep, the_wire) {
        if an_edge.is_same(the_edge) {
            edge_to_reverse = Some(an_edge);
            break;
        }
    }

    if let Some(mut edge) = edge_to_reverse {
        let mut bb = BRepBuilder::new();
        edge.orientation = match edge.orientation {
            rcad_kernel::topo::topods::Orientation::Forward => {
                rcad_kernel::topo::topods::Orientation::Reversed
            }
            rcad_kernel::topo::topods::Orientation::Reversed => {
                rcad_kernel::topo::topods::Orientation::Forward
            }
            o => o,
        };
        bb.remove_from_wire(brep, the_wire.clone(), edge.clone());
        bb.add_to_wire(brep, the_wire.clone(), edge);
    }
}

/// OCCT static ReverseModifiedEdges (cxx L1902-1936).
pub(super) fn reverse_modified_edges(brep: &mut BRep, the_wire: &mut Shape, the_emap: &ShapeSet) {
    if the_emap.is_empty() {
        return;
    }

    let mut bb = BRepBuilder::new();
    let ledges = crate::brep_fill::compatible_wires::wire_edges(brep, the_wire);
    for e in &ledges {
        bb.remove_from_wire(brep, the_wire.clone(), e.clone());
    }

    for e in ledges {
        let mut an_edge = e;
        if the_emap.contains(&ShapeKey(an_edge.ptr_id())) {
            an_edge = shape_reversed(&an_edge);
        }
        bb.add_to_wire(brep, the_wire.clone(), an_edge);
    }
}

// ---------------------------------------------------------------------------
// BRepFill_Sweep
// ---------------------------------------------------------------------------

/// OCCT BRepFill_Sweep (hxx L43-189): Topological Sweep Algorithm.
pub struct BRepFillSweep {
    /// OCCT bool isDone.
    pub is_done: bool,
    /// OCCT bool KPart.
    pub kpart: bool,
    /// OCCT double myTol3d.
    pub my_tol3d: f64,
    /// OCCT double myBoundTol.
    pub my_bound_tol: f64,
    /// OCCT double myTol2d.
    pub my_tol2d: f64,
    /// OCCT double myTolAngular.
    pub my_tol_angular: f64,
    /// OCCT double myAngMin.
    pub my_ang_min: f64,
    /// OCCT double myAngMax.
    pub my_ang_max: f64,
    /// OCCT GeomFill_ApproxStyle myApproxStyle.
    pub my_approx_style: GeomFillApproxStyle,
    /// OCCT GeomAbs_Shape myContinuity.
    pub my_continuity: GeomAbsShape,
    /// OCCT int myDegmax.
    pub my_degmax: i32,
    /// OCCT int mySegmax.
    pub my_segmax: i32,
    /// OCCT bool myForceApproxC1.
    pub my_force_approx_c1: bool,
    /// OCCT TopoDS_Shape myShape.
    pub my_shape: Shape,
    /// OCCT handle(BRepFill_LocationLaw) myLoc.
    pub my_loc: Rc<RefCell<dyn BRepFillLocationLawOps>>,
    /// OCCT handle(BRepFill_SectionLaw) mySec.
    pub my_sec: Rc<RefCell<dyn crate::brep_fill::brep_fill_section_law::BRepFillSectionLawOps>>,
    /// OCCT handle(NCollection_HArray2<TopoDS_Shape>) myUEdges.
    pub my_u_edges: Option<ShapeHArray2>,
    /// OCCT handle(NCollection_HArray2<TopoDS_Shape>) myVEdges.
    pub my_v_edges: Option<ShapeHArray2>,
    /// OCCT NCollection_DataMap<TopoDS_Shape, TopoDS_Shape> myVEdgesModified.
    pub my_v_edges_modified: std::collections::HashMap<ShapeKey, Shape>,
    /// OCCT handle(NCollection_HArray2<TopoDS_Shape>) myFaces.
    pub my_faces: Option<ShapeHArray2>,
    /// OCCT NCollection_List<TopoDS_Shape> myAuxShape.
    pub my_aux_shape: Vec<Shape>,
    /// OCCT handle(NCollection_HArray1<TopoDS_Shape>) myTapes.
    pub my_tapes: Option<Vec<Shape>>,
    /// OCCT double Error.
    pub error: f64,
    /// OCCT TopoDS_Wire FirstShape.
    pub first_shape: Shape,
    /// OCCT TopoDS_Wire LastShape.
    pub last_shape: Shape,
    /// The BuildWire v-close flag (the wire of cxx L2094-2096 — the
    /// wire under construction, carried for the pool-based builder).
    pub(crate) brep_wire: Option<Shape>,
}

impl BRepFillSweep {
    /// OCCT BRepFill_Sweep(Section, Location, WithKPart) (cxx L1940-1958).
    pub fn new(
        brep: &mut BRep,
        section: Rc<RefCell<dyn crate::brep_fill::brep_fill_section_law::BRepFillSectionLawOps>>,
        location: Rc<RefCell<dyn BRepFillLocationLawOps>>,
        with_kpart: bool,
    ) -> Self {
        let mut this = BRepFillSweep {
            is_done: false,
            kpart: with_kpart,
            my_tol3d: 0.0,
            my_bound_tol: 0.0,
            my_tol2d: 0.0,
            my_tol_angular: 0.0,
            my_ang_min: 0.0,
            my_ang_max: 0.0,
            my_approx_style: GeomFillApproxStyle::GeomFill_Location,
            my_continuity: GeomAbsShape::C2,
            my_degmax: 11,
            my_segmax: 30,
            my_force_approx_c1: false,
            my_shape: Shape::null(),
            my_loc: location,
            my_sec: section,
            my_u_edges: None,
            my_v_edges: None,
            my_v_edges_modified: std::collections::HashMap::new(),
            my_faces: None,
            my_aux_shape: Vec::new(),
            my_tapes: None,
            error: 0.0,
            first_shape: Shape::null(),
            last_shape: Shape::null(),
            brep_wire: None,
        };

        this.set_tolerance(brep, 1.0e-4, 1.0, 1.0e-5, 1.0e-2);
        this.set_angular_control(0.01, 6.0);
        // myAuxShape.Clear();
        this.my_aux_shape.clear();

        this.my_approx_style = GeomFillApproxStyle::GeomFill_Location;
        this.my_continuity = GeomAbsShape::C2;
        this.my_degmax = 11;
        this.my_segmax = 30;
        this.my_force_approx_c1 = false;
        this
    }

    /// OCCT SetBounds (cxx L1962-2007) — It is necessary to check the
    /// SameRange on its (PRO13551).  `BRepLib::CheckSameRange(E)` maps to
    /// the stored `TEdgeData::same_range` flag.
    pub fn set_bounds(&mut self, brep: &mut BRep, first: &Shape, last: &Shape) {
        self.first_shape = first.clone();
        self.last_shape = last.clone();

        let mut b = BRepBuilder::new();
        if !self.first_shape.is_null() {
            for e in crate::brep_fill::compatible_wires::wire_edges(brep, &self.first_shape) {
                let same_range = brep.edge(e.clone()).same_range;
                if !same_range {
                    b.set_edge_same_range(brep, e.clone(), false);
                    b.set_edge_same_parameter(brep, e.clone(), false);
                }
            }
        }

        if !self.last_shape.is_null() {
            for e in crate::brep_fill::compatible_wires::wire_edges(brep, &self.last_shape) {
                let same_range = brep.edge(e.clone()).same_range;
                if !same_range {
                    b.set_edge_same_range(brep, e.clone(), false);
                    b.set_edge_same_parameter(brep, e.clone(), false);
                }
            }
        }
    }

    /// OCCT SetTolerance (cxx L2011-2020).
    pub fn set_tolerance(
        &mut self,
        _brep: &mut BRep,
        tol3d: f64,
        bound_tol: f64,
        tol2d: f64,
        tol_angular: f64,
    ) {
        self.my_tol3d = tol3d;
        self.my_bound_tol = bound_tol;
        self.my_tol2d = tol2d;
        self.my_tol_angular = tol_angular;
    }

    /// OCCT SetAngularControl (cxx L2024-2028).
    pub fn set_angular_control(&mut self, min_angle: f64, max_angle: f64) {
        self.my_ang_min = min_angle.max(rcad_kernel::core::precision::ANGULAR);
        self.my_ang_max = max_angle.min(6.28);
    }

    /// OCCT SetForceApproxC1 (cxx L2036-2039).
    pub fn set_force_approx_c1(&mut self, force_approx_c1: bool) {
        self.my_force_approx_c1 = force_approx_c1;
    }

    /// OCCT BuildWire (cxx L2078-2209) — the vertex-section sweep: one
    /// GeomFill_Sweep per path edge, the iso edge built over each surface.
    pub fn build_wire(&mut self, brep: &mut BRep, _transition: BRepFillTransitionStyle) -> bool {
        let mut isec = 1i32;
        let mut p1;

        let mut b = BRepBuilder::new();
        let nb_path = self.my_loc.borrow().base().nb_law();
        // vclose = (myLoc->IsClosed() && (myLoc->IsG1(0, myTol3d) >= 0));
        let vclose = self.my_loc.borrow().base().is_closed(brep)
            && self.my_loc.borrow().base().is_g1(brep, 0, self.my_tol3d, 1.0e-4) >= 0;
        self.error = 0.0;

        // TopoDS_Wire wire; B.MakeWire(wire);
        let wire = b.make_wire(brep);

        // (1.1) Construction of Tables
        // myFaces  = new HArray2(1, 1, 1, NbPath);
        // myUEdges = new HArray2(1, 2, 1, NbPath);
        // myVEdges = new HArray2(1, 1, 1, NbPath + 1);
        self.my_faces = Some(vec![vec![Shape::null(); nb_path as usize]]);
        self.my_u_edges = Some(vec![vec![Shape::null(); nb_path as usize]; 2]);
        self.my_v_edges = Some(vec![vec![Shape::null(); nb_path as usize + 1]]);

        // (1.2) Calculate curves / vertex / edge
        for ipath in 1..=nb_path {
            // GeomFill_Sweep Sweep(myLoc->Law(ipath), KPart);
            let loc_law = self.my_loc.borrow().base().law(ipath);
            let mut sweep =
                crate::geomalgo::geomfill::sweep::Sweep::new(loc_law, self.kpart);
            sweep.set_tolerance_sweep(
                self.my_tol3d,
                self.my_bound_tol,
                self.my_tol2d,
                self.my_tol_angular,
            );
            sweep.set_force_approx_c1(self.my_force_approx_c1);
            let sec_law = self.my_sec.borrow().base().law(isec);
            sweep.build(
                sec_law,
                self.my_approx_style,
                self.my_continuity,
                self.my_degmax,
                self.my_segmax,
            );
            if !sweep.is_done() {
                return false;
            }
            let surf = sweep.surface().expect("null surface").clone();
            // S->Bounds(...) and the iso parameter (cxx L2119-2141).
            let (_basis, bounds) = surface_adaptor_basis_and_bounds(&surf);
            let (umin, umax, vmin, vmax) = (bounds[0], bounds[1], bounds[2], bounds[3]);
            let (mut first, mut last, val);
            let iso: Curve3;
            if sweep.exchange_uv() {
                if sweep.u_reversed() {
                    first = umin;
                    last = umax;
                    val = vmax;
                } else {
                    first = umin;
                    last = umax;
                    val = vmin;
                }
                iso = surface_viso(&surf, val);
            } else {
                if sweep.u_reversed() {
                    val = umax;
                    first = vmin;
                    last = vmax;
                } else {
                    val = umin;
                    first = vmin;
                    last = vmax;
                }
                iso = surface_uiso(&surf, val);
            }

            // Vertex by position (cxx L2143-2189)
            let my_v_edges = self.my_v_edges.as_mut().expect("myVEdges");
            if ipath < nb_path {
                let mut v = my_v_edges[0][ipath as usize].clone();
                build_vertex(brep, &iso, false, first, last, &mut v);
                my_v_edges[0][ipath as usize] = v;
            } else if vclose {
                let mut v = my_v_edges[0][0].clone();
                my_v_edges[0][ipath as usize] = v.clone();
                // Iso->D0(Last, P1); Tol = P1.Distance(BRep_Tool::Pnt(V));
                p1 = iso.point_at(last);
                let tol = p1.distance(brep.vertex(v.clone()).point);
                b.update_vertex_tolerance(brep, v, tol);
            } else if !self.first_shape.is_null() {
                // myVEdges->SetValue(1, NbPath, FirstShape);
                my_v_edges[0][nb_path as usize - 1] = self.first_shape.clone();
            } else {
                let mut v = my_v_edges[0][nb_path as usize].clone();
                build_vertex(brep, &iso, false, first, last, &mut v);
                my_v_edges[0][nb_path as usize] = v;
            }

            if ipath > 1 {
                p1 = iso.point_at(first);
                let v = my_v_edges[0][(ipath - 1) as usize].clone();
                let tol = p1.distance(brep.vertex(v.clone()).point);
                b.update_vertex_tolerance(brep, v, tol);
            }
            if ipath == 1 {
                if !self.first_shape.is_null() {
                    my_v_edges[0][0] = self.first_shape.clone();
                } else {
                    let mut v = my_v_edges[0][0].clone();
                    build_vertex(brep, &iso, true, first, last, &mut v);
                    my_v_edges[0][0] = v;
                }
            }

            // Construction of the edge (cxx L2192-2205)
            // BRepLib_MakeEdge MkE; MkE.Init(Iso, V(1, ipath), V(1, ipath+1),
            //                                Iso->FirstParameter(),
            //                                Iso->LastParameter());
            let v_edges = self.my_v_edges.as_ref().expect("myVEdges");
            let v1 = v_edges[0][(ipath - 1) as usize].clone();
            let v2 = v_edges[0][ipath as usize].clone();
            let e = b.add_edge(
                brep,
                Some(iso.clone()),
                v1,
                v2,
                [
                    curve_first_parameter_of(&iso),
                    curve_last_parameter_of(&iso),
                ],
            );
            // if (!MkE.IsDone()) return false; — the rcad add_edge is
            // infallible at this layer (the OCCT error arm is empty too).
            b.update_edge_tolerance(brep, e.clone(), sweep.error_on_surface());
            b.add_to_wire(brep, wire.clone(), e.clone());
            self.my_faces.as_mut().expect("myFaces")[0][(ipath - 1) as usize] = e;
        }
        self.my_shape = wire;
        true
    }

    /// OCCT MergeVertex (cxx L3984-4004) — Make V2 = V1 if V2 is too close
    /// to V1.
    pub fn merge_vertex(&self, brep: &BRep, v1: &Shape, v2: &mut Shape) -> bool {
        let mut tol = brep.vertex(v1.clone()).tolerance.max(brep.vertex(v2.clone()).tolerance);
        if tol < self.my_tol3d {
            tol = self.my_tol3d;
        }
        let p1 = brep.vertex(v1.clone()).point;
        let p2 = brep.vertex(v2.clone()).point;
        if p1.distance(p2) <= tol {
            *v2 = v1.clone();
            return true;
        }
        false
    }

    /// OCCT UpdateVertex (cxx L4010-4040) — Update the Tolerance of
    /// Vertices depending on Laws.
    pub fn update_vertex(
        &self,
        brep: &mut BRep,
        ipath: i32,
        isec: i32,
        err_app: f64,
        param: f64,
        v: &mut Shape,
    ) {
        let mut b = BRepBuilder::new();
        let the_v = v.clone();
        let mut vv = Shape::null();
        let (sec_vertex, sec_tol) = {
            let sec = self.my_sec.borrow();
            let v = sec.vertex(brep, isec, param);
            let t = sec.vertex_tol(brep, isec - 1, param);
            (v, t)
        };
        self.my_loc
            .borrow()
            .base()
            .perform_vertex(brep, ipath, &sec_vertex, err_app + sec_tol, &mut vv, 0);
        let p1 = brep.vertex(vv.clone()).point;
        let p2 = brep.vertex(the_v.clone()).point;

        let mut tol = brep.vertex(vv.clone()).tolerance;
        tol += p1.distance(p2);

        if tol > brep.vertex(the_v.clone()).tolerance {
            b.update_vertex_tolerance(brep, the_v, tol);
        }
    }

    /// OCCT EvalExtrapol (cxx L3910-3978).  (the &mut receiver carries the
    /// CurvilinearBounds myLength mutation through the OCCT array handle —
    /// the const-mutation pattern documented in brep_fill_location_law.rs)
    pub fn eval_extrapol(
        &mut self,
        brep: &mut BRep,
        index: i32,
        transition: BRepFillTransitionStyle,
    ) -> f64 {
        let mut extrap = 0.0f64;
        if transition == BRepFillTransitionStyle::RightCorner {
            let i1;
            let i2;
            if (index == 1) || (index == self.my_loc.borrow().base().nb_law() + 1) {
                if !self.my_loc.borrow().base().is_closed(brep)
                    || !self.my_sec.borrow().base().is_vclosed()
                {
                    return extrap;
                }
                i1 = self.my_loc.borrow().base().nb_law();
                i2 = 1;
            } else {
                i1 = index - 1;
                i2 = index;
            }

            let mut v1 = DVec3::ZERO;
            let mut v2 = DVec3::ZERO;
            let mut m1 = crate::geomalgo::geomfill::gp_mat::GpMat::identity();
            let mut m2 = crate::geomalgo::geomfill::gp_mat::GpMat::identity();

            let mut f = 0.0;
            let mut l = 0.0;
            self.my_loc.borrow().base().law(i1).borrow().get_domain(&mut f, &mut l);
            self.my_loc.borrow().base().law(i1).borrow().d0(l, &mut m1, &mut v1);
            let t1 = crate::brep_fill::brep_fill_location_law::gp_mat_column(&m1, 3);
            self.my_loc.borrow().base().law(i2).borrow().get_domain(&mut f, &mut l);
            self.my_loc.borrow().base().law(i2).borrow().d0(f, &mut m2, &mut v2);
            let t2 = crate::brep_fill::brep_fill_location_law::gp_mat_column(&m2, 3);

            let alpha = gv_angle(t1, t2);
            if (alpha > self.my_ang_max) || (alpha < self.my_ang_min) {
                // Angle too great => No "straight" connection
                // Angle too small => No connection
                return extrap; // = 0.0
            }

            // Sec = mySec->ConcatenedLaw();
            let sec = self.my_sec.borrow().concatened_law(brep);
            let Some(sec) = sec else { return extrap };

            // Calculating parameter U (cxx L3954-3960)
            let (mut lf, mut length) = (0.0f64, 0.0f64);
            self.my_loc.borrow_mut().base_mut().curvilinear_bounds(
                self.my_loc.borrow().base().nb_law(),
                &mut lf,
                &mut length,
            );
            let (mut sec_first, mut sec_len) = (0.0f64, 0.0f64);
            self.my_sec
                .borrow()
                .base()
                .law(1)
                .borrow()
                .get_domain(&mut sec_first, &mut sec_len);
            sec_len -= sec_first;
            let (mut lf2, mut ll) = (0.0f64, 0.0f64);
            self.my_loc
                .borrow_mut()
                .base_mut()
                .curvilinear_bounds(i1, &mut lf2, &mut ll);
            let u = sec_first + (ll / length) * sec_len;

            // Box(Sec, U, box); box.Get(Xmin, ..., Zmax);
            let (min, max) = box_of_section(&sec, u);

            let mut r = min[0]
                .abs()
                .max(max[0].abs())
                .max(min[1].abs().max(max[1].abs()));
            // double coef = 2.;
            let coef = 2.0;
            r *= coef;
            extrap = min[2].abs().max(max[2].abs()) + 100.0 * self.my_tol3d;
            extrap += r * (alpha / 2.0).tan();
        }
        extrap
    }

    // -------------------------------------------------------------------
    // Accessors (cxx L3534-3583)
    // -------------------------------------------------------------------

    /// OCCT IsDone (cxx L3534-3537).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// OCCT Shape (cxx L3541-3544).
    pub fn shape(&self) -> Shape {
        self.my_shape.clone()
    }

    /// OCCT ErrorOnSurface (cxx L3548-3551).
    pub fn error_on_surface(&self) -> f64 {
        self.error
    }

    /// OCCT SubShape (cxx L3555-3558) — myFaces.
    pub fn sub_shape(&self) -> Option<ShapeHArray2> {
        self.my_faces.clone()
    }

    /// OCCT InterFaces (cxx L3562-3565) — myUEdges.
    pub fn inter_faces(&self) -> Option<ShapeHArray2> {
        self.my_u_edges.clone()
    }

    /// OCCT Sections (cxx L3571-3574) — myVEdges.
    pub fn sections(&self) -> Option<ShapeHArray2> {
        self.my_v_edges.clone()
    }

    /// OCCT Tape (cxx L3580-3583).
    pub fn tape(&self, index: i32) -> Shape {
        self.my_tapes
            .as_ref()
            .map(|t| t[(index - 1) as usize].clone())
            .unwrap_or_else(Shape::null)
    }
}

/// OCCT gp_Vec::Angle (gp_XYZ::Angle).
fn gv_angle(v1: DVec3, v2: DVec3) -> f64 {
    let an_norm = v1.length();
    let a_no_norm = v2.length();
    let mut value = v1.dot(v2) / (an_norm * a_no_norm);
    if value > 1.0 {
        value = 1.0;
    } else if value < -1.0 {
        value = -1.0;
    }
    value.acos()
}

/// OCCT Geom_Curve::FirstParameter over the rcad Curve3.
fn curve_first_parameter_of(c: &Curve3) -> f64 {
    crate::geomalgo::geomfill::trihedron_law::curve_first_parameter(c)
}

/// OCCT Geom_Curve::LastParameter over the rcad Curve3.
fn curve_last_parameter_of(c: &Curve3) -> f64 {
    crate::geomalgo::geomfill::trihedron_law::curve_last_parameter(c)
}

// The gp::Resolution threshold consumed by the re-hosts.
const _: f64 = GP_RESOLUTION;

// The unused-import guards for the part-B symbols consumed from this module.
#[allow(unused_imports)]
use crate::brep_fill::brep_fill_pipe_shell_b::ShapeToArray2Map as _PartBShapeToArray2Map;
#[allow(dead_code)]
fn _part_b_type_guards(_t: Option<TrimmedCurve3>, _s: Option<HashSet<ShapeKey>>) {}

#[cfg(test)]
mod bezier_surface_iso_tests {
    //! Regression guard for the Geom_BezierSurface::UIso / VIso arms: the iso
    //! curve must reproduce the surface evaluation at the fixed parameter.

    use super::*;
    use rcad_kernel::geom::SurfaceEval;

    /// A 3 x 3 (degree 2 x 2) Bezier surface; `weights` selects the arm.
    fn sample_surface(weights: Vec<Vec<f64>>) -> BezierSurface {
        BezierSurface {
            control_points: vec![
                vec![
                    DVec3::new(0.0, 0.0, 0.0),
                    DVec3::new(0.0, 1.0, 0.5),
                    DVec3::new(0.0, 2.0, 0.0),
                ],
                vec![
                    DVec3::new(1.0, 0.0, 1.0),
                    DVec3::new(1.0, 1.0, 1.5),
                    DVec3::new(1.0, 2.0, 1.0),
                ],
                vec![
                    DVec3::new(2.0, 0.0, 0.0),
                    DVec3::new(2.0, 1.0, 0.5),
                    DVec3::new(2.0, 2.0, 0.0),
                ],
            ],
            weights,
        }
    }

    fn check_iso(iso: &Curve3, surf: &Surface3, fixed: f64, is_u: bool) {
        for jj in 0..=8 {
            let t = jj as f64 / 8.0;
            let p_iso = iso.point_at(t);
            let p_surf = if is_u {
                surf.point_at(fixed, t)
            } else {
                surf.point_at(t, fixed)
            };
            assert!(
                (p_iso - p_surf).length() < 1e-12,
                "t={t} iso={p_iso:?} surf={p_surf:?}"
            );
        }
    }

    /// The degree-2 Bernstein basis.
    fn bern2(t: f64) -> [f64; 3] {
        [(1.0 - t) * (1.0 - t), 2.0 * t * (1.0 - t), t * t]
    }

    /// An independent rational tensor-product Bezier reference,
    /// `sum_ij B_i(u) B_j(v) w_ij P_ij / sum_ij B_i(u) B_j(v) w_ij` — pole and
    /// weight read at the same indices, as OCCT's `Poles(i, j)` /
    /// `Weights(i, j)` in `BSplSLib::Iso` (BSplSLib.cxx L1678-1681).
    fn rational_reference(surf: &BezierSurface, u: f64, v: f64) -> DVec3 {
        let bu = bern2(u);
        let bv = bern2(v);
        let mut acc = DVec3::ZERO;
        let mut den = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                let w = bu[i] * bv[j] * surf.weights[i][j];
                acc += w * surf.control_points[i][j];
                den += w;
            }
        }
        acc / den
    }

    /// Check an iso curve against the rational reference along the varying
    /// direction (`is_u` = the UIso arm, fixed `u`).
    fn check_iso_against(iso: &Curve3, surf: &BezierSurface, fixed: f64, is_u: bool) {
        for jj in 0..=8 {
            let t = jj as f64 / 8.0;
            let p_iso = iso.point_at(t);
            let p_ref = if is_u {
                rational_reference(surf, fixed, t)
            } else {
                rational_reference(surf, t, fixed)
            };
            assert!(
                (p_iso - p_ref).length() < 1e-12,
                "t={t} iso={p_iso:?} ref={p_ref:?}"
            );
        }
    }

    /// The elementary-surface arms must follow `Geom_*Surface::UIso/VIso`,
    /// which are thin wrappers over the ElSLib iso constructors.  The
    /// surface kind of each iso also matters: a cylinder's U-iso is a LINE
    /// (the ruling) and its V-iso is a CIRCLE, and a sphere's U-iso is the
    /// meridian circle wrapped in a `Geom_TrimmedCurve(-M_PI/2, M_PI/2)`
    /// (Geom_SphericalSurface.cxx L292-297).
    #[test]
    fn elementary_surface_isos_follow_the_occt_wrappers() {
        use rcad_kernel::geom::{ConicalSurface, CylindricalSurface, Plane, SphericalSurface};
        use rcad_kernel::geom::ToroidalSurface;

        let axis = DVec3::Z;
        let x = DVec3::X;

        let cases: Vec<(Surface3, f64, f64)> = vec![
            (
                Surface3::Plane(Plane {
                    origin: DVec3::new(1.0, 2.0, 3.0),
                    normal: axis,
                    u_dir: x,
                    v_dir: axis.cross(x),
                }),
                0.4,
                0.7,
            ),
            (
                Surface3::Cylinder(CylindricalSurface {
                    origin: DVec3::new(1.0, 2.0, 3.0),
                    axis,
                    radius: 2.5,
                    ref_dir: x,
                    y_dir: None,
                }),
                0.4,
                0.7,
            ),
            (
                Surface3::Sphere(SphericalSurface {
                    center: DVec3::new(1.0, 2.0, 3.0),
                    axis,
                    radius: 2.5,
                    ref_dir: x,
                }),
                0.4,
                0.5,
            ),
            (
                Surface3::Cone(ConicalSurface::new_with_ref_dir(
                    DVec3::new(1.0, 2.0, 3.0),
                    axis,
                    2.5,
                    0.4,
                    x,
                )),
                0.4,
                0.5,
            ),
            (
                Surface3::Torus(ToroidalSurface {
                    center: DVec3::new(1.0, 2.0, 3.0),
                    axis,
                    ref_dir: x,
                    major_radius: 5.0,
                    minor_radius: 1.5,
                }),
                0.4,
                0.5,
            ),
        ];

        for (surf, u, v) in cases {
            let uiso = surface_uiso(&surf, u);
            check_iso(&uiso, &surf, u, true);
            let viso = surface_viso(&surf, v);
            check_iso(&viso, &surf, v, false);
        }

        // The cylinder U-iso is the ruling LINE and the V-iso is a CIRCLE
        // (Geom_CylindricalSurface.cxx L294-306); the two used to be swapped.
        let cy = Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::new(1.0, 2.0, 3.0),
            axis,
            radius: 2.5,
            ref_dir: x,
            y_dir: None,
        });
        assert!(matches!(surface_uiso(&cy, 0.4), Curve3::Line(_)));
        assert!(matches!(surface_viso(&cy, 0.7), Curve3::Circle(_)));

        // The sphere U-iso carries the OCCT TrimmedCurve latitude clamp.
        let sp = Surface3::Sphere(SphericalSurface {
            center: DVec3::new(1.0, 2.0, 3.0),
            axis,
            radius: 2.5,
            ref_dir: x,
        });
        let Curve3::Trimmed(t) = surface_uiso(&sp, 0.4) else {
            panic!("Geom_SphericalSurface::UIso is a Geom_TrimmedCurve")
        };
        assert!((t.first + 0.5 * std::f64::consts::PI).abs() < 1e-12);
        assert!((t.last - 0.5 * std::f64::consts::PI).abs() < 1e-12);
    }

    #[test]
    fn non_rational_uiso_and_viso_follow_the_surface() {
        let surf = Surface3::Bezier(sample_surface(vec![vec![1.0; 3]; 3]));
        let iso = surface_uiso(&surf, 0.3);
        let Curve3::Bezier(curve) = &iso else {
            panic!("UIso of a Bezier surface must be a Geom_BezierCurve")
        };
        // myPoles.RowLength() = the number of V poles.
        assert_eq!(curve.control_points.len(), 3);
        check_iso(&iso, &surf, 0.3, true);

        let iso = surface_viso(&surf, 0.7);
        let Curve3::Bezier(curve) = &iso else {
            panic!("VIso of a Bezier surface must be a Geom_BezierCurve")
        };
        // myPoles.ColLength() = the number of U poles.
        assert_eq!(curve.control_points.len(), 3);
        check_iso(&iso, &surf, 0.7, false);
    }

    #[test]
    fn u_rational_uiso_is_a_rational_bezier_curve() {
        // myURational: the weights vary along the second (V) pole index.
        let weights: Vec<Vec<f64>> = (0..3).map(|_| vec![1.0, 1.5, 3.0]).collect();
        let surf = sample_surface(weights);
        let iso = surface_uiso(&Surface3::Bezier(surf.clone()), 0.25);
        let Curve3::Bezier(curve) = &iso else {
            panic!("UIso of a Bezier surface must be a Geom_BezierCurve")
        };
        assert_eq!(curve.control_points.len(), 3);
        assert!(curve.weights.iter().any(|w| *w != 1.0));
        check_iso_against(&iso, &surf, 0.25, true);
    }

    /// The VIso rational arm pairs each pole with its weight at identical
    /// indices, as OCCT BSplSLib::Iso does: `Poles(j, index)` and
    /// `(*Weights)(j, index)` (BSplSLib.cxx L1678-1681).  The shared kernel
    /// helper `bspl_slib_iso` (rcad-kernel/src/math/bspl_lib.rs) originally
    /// read the weight grid transposed in its `is_u == false` branch; that
    /// was fixed to the plain `Weights(i, j)` accessor, so this test now
    /// guards the pairing.
    #[test]
    fn v_rational_viso_is_a_rational_bezier_curve() {
        // myVRational: the weights vary along the first (U) pole index.
        let weights: Vec<Vec<f64>> = vec![vec![1.0; 3], vec![1.5; 3], vec![3.0; 3]];
        let surf = sample_surface(weights);
        let iso = surface_viso(&Surface3::Bezier(surf.clone()), 0.4);
        let Curve3::Bezier(curve) = &iso else {
            panic!("VIso of a Bezier surface must be a Geom_BezierCurve")
        };
        assert_eq!(curve.control_points.len(), 3);
        assert!(curve.weights.iter().any(|w| *w != 1.0));
        check_iso_against(&iso, &surf, 0.4, false);
    }
}

#[cfg(test)]
#[path = "brep_fill_sweep_offset_iso_tests.rs"]
mod offset_surface_iso_tests;

#[cfg(test)]
#[path = "brep_fill_sweep_iso_tests.rs"]
mod bspline_surface_iso_tests;
