//! OCCT TopOpeBRepTool_CurveTool.cxx (TKBool/TopOpeBRepTool) — the pcurve
//! projection entry points consumed by the FC2D family
//! (`fillet/topopebrep_tool_2d.rs`).
//!
//! Translated here:
//! - the static `::MakePCurve(PC)` (TopOpeBRepTool_CurveTool.cxx L111-140)
//! - `TopOpeBRepTool_CurveTool::MakePCurveOnFace` (L1040-1155)
//!
//! Supporting dependencies:
//! - `BRepTools::UVBounds(F, ...)` / `AddUVBounds(F, B)` /
//!   `AddUVBounds(F, E, B)` (TKBRep/BRepTools/BRepTools.cxx L64-77,
//!   L123-158, L170-362) — the kernel canonical 1:1 bodies in
//!   `rcad-kernel/src/base/proj_lib/brep_adaptor.rs` (the former private
//!   copy of this file was deleted; the dedup rule keeps the complete
//!   aligned body only).
//! - `BRepAdaptor_Surface(F, R)` — the consumed subset of
//!   BRepAdaptor_Surface::Initialize (TKBRep/BRepAdaptor/
//!   BRepAdaptor_Surface.cxx L65-92).
//! - `BRepAdaptor_Curve(E, F)` — the consumed
//!   Adaptor3d_CurveOnSurface form of BRepAdaptor_Curve::Initialize
//!   (TKBRep/BRepAdaptor/BRepAdaptor_Curve.cxx L117-137), plus its Trim
//!   chain (BRepAdaptor_Curve.cxx L154-180 -> Adaptor3d_CurveOnSurface::
//!   Trim) consumed by ProjLib_ProjectedCurve::Perform TrimC3d.
//!
//! The ProjLib side of the chain is NOT re-translated here: the
//! ProjLib_ProjectedCurve class + the analytic ProjLib_Plane / Cylinder /
//! Cone / Sphere / Torus members + the Perform dispatch already live in
//! `rcad-kernel/src/base/proj_lib/` (proj_lib_projected_curve[_b].rs).
//! Their GAP leaves (ProjLib_ComputeApproxOnPolarSurface, the
//! ProjLib_HCompProjectedCurve + Approx_CurveOnSurface branch, the
//! ProjLib_ComputeApprox payload) keep the OCCT failure paths documented
//! there.
//!
//! GAP carriers (outside-scope deps, each keeps the OCCT failure path):
//! - none in this file.  `::MakePCurve` throws Standard_NotImplemented for
//!   a non-analytic projected type (CurveTool.cxx L134-136); the panic
//!   below is that failure path, not a carrier.

#![allow(dead_code)]

use std::sync::Arc;

use glam::DVec2;

use rcad_kernel::base::proj_lib::proj_lib_projected_curve::{
    Adaptor3dCurveGeom, GeomCurveAdaptor, GeomSurfaceAdaptor,
};
use rcad_kernel::base::proj_lib::proj_lib_projected_curve_b::{
    GeomCurveHandle, GeomSurfaceHandle, ProjLibProjectedCurve,
};
use rcad_kernel::base::proj_lib::{
    Adaptor3dCurve, Adaptor3dSurface, CurveType, CurveOnSurface, Geom2dCurveAdaptor,
};
use rcad_kernel::geom::{Curve2d, Curve2dEval as _, Plane, Surface3, TrimmedCurve2};
use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::topo::topods::{BRep, BRepTool as _, Shape};

// =========================================================================
// OCCT BRepAdaptor_Surface::Initialize (TKBRep/BRepAdaptor/
// BRepAdaptor_Surface.cxx L65-92) — the consumed constructor form.
//
// Architecture difference: the OCCT adaptor composes the face location
// through GeomAdaptor_TransformedSurface; rcad reads the world (located)
// face surface through BRep::face_surface_world, and the natural-bounds /
// UV-window domain through GeomSurfaceAdaptor (the base::proj_lib
// GeomAdaptor_Surface encoding consumed by ProjLib_ProjectedCurve).
// =========================================================================

/// OCCT BRepAdaptor_Surface(F, R).
pub(crate) fn brep_adaptor_surface(brep: &BRep, f: &Shape, restriction: bool) -> GeomSurfaceAdaptor {
    // OCCT L68-71: null face -> an undefined adaptor (rcad keeps the
    // natural-domain plane default of GeomSurfaceAdaptor::new — the call
    // sites never pass a null face).
    // OCCT L74: aSurface = BRep_Tool::Surface(F, L); null -> return.
    let Some(world) = brep.face_surface_world(f) else {
        return GeomSurfaceAdaptor::new(Surface3::Plane(Plane {
            origin: Default::default(),
            normal: glam::DVec3::Z,
            u_dir: glam::DVec3::X,
            v_dir: glam::DVec3::Y,
        }));
    };
    if restriction {
        // OCCT L77-80: UVBounds(F, ...) then Load(S, U1, U2, V1, V2, Trsf)
        // — the kernel canonical BRepTools::UVBounds
        // (base::proj_lib::brep_adaptor.rs, BRepTools.cxx L64-77).
        let (umin, umax, vmin, vmax) = rcad_kernel::base::proj_lib::brep_tools_uv_bounds(brep, f);
        let mut bas = GeomSurfaceAdaptor::new(world);
        bas.load_with_window(bas.surface.clone(), umin, umax, vmin, vmax);
        bas
    } else {
        // OCCT L82: Load(aSurface, L.Transformation()) — the natural bounds.
        GeomSurfaceAdaptor::new(world)
    }
}

// =========================================================================
// OCCT BRepAdaptor_Curve(E, F) (TKBRep/BRepAdaptor/BRepAdaptor_Curve.cxx
// L117-137) — the Adaptor3d_CurveOnSurface form (a 3D curve defined as
// the image of the E-on-F pcurve on the F surface).
// =========================================================================

/// OCCT BRepAdaptor_Curve — the (E, F) form whose payload is an
/// Adaptor3d_CurveOnSurface (myConSurf).  The 3D-curve form is not
/// consumed by this chain (the GeomCurveAdaptor covers it).
pub(crate) struct BRepAdaptorCurveOnFace {
    /// OCCT: Adaptor3d_CurveOnSurface myConSurf.
    pub my_con_surf: CurveOnSurface,
}

impl BRepAdaptorCurveOnFace {
    /// OCCT BRepAdaptor_Curve::Initialize(E, F) (L117-137).
    pub fn new(brep: &BRep, e: &Shape, f: &Shape) -> Self {
        // OCCT L123: S = BRep_Tool::Surface(F, L) — the LOCAL surface.
        let s = brep
            .face_surface(f)
            .expect("BRepAdaptor_Curve::Initialize: face without surface")
            .clone();
        // OCCT L124: PC = BRep_Tool::CurveOnSurface(E, F, pf, pl) — a null
        // PC throws Standard_NullObject further down (the OCCT null-handle
        // load); rcad keeps that failure as a panic.
        let (pc, pf, pl) = brep
            .curve_on_surface(e, f)
            .expect("BRepAdaptor_Curve::Initialize: no pcurve on the face");
        // OCCT L126-127: HS->Load(S); HC->Load(PC, pf, pl) — the restricted
        // Geom2dAdaptor domain is encoded with the Trimmed wrapper (the
        // adaptor.rs restriction precedent).
        let hs: GeomSurfaceHandle = Arc::new(GeomSurfaceAdaptor::new(s));
        let hc = Arc::new(Geom2dCurveAdaptor::new(Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(pc),
            t_min: pf,
            t_max: pl,
        })));
        // OCCT L128-129: myConSurf = new Adaptor3d_CurveOnSurface();
        // myConSurf->Load(HC, HS).
        BRepAdaptorCurveOnFace {
            my_con_surf: CurveOnSurface::new(hc, hs),
        }
    }

    /// OCCT Adaptor3d_CurveOnSurface::Trim(First, Last, Tol) — the 2D
    /// adaptor is trimmed, the surface is kept (Adaptor3d_CurveOnSurface.cxx
    /// Trim); consumed through BRepAdaptor_Curve::Trim
    /// (BRepAdaptor_Curve.cxx L154-180) by the TrimC3d of Perform.
    fn trim_impl(
        &self,
        first: f64,
        last: f64,
        tol: f64,
    ) -> BRepAdaptorCurveOnFace {
        BRepAdaptorCurveOnFace {
            my_con_surf: CurveOnSurface::new(
                self.my_con_surf.get_curve().trim(first, last, tol),
                self.my_con_surf.get_surface().clone(),
            ),
        }
    }
}

impl rcad_kernel::base::proj_lib::Adaptor3dCurve for BRepAdaptorCurveOnFace {
    fn first_parameter(&self) -> f64 {
        self.my_con_surf.first_parameter()
    }

    fn last_parameter(&self) -> f64 {
        self.my_con_surf.last_parameter()
    }

    fn value(&self, u: f64) -> glam::DVec3 {
        // OCCT Adaptor3d_CurveOnSurface::Value(U).
        self.my_con_surf.value(u)
    }

    fn d1(&self, u: f64) -> (glam::DVec3, glam::DVec3) {
        self.my_con_surf.d1(u)
    }

    fn d2(&self, u: f64) -> (glam::DVec3, glam::DVec3, glam::DVec3) {
        self.my_con_surf.d2(u)
    }

    fn continuity(&self) -> GeomAbsShape {
        self.my_con_surf.continuity()
    }

    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        self.my_con_surf.nb_intervals(s)
    }

    fn intervals(&self, s: GeomAbsShape) -> Vec<f64> {
        self.my_con_surf.intervals(s)
    }

    fn trim(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor3dCurve> {
        Arc::new(self.trim_impl(first, last, tol))
    }
}

impl Adaptor3dCurveGeom for BRepAdaptorCurveOnFace {
    /// OCCT Adaptor3d_CurveOnSurface::GetType() — GeomAbs_OtherCurve.
    fn get_type(&self) -> CurveType {
        CurveType::Other
    }

    // OCCT: the geometry accessors of a curve-on-surface composite are
    // invalid (no elementary kind); the Project dispatch checks GetType()
    // first and never reaches them.
    fn line(&self) -> rcad_kernel::geom::Line3 {
        panic!("BRepAdaptor_Curve(E, F): Line() on a curve-on-surface composite")
    }

    fn circle(&self) -> rcad_kernel::geom::Circle3 {
        panic!("BRepAdaptor_Curve(E, F): Circle() on a curve-on-surface composite")
    }

    fn ellipse(&self) -> rcad_kernel::geom::Ellipse3 {
        panic!("BRepAdaptor_Curve(E, F): Ellipse() on a curve-on-surface composite")
    }

    fn parabola(&self) -> rcad_kernel::geom::Parabola3 {
        panic!("BRepAdaptor_Curve(E, F): Parabola() on a curve-on-surface composite")
    }

    fn hyperbola(&self) -> rcad_kernel::geom::Hyperbola3 {
        panic!("BRepAdaptor_Curve(E, F): Hyperbola() on a curve-on-surface composite")
    }

    fn trim_geom(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor3dCurveGeom> {
        Arc::new(self.trim_impl(first, last, tol))
    }
}

// =========================================================================
// OCCT gp_Trsf2d::SetMirror(gp_Ax2d) + Geom2d_Curve::Transform — the
// exact per-type composition consumed by the MakePCurveOnFace sphere
// quasiperiodic block.
// =========================================================================

/// OCCT gp_Trsf2d mirror about the 2D axis (axis_loc, axis_dir): the point
/// transform p' = A + 2*((p-A).D)D - (p-A).
fn trsf2d_mirror_point(axis_loc: DVec2, axis_dir: DVec2, p: DVec2) -> DVec2 {
    let d = p - axis_loc;
    axis_loc + axis_dir * (2.0 * d.dot(axis_dir)) - d
}

/// OCCT gp_Trsf2d mirror applied to a direction (no translation part).
fn trsf2d_mirror_dir(axis_dir: DVec2, v: DVec2) -> DVec2 {
    axis_dir * (2.0 * v.dot(axis_dir)) - v
}

/// OCCT Geom2d_Curve::Transform(gp_Trsf2d) for the mirror of the sphere
/// block, per dynamic type.  For the elementary conics the OCCT transform
/// maps the frame (gp_Ax22d::Transform maps location + XDir + YDir as
/// vectors), so mirroring the defining frame is pointwise exact.
fn curve2d_transform_mirror(c: &Curve2d, axis_loc: DVec2, axis_dir: DVec2) -> Curve2d {
    match c {
        // OCCT Geom2d_Line::Transform: the position is transformed.
        Curve2d::Line(l) => {
            let origin = trsf2d_mirror_point(axis_loc, axis_dir, l.origin);
            let dir = trsf2d_mirror_dir(axis_dir, l.direction).normalize_or_zero();
            Curve2d::Line(rcad_kernel::geom::Line2d { origin, direction: dir })
        }
        // OCCT Geom2d_Circle::Transform: radius * |ScaleFactor| (1 for the
        // mirror) + the transformed position frame.
        Curve2d::Circle(cir) => Curve2d::Circle(rcad_kernel::geom::Circle2d {
            center: trsf2d_mirror_point(axis_loc, axis_dir, cir.center),
            x_dir: trsf2d_mirror_dir(axis_dir, cir.x_dir).normalize_or_zero(),
            y_dir: trsf2d_mirror_dir(axis_dir, cir.y_dir).normalize_or_zero(),
            radius: cir.radius,
        }),
        // OCCT Geom2d_Ellipse::Transform — same frame mapping.
        Curve2d::Ellipse(el) => Curve2d::Ellipse(rcad_kernel::geom::Ellipse2d {
            center: trsf2d_mirror_point(axis_loc, axis_dir, el.center),
            major_dir: trsf2d_mirror_dir(axis_dir, el.major_dir).normalize_or_zero(),
            minor_dir: trsf2d_mirror_dir(axis_dir, el.minor_dir).normalize_or_zero(),
            major_radius: el.major_radius,
            minor_radius: el.minor_radius,
        }),
        // OCCT Geom2d_Parabola / Geom2d_Hyperbola::Transform keep a full
        // (XDir, YDir) frame; the rcad structs imply YDir from XDir with a
        // fixed handedness, so a handedness-flipping mirror is not
        // representable.  Unreachable from the sphere block: the analytic
        // ProjLib_Sphere result is a Line, and the approximation fallback
        // (a BSpline) is a GAP'd dependency that fails before this block.
        Curve2d::Parabola(_) | Curve2d::Hyperbola(_) => {
            panic!("Geom2d_Curve::Transform(mirror): the rcad Parabola2d/Hyperbola2d cannot represent a handedness-flipped frame")
        }
        // OCCT Geom2d_BSplineCurve::Transform: the poles are transformed
        // (the mirror is affine, so this is pointwise exact).
        Curve2d::BSpline(b) => {
            let mut b = b.clone();
            for p in b.control_points.iter_mut() {
                *p = trsf2d_mirror_point(axis_loc, axis_dir, *p);
            }
            Curve2d::BSpline(b)
        }
        // OCCT Geom2d_BezierCurve::Transform: the poles are transformed.
        Curve2d::Bezier(b) => {
            let mut b = b.clone();
            for p in b.control_points.iter_mut() {
                *p = trsf2d_mirror_point(axis_loc, axis_dir, *p);
            }
            Curve2d::Bezier(b)
        }
        _ => panic!("Geom2d_Curve::Transform(mirror): curve kind not translated"),
    }
}

// =========================================================================
// OCCT TopOpeBRepTool_CurveTool.cxx L111-140 — static ::MakePCurve
// =========================================================================

/// OCCT static MakePCurve(const ProjLib_ProjectedCurve& PC)
/// (TopOpeBRepTool_CurveTool.cxx L111-140) — the ProjLib_ProjectedCurve ->
/// Geom2d_Curve conversion.  The OCCT Standard_NotImplemented throw of the
/// `default` arm (L134-136) is the panic below.
pub(crate) fn curve_tool_make_pcurve(pc: &ProjLibProjectedCurve) -> Option<Curve2d> {
    // OCCT L113: occ::handle<Geom2d_Curve> C2D;
    let c2d: Option<Curve2d>;
    // OCCT L114-137: the type switch.
    match pc.get_type() {
        CurveType::Line => c2d = Some(Curve2d::Line(pc.line())),
        CurveType::Circle => c2d = Some(Curve2d::Circle(pc.circle())),
        CurveType::Ellipse => c2d = Some(Curve2d::Ellipse(pc.ellipse())),
        CurveType::Parabola => c2d = Some(Curve2d::Parabola(pc.parabola())),
        CurveType::Hyperbola => c2d = Some(Curve2d::Hyperbola(pc.hyperbola())),
        // OCCT L131-133: C2D = PC.BSpline() — a null handle stays null.
        CurveType::BSpline => c2d = pc.bspline().map(Curve2d::BSpline),
        // OCCT L134-136: default -> throw Standard_NotImplemented.
        _ => panic!("Standard_NotImplemented: CurveTool::MakePCurve"),
    }
    // OCCT L138: return C2D;
    c2d
}

// =========================================================================
// OCCT TopOpeBRepTool_CurveTool.cxx L1040-1155 — MakePCurveOnFace
// =========================================================================

/// OCCT TopOpeBRepTool_CurveTool::MakePCurveOnFace(S, C3D, TolReached2d,
/// first, last) (L1040-1155).
pub(crate) fn curve_tool_make_pcurve_on_face(
    brep: &BRep,
    s: &Shape,
    c3d: &rcad_kernel::geom::Curve3,
    tol_reached2d: &mut f64,
    first: f64,
    last: f64,
) -> Option<Curve2d> {
    // OCCT L1048-1052: bool trim = (first < last).
    let trim = first < last;

    // OCCT L1054-1056: F = TopoDS::Face(S); BAS(F, compminmaxUV = false).
    let f = s;
    let compminmax_uv = false;
    let bas = brep_adaptor_surface(brep, f, compminmax_uv);
    // OCCT L1057-1065: GAC.Load(C3D[, first, last]).
    let gac = if trim {
        GeomCurveAdaptor::with_range(c3d.clone(), first, last)
    } else {
        GeomCurveAdaptor::new(c3d.clone())
    };
    // OCCT L1066-1067: BAHS = new BRepAdaptor_Surface(BAS);
    // BAHC = new GeomAdaptor_Curve(GAC).
    let bahs: GeomSurfaceHandle = Arc::new(bas.clone());
    let bahc: GeomCurveHandle = Arc::new(gac.clone());
    // OCCT L1068: ProjLib_ProjectedCurve projcurv(BAHS, BAHC).
    let projcurv = ProjLibProjectedCurve::with_surface_curve(bahs.clone(), bahc);
    // OCCT L1069-1070: C2D = ::MakePCurve(projcurv);
    // TolReached2d = projcurv.GetTolerance().
    let mut c2d = curve_tool_make_pcurve(&projcurv);
    *tol_reached2d = projcurv.get_tolerance();

    // OCCT L1072-1073: BRepTools::UVBounds(F, UMin, UMax, VMin, VMax) —
    // the kernel canonical body (base::proj_lib::brep_adaptor.rs,
    // BRepTools.cxx L64-77).
    let (u_min, u_max, v_min, v_max) = rcad_kernel::base::proj_lib::brep_tools_uv_bounds(brep, f);

    // OCCT L1075-1081: the mid-parameter test point.
    let f_par = gac.first_parameter();
    let l_par = gac.last_parameter();
    let t = (f_par + l_par) * 0.5;
    let (mut u2, mut v2) = match &c2d {
        Some(c) => {
            let p = c.point_at(t);
            (p.x, p.y)
        }
        // OCCT dereferences the null C2D here; rcad panics the same way.
        None => panic!("TopOpeBRepTool_CurveTool::MakePCurveOnFace: null pcurve"),
    };

    // OCCT L1083-1116: the sphere quasiperiodic shift (MSV).
    if bas.get_type() == rcad_kernel::base::proj_lib::GeomAbsSurfaceType::Sphere {
        // OCCT L1086-1090.
        let v_first = bas.first_v_parameter();
        let v_last = bas.last_v_parameter();
        let mincond = v2 < v_first;
        let maxcond = v2 > v_last;
        if mincond || maxcond {
            // OCCT L1092-1101: the copy + the mirror about the apex isoline
            // -PI/2 or PI/2.
            let pc_t = match &c2d {
                Some(c) => curve2d_transform_mirror(c, DVec2::new(0.0, {
                    // OCCT L1095-1099: po(0, -PI/2); maxcond -> po.SetY(PI/2).
                    if maxcond {
                        std::f64::consts::FRAC_PI_2
                    } else {
                        -std::f64::consts::FRAC_PI_2
                    }
                }), glam::DVec2::X),
                None => unreachable!(),
            };
            // OCCT L1102-1109: the PI translation along U, reversed when
            // the test point already passed the +PI half.
            let mut vec = DVec2::new(std::f64::consts::PI, 0.0);
            let u_first = bas.first_u_parameter();
            if u2 - u_first - std::f64::consts::PI > -1.0e-7 {
                vec = -vec;
            }
            let mut c2d_t = Some(pc_t);
            if let Some(curv) = c2d_t.as_mut() {
                *curv = rcad_kernel::geom::translate_curve2d(curv, vec);
            }
            // OCCT L1110-1114: C2D = PCT; recompute the test point.
            c2d = c2d_t;
            let p = c2d.as_ref().unwrap().point_at(t);
            u2 = p.x;
            v2 = p.y;
        }
    }

    // OCCT L1118-1131: the U decal for a U-periodic surface.
    let mut du = 0.0_f64;
    if bahs.is_u_periodic() {
        // OCCT L1122-1127.
        let mincond = u_min - u2 > 1.0e-7;
        let maxcond = u2 - u_max > 1.0e-7;
        let decalu = mincond || maxcond;
        if decalu {
            du = if mincond { bahs.u_period() } else { -bahs.u_period() };
        }
    }
    // OCCT L1132-1145: the V decal for a V-periodic surface.
    let mut dv = 0.0_f64;
    if bahs.is_v_periodic() {
        // OCCT L1136-1141.
        let mincond = v_min - v2 > 1.0e-7;
        let maxcond = v2 - v_max > 1.0e-7;
        let decalv = mincond || maxcond;
        if decalv {
            dv = if mincond { bahs.v_period() } else { -bahs.v_period() };
        }
    }

    // OCCT L1147-1152.
    if du != 0.0 || dv != 0.0 {
        if let Some(curv) = c2d.as_mut() {
            *curv = rcad_kernel::geom::translate_curve2d(curv, DVec2::new(du, dv));
        }
    }

    // OCCT L1154: return C2D.
    c2d
}
