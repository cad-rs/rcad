//! OCCT BRepFill_Sweep.cxx, part B (file statics, cxx L96-1876) — the
//! companion of [`super::brep_fill_sweep`] (part A).  Translated units:
//! NumberOfPoles, Couture, CheckSameParameter, CheckSameParameterExact,
//! SameParameter, CorrectSameParameter, Oriente, UpdateEdgeOnPlane,
//! BuildFace, BuildEdge (3d form), Filling, Substitute, KeepEdge, BuildEdge
//! (iso form), UpdateEdge, IsDegen.  (Translate / Box / BuildVertex /
//! NullEdge / HasPCurves / ReverseEdgeInFirstOrLastWire /
//! ReverseModifiedEdges are in part A.)
//!
//! Architecture differences:
//! - `handle(Geom_Surface)` / `handle(Geom_Curve)` / `handle(Geom2d_Curve)`
//!   map to the rcad value types `Surface3` / `Curve3` / `Curve2d`; the
//!   OCCT null handle maps to `Option<...>`.
//! - `GeomAdaptor_Curve(C, f, l)` / `BRepAdaptor_Curve(E)` map to the
//!   (curve value, parameter range) pair read through the pool.
//! - OCCT stores a curve-on-surface representation as (surface, location)
//!   without a `TopoDS_Face`; the rcad pool keys pcurves by face identity.
//!   The sweep mints an anonymous face TShape per surface pair
//!   ([`intern_surface_face`]) — the (S, L) carrier — and every
//!   face-resolved read resolves through the same surface identity.
//! - `TopoDS_Iterator` maps to [`topods_iterator`] (the direct sub-shapes).
//! - `BRep_Builder` maps to the pool-leading `BRepBuilder` convention
//!   (part A).

use std::collections::HashMap;
use std::sync::Arc;

use glam::{DVec2, DVec3};

use rcad_kernel::core::precision::{p_confusion, CONFUSION};
use rcad_kernel::geom::{
    Curve2d, Curve2dEval, Curve3, CurveEval, Line2d, Line3, Plane, Surface3, SurfaceEval,
    TrimmedCurve2, TrimmedCurve3,
};
use rcad_kernel::base::gcpnts::abscissa_point::arc_length;
use rcad_kernel::topo::topods::{
    surface_adaptor_basis_and_bounds, surface_same, BRep, BRepBuilder, Orientation, Shape, TShape,
};

use crate::brep_fill::brep_fill_pipe_shell_b::ShapeHArray2;
use crate::brep_fill::brep_fill_sweep::{has_pcurves, null_edge};
use crate::brep_fill::generator::ShapeKey;

// ---------------------------------------------------------------------------
// Shared tiny helpers (gp / TopoDS pure-math re-hosts)
// ---------------------------------------------------------------------------

/// OCCT TopoDS_Iterator — the direct sub-shapes of a shape.
pub(super) fn topods_iterator(s: &Shape) -> Vec<Shape> {
    if s.is_null() {
        return Vec::new();
    }
    match s.data.as_ref() {
        TShape::Vertex(_) => Vec::new(),
        TShape::Edge(ed) => {
            let mut out = Vec::new();
            if !ed.first.is_null() {
                out.push(ed.first.clone());
            }
            if !ed.last.is_null() {
                out.push(ed.last.clone());
            }
            out
        }
        TShape::Wire(wd) => wd.edges.clone(),
        TShape::Face(fd) => {
            let mut out = Vec::new();
            if !fd.outer_wire.is_null() {
                out.push(fd.outer_wire.clone());
            }
            out.extend(fd.inner_wires.iter().cloned());
            out
        }
        TShape::Shell(sd) => sd.faces.clone(),
        TShape::Solid(sd) => sd.shells.clone(),
        TShape::CompSolid(sd) => sd.clone(),
        TShape::Compound(cd) => cd.clone(),
    }
}

/// OCCT gp_Ax2 — Location + X/Y/Direction.
#[derive(Debug, Clone, Copy)]
pub(super) struct GpAx2 {
    /// OCCT Location.
    pub location: DVec3,
    /// OCCT Direction (the main direction).
    pub direction: DVec3,
    /// OCCT XDirection.
    pub x_direction: DVec3,
    /// OCCT YDirection.
    pub y_direction: DVec3,
}

impl GpAx2 {
    /// OCCT gp_Ax2(P, N, Dx) — YDirection is built N ^ Dx (gp_Ax2.cxx).
    pub fn new(p: DVec3, n: DVec3, dx: DVec3) -> Self {
        let n = n.normalize_or_zero();
        let dx = dx.normalize_or_zero();
        GpAx2 {
            location: p,
            direction: n,
            x_direction: dx,
            y_direction: n.cross(dx).normalize_or_zero(),
        }
    }

    /// OCCT gp_Ax2(P, N) — the X direction is solver-picked; the rcad
    /// default keeps a deterministic orthogonal.
    pub fn from_normal(p: DVec3, n: DVec3) -> Self {
        let n = n.normalize_or_zero();
        let dx = if n.x.abs() < 0.9 {
            DVec3::X.cross(n)
        } else {
            DVec3::Y.cross(n)
        }
        .normalize_or_zero();
        GpAx2::new(p, n, dx)
    }
}

/// OCCT gp_Trsf::SetTransformation(gp_Ax2) — the transform of a world point
/// into the Ax2 coordinate frame (pure math).
pub(super) fn point_in_ax2_frame(axe: &GpAx2, p: DVec3) -> DVec3 {
    let d = p - axe.location;
    DVec3::new(
        d.dot(axe.x_direction),
        d.dot(axe.y_direction),
        d.dot(axe.direction),
    )
}

/// OCCT gp_Vec::Angle — the unsigned angle [0, PI] between two vectors
/// (gp_XYZ::Angle, pure math).
pub(super) fn gp_vec_angle(v1: DVec3, v2: DVec3) -> f64 {
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

/// OCCT gp_Vec2d::Angle — the signed angle (-PI, PI] between two 2d vectors
/// (pure math).
pub(super) fn gp_vec2d_angle(v1: DVec2, v2: DVec2) -> f64 {
    v2.y.atan2(v2.x) - v1.y.atan2(v1.x)
}

/// OCCT gp_Vec2d::IsParallel(TheOther, AngularTolerance) — collinear in
/// either direction within the angular tolerance (pure math; expressed
/// through the angle test of gp_Vec2d::Angle).
pub(super) fn gp_vec2d_is_parallel(v1: DVec2, v2: DVec2, angular_tolerance: f64) -> bool {
    let ang = gp_vec2d_angle(v1, v2).abs();
    ang <= angular_tolerance || (std::f64::consts::PI - ang).abs() <= angular_tolerance
}

/// OCCT gp_Vec2d::IsOpposite(TheOther, AngularTolerance) — opposite
/// directions within the angular tolerance (pure math).
pub(super) fn gp_vec2d_is_opposite(v1: DVec2, v2: DVec2, angular_tolerance: f64) -> bool {
    let ang = gp_vec2d_angle(v1, v2).abs();
    (std::f64::consts::PI - ang).abs() <= angular_tolerance
}

/// OCCT IsEqual(v1, v2) over Standard_Real — |v1 - v2| < RealSmall()
/// (Standard_Real.hxx L148-151).
pub(super) fn is_equal_real(v1: f64, v2: f64) -> bool {
    (v1 - v2).abs() < f64::MIN_POSITIVE
}

/// The 2d parameter bounds stored with a pcurve (the OCCT GCurve
/// First/Last of the representation). The unbounded-curve fallback mirrors
/// the OCCT Geom2d_Line parameters (Geom2d_Line.cxx L142-151:
/// -/+ Precision::Infinite()).
pub(super) fn curve2d_param_bounds(pc: &Curve2d) -> (f64, f64) {
    match pc {
        Curve2d::Trimmed(t) => (t.t_min, t.t_max),
        _ => (
            -rcad_kernel::core::precision::INFINITE_VALUE,
            rcad_kernel::core::precision::INFINITE_VALUE,
        ),
    }
}

/// The (S, L) carrier: OCCT BRep_CurveOnSurface stores (surface, location)
/// without a face; the rcad pool keys pcurves by face identity, so the
/// sweep interns an anonymous face TShape per surface (surface_same) and
/// keys by its identity.
pub(super) fn intern_surface_face(brep: &mut BRep, s: &Surface3) -> Shape {
    for (idx, ts) in brep.tshapes.iter().enumerate() {
        if let TShape::Face(fd) = ts.as_ref() {
            if let Some(fs) = &fd.surface {
                if surface_same(fs, s) {
                    return brep.shape_at(idx);
                }
            }
        }
    }
    let tshape = Arc::new(TShape::Face(rcad_kernel::topo::topods::TFaceData {
        my_shapes: Vec::new(),
        flags: rcad_kernel::topo::topods::tshape_flags::DEFAULT,
        surface: Some(s.clone()),
        surface_location: 0,
        outer_wire: Shape::null(),
        inner_wires: Vec::new(),
        sample_point: None,
        uv_domain: None,
        internal_vertices: Vec::new(),
        tolerance: 0.0,
        natural_restriction: false,
    }));
    let index = brep.tshapes.len();
    brep.tshapes.push(tshape);
    brep.shape_at(index)
}

/// OCCT BRep_Builder::UpdateEdge(E, C2d, S, L, Tol) over a bare surface
/// (BRep_Builder.cxx L660-700) — the curve-on-surface representation is
/// appended (the anonymous (S, L) face carries the key).
pub(super) fn update_edge_pcurve_on_surf(
    brep: &mut BRep,
    e: &Shape,
    c2d: &Curve2d,
    s: &Surface3,
    loc: u32,
    tol: f64,
) {
    let f = intern_surface_face(brep, s);
    let key = (f.ptr_id(), loc);
    let (ta, tb) = curve2d_param_bounds(c2d);
    let ed = brep.edge_mut_inplace(e.clone());
    ed.pcurves.insert(key, (c2d.clone(), ta, tb));
    ed.representations
        .push(rcad_kernel::topo::topods::CurveRepresentation::CurveOnSurface {
            face: key,
            pcurve: c2d.clone(),
            range: [ta, tb],
        });
    ed.tolerance = ed.tolerance.max(tol);
}

/// OCCT BRep_Builder::UpdateEdge(E, C2d1, C2d2, S, L, Tol) — the seam
/// (closed-surface) two-pcurve form (BRep_Builder.cxx L704-741).
pub(super) fn update_edge_pcurve2_on_surf(
    brep: &mut BRep,
    e: &Shape,
    c2d1: &Curve2d,
    c2d2: &Curve2d,
    s: &Surface3,
    loc: u32,
    tol: f64,
) {
    let f = intern_surface_face(brep, s);
    let key = (f.ptr_id(), loc);
    let (ta, tb) = curve2d_param_bounds(c2d1);
    let ed = brep.edge_mut_inplace(e.clone());
    ed.pcurves.insert(key, (c2d1.clone(), ta, tb));
    ed.representations
        .push(rcad_kernel::topo::topods::CurveRepresentation::CurveOnClosedSurface {
            face: key,
            pcurve1: c2d1.clone(),
            pcurve2: c2d2.clone(),
            range: [ta, tb],
        });
    ed.tolerance = ed.tolerance.max(tol);
}

/// OCCT BRep_Tool::CurveOnSurface(E, S, L, f, l) over a bare surface — the
/// pcurve of E on the surface S (identity by surface_same + location).
pub(super) fn curve_on_surface_by_surf(
    brep: &BRep,
    e: &Shape,
    s: &Surface3,
    loc: u32,
) -> Option<(Curve2d, f64, f64)> {
    let ed = brep.edge(e.clone());
    for (key, (pc, tf, tl)) in ed.pcurves.iter() {
        if key.1 != loc {
            continue;
        }
        if let Some(idx) = brep.index_by_ptr(key.0) {
            if let TShape::Face(fd) = brep.tshapes[idx].as_ref() {
                if let Some(fs) = &fd.surface {
                    if surface_same(fs, s) {
                        return Some((pc.clone(), *tf, *tl));
                    }
                }
            }
        }
    }
    None
}

/// OCCT BRep_TVertex::Tolerance(theTol) — the SETTER (not the max-update):
/// BRep_Builder keeps the max, the sweep lowers vertex tolerances explicitly
/// (Filling cxx L1203).
pub(super) fn set_vertex_tolerance(brep: &mut BRep, v: &Shape, tol: f64) {
    brep.vertex_mut(v.clone()).tolerance = tol;
}

/// OCCT BRep_TEdge::Tolerance(theTol) — the SETTER (RebuildTopOrBottomEdge
/// cxx L4105).
pub(super) fn set_edge_tolerance(brep: &mut BRep, e: &Shape, tol: f64) {
    brep.edge_mut_inplace(e.clone()).tolerance = tol;
}

/// OCCT GeomAdaptor default domain over the rcad Curve3.
pub(super) fn curve3_domain(c: &Curve3) -> (f64, f64) {
    let d = c.default_domain();
    (d[0], d[1])
}

// ---------------------------------------------------------------------------
// NumberOfPoles (cxx L96-152)
// ---------------------------------------------------------------------------/// OCCT static NumberOfPoles (cxx L96-152).
pub(super) fn number_of_poles(brep: &BRep, w: &Shape) -> i32 {
    let mut nb_points = 0i32;

    // TopoDS_Iterator iter(W);
    for it in topods_iterator(w) {
        // BRepAdaptor_Curve c(TopoDS::Edge(iter.Value()));
        let ed = brep.edge(it.clone());
        let c = ed
            .curve
            .clone()
            .expect("BRepAdaptor_Curve: edge without 3d curve");
        let df_uf = ed.range[0];
        let df_ul = ed.range[1];
        if is_equal_real(df_uf, df_ul) {
            // Degenerate
            continue;
        }

        match &c {
            // case GeomAbs_BezierCurve: put all poles for bezier
            Curve3::Bezier(gc) => {
                let i_nb_pol = gc.control_points.len() as i32;
                if i_nb_pol >= 2 {
                    nb_points += i_nb_pol;
                }
            }
            // case GeomAbs_BSplineCurve: put all poles for bspline
            Curve3::BSpline(gc) => {
                let i_nb_pol = gc.control_points.len() as i32;
                if i_nb_pol >= 2 {
                    nb_points += i_nb_pol;
                }
            }
            // case GeomAbs_Line
            Curve3::Line(_) => {
                nb_points += 2;
            }
            // case GeomAbs_Circle / Ellipse / Hyperbola / Parabola
            Curve3::Circle(_) | Curve3::Ellipse(_) | Curve3::Hyperbola(_) | Curve3::Parabola(_) => {
                nb_points += 4;
            }
            // default: 15 + c.NbIntervals(GeomAbs_C3) — the BRepAdaptor
            // interval computation is a GAP carrier (plan section 0.6).
            _ => panic!(
                "GAP: BRepAdaptor_Curve::NbIntervals(GeomAbs_C3) (TKBRep/\
                 BRepAdaptor) — BRepFill_Sweep::NumberOfPoles default branch"
            ),
        }
    }

    nb_points
}

// ---------------------------------------------------------------------------
// Couture (cxx L213-244)
// ---------------------------------------------------------------------------

/// OCCT static Couture (cxx L213-244) — check if E is an edge of sewing on S
/// and make the representation HadHoc.
pub(super) fn couture(brep: &BRep, e: &Shape, s: &Surface3, l: u32) -> Option<Curve2d> {
    // bool Eisreversed = (E.Orientation() == TopAbs_REVERSED);
    let e_is_reversed = e.orientation == Orientation::Reversed;

    // find the representation (the ChangeCurves list walk)
    let ed = brep.edge(e.clone());
    for cr in &ed.representations {
        match cr {
            rcad_kernel::topo::topods::CurveRepresentation::CurveOnSurface {
                face,
                pcurve,
                ..
            } => {
                if representation_on_surface(brep, *face, s, l) {
                    return Some(pcurve.clone());
                }
            }
            rcad_kernel::topo::topods::CurveRepresentation::CurveOnClosedSurface {
                face,
                pcurve1,
                pcurve2,
                ..
            } => {
                if representation_on_surface(brep, *face, s, l) {
                    // if (GC->IsCurveOnClosedSurface() && Eisreversed)
                    if e_is_reversed {
                        return Some(pcurve2.clone());
                    } else {
                        return Some(pcurve1.clone());
                    }
                }
            }
            _ => {}
        }
    }
    None // pc.Nullify()
}

/// OCCT BRep_CurveRepresentation::IsCurveOnSurface(S, L) — the
/// (surface, location) identity test over the rcad representation key.
fn representation_on_surface(brep: &BRep, face_key: (u64, u32), s: &Surface3, l: u32) -> bool {
    if face_key.1 != l {
        return false;
    }
    if let Some(idx) = brep.index_by_ptr(face_key.0) {
        if let TShape::Face(fd) = brep.tshapes[idx].as_ref() {
            if let Some(fs) = &fd.surface {
                return surface_same(fs, s);
            }
        }
    }
    false
}

// ---------------------------------------------------------------------------
// CheckSameParameter (cxx L251-282)
// ---------------------------------------------------------------------------

/// OCCT static CheckSameParameter (cxx L251-282) — check a posteriori that
/// sameparameter has worked correctly.  (The Adaptor3d_Curve /
/// Adaptor3d_Surface handles map to the rcad (curve, range) pair and the
/// Surface3 value.)
pub(super) fn check_same_parameter(
    c3d: &Curve3,
    c3d_first: f64,
    c3d_last: f64,
    pcurv: &Curve2d,
    s: &Surface3,
    tol3d: f64,
    tolreached: &mut f64,
) -> bool {
    *tolreached = 0.0;
    let f = c3d_first;
    let l = c3d_last;
    let nbp = 45i32;
    let step = 1.0 / (nbp as f64 - 1.0);
    for i in 0..nbp {
        let mut t;
        t = step * i as f64;
        t = (1.0 - t) * f + t * l;
        let (u, v) = {
            let uv = pcurv.point_at(t);
            (uv.x, uv.y)
        };
        let p_s = s.point_at(u, v);
        let p_c = c3d.point_at(t);
        let d2 = p_s.distance_squared(p_c);
        *tolreached = tolreached.max(d2);
    }
    *tolreached = tolreached.sqrt();
    if *tolreached > tol3d {
        *tolreached *= 2.0;
        return false;
    }
    *tolreached *= 2.0;
    *tolreached = tolreached.max(CONFUSION);
    true
}

// ---------------------------------------------------------------------------
// CheckSameParameterExact (cxx L289-310)
// ---------------------------------------------------------------------------

/// OCCT static CheckSameParameterExact (cxx L289-310) — the exact
/// calculation method of edge tolerance.
pub(super) fn check_same_parameter_exact(
    c3d: &Curve3,
    c3d_first: f64,
    c3d_last: f64,
    curve_on_surface: &(Curve2d, Surface3),
    tol3d: f64,
    tolreached: &mut f64,
) -> bool {
    let _ = (c3d_first, c3d_last);
    // GeomLib_CheckCurveOnSurface aCheckCurveOnSurface(C3d);
    let mut a_check_curve_on_surface =
        crate::geomalgo::geom_lib_check_curve_on_surface::GeomLibCheckCurveOnSurface::new(c3d);
    // aCheckCurveOnSurface.SetParallel(false);
    a_check_curve_on_surface.set_parallel(false);
    // aCheckCurveOnSurface.Perform(curveOnSurface);
    a_check_curve_on_surface.perform(curve_on_surface);

    // tolreached = aCheckCurveOnSurface.MaxDistance();
    *tolreached = a_check_curve_on_surface.max_distance();

    if *tolreached > tol3d {
        return false;
    } else {
        *tolreached = tolreached.max(CONFUSION);
        *tolreached *= 1.05;
    }
    true
}

// ---------------------------------------------------------------------------
// SameParameter (cxx L319-392)
// ---------------------------------------------------------------------------

/// OCCT static SameParameter (cxx L319-392) — encapsulation of Sameparameter.
/// Boolean informs if the pcurve was computed or not...  The tolerance is
/// always OK.  (Pcurv is the in-out handle.)
pub(super) fn same_parameter(
    brep: &mut BRep,
    e: &Shape,
    pcurv: &mut Curve2d,
    surf: &Surface3,
    tol3d: f64,
    tolreached: &mut f64,
) -> bool {
    // occ::handle<Geom_Curve> C3d = BRep_Tool::Curve(E, f, l);
    let ed = brep.edge(e.clone());
    let mut c3d = ed
        .curve
        .clone()
        .expect("SameParameter: edge without 3d curve");
    let (f, l) = (ed.range[0], ed.range[1]);
    // GeomAdaptor_Curve GAC3d(C3d, f, l); HC3d = new(GAC3d) — the rcad
    // (curve, range) pair carrier.
    let (hc3d_first, hc3d_last) = (f, l);

    let mut res_tol = 0.0f64;

    // if (CheckSameParameter(HC3d, Pcurv, S, tol3d, tolreached))
    if check_same_parameter(&c3d, hc3d_first, hc3d_last, pcurv, surf, tol3d, tolreached) {
        return true;
    }

    if !has_pcurves(brep, e) {
        // occ::handle<Geom2dAdaptor_Curve> HC2d = new Geom2dAdaptor_Curve(Pcurv);
        let hc2d: rcad_kernel::base::proj_lib::adaptor::Curve2dHandle = Arc::new(
            rcad_kernel::base::proj_lib::adaptor::Geom2dCurveAdaptor::new(pcurv.clone()),
        );
        // occ::handle<GeomAdaptor_Surface> S = new (GeomAdaptor_Surface)(Surf);
        let s_handle: rcad_kernel::base::proj_lib::adaptor::SurfaceHandle = Arc::new(
            rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomSurfaceAdaptor::new(
                surf.clone(),
            ),
        );
        // Approx_CurveOnSurface AppCurve(HC2d, S, HC2d->FirstParameter(),
        //                                HC2d->LastParameter(),
        //                                Precision::Confusion());
        let (pc_first, pc_last) = curve2d_param_bounds(pcurv);
        let mut app_curve = crate::geomalgo::approx_curve_on_surface::ApproxCurveOnSurface::new(
            hc2d,
            s_handle,
            pc_first,
            pc_last,
            CONFUSION,
        );
        // AppCurve.Perform(10, 10, GeomAbs_C1, true); (Only2d = false)
        app_curve.perform(10, 10, rcad_kernel::math::GeomAbsShape::C1, true, false);
        if app_curve.is_done() && app_curve.has_result() {
            // C3d = AppCurve.Curve3d();
            c3d = Curve3::BSpline(app_curve.curve3d().expect("AppCurve.Curve3d"));
            // tolreached = AppCurve.MaxError3d();
            *tolreached = app_curve.max_error3d();
            // B.UpdateEdge(E, C3d, tolreached);
            let ed = brep.edge_mut_inplace(e.clone());
            ed.curve = Some(c3d.clone());
            ed.tolerance = ed.tolerance.max(*tolreached);
            return true;
        }
    }

    // const occ::handle<Adaptor3d_Curve>& aHCurve = HC3d;
    // Approx_SameParameter sp(aHCurve, Pcurv, S, tol3d);
    let sp = crate::geomalgo::approx_same_parameter::ApproxSameParameter::new(
        &c3d,
        hc3d_first,
        hc3d_last,
        pcurv,
        surf,
        tol3d,
    );
    if sp.is_done() && !sp.is_same_parameter() {
        *pcurv = sp.curve2d();
    } else if !sp.is_done() && !sp.is_same_parameter() {
        // "echec SameParameter"
        return false;
    }

    // occ::handle<Adaptor3d_Curve> curve3d = sp.Curve3d();
    // occ::handle<Adaptor3d_CurveOnSurface> curveOnSurface = sp.CurveOnSurface();
    let curve3d = sp.curve3d();
    let (c3d_first, c3d_last) = curve3_domain(&curve3d);
    let curve_on_surface = sp.curve_on_surface();

    if !check_same_parameter_exact(&curve3d, c3d_first, c3d_last, &curve_on_surface, tol3d, &mut res_tol)
        && res_tol > *tolreached
    {
        // "SameParameter : Tolerance not reached!"
        return false;
    } else {
        *tolreached = 1.1 * res_tol;
        if sp.is_done() && !sp.is_same_parameter() {
            *pcurv = sp.curve2d();
        }
    }
    true
}

// ---------------------------------------------------------------------------
// CorrectSameParameter (cxx L394-449)
// ---------------------------------------------------------------------------

/// OCCT static CorrectSameParameter (cxx L394-449).
pub(super) fn correct_same_parameter(brep: &mut BRep, the_edge: &Shape, the_face1: &Shape, the_face2: &Shape) {
    // if (BRep_Tool::Degenerated(theEdge)) return;
    if brep.edge(the_edge.clone()).degenerated {
        return;
    }

    // occ::handle<Geom_Curve> aCurve = BRep_Tool::Curve(theEdge, fpar, lpar);
    let ed = brep.edge(the_edge.clone());
    let a_curve = ed
        .curve
        .clone()
        .expect("CorrectSameParameter: edge without 3d curve");
    let (fpar, lpar) = (ed.range[0], ed.range[1]);

    // bool PCurveExists[2] = {false, false}; BRepAdaptor_Curve BAcurve[2];
    let mut pcurve_exists = [false, false];
    let mut ba_curve: [Option<(Curve2d, Surface3)>; 2] = [None, None];

    if !the_face1.is_null() {
        pcurve_exists[0] = true;
        // BAcurve[0].Initialize(theEdge, theFace1);
        ba_curve[0] = brep_adaptor_curve_on_surface(brep, the_edge, the_face1);
    }
    if !the_face1.is_null() && the_face1.is_same(the_face2) {
        // theEdge.Reverse();
        let mut e = the_edge.clone();
        e.orientation = match e.orientation {
            Orientation::Forward => Orientation::Reversed,
            Orientation::Reversed => Orientation::Forward,
            o => o,
        };
        let _ = e;
    }
    if !the_face2.is_null() {
        pcurve_exists[1] = true;
        // BAcurve[1].Initialize(theEdge, theFace2);
        ba_curve[1] = brep_adaptor_curve_on_surface(brep, the_edge, the_face2);
    }

    let mut max_sq_dist = 0.0f64;
    let ncontrol = 23i32;
    let delta = (lpar - fpar) / ncontrol as f64;

    for i in 0..=ncontrol {
        let a_param = fpar + i as f64 * delta;
        let a_pnt = a_curve.point_at(a_param);
        for j in 0..2 {
            if pcurve_exists[j] {
                // gp_Pnt aPntFromFace = BAcurve[j].Value(aParam);
                let (pc, surf) = ba_curve[j].as_ref().expect("BAcurve");
                let a_pnt_from_face = {
                    let uv = pc.point_at(a_param);
                    surf.point_at(uv.x, uv.y)
                };
                let a_sq_dist = a_pnt.distance_squared(a_pnt_from_face);
                if a_sq_dist > max_sq_dist {
                    max_sq_dist = a_sq_dist;
                }
            }
        }
    }

    let a_tol = max_sq_dist.sqrt();
    // BB.UpdateEdge(theEdge, aTol);
    BRepBuilder::new().update_edge_tolerance(brep, the_edge.clone(), a_tol);
}

/// OCCT BRepAdaptor_Curve(E, F) — the (pcurve, surface) carrier read.
fn brep_adaptor_curve_on_surface(brep: &BRep, e: &Shape, f: &Shape) -> Option<(Curve2d, Surface3)> {
    let (pc, _tf, _tl) = crate::brep_algo::tool::brep_tool_curve_on_surface(e, f)?;
    let surf = rcad_kernel::topo::topods::face_surface_value(brep, f)?;
    Some((pc, surf))
}

// ---------------------------------------------------------------------------
// Oriente (cxx L455-491)
// ---------------------------------------------------------------------------

/// OCCT static Oriente (cxx L455-491) — orientate an edge of natural
/// restriction.
pub(super) fn oriente(brep: &BRep, s: &Surface3, e: &mut Shape) {
    let u_ref = DVec2::new(1.0, 0.0);
    let v_ref = DVec2::new(0.0, 1.0);
    // S->Bounds(UFirst, ULast, VFirst, VLast);
    let (_basis, bounds) = surface_adaptor_basis_and_bounds(s);
    let (u_first, _u_last, v_first, _v_last) = (bounds[0], bounds[1], bounds[2], bounds[3]);

    // C = BRep_Tool::CurveOnSurface(E, S, bid, f, l);
    let (c, f, l) =
        curve_on_surface_by_surf(brep, e, s, 0).expect("Oriente: no pcurve on surface");
    // C->D1((f + l) / 2, P, D);
    let mid = (f + l) / 2.0;
    let p = c.point_at(mid);
    let d = {
        // the OCCT D1 first derivative via the rcad central difference
        let h = 1.0e-6;
        (c.point_at(mid + h) - c.point_at(mid - h)) / (2.0 * h)
    };

    let is_uiso = gp_vec2d_is_parallel(d, v_ref, 0.1);

    let is_first;
    let is_opposite;
    if is_uiso {
        is_first = (p.x - u_first).abs() < CONFUSION;
        is_opposite = gp_vec2d_is_opposite(d, v_ref, 0.1);
        e.orientation = Orientation::Reversed;
    } else {
        is_first = (p.y - v_first).abs() < CONFUSION;
        is_opposite = gp_vec2d_is_opposite(d, u_ref, 0.1);
        e.orientation = Orientation::Forward;
    }

    if !is_first {
        // E.Reverse();
        e.orientation = match e.orientation {
            Orientation::Forward => Orientation::Reversed,
            Orientation::Reversed => Orientation::Forward,
            o => o,
        };
    }
    if is_opposite {
        // E.Reverse();
        e.orientation = match e.orientation {
            Orientation::Forward => Orientation::Reversed,
            Orientation::Reversed => Orientation::Forward,
            o => o,
        };
    }
}

// ---------------------------------------------------------------------------
// UpdateEdgeOnPlane (cxx L494-517)
// ---------------------------------------------------------------------------

/// OCCT static UpdateEdgeOnPlane (cxx L494-517) — the OCC500(apo) block.
pub(super) fn update_edge_on_plane(brep: &mut BRep, f: &Shape, e: &Shape, _bb: &mut BRepBuilder) {
    // occ::handle<Geom2d_Curve> C2d = BRep_Tool::CurveOnSurface(E, F, f, l);
    let (c2d, _tf, _tl) = crate::brep_algo::tool::brep_tool_curve_on_surface(e, f)
        .expect("UpdateEdgeOnPlane: no pcurve on face");
    // occ::handle<Geom_Surface> S = BRep_Tool::Surface(F);
    let s = rcad_kernel::topo::topods::face_surface_value(brep, f)
        .expect("UpdateEdgeOnPlane: face without surface");
    // double Tol = BRep_Tool::Tolerance(E);
    let mut tol = brep.edge(e.clone()).tolerance;
    // BB.UpdateEdge(E, C2d, S, Loc, Tol);
    update_edge_pcurve_on_surf(brep, e, &c2d, &s, 0, tol);
    // BRepCheck_Edge Check(E); Tol = std::max(Tol, Check.Tolerance());
    tol = tol.max(brep_check_edge_tolerance(brep, e, &c2d, &s));
    // BB.UpdateEdge(E, Tol);
    BRepBuilder::new().update_edge_tolerance(brep, e.clone(), tol);
    // TopoDS_Vertex V; Tol *= 1.01;
    tol *= 1.01;
    // V = TopExp::FirstVertex(E);
    let (vf, vl) = crate::brep_fill::brep_fill_pipe::top_exp_vertices(e);
    let v = vf;
    if brep.vertex(v.clone()).tolerance < tol {
        // BB.UpdateVertex(V, Tol);
        BRepBuilder::new().update_vertex_tolerance(brep, v.clone(), tol);
    }
    // V = TopExp::LastVertex(E);
    let v = vl;
    if brep.vertex(v.clone()).tolerance < tol {
        BRepBuilder::new().update_vertex_tolerance(brep, v.clone(), tol);
    }
}

/// OCCT BRepCheck_Edge::Tolerance (BRepCheck_Edge.cxx) — the max deviation
/// between the 3d curve and the pcurve representation; carried through the
/// kernel sampling carrier of GeomLib_CheckCurveOnSurface (plan section 0.6).
fn brep_check_edge_tolerance(brep: &BRep, e: &Shape, pc: &Curve2d, s: &Surface3) -> f64 {
    let _ = brep;
    let ed = brep.edge(e.clone());
    if let Some(c3d) = &ed.curve {
        let (tf, tl) = curve2d_param_bounds(pc);
        let mut check = rcad_kernel::base::geom_lib::CheckCurveOnSurface::with_curve(c3d, 0.0);
        check.perform(c3d, pc, s);
        if check.is_done() {
            let _ = (tf, tl);
            return check.max_distance();
        }
    }
    0.0
}

// ---------------------------------------------------------------------------
// BuildFace (cxx L527-762)
// ---------------------------------------------------------------------------

/// OCCT static BuildFace (cxx L527-762) — construct a Face via a surface and
/// 4 Edges (natural borders).
#[allow(clippy::too_many_arguments)]
pub(super) fn build_face(
    brep: &mut BRep,
    s: &Surface3,
    e1: &Shape,
    e2: &Shape,
    e3: &Shape,
    e4: &Shape,
    eemap: &mut HashMap<ShapeKey, Shape>,
    exch_uv: bool,
    u_reverse: bool,
    f: &mut Shape,
) {
    let mut b = BRepBuilder::new();

    // Is the surface planar ?
    // Tol1..Tol4 = BRep_Tool::Tolerance(E1..E4);
    let tol1 = brep.edge(e1.clone()).tolerance;
    let tol2 = brep.edge(e2.clone()).tolerance;
    let tol3 = brep.edge(e3.clone()).tolerance;
    let tol4 = brep.edge(e4.clone()).tolerance;
    let mut tol = tol1.min(tol2);
    tol = tol.min(tol3.min(tol4));
    let mut is_plan = false;
    let mut the_plane: Option<Surface3> = None;

    if !e1.is_same(e3) && !e2.is_same(e4) {
        // exclude cases with seam edges: they are not planar
        // GeomLib_IsPlanarSurface IsP(S, Tol);
        let is_p = rcad_kernel::base::geom_lib::IsPlanarSurface::new(s, tol);
        if is_p.is_planar() {
            is_plan = true;
            the_plane = Some(Surface3::Plane(is_p.plan().clone()));
        } else {
            // TE1..TE4->Tolerance(Precision::Confusion());
            set_edge_tolerance(brep, e1, CONFUSION);
            set_edge_tolerance(brep, e2, CONFUSION);
            set_edge_tolerance(brep, e3, CONFUSION);
            set_edge_tolerance(brep, e4, CONFUSION);

            // TopoDS_Wire theWire = BRepLib_MakeWire(E1, E2, E3, E4);
            let the_wire = b.build_wire(
                brep,
                vec![e1.clone(), e2.clone(), e3.clone(), e4.clone()],
            );
            // int NbPoints = NumberOfPoles(theWire);
            let nb_points = number_of_poles(brep, &the_wire);
            if nb_points <= 100 {
                // limitation for CPU
                // BRepLib_FindSurface FS(theWire, -1, true);
                let fs = crate::topalgo::brep_lib_find_surface::BRepLibFindSurface::new(
                    brep, &the_wire, -1.0, true,
                );
                if fs.found() {
                    is_plan = true;
                    the_plane = fs.surface();
                }
            }
            // BB.UpdateEdge(E1, Tol1); ... (restore the tolerances)
            BRepBuilder::new().update_edge_tolerance(brep, e1.clone(), tol1);
            BRepBuilder::new().update_edge_tolerance(brep, e2.clone(), tol2);
            BRepBuilder::new().update_edge_tolerance(brep, e3.clone(), tol3);
            BRepBuilder::new().update_edge_tolerance(brep, e4.clone(), tol4);
        }
    }

    // Construction of the wire
    // e1 = E1; Oriente(S, e1);
    let mut e1l = e1.clone();
    oriente(brep, s, &mut e1l);
    let mut ww = b.make_wire(brep);
    if !is_plan || !brep.edge(e1l.clone()).degenerated {
        // B.Add(e1);
        b.add_to_wire(brep, ww.clone(), e1l.clone());
    }

    let mut e2l = e2.clone();
    oriente(brep, s, &mut e2l);
    if !is_plan || !brep.edge(e2l.clone()).degenerated {
        b.add_to_wire(brep, ww.clone(), e2l.clone());
        if !brep.edge(e2l.clone()).degenerated {
            // WW = B.Wire();
            ww = wire_of(brep, &ww);
            let mut new_edge = Shape::null();
            // take the last edge added to WW
            for iter in topods_iterator(&ww) {
                new_edge = iter;
            }
            if !e2l.is_same(&new_edge) {
                // EEmap.Bind(e2, NewEdge);
                eemap.insert(ShapeKey(e2l.ptr_id()), new_edge);
            }
        }
    }

    let e: Shape;
    if e3.is_same(e1) {
        // E = e1; E.Reverse();
        let mut ee = e1l.clone();
        ee.orientation = match ee.orientation {
            Orientation::Forward => Orientation::Reversed,
            Orientation::Reversed => Orientation::Forward,
            o => o,
        };
        e = ee;
    } else {
        let mut ee = e3.clone();
        oriente(brep, s, &mut ee);
        e = ee;
    }
    if !is_plan || !brep.edge(e.clone()).degenerated {
        b.add_to_wire(brep, ww.clone(), e.clone());
        if !brep.edge(e.clone()).degenerated {
            ww = wire_of(brep, &ww);
            let mut new_edge = Shape::null();
            for iter in topods_iterator(&ww) {
                new_edge = iter;
            }
            if !e.is_same(&new_edge) {
                eemap.insert(ShapeKey(e.ptr_id()), new_edge);
            }
        }
    }

    let e: Shape;
    if e4.is_same(e2) {
        let mut ee = e2l.clone();
        ee.orientation = match ee.orientation {
            Orientation::Forward => Orientation::Reversed,
            Orientation::Reversed => Orientation::Forward,
            o => o,
        };
        e = ee;
    } else {
        let mut ee = e4.clone();
        oriente(brep, s, &mut ee);
        e = ee;
    }
    if !is_plan || !brep.edge(e.clone()).degenerated {
        b.add_to_wire(brep, ww.clone(), e.clone());
        if !brep.edge(e.clone()).degenerated {
            ww = wire_of(brep, &ww);
            let mut new_edge = Shape::null();
            for iter in topods_iterator(&ww) {
                new_edge = iter;
            }
            if !e.is_same(&new_edge) {
                eemap.insert(ShapeKey(e.ptr_id()), new_edge);
            }
        }
    }

    // WW = B.Wire();
    ww = wire_of(brep, &ww);

    // Construction of the face.
    if is_plan {
        // Suspend representation 2d and construct face Plane
        let the_plane = the_plane.expect("BuildFace: planar without plane");
        let mut the_plane = the_plane;
        // gp_Pnt aPnt; gp_Vec DU, DV, NS, NP;
        let (_basis, bounds) = surface_adaptor_basis_and_bounds(s);
        let (ufirst, ulast, vfirst, vlast) = (bounds[0], bounds[1], bounds[2], bounds[3]);
        let (p, du, dv) = s.derivatives((ufirst + ulast) / 2.0, (vfirst + vlast) / 2.0);
        let _ = p;
        // NS = DU ^ DV;
        let ns = du.cross(dv);
        // NP = thePlane->Pln().Axis().Direction();
        let np = match &the_plane {
            Surface3::Plane(pl) => pl.normal,
            _ => unreachable!("BuildFace: thePlane"),
        };
        if ns.dot(np) < 0.0 {
            // thePlane->UReverse(); (gp_Ax3::UReverse — the X direction flips)
            if let Surface3::Plane(pl) = &mut the_plane {
                pl.u_dir = -pl.u_dir;
            }
        }
        // BRepLib_MakeFace MkF(thePlane, WW);
        let mkf_face = b.make_face(brep, Some(the_plane.clone()), ww.clone());
        // if (MkF.Error() != BRepLib_FaceDone) { debug print } else { ... }
        // (the OCCT error arm is OCCT_DEBUG-only; the rcad make_face is
        // infallible at this layer)
        {
            // occ::handle<Geom2d_Curve> NullC2d; TopLoc_Location Loc;
            // BB.UpdateEdge(E1, NullC2d, S, Loc, Tol1); ... (E2/E3/E4)
            for (ee, t) in [(e1, tol1), (e2, tol2), (e3, tol3), (e4, tol4)] {
                remove_pcurves_on_surf(brep, ee, s);
                BRepBuilder::new().update_edge_tolerance(brep, ee.clone(), t);
            }

            // F = MkF.Face();
            *f = mkf_face;
            // UpdateEdgeOnPlane(F, E1, BB); ... (E2/E3/E4)
            let mut bb = BRepBuilder::new();
            update_edge_on_plane(brep, f, e1, &mut bb);
            update_edge_on_plane(brep, f, e2, &mut bb);
            update_edge_on_plane(brep, f, e3, &mut bb);
            update_edge_on_plane(brep, f, e4, &mut bb);
        }
    }

    if !is_plan {
        // Cas Standard : Ajout
        // BB.MakeFace(F, S, Precision::Confusion()); BB.Add(F, WW);
        *f = b.make_face(brep, Some(s.clone()), ww.clone());
    }

    // Reorientation
    if exch_uv {
        // F.Reverse();
        f.orientation = match f.orientation {
            Orientation::Forward => Orientation::Reversed,
            Orientation::Reversed => Orientation::Forward,
            o => o,
        };
    }
    if u_reverse {
        f.orientation = match f.orientation {
            Orientation::Forward => Orientation::Reversed,
            Orientation::Reversed => Orientation::Forward,
            o => o,
        };
    }
}

/// OCCT BRepBuilderAPI_MakeWire::Wire() — the wire under construction (the
/// pool wire carries its edges in TWireData).
fn wire_of(_brep: &BRep, w: &Shape) -> Shape {
    w.clone()
}

/// OCCT BRep_Builder::UpdateEdge(E, NullC2d, S, Loc, Tol) — the pcurve
/// representation on (S, Loc) is dropped.
fn remove_pcurves_on_surf(brep: &mut BRep, e: &Shape, s: &Surface3) {
    let f = intern_surface_face(brep, s);
    let key = (f.ptr_id(), 0u32);
    let ed = brep.edge_mut_inplace(e.clone());
    ed.pcurves.shift_remove(&key);
    ed.representations
        .retain(|cr| match cr {
            rcad_kernel::topo::topods::CurveRepresentation::CurveOnSurface { face, .. } => {
                *face != key
            }
            rcad_kernel::topo::topods::CurveRepresentation::CurveOnClosedSurface {
                face,
                ..
            } => *face != key,
            _ => true,
        });
}

// ---------------------------------------------------------------------------
// BuildEdge — the 3d-curve form (cxx L766-861)
// ---------------------------------------------------------------------------

/// OCCT static BuildEdge (cxx L766-861) — the explicit 3d + 2d form.
pub(super) fn build_edge_c3d(
    brep: &mut BRep,
    c3d: &Option<Curve3>,
    c2d: &Curve2d,
    s: &Surface3,
    vf: &Shape,
    vl: &Shape,
    f: f64,
    l: f64,
    tol3d: f64,
) -> Shape {
    let mut b = BRepBuilder::new();
    let mut e: Shape;

    let p1 = brep.vertex(vf.clone()).point;
    let tol1 = brep.vertex(vf.clone()).tolerance;
    let p2 = brep.vertex(vl.clone()).point;
    let tol2 = brep.vertex(vl.clone()).tolerance;
    let mut tol = tol1.max(tol2);

    if vf.is_same(vl) || p1.distance(p2) < tol {
        // Degenerated case
        let mut d;
        let p2d = c2d.point_at(f);
        let p = s.point_at(p2d.x, p2d.y);
        d = p1.distance(p);
        if d > tol {
            tol = d;
        }
        let p2d = c2d.point_at(l);
        let p = s.point_at(p2d.x, p2d.y);
        d = p2.distance(p);
        if d > tol {
            tol = d;
        }

        b.update_vertex_tolerance(brep, vf.clone(), tol);
        b.update_vertex_tolerance(brep, vl.clone(), tol);

        // B.MakeEdge(E); B.UpdateEdge(E, C2d, S, TopLoc_Location(), Tol);
        e = b.add_edge(brep, None, vf.clone(), vl.clone(), [f, l]);
        update_edge_pcurve_on_surf(brep, &e, c2d, s, 0, tol);
        // B.Add(E, VF); B.Add(E, VL);
        b.add_to_edge(brep, e.clone(), vf.clone());
        b.add_to_edge(brep, e.clone(), vl.clone());
        // B.Range(E, f, l); B.Degenerated(E, true);
        b.set_edge_range(brep, e.clone(), f, l);
        b.set_edge_degenerated(brep, e.clone(), true);

        return e;
    }

    if let Some(c3d) = c3d {
        // C3d->D0(f, P); d = P1.Distance(P);
        let p = c3d.point_at(f);
        let d = p1.distance(p);
        if d > tol1 {
            b.update_vertex_tolerance(brep, vf.clone(), d);
        }

        let p = c3d.point_at(l);
        let d = p2.distance(p);
        if d > tol2 {
            b.update_vertex_tolerance(brep, vl.clone(), d);
        }

        // BRepLib_MakeEdge MkE(C3d, VF, VL, f, l);
        e = b.add_edge(brep, Some(c3d.clone()), vf.clone(), vl.clone(), [f, l]);
        // if (!MkE.IsDone()) { } — the OCCT error arm is empty (the rcad
        // add_edge is infallible at this layer)
    } else {
        panic!("BuildEdge: null C3d (the OCCT null handle dereference)");
    }
    // TopLoc_Location Loc; B.UpdateEdge(E, C2d, S, Loc, Tol3d);
    update_edge_pcurve_on_surf(brep, &e, c2d, s, 0, tol3d);

    // if (BOPTools_AlgoTools::IsMicroEdge(E, aNullCtx)) — GAP carrier (the
    // bop/ module owns BOPTools).
    let is_micro = is_micro_edge(brep, &e);
    if is_micro {
        let a_v = vf.clone();
        b.update_vertex_tolerance(brep, a_v.clone(), p1.distance(p2));
        // B.MakeEdge(E); B.UpdateEdge(E, C2d, S, TopLoc_Location(), Tol);
        e = b.add_edge(brep, None, a_v.clone(), a_v.clone(), [f, l]);
        update_edge_pcurve_on_surf(brep, &e, c2d, s, 0, tol);
        // B.Add(E, TopoDS::Vertex(aV.Oriented(TopAbs_FORWARD)));
        let mut vfwd = a_v.clone();
        vfwd.orientation = Orientation::Forward;
        b.add_to_edge(brep, e.clone(), vfwd);
        // B.Add(E, TopoDS::Vertex(aV.Oriented(TopAbs_REVERSED)));
        let mut vrev = a_v.clone();
        vrev.orientation = Orientation::Reversed;
        b.add_to_edge(brep, e.clone(), vrev);
        b.set_edge_range(brep, e.clone(), f, l);
        b.set_edge_degenerated(brep, e.clone(), true);
    }

    e
}

/// OCCT BOPTools_AlgoTools::IsMicroEdge (BOPTools_AlgoTools.cxx) — GAP
/// carrier (the bop/ module owns BOPTools; plan section 0.6).
fn is_micro_edge(brep: &BRep, e: &Shape) -> bool {
    let _ = (brep, e);
    panic!(
        "GAP: BOPTools_AlgoTools::IsMicroEdge (TKBool/BOPTools) is not \
         translated — see BRepFill_Sweep.cxx L847-858 (plan section 0.6)"
    );
}

// ---------------------------------------------------------------------------
// Substitute (cxx L1341-1368)
// ---------------------------------------------------------------------------

/// OCCT static Substitute (cxx L1341-1368).
pub(super) fn substitute(
    _brep: &mut BRep,
    a_substitute: &mut crate::topalgo::brep_tools_substitution::BRepToolsSubstitution,
    old: &Shape,
    new: &Shape,
) {
    let mut list_shape: Vec<Shape> = Vec::new();

    // TopExp::Vertices(Old, OldV1, OldV2); TopExp::Vertices(New, NewV1, NewV2);
    let (old_v1, old_v2) = crate::brep_fill::brep_fill_pipe::top_exp_vertices(old);
    let (new_v1, new_v2) = crate::brep_fill::brep_fill_pipe::top_exp_vertices(new);

    if !a_substitute.is_copied(&old_v1) {
        let mut fwd = new_v1.clone();
        fwd.orientation = Orientation::Forward;
        list_shape.push(fwd);
        a_substitute.substitute(&old_v1, list_shape.clone());
        list_shape.clear();
    }
    if !a_substitute.is_copied(&old_v2) {
        let mut fwd = new_v2.clone();
        fwd.orientation = Orientation::Forward;
        list_shape.push(fwd);
        a_substitute.substitute(&old_v2, list_shape.clone());
        list_shape.clear();
    }
    if !a_substitute.is_copied(old) {
        let mut fwd = new.clone();
        fwd.orientation = Orientation::Forward;
        list_shape.push(fwd);
        a_substitute.substitute(old, list_shape);
    }
}

// ---------------------------------------------------------------------------
// KeepEdge (cxx L1413-1436)
// ---------------------------------------------------------------------------

/// OCCT static KeepEdge (cxx L1413-1436) — find edges of the face supported
/// by the same Curve.
pub(super) fn keep_edge(brep: &BRep, face: &Shape, edge: &Shape) -> Vec<Shape> {
    let mut list: Vec<Shape> = Vec::new();
    let _f = 0.0;
    let _l = 0.0;
    // TopExp_Explorer Exp(Face, TopAbs_EDGE);
    let exp = crate::brep_algo::tool::explorer(
        face,
        rcad_kernel::topo::topods::ShapeType::Edge,
        rcad_kernel::topo::topods::ShapeType::Compound,
    );
    // Cref = BRep_Tool::Curve(TopoDS::Edge(Edge), Lref, f, l);
    let ed_ref = brep.edge(edge.clone());
    let cref = ed_ref.curve.clone();
    let lref = edge.location;

    for current in exp {
        // C = BRep_Tool::Curve(TopoDS::Edge(Exp.Current()), L, f, l);
        let ed_cur = brep.edge(current.clone());
        let c = ed_cur.curve.clone();
        let loc = current.location;
        // if ((Cref == C) && (Lref == L))
        if geom_curve_handle_same(&cref, &c) && lref == loc {
            // List.Append(Exp.Current());
            list.push(current);
        }
    }
    list
}

/// OCCT `Cref == C` — the Geom_Curve handle identity.  The rcad curve
/// travels as a value, so the identity stand-in is the edge TShape identity
/// plus a structural comparison of the curve payloads.
fn geom_curve_handle_same(a: &Option<Curve3>, b: &Option<Curve3>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(ca), Some(cb)) => {
            std::mem::discriminant(ca) == std::mem::discriminant(cb)
                && format!("{:?}", ca) == format!("{:?}", cb)
        }
        _ => false,
    }
}

// ---------------------------------------------------------------------------
// BuildEdge — the iso form (cxx L1476-1655)
// ---------------------------------------------------------------------------

/// OCCT static BuildEdge (cxx L1476-1655) — the iso form.
pub(super) fn build_edge_iso(
    brep: &mut BRep,
    s: &Surface3,
    is_uiso: bool,
    val_iso: f64,
    v_first: &Shape,
    v_last: &Shape,
    tol: f64,
) -> Shape {
    let mut b = BRepBuilder::new();
    let mut e: Shape;
    // occ::handle<Geom_Curve> Iso;
    let iso = if is_uiso {
        crate::brep_fill::brep_fill_sweep::surface_uiso(s, val_iso)
    } else {
        crate::brep_fill::brep_fill_sweep::surface_viso(s, val_iso)
    };
    let mut sing = false;

    let (iso_first, iso_last) = curve3_domain(&iso);

    if v_first.is_same(v_last) {
        // Singular case ?
        let v = v_first;
        let mut t_tol = brep.vertex(v.clone()).tolerance;
        if tol > t_tol {
            t_tol = tol;
        }
        let p = iso.point_at((iso_first + iso_last) / 2.0);
        if p.distance(brep.vertex(v.clone()).point) < t_tol {
            // GeomAdaptor_Curve AC(Iso);
            // sing = GCPnts_AbscissaPoint::Length(AC, tol / 4) < tol;
            let length = arc_length(&iso, iso_first, iso_last);
            sing = length < t_tol;
        }
    }

    if sing {
        // Singular case
        let mut v = v_first.clone();
        e = null_edge(brep, &mut v);
        // B.Degenerated(E, true);
        b.set_edge_degenerated(brep, e.clone(), true);
    } else {
        // Construction Via 3d
        let mut fwd = false;
        let p1 = iso.point_at(iso_first);
        let p2 = iso.point_at(iso_last);

        let t1 = brep.vertex(v_first.clone()).tolerance;
        let t2 = brep.vertex(v_last.clone()).tolerance;

        let mut bb = BRepBuilder::new();

        let p11 = p1.distance(brep.vertex(v_first.clone()).point);
        let p22 = p2.distance(brep.vertex(v_last.clone()).point);
        let p12 = p1.distance(brep.vertex(v_last.clone()).point);
        let p21 = p2.distance(brep.vertex(v_first.clone()).point);

        if p11 < p12 && p22 < p21 {
            fwd = true;
        }

        if fwd {
            // OCC500(apo)
            if p11 >= t1 {
                bb.update_vertex_tolerance(brep, v_first.clone(), 1.01 * p11);
            }
            if p22 >= t2 {
                bb.update_vertex_tolerance(brep, v_last.clone(), 1.01 * p22);
            }
        } else {
            if p12 >= t2 {
                bb.update_vertex_tolerance(brep, v_last.clone(), 1.01 * p12);
            }
            if p21 >= t1 {
                bb.update_vertex_tolerance(brep, v_first.clone(), 1.01 * p21);
            }
        }

        // BRepLib_MakeEdge MkE;
        let mut mk_e = BRepBuilder::new();
        // if (fwd) MkE.Init(Iso, VFirst, VLast, First, Last);
        let mk_e_edge = if fwd {
            mk_e.add_edge(
                brep,
                Some(iso.clone()),
                v_first.clone(),
                v_last.clone(),
                [iso_first, iso_last],
            )
        } else {
            // MkE.Init(Iso, VLast, VFirst, First, Last);
            mk_e.add_edge(
                brep,
                Some(iso.clone()),
                v_last.clone(),
                v_first.clone(),
                [iso_first, iso_last],
            )
        };

        // if (!MkE.IsDone()) throw Standard_ConstructionError — the rcad
        // add_edge is infallible at this layer.
        e = mk_e_edge;
    }

    // Associate 2d
    // occ::handle<Geom2d_Line> L; TopLoc_Location Loc;
    let (_u_min, _u_max, v_min, _v_max) = {
        let (_basis, bounds) = surface_adaptor_basis_and_bounds(s);
        (bounds[0], bounds[1], bounds[2], bounds[3])
    };
    let (u_min, _u_max, _v_min, _v_max) = {
        let (_basis, bounds) = surface_adaptor_basis_and_bounds(s);
        (bounds[0], bounds[1], bounds[2], bounds[3])
    };
    let line2d: Curve2d = if is_uiso {
        // gp_Pnt2d P(ValIso, Vmin - Iso->FirstParameter()); gp_Vec2d V(0., 1.);
        Curve2d::Line(Line2d::new(
            DVec2::new(val_iso, v_min - iso_first),
            DVec2::new(0.0, 1.0),
        ))
    } else {
        // gp_Pnt2d P(Umin - Iso->FirstParameter(), ValIso); gp_Vec2d V(1., 0.);
        Curve2d::Line(Line2d::new(
            DVec2::new(u_min - iso_first, val_iso),
            DVec2::new(1.0, 0.0),
        ))
    };

    // B.UpdateEdge(E, L, S, Loc, Precision::Confusion());
    update_edge_pcurve_on_surf(brep, &e, &line2d, s, 0, CONFUSION);
    if sing {
        // B.Range(E, S, Loc, Iso->FirstParameter(), Iso->LastParameter());
        set_pcurve_range_on_surf(brep, &e, s, 0, iso_first, iso_last);
    }

    // double MaxTol = 1.e-4; double theTol;
    let max_tol = 1.0e-4;
    let mut the_tol = 0.0f64;
    // GeomAdaptor_Curve GAiso(Iso); GAHiso = new(GAiso);
    // GeomAdaptor_Surface GAsurf(S); GAHsurf = new(GAsurf);
    // CheckSameParameter(GAHiso, L, GAHsurf, MaxTol, theTol);
    check_same_parameter(&iso, iso_first, iso_last, &line2d, s, max_tol, &mut the_tol);
    // B.UpdateEdge(E, theTol);
    BRepBuilder::new().update_edge_tolerance(brep, e.clone(), the_tol);

    e
}

/// OCCT BRep_Builder::Range(E, S, L, First, Last) — the pcurve range
/// restriction (BRep_Builder.cxx L560-600).
fn set_pcurve_range_on_surf(brep: &mut BRep, e: &Shape, s: &Surface3, loc: u32, first: f64, last: f64) {
    let f = intern_surface_face(brep, s);
    let key = (f.ptr_id(), loc);
    let ed = brep.edge_mut_inplace(e.clone());
    if let Some(entry) = ed.pcurves.get_mut(&key) {
        entry.1 = first;
        entry.2 = last;
    }
}

// ---------------------------------------------------------------------------
// UpdateEdge (cxx L1659-1820)
// ---------------------------------------------------------------------------

/// OCCT static UpdateEdge (cxx L1659-1820).
pub(super) fn update_edge_iso(brep: &mut BRep, e: &Shape, s: &Surface3, is_uiso: bool, val_iso: f64) {
    let mut b = BRepBuilder::new();
    // occ::handle<Geom2d_Line> L; occ::handle<Geom2d_Curve> PCurve, CL;
    let mut f2d;
    let mut l2d;
    let (_u_first, _u_last, _v_first, _v_last) = {
        let (_basis, bounds) = surface_adaptor_basis_and_bounds(s);
        (bounds[0], bounds[1], bounds[2], bounds[3])
    };
    let (u_first, u_last, v_first, v_last) = {
        let (_basis, bounds) = surface_adaptor_basis_and_bounds(s);
        (bounds[0], bounds[1], bounds[2], bounds[3])
    };

    // bool sing = false; occ::handle<Geom_Curve> Iso;
    let mut sing = false;
    let iso = if is_uiso {
        crate::brep_fill::brep_fill_sweep::surface_uiso(s, val_iso)
    } else {
        crate::brep_fill::brep_fill_sweep::surface_viso(s, val_iso)
    };
    let (iso_first, iso_last) = curve3_domain(&iso);

    // TopExp::Vertices(E, Vf, Vl);
    let (vf, vl) = crate::brep_fill::brep_fill_pipe::top_exp_vertices(e);
    if vf.is_same(&vl) {
        // Singular case ?
        // double tol = BRep_Tool::Tolerance(Vf);
        let t_tol = brep.vertex(vf.clone()).tolerance;
        let pmid = iso.point_at((iso_first + iso_last) / 2.0);
        if pmid.distance(brep.vertex(vf.clone()).point) < t_tol {
            let length = arc_length(&iso, iso_first, iso_last);
            sing = length < t_tol;
        }
    }

    let line: Curve2d = if is_uiso {
        // gp_Pnt2d P(ValIso, 0); gp_Vec2d V(0., 1.); F2d = VFirst; L2d = VLast;
        f2d = v_first;
        l2d = v_last;
        Curve2d::Line(Line2d::new(DVec2::new(val_iso, 0.0), DVec2::new(0.0, 1.0)))
    } else {
        // gp_Pnt2d P(0., ValIso); gp_Vec2d V(1., 0.); F2d = UFirst; L2d = ULast;
        f2d = u_first;
        l2d = u_last;
        Curve2d::Line(Line2d::new(DVec2::new(0.0, val_iso), DVec2::new(1.0, 0.0)))
    };
    // CL = new (Geom2d_TrimmedCurve)(L, F2d, L2d);
    let mut cl = Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(line.clone()),
        t_min: f2d,
        t_max: l2d,
    });

    // Control direction & Range — double R, First, Last, Tol = 1.e-4;
    let mut first;
    let mut last;
    let tol = 1.0e-4;
    let mut r = 0.0f64;
    let mut reverse = false;

    // BRep_Tool::Range(E, First, Last);
    let ed = brep.edge(e.clone());
    first = ed.range[0];
    last = ed.range[1];

    if !vf.is_same(&vl) {
        // Test distances between "FirstPoint" and "Vertex"
        // P2d = CL->Value(F2d); POnS = S->Value(P2d.X(), P2d.Y());
        let p2d = cl.point_at(f2d);
        let p_on_s = s.point_at(p2d.x, p2d.y);
        // reverse = POnS.Distance(BRep_Tool::Pnt(Vl)) < POnS.Distance(BRep_Tool::Pnt(Vf));
        reverse = p_on_s.distance(brep.vertex(vl.clone()).point)
            < p_on_s.distance(brep.vertex(vf.clone()).point);
    } else if !sing {
        // Test angle between "First Tangente"
        // BRepAdaptor_Curve C3d(E); C3d.D1(First, POnS, dC3d);
        let ed = brep.edge(e.clone());
        let c3d = ed.curve.clone().expect("UpdateEdge: edge without 3d curve");
        let p_on_s = c3d.point_at(first);
        let d_c3d = c3d.derivative_at(first);
        // CL->D1(F2d, P2d, V2d);
        let p2d = cl.point_at(f2d);
        let v2d = {
            let h = 1.0e-6;
            (cl.point_at(f2d + h) - cl.point_at(f2d - h)) / (2.0 * h)
        };
        // S->D1(P2d.X(), P2d.Y(), POnS, du, dv);
        let (_p, du, dv) = s.derivatives(p2d.x, p2d.y);
        let _ = p_on_s;
        // V3d.SetLinearForm(V2d.X(), du, V2d.Y(), dv);
        let v3d = du * v2d.x + dv * v2d.y;
        // reverse = (dC3d.Angle(V3d) > Tol);
        reverse = gp_vec_angle(d_c3d, v3d) > tol;
    }
    if reverse {
        // Return curve 2d
        // CL = new (Geom2d_TrimmedCurve)(L, F2d, L2d); CL->Reverse();
        cl = reversed_trimmed_line(&line, f2d, l2d);
        let (f2d_r, l2d_r) = curve2d_param_bounds(&cl);
        f2d = f2d_r;
        l2d = l2d_r;
    }

    if sing {
        // occ::handle<Geom_Curve> NullCurve; B.UpdateEdge(E, NullCurve, 0.);
        brep.edge_mut_inplace(e.clone()).curve = None;
        // B.Degenerated(E, true);
        b.set_edge_degenerated(brep, e.clone(), true);
        // B.Range(E, F2d, L2d);
        b.set_edge_range(brep, e.clone(), f2d, l2d);
        first = f2d;
        last = l2d;
    }

    if first != f2d || last != l2d {
        // occ::handle<Geom2d_Curve> C2d;
        // GeomLib::SameRange(PConfusion(), CL, F2d, L2d, First, Last, C2d);
        let c2d = crate::geomalgo::geom_lib_same_range::same_range(
            p_confusion(),
            &cl,
            f2d,
            l2d,
            first,
            last,
        );
        // CL = new (Geom2d_TrimmedCurve)(C2d, First, Last);
        cl = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(c2d),
            t_min: first,
            t_max: last,
        });
    }

    // Update des Vertex
    // P2d = CL->Value(First); POnS = S->Value(P2d.X(), P2d.Y());
    let p2d = cl.point_at(first);
    let p_on_s = s.point_at(p2d.x, p2d.y);
    // V = TopExp::FirstVertex(E); R = POnS.Distance(BRep_Tool::Pnt(V));
    let (v_first_v, v_last_v) = crate::brep_fill::brep_fill_pipe::top_exp_vertices(e);
    let v = v_first_v;
    r = p_on_s.distance(brep.vertex(v.clone()).point);
    b.update_vertex_tolerance(brep, v.clone(), r);

    let p2d = cl.point_at(last);
    let p_on_s = s.point_at(p2d.x, p2d.y);
    let v = v_last_v;
    r = p_on_s.distance(brep.vertex(v.clone()).point);
    b.update_vertex_tolerance(brep, v.clone(), r);

    // Update Edge
    // if (!sing && SameParameter(E, CL, S, Tol, R)) B.UpdateEdge(E, R);
    if !sing && same_parameter(brep, e, &mut cl, s, tol, &mut r) {
        b.update_edge_tolerance(brep, e.clone(), r);
    }

    // PCurve = Couture(E, S, Loc);
    let pcurve = couture(brep, e, s, 0);
    if pcurve.is_none() {
        // B.UpdateEdge(E, CL, S, Loc, Precision::Confusion());
        update_edge_pcurve_on_surf(brep, e, &cl, s, 0, CONFUSION);
    } else {
        // Sewing edge
        let pcurve = pcurve.expect("PCurve");
        let mut ee = e.clone();
        oriente(brep, s, &mut ee);
        if ee.orientation == Orientation::Reversed {
            // B.UpdateEdge(E, CL, PCurve, S, Loc, Precision::Confusion());
            update_edge_pcurve2_on_surf(brep, e, &cl, &pcurve, s, 0, CONFUSION);
        } else {
            // B.UpdateEdge(E, PCurve, CL, S, Loc, Precision::Confusion());
            update_edge_pcurve2_on_surf(brep, e, &pcurve, &cl, s, 0, CONFUSION);
        }
    }

    // Attention to case not SameRange on its shapes (PRO13551)
    // if (!BRep_Tool::SameRange(E)) B.Range(E, S, Loc, First, Last);
    if !brep.edge(e.clone()).same_range {
        set_pcurve_range_on_surf(brep, e, s, 0, first, last);
    }
}

/// OCCT Geom2d_TrimmedCurve(L, F, L) + Reverse() — the basis line direction
/// is reversed and the bounds swapped (Geom2d_TrimmedCurve.cxx Reverse /
/// Geom2d_Line::Reversed).
fn reversed_trimmed_line(line: &Curve2d, f2d: f64, l2d: f64) -> Curve2d {
    let reversed_basis = match line {
        Curve2d::Line(l) => Curve2d::Line(Line2d {
            origin: l.origin,
            direction: -l.direction,
        }),
        _ => panic!(
            "GAP: Geom2d_Line::Reversed (the rcad re-host covers the line \
             basis) — BRepFill_Sweep::UpdateEdge"
        ),
    };
    Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(reversed_basis),
        t_min: l2d,
        t_max: f2d,
    })
}

// ---------------------------------------------------------------------------
// IsDegen (cxx L1825-1875)
// ---------------------------------------------------------------------------

/// OCCT static IsDegen (cxx L1825-1875) — check if a surface is degenerated.
pub(super) fn is_degen(s: &Surface3, tol: f64) -> bool {
    let nb = 5i32;
    let mut b;
    let (_basis, bounds) = surface_adaptor_basis_and_bounds(s);
    let (u_min, u_max, v_min, v_max) = (bounds[0], bounds[1], bounds[2], bounds[3]);
    let mut t;
    let mut dt;
    let mut l;

    // Check the length of Iso-U
    t = (u_min + u_max) / 2.0;
    let p1 = s.point_at(t, v_min);
    let p2 = s.point_at(t, (v_min + v_max) / 2.0);
    let p3 = s.point_at(t, v_max);
    b = (p1.distance(p2) + p2.distance(p3)) < tol;

    dt = (u_max - u_min) / (nb as f64 + 1.0);
    let mut ii = 1i32;
    while b && ii <= nb {
        t = u_min + ii as f64 * dt;
        let iso = crate::brep_fill::brep_fill_sweep::surface_uiso(s, t);
        let (iso_first, iso_last) = curve3_domain(&iso);
        // GeomAdaptor_Curve AC(Iso);
        // l = GCPnts_AbscissaPoint::Length(AC, Tol / 4);
        l = arc_length(&iso, iso_first, iso_last);
        b = l <= tol;
        ii += 1;
    }

    if b {
        return true;
    }

    // Check the length of Iso-V
    t = (v_min + v_max) / 2.0;
    let p1 = s.point_at(u_min, t);
    let p2 = s.point_at((u_min + u_max) / 2.0, t);
    let p3 = s.point_at(u_max, t);
    b = (p1.distance(p2) + p2.distance(p3)) < tol;

    dt = (v_max - v_min) / (nb as f64 + 1.0);
    let mut ii = 1i32;
    while b && ii <= nb {
        t = v_min + ii as f64 * dt;
        let iso = crate::brep_fill::brep_fill_sweep::surface_viso(s, t);
        let (iso_first, iso_last) = curve3_domain(&iso);
        l = arc_length(&iso, iso_first, iso_last);
        b = l <= tol;
        ii += 1;
    }

    b
}


// Silence the unused warnings for the params carried for form fidelity.
#[allow(unused_variables)]
fn _form_guards(
    _h2: Option<ShapeHArray2>,
    _c: Option<Curve3>,
    _t: Option<TrimmedCurve3>,
    _l: Option<Line3>,
    _p: Option<Plane>,
) {
}
