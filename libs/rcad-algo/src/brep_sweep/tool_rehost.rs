//! Shared leaf re-hosts for the BRepSweep translation (pool-free style,
//! brep_algo/tool.rs precedent).
//!
//! Each function carries its OCCT anchor. These are TKBRep / TKG2d / TKG3d /
//! TKMath leaf helpers consumed by the TKPrim BRepSweep package; they are
//! re-hosted here (not GAP carriers) so the sweep engine runs end-to-end.
//! The owning-package batches can re-home them later.

use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor;
use rcad_kernel::geom::{
    transform_curve, Circle2d, Circle3, ConicalSurface, Curve2d, Curve3, CurveEval,
    CylindricalSurface, Line2d, Line3, Plane, RevolutionSurface, SphericalSurface, Surface3,
    ToroidalSurface, TrimmedCurve3,
};
use rcad_kernel::precision::{ANGULAR, CONFUSION, is_infinite_value, is_negative_infinite_value, is_positive_infinite_value};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{tshape_flags, Orientation, ShapeType, TShape};
use std::collections::HashSet;
use std::sync::Arc;

use crate::brep_algo::tool::{brep_tool_uv_points, explorer, shape_key, top_exp_vertices_raw};

// ---------------------------------------------------------------------------
// gp_Dir / gp_Vec comparison semantics (gp_Dir.hxx L164-205, gp_Vec.hxx
// IsEqual) — used by the GeomAdaptor classification and BRepLProp below.
// ---------------------------------------------------------------------------

/// OCCT gp_Dir::IsEqual(theOther, theAngularTolerance) (gp_Dir.hxx L164-167):
/// the angle between the two directions is lower or equal to the tolerance.
pub fn gp_dir_is_equal(a: glam::DVec3, b: glam::DVec3, ang_tol: f64) -> bool {
    dir_angle(a, b) <= ang_tol
}

/// OCCT gp_Dir::IsNormal(theOther, theAngularTolerance) (gp_Dir.hxx L171-180):
/// the angle is equal to Pi/2 within the tolerance.
pub fn gp_dir_is_normal(a: glam::DVec3, b: glam::DVec3, ang_tol: f64) -> bool {
    let mut an_ang = std::f64::consts::FRAC_PI_2 - dir_angle(a, b);
    if an_ang < 0.0 {
        an_ang = -an_ang;
    }
    an_ang <= ang_tol
}

/// OCCT gp_Dir::IsParallel(theOther, theAngularTolerance) (gp_Dir.hxx
/// L191-199): the angle is equal to 0 or to Pi within the tolerance.
pub fn gp_dir_is_parallel(a: glam::DVec3, b: glam::DVec3, ang_tol: f64) -> bool {
    let an_ang = dir_angle(a, b);
    an_ang <= ang_tol || std::f64::consts::PI - an_ang <= ang_tol
}

/// OCCT gp_Dir::Angle(theOther) — the angular value in [0, PI] between two
/// directions (gp_Dir.cxx Angle: acos of the dot product clamped to [-1, 1]).
fn dir_angle(a: glam::DVec3, b: glam::DVec3) -> f64 {
    let d = a.dot(b).clamp(-1.0, 1.0);
    d.acos()
}

/// OCCT gp_Vec::IsEqual(Other, LinearTolerance, AngularTolerance)
/// (gp_Vec.hxx IsEqual): both magnitudes below the linear tolerance are
/// equal; otherwise the magnitude difference must be within the linear
/// tolerance and the directions equal within the angular tolerance.
pub fn gp_vec_is_equal(a: glam::DVec3, b: glam::DVec3, lin_tol: f64, ang_tol: f64) -> bool {
    let the_norm = a.length();
    let the_other_norm = b.length();
    if the_norm <= lin_tol && the_other_norm <= lin_tol {
        return true;
    }
    let the_diff = the_norm - the_other_norm;
    if the_diff > lin_tol || -the_diff > lin_tol {
        return false;
    }
    gp_dir_is_equal(a, b, ang_tol)
}

// ---------------------------------------------------------------------------
// BRep_Tool::IsClosed(theShape) (BRep_Tool.cxx L1707-1763).
// ---------------------------------------------------------------------------

/// OCCT TopoDS_Shape::Oriented(theOr) — a copy carrying the orientation.
pub fn oriented(s: &Shape, the_or: Orientation) -> Shape {
    let mut c = s.clone();
    c.orientation = the_or;
    c
}

/// OCCT BRep_Tool::Degenerated(E) (BRep_Tool.cxx L1162-1171) — the edge
/// degenerated flag read.
pub fn brep_tool_degenerated(e: &Shape) -> bool {
    match e.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT BRep_Tool::IsClosed(theShape) (BRep_Tool.cxx L1707-1763): SHELL walks
/// its edges (add/remove on the occurrence map), WIRE walks its vertices,
/// EDGE compares the extremities, other types read the Closed flag.
pub fn brep_tool_is_closed_shape(the_shape: &Shape) -> bool {
    if the_shape.shape_type() == ShapeType::Shell {
        let mut a_map: HashSet<(u64, u32)> = HashSet::new();
        let mut has_bound = false;
        for e in explorer(
            &oriented(the_shape, Orientation::Forward),
            ShapeType::Edge,
            ShapeType::Shape,
        ) {
            let degenerated = match e.data.as_ref() {
                TShape::Edge(ed) => ed.degenerated,
                _ => false,
            };
            if degenerated
                || e.orientation == Orientation::Internal
                || e.orientation == Orientation::External
            {
                continue;
            }
            has_bound = true;
            if !a_map.insert(shape_key(&e)) {
                a_map.remove(&shape_key(&e));
            }
        }
        return has_bound && a_map.is_empty();
    } else if the_shape.shape_type() == ShapeType::Wire {
        let mut a_map: HashSet<(u64, u32)> = HashSet::new();
        let mut has_bound = false;
        for v in explorer(
            &oriented(the_shape, Orientation::Forward),
            ShapeType::Vertex,
            ShapeType::Shape,
        ) {
            if v.orientation == Orientation::Internal || v.orientation == Orientation::External {
                continue;
            }
            has_bound = true;
            if !a_map.insert(shape_key(&v)) {
                a_map.remove(&shape_key(&v));
            }
        }
        return has_bound && a_map.is_empty();
    } else if the_shape.shape_type() == ShapeType::Edge {
        let (a_v_first, a_v_last) = top_exp_vertices_raw(the_shape);
        return match (&a_v_first, &a_v_last) {
            (Some(vf), Some(vl)) => vf.is_same(vl),
            _ => false,
        };
    }
    crate::brep_algo::tool::shape_is_closed(the_shape)
}

// ---------------------------------------------------------------------------
// BRep_Tool::IsClosed(E, F) (BRep_Tool.cxx L795-812) — the edge occurrence
// on the face; the Triangulation branch has no rcad equivalent (the sweep
// pipeline carries no triangulations — architecture difference).
// ---------------------------------------------------------------------------

pub fn brep_tool_is_closed_edge_face(e: &Shape, f: &Shape) -> bool {
    // OCCT L797-802: IsClosed(E, S, L) (BRep_Tool.cxx L814-841) — the edge
    // carries two pcurves on the face surface (a CurveOnClosedSurface
    // representation matching the face surface).
    let ed = match e.data.as_ref() {
        TShape::Edge(ed) => ed,
        _ => return false,
    };
    let fkey = shape_key(f);
    ed.representations
        .iter()
        .any(|cr| matches!(cr, rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface { face, .. } if *face == fkey))
}

/// OCCT BRep_Tool::CurveOnSurface(E, F, f, l) (BRep_Tool.cxx L339-368): a
/// seam edge occurrence with REVERSED orientation selects the second pcurve
/// of its BRep_CurveOnClosedSurface pair; otherwise the pcurve keyed by the
/// face (the L.Predivided(E.Location()) form).
pub fn brep_tool_curve_on_surface_seam(e: &Shape, f: &Shape) -> Option<(Curve2d, f64, f64)> {
    let ed = match e.data.as_ref() {
        TShape::Edge(ed) => ed,
        _ => return None,
    };
    let fkey = shape_key(f);
    for r in &ed.representations {
        match r {
            rcad_kernel::topods::CurveRepresentation::CurveOnSurface { face, pcurve, range } => {
                if *face == fkey {
                    return Some((pcurve.clone(), range[0], range[1]));
                }
            }
            rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                face,
                pcurve1,
                pcurve2,
                range,
            } => {
                if *face == fkey {
                    if e.orientation == Orientation::Reversed {
                        return Some((pcurve2.clone(), range[0], range[1]));
                    }
                    return Some((pcurve1.clone(), range[0], range[1]));
                }
            }
            _ => {}
        }
    }
    ed.pcurves
        .get(&fkey)
        .map(|(c, f0, l0)| (c.clone(), *f0, *l0))
}

// ---------------------------------------------------------------------------
// BRep_Tool::Parameters(V, F) (BRep_Tool.cxx L1661-1703).
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Parameters(V, F) (BRep_Tool.cxx L1661-1703): first the
/// vertex PointOnSurface representations are consulted, then (PMN 4/06/97)
/// the face edges are walked for an extremity match and the UVPoints of the
/// matching end are returned.
pub fn brep_tool_parameters(v: &Shape, f: &Shape) -> glam::DVec2 {
    // OCCT L1667-1677: the PointRepresentation branch.  The rcad
    // PointRepresentation::PointOnSurface carries a pool index (architecture
    // difference: no surface value / location pair), and the pool-free sweep
    // pipeline never creates vertex surface points, so the branch cannot
    // match here — the edge walk below is the producing path.
    // OCCT L1686-1700: the edge walk.
    for e in explorer(f, ShapeType::Edge, ShapeType::Shape) {
        let (v_f, v_l) = top_exp_vertices_raw(&e);
        let same_f = v_f.as_ref().map_or(false, |vf| v.is_same(vf));
        let same_l = v_l.as_ref().map_or(false, |vl| v.is_same(vl));
        if same_f || same_l {
            let (p_f, p_l) = brep_tool_uv_points(&e, f);
            if same_f {
                return p_f;
            } else {
                // Ambiguity (natural) for degenerated edges.
                return p_l;
            }
        }
    }
    // OCCT: throw Standard_NoSuchObject("BRep_Tool:: no parameters on surface").
    panic!("Standard_NoSuchObject: BRep_Tool:: no parameters on surface");
}

// ---------------------------------------------------------------------------
// BRepTools::IsReallyClosed(E, F) (BRepTools.cxx L1204-1220).
// ---------------------------------------------------------------------------

/// OCCT BRepTools::IsReallyClosed(E, F) (BRepTools.cxx L1204-1220): the edge
/// is a seam (two pcurves on the face) and occurs exactly twice among the
/// face edges.
pub fn brep_tools_is_really_closed(e: &Shape, f: &Shape) -> bool {
    if !brep_tool_is_closed_edge_face(e, f) {
        return false;
    }
    let mut nbocc = 0;
    for cur in explorer(f, ShapeType::Edge, ShapeType::Shape) {
        if cur.is_same(e) {
            nbocc += 1;
        }
    }
    nbocc == 2
}

// ---------------------------------------------------------------------------
// ElCLib::AdjustPeriodic (ElCLib.cxx).
// ---------------------------------------------------------------------------

/// OCCT ElCLib::AdjustPeriodic(UFirst, ULast, Preci, U1, U2)
/// (ElCLib.cxx AdjustPeriodic): both values are folded into the periodic
/// range, U1 before U2, preserving their relative order.
pub fn elclib_adjust_periodic(u_first: f64, u_last: f64, preci: f64, u1: &mut f64, u2: &mut f64) {
    // OCCT ElCLib.cxx L122-128: Precision::IsInfinite(UFirst/ULast) guard
    // (Precision.hxx L350-353).
    if is_infinite_value(u_first) || is_infinite_value(u_last) {
        *u1 = u_first;
        *u2 = u_last;
        return;
    }
    let a_period = u_last - u_first;
    if a_period < f64::EPSILON * u_last.abs() {
        // In order to avoid FLT_Overflow exception (test bugs moddata_1
        // bug22757).
        *u1 = u_first;
        *u2 = u_last;
        return;
    }
    *u1 -= ((*u1 - u_first) / a_period).floor() * a_period;
    if u_last - *u1 < preci {
        *u1 -= a_period;
    }
    *u2 -= ((*u2 - *u1) / a_period).floor() * a_period;
    if *u2 - *u1 < preci {
        *u2 += a_period;
    }
}

// ---------------------------------------------------------------------------
// GeomAdaptor_Surface (the generic adaptor over a concrete surface value) —
// the GetType passthrough to the surface's own type.
// ---------------------------------------------------------------------------

/// The GeomAbs_SurfaceType encoding (GeomAbs_SurfaceType.hxx).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomAbsSurfaceType {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    BezierSurface,
    BSplineSurface,
    SurfaceOfRevolution,
    SurfaceOfExtrusion,
    OffsetSurface,
    OtherSurface,
}

/// OCCT GeomAdaptor_Surface::GetType() over a concrete Geom_Surface — the
/// passthrough to the surface's own type.
pub fn geom_adaptor_surface_get_type(s: &Surface3) -> GeomAbsSurfaceType {
    match s {
        Surface3::Plane(_) => GeomAbsSurfaceType::Plane,
        Surface3::Cylinder(_) => GeomAbsSurfaceType::Cylinder,
        Surface3::Cone(_) => GeomAbsSurfaceType::Cone,
        Surface3::Sphere(_) => GeomAbsSurfaceType::Sphere,
        Surface3::Torus(_) => GeomAbsSurfaceType::Torus,
        Surface3::Bezier(_) => GeomAbsSurfaceType::BezierSurface,
        Surface3::BSpline(_) => GeomAbsSurfaceType::BSplineSurface,
        Surface3::Revolution(_) => GeomAbsSurfaceType::SurfaceOfRevolution,
        Surface3::LinearExtrusion(_) => GeomAbsSurfaceType::SurfaceOfExtrusion,
        Surface3::Offset(_) => GeomAbsSurfaceType::OffsetSurface,
        _ => GeomAbsSurfaceType::OtherSurface,
    }
}

// ---------------------------------------------------------------------------
// The GeomAbs_CurveType of the adaptor basis curves (GeomAdaptor_Curve /
// Geom_Curve::GetType).  A trimmed curve reports its basis type
// (Geom_TrimmedCurve::GetType).
// ---------------------------------------------------------------------------

// Canonical GeomAbs_CurveType lives in rcad_kernel::math (OCCT
// GeomAbs_CurveType.hxx L23-31, nine enumerators in OCCT order).  Re-exported
// here so the historic import path super::tool_rehost::GeomAbsCurveType
// (brep_sweep/rotation.rs) keeps resolving; the local copy with the shortened
// variant names (Bezier/BSpline/Offset/Other) was deleted (Rule 4).
pub use rcad_kernel::math::GeomAbsCurveType;

/// The GeomAbs_CurveType of an rcad curve value.
pub fn curve_type(c: &Curve3) -> GeomAbsCurveType {
    match c {
        Curve3::Line(_) => GeomAbsCurveType::Line,
        Curve3::Circle(_) => GeomAbsCurveType::Circle,
        Curve3::Ellipse(_) => GeomAbsCurveType::Ellipse,
        Curve3::Hyperbola(_) => GeomAbsCurveType::Hyperbola,
        Curve3::Parabola(_) => GeomAbsCurveType::Parabola,
        Curve3::Bezier(_) => GeomAbsCurveType::BezierCurve,
        Curve3::BSpline(_) => GeomAbsCurveType::BSplineCurve,
        Curve3::Offset(_) => GeomAbsCurveType::OffsetCurve,
        Curve3::Trimmed(t) => curve_type(&t.curve),
        _ => GeomAbsCurveType::OtherCurve,
    }
}

/// Distance of a point from a gp_Lin (Location, Direction).
fn line_point_distance(loc: glam::DVec3, dir: glam::DVec3, p: glam::DVec3) -> f64 {
    let v = p - loc;
    (v - v.dot(dir) * dir).length()
}

// ---------------------------------------------------------------------------
// GeomAdaptor_SurfaceOfLinearExtrusion re-host (TKG3d/GeomAdaptor,
// GeomAdaptor_SurfaceOfLinearExtrusion.cxx) — Load (L90-113), GetType
// (L269-328), Plane (L332-371), Cylinder (L375-387).
// ---------------------------------------------------------------------------

/// OCCT GeomAdaptor_SurfaceOfLinearExtrusion (the consumed subset).  OCCT
/// `Handle(Adaptor3d_Curve) myBasisCurve` is a loaded `GeomAdaptor_Curve`
/// (the callers construct it as `GeomAdaptor_Curve(C, First, Last)` —
/// BRepSweep_Translation.cxx L243) whose `load` unwraps a Geom_TrimmedCurve
/// to its basis (GeomAdaptor_Curve.cxx L252-255); the rcad equivalent is the
/// kernel GeomCurveAdaptor carrying the same unwrap and the loaded range.
pub struct GeomAdaptorSurfaceOfLinearExtrusion {
    /// OCCT myBasisCurve.
    my_basis_curve: GeomCurveAdaptor,
    /// OCCT myDirection — the (normalized) extrusion direction.
    my_direction: glam::DVec3,
}

impl GeomAdaptorSurfaceOfLinearExtrusion {
    /// OCCT GeomAdaptor_SurfaceOfLinearExtrusion(C, V) — Load(C); Load(V)
    /// (cxx L90-113).
    pub fn load(the_curve: &GeomCurveAdaptor, the_dir: glam::DVec3) -> Self {
        GeomAdaptorSurfaceOfLinearExtrusion {
            my_basis_curve: the_curve.clone(),
            my_direction: the_dir.normalize_or_zero(),
        }
    }

    /// OCCT GeomAdaptor_SurfaceOfLinearExtrusion::GetType()
    /// (GeomAdaptor_SurfaceOfLinearExtrusion.cxx L269-328): the basis-curve
    /// type decides the canonization.
    pub fn get_type(&self) -> GeomAbsSurfaceType {
        match self.my_basis_curve.my_type_curve {
            GeomAbsCurveType::Line => {
                // OCCT L274-281: D = myBasisCurve->Line().Direction().
                let d = self.my_basis_curve.line().direction;
                if !gp_dir_is_parallel(self.my_direction, d, ANGULAR) {
                    return GeomAbsSurfaceType::Plane;
                }
            }
            GeomAbsCurveType::Circle => {
                // OCCT L283-294: D = myBasisCurve->Circle().Axis().Direction().
                let d = self.my_basis_curve.circle().normal;
                if gp_dir_is_parallel(self.my_direction, d, ANGULAR) {
                    return GeomAbsSurfaceType::Cylinder;
                } else if gp_dir_is_normal(self.my_direction, d, ANGULAR) {
                    return GeomAbsSurfaceType::Plane;
                }
            }
            GeomAbsCurveType::Ellipse => {
                // OCCT L296-303: D = myBasisCurve->Ellipse().Axis().Direction().
                let d = self.my_basis_curve.ellipse().normal;
                if gp_dir_is_normal(self.my_direction, d, ANGULAR) {
                    return GeomAbsSurfaceType::Plane;
                }
            }
            GeomAbsCurveType::Parabola => {
                // OCCT L305-312: D = myBasisCurve->Parabola().Axis().Direction().
                let d = self.my_basis_curve.parabola().normal;
                if gp_dir_is_normal(self.my_direction, d, ANGULAR) {
                    return GeomAbsSurfaceType::Plane;
                }
            }
            GeomAbsCurveType::Hyperbola => {
                // OCCT L314-321: D = myBasisCurve->Hyperbola().Axis().Direction().
                let d = self.my_basis_curve.hyperbola().normal;
                if gp_dir_is_normal(self.my_direction, d, ANGULAR) {
                    return GeomAbsSurfaceType::Plane;
                }
            }
            _ => {}
        }
        GeomAbsSurfaceType::SurfaceOfExtrusion
    }

    /// OCCT GeomAdaptor_SurfaceOfLinearExtrusion::Plane() (L332-371): the
    /// plane frame built from a sampled D1 point/normal search over the
    /// loaded basis-curve parameter range (FirstParameter/LastParameter).
    pub fn plane(&self) -> Plane {
        let mut p = glam::DVec3::ZERO;
        let mut d1u = glam::DVec3::ZERO;
        let mut new_z = glam::DVec3::ZERO;
        // OCCT L339-340: UFirst = myBasisCurve->FirstParameter();
        // ULast = myBasisCurve->LastParameter() — the loaded range.
        let mut u_first = self.my_basis_curve.first;
        let mut u_last = self.my_basis_curve.last;
        // OCCT L341-352: Precision::IsNegativeInfinite(UFirst) /
        // Precision::IsPositiveInfinite(ULast) gates (Precision.hxx L357-367).
        if is_negative_infinite_value(u_first) && is_positive_infinite_value(u_last) {
            u_first = -100.0;
            u_last = 100.0;
        } else if is_negative_infinite_value(u_first) {
            u_first = u_last - 200.0;
        } else if is_positive_infinite_value(u_last) {
            u_last = u_first + 200.0;
        }
        let deltau = (u_last - u_first) / 20.0;
        for i in 1..=21 {
            let prm = u_first + (i - 1) as f64 * deltau;
            // OCCT L358: myBasisCurve->D1(prm, P, D1u).
            let (pp, dd) = self.my_basis_curve.d1_at(prm);
            p = pp;
            d1u = dd;
            new_z = d1u.normalize_or_zero().cross(self.my_direction);
            if new_z.length() > 1.0e-12 {
                break;
            }
        }
        // OCCT: gp_Ax3 Ax3(P, gp_Dir(newZ), gp_Dir(D1u)) — the X direction is
        // D1u projected perpendicular to newZ; YDirection = newZ x XDirection;
        // if myDirection dot YDirection < 0 the Y direction is reversed.
        let normal = new_z.normalize_or_zero();
        let mut x_dir = d1u - normal * d1u.dot(normal);
        if x_dir.length_squared() < 1e-24 {
            x_dir = any_perpendicular(normal);
        }
        let x_dir = x_dir.normalize();
        let mut y_dir = normal.cross(x_dir);
        if self.my_direction.dot(y_dir) < 0.0 {
            // OCCT: Ax3.YReverse() (gp_Ax3.hxx L128 — vydir.Reverse()).
            y_dir = -y_dir;
        }
        Plane {
            origin: p,
            normal,
            u_dir: x_dir,
            v_dir: y_dir,
        }
    }

    /// OCCT GeomAdaptor_SurfaceOfLinearExtrusion::Cylinder() (L375-387): the
    /// basis circle frame with the axis sign aligned to the extrusion
    /// direction.
    pub fn cylinder(&self) -> CylindricalSurface {
        // OCCT L380: gp_Circ C = myBasisCurve->Circle().
        let c = self.my_basis_curve.circle();
        let (center, normal, x_dir, radius) = (c.center, c.normal, c.x_dir, c.radius);
        // OCCT gp_Cylinder(Ax3, R): the frame Y = axis x X at construction;
        // the ZReverse below flips ONLY the axis (the Y is preserved, the
        // left-handed swept-lateral frame).
        let y_dir = normal.cross(x_dir);
        let mut axis = normal;
        let mut ref_dir = x_dir;
        if self.my_direction.dot(normal) < 0.0 {
            // OCCT: Ax3.ZReverse() (gp_Ax3.hxx L131 — axis.Reverse(); the X
            // and Y directions are kept).
            axis = -axis;
            ref_dir = ref_dir - axis * ref_dir.dot(axis);
            if ref_dir.length_squared() < 1e-24 {
                ref_dir = any_perpendicular(axis);
            }
            ref_dir = ref_dir.normalize();
        }
        CylindricalSurface {
            origin: center,
            axis,
            radius,
            ref_dir,
            y_dir: Some(y_dir),
        }
    }
}

// ---------------------------------------------------------------------------
// GeomAdaptor_SurfaceOfRevolution re-host (TKG3d/GeomAdaptor,
// GeomAdaptor_SurfaceOfRevolution.cxx) — Load (L89-199), GetType (L370-468),
// Plane (L472-489), Cylinder (L493-501), Cone (L505-527), Sphere (L531-540),
// Torus (L544-552).
// ---------------------------------------------------------------------------

/// OCCT GeomAdaptor_SurfaceOfRevolution (the consumed subset).  OCCT
/// `Handle(Adaptor3d_Curve) myBasisCurve` is a loaded `GeomAdaptor_Curve`
/// (the callers construct it as `GeomAdaptor_Curve(); HC->Load(C, First,
/// Last)` — BRepSweep_Rotation.cxx L333-334) whose `load` unwraps a
/// Geom_TrimmedCurve to its basis (GeomAdaptor_Curve.cxx L252-255); the rcad
/// equivalent is the kernel GeomCurveAdaptor carrying the same unwrap and
/// the loaded range.  myAxeRev is the gp_Ax3 frame (Location, Direction,
/// XDirection, YDirection).
pub struct GeomAdaptorSurfaceOfRevolution {
    /// OCCT myBasisCurve.
    my_basis_curve: GeomCurveAdaptor,
    /// OCCT myAxis — the loaded revolution axis (gp_Ax1: Location, Direction).
    my_axis: (glam::DVec3, glam::DVec3),
    /// OCCT myAxeRev — the gp_Ax3 frame (Location, Direction, XDir, YDir).
    my_axe_rev: (glam::DVec3, glam::DVec3, glam::DVec3, glam::DVec3),
}

impl GeomAdaptorSurfaceOfRevolution {
    /// OCCT GeomAdaptor_SurfaceOfRevolution::Load(C, V) (cxx L89-199): the
    /// meridian curve, the axis, and the myAxeRev frame determination.
    pub fn load(the_curve: &GeomCurveAdaptor, axis_loc: glam::DVec3, axis_dir: glam::DVec3) -> Self {
        // The OCCT fresh adaptor carries the default gp_Ax3 frame; Load
        // recomputes it.  TheGetType call inside Load (the Cone case) reads
        // the frame as-is at that point — the same literal order.
        let mut me = GeomAdaptorSurfaceOfRevolution {
            my_basis_curve: the_curve.clone(),
            my_axis: (axis_loc, axis_dir.normalize_or_zero()),
            my_axe_rev: (
                glam::DVec3::ZERO,
                glam::DVec3::Z,
                glam::DVec3::X,
                glam::DVec3::Y,
            ),
        };
        me.load_axis();
        me
    }

    /// OCCT Value(u, v) — the meridian point at parameter v rotated by the
    /// angle u about the axis; only the u = 0 form is consumed (the unrotated
    /// meridian point).  The basis-curve evaluation reads through the loaded
    /// adaptor (`myBasisCurve->Value(v)`).
    fn value00(&self, v: f64) -> glam::DVec3 {
        self.my_basis_curve.value_at(v)
    }

    /// OCCT Load(const gp_Ax1& V) (cxx L100-199) — the myAxeRev frame.
    fn load_axis(&mut self) {
        let mut o = self.my_axis.0;
        let mut oz = self.my_axis.1;
        let mut yrev = false;
        // OCCT L119-126: if (myBasisCurve->GetType() == GeomAbs_Line)
        //   if ((myBasisCurve->Line().Direction()).Dot(Oz) < 0.) { yrev; Oz.Reverse(); }
        if self.my_basis_curve.my_type_curve == GeomAbsCurveType::Line {
            if self.my_basis_curve.line().direction.dot(oz) < 0.0 {
                yrev = true;
                oz = -oz;
            }
        }

        let p: glam::DVec3;
        let q: glam::DVec3;
        // OCCT L128-155.
        if self.my_basis_curve.my_type_curve == GeomAbsCurveType::Circle {
            // OCCT L130: Q = P = (myBasisCurve->Circle()).Location().
            p = self.my_basis_curve.circle().center;
            q = p;
        } else {
            // OCCT L134: First = myBasisCurve->FirstParameter() — the loaded
            // range bound of the adaptor.
            let first = self.my_basis_curve.first;
            p = self.value00(0.0); // which does not mean much
            if self.get_type() == GeomAbsSurfaceType::Cone {
                if line_point_distance(self.my_axis.0, self.my_axis.1, p) <= CONFUSION {
                    // OCCT L140: Q = ElCLib::Value(1., myBasisCurve->Line()).
                    let l = self.my_basis_curve.line();
                    q = l.origin + l.direction;
                } else {
                    q = p;
                }
            // OCCT GeomAdaptor_SurfaceOfRevolution.cxx L147:
            // Precision::IsInfinite(First).
            } else if is_infinite_value(first) {
                q = p;
            } else {
                // OCCT L153: Q = Value(0., First).
                q = self.value00(first);
            }
        }

        // OCCT: gp_Dir DZ = myAxis.Direction() — the ORIGINAL axis direction.
        let dz = self.my_axis.1;
        o = o + (p - o).dot(dz) * dz;
        let ox: glam::DVec3;
        if line_point_distance(self.my_axis.0, dz, q) > CONFUSION {
            ox = (q - o).normalize_or_zero();
        } else {
            // OCCT L165-166: First/Last = myBasisCurve->FirstParameter() /
            // LastParameter() — the loaded range.
            let first = self.my_basis_curve.first;
            let last = self.my_basis_curve.last;
            let mut ratio = 1;
            let mut dist;
            let mut pp;
            loop {
                pp = self.my_basis_curve.value_at(first + (last - first) / ratio as f64);
                dist = line_point_distance(self.my_axis.0, dz, pp);
                ratio += 1;
                if !(dist < CONFUSION && ratio < 100) {
                    break;
                }
            }
            if ratio >= 100 {
                // OCCT: throw Standard_ConstructionError(
                //   "Adaptor3d_SurfaceOfRevolution : Axe and meridian are confused").
                panic!("Standard_ConstructionError: Adaptor3d_SurfaceOfRevolution : Axe and meridian are confused");
            }
            // OCCT: Ox = ((Oz ^ gp_Vec(PP.XYZ() - O.XYZ())) ^ Oz).
            ox = oz.cross(pp - o).cross(oz).normalize_or_zero();
        }

        // OCCT: myAxeRev = gp_Ax3(O, Oz, Ox) — YDirection = Oz x Ox.
        let mut vy = oz.cross(ox);
        self.my_axe_rev = (o, oz, ox, vy);

        if yrev {
            // OCCT: myAxeRev.YReverse() — the Y direction flips.
            vy = -vy;
            self.my_axe_rev.3 = vy;
        } else if self.my_basis_curve.my_type_curve == GeomAbsCurveType::Circle {
            // OCCT L193: DC = (myBasisCurve->Circle()).Axis().Direction().
            let dc = self.my_basis_curve.circle().normal;
            if ox.cross(oz).dot(dc) < 0.0 {
                // OCCT: myAxeRev.ZReverse() — the direction flips; the X and
                // Y directions are kept (gp_Ax3.hxx L131).
                self.my_axe_rev.1 = -self.my_axe_rev.1;
            }
        }
    }

    /// OCCT GeomAdaptor_SurfaceOfRevolution::GetType() (cxx L370-468).
    pub fn get_type(&self) -> GeomAbsSurfaceType {
        let tol_conf = CONFUSION;
        let tol_ang = ANGULAR;
        let tol_cone_semi_ang = CONFUSION;

        match self.my_basis_curve.my_type_curve {
            GeomAbsCurveType::Line => {
                // OCCT L379: gp_Ax1 Axe = myBasisCurve->Line().Position().
                let l = self.my_basis_curve.line();
                let (axe_loc, axe_dir) = (l.origin, l.direction);
                if gp_dir_is_parallel(self.my_axis.1, axe_dir, tol_ang) {
                    let p = self.value00(0.0);
                    let axe_rev = self.my_axe_rev;
                    let r = (p - axe_rev.0).dot(axe_rev.2);
                    if r > tol_conf {
                        return GeomAbsSurfaceType::Cylinder;
                    }
                } else if gp_dir_is_normal(self.my_axis.1, axe_dir, tol_ang) {
                    return GeomAbsSurfaceType::Plane;
                } else {
                    // OCCT L396-398: uf/ul = myBasisCurve->FirstParameter() /
                    // LastParameter() — the loaded range.
                    let uf = self.my_basis_curve.first;
                    let ul = self.my_basis_curve.last;
                    // OCCT GeomAdaptor_SurfaceOfRevolution.cxx L398:
                    // istrim = !Precision::IsInfinite(uf) && !Precision::IsInfinite(ul).
                    let istrim = !is_infinite_value(uf) && !is_infinite_value(ul);
                    if istrim {
                        let pf = self.my_basis_curve.value_at(uf);
                        let pl = self.my_basis_curve.value_at(ul);
                        let len = pf.distance(pl);
                        // Compute the distance projected onto the axis.
                        let vlin = pl - pf;
                        let vaxe = self.my_axis.1;
                        let projlen = vaxe.dot(vlin).abs();
                        if len - projlen <= tol_conf {
                            let p = self.value00(0.0);
                            let axe_rev = self.my_axe_rev;
                            let r = (p - axe_rev.0).dot(axe_rev.2);
                            if r > tol_conf {
                                return GeomAbsSurfaceType::Cylinder;
                            }
                        } else if projlen <= tol_conf {
                            return GeomAbsSurfaceType::Plane;
                        }
                    }
                    // OCCT L422: gp_Vec V(myAxis.Location(),
                    // myBasisCurve->Line().Location()).
                    let v = axe_loc - self.my_axis.0;
                    let w = axe_dir;
                    let axis_dir = self.my_axis.1;
                    let proj = w.dot(axis_dir).abs();
                    if v.dot(axis_dir.cross(w)).abs() <= tol_conf
                        && proj >= tol_cone_semi_ang
                        && proj <= 1.0 - tol_cone_semi_ang
                    {
                        return GeomAbsSurfaceType::Cone;
                    }
                }
            }
            GeomAbsCurveType::Circle => {
                // OCCT L439-441: gp_Circ C = myBasisCurve->Circle().
                let c = self.my_basis_curve.circle();
                let (center, normal, _x_dir, radius) = (c.center, c.normal, c.x_dir, c.radius);
                let a_r = radius;
                // OCCT: C.Position().IsCoplanar(myAxis, TolConf, TolAng)
                // (gp_Ax2.hxx L593-607): |circleNormal . (axisLoc - center)|
                // <= TolConf and circleNormal IsNormal(axisDir, TolAng).
                let d1 = normal.dot(self.my_axis.0 - center).abs();
                let coplanar =
                    d1 <= tol_conf && gp_dir_is_normal(normal, self.my_axis.1, tol_ang);
                if !coplanar {
                    return GeomAbsSurfaceType::SurfaceOfRevolution;
                } else if line_point_distance(self.my_axis.0, self.my_axis.1, center) <= tol_conf {
                    return GeomAbsSurfaceType::Sphere;
                } else {
                    let major_radius =
                        line_point_distance(self.my_axis.0, self.my_axis.1, center);
                    if major_radius > a_r {
                        return GeomAbsSurfaceType::Torus;
                    }
                }
            }
            _ => {}
        }
        GeomAbsSurfaceType::SurfaceOfRevolution
    }

    /// OCCT GeomAdaptor_SurfaceOfRevolution::Plane() (cxx L472-489).
    pub fn plane(&self) -> Plane {
        let (loc, dir, mut x_dir, y_dir) = self.my_axe_rev;
        let a_pon_curve = self.value00(0.0);
        let a_dot = (a_pon_curve - self.my_axis.0).dot(self.my_axis.1);
        let p = self.my_axis.0 + a_dot * self.my_axis.1;
        // OCCT: Axe.SetLocation(P).
        let _ = loc;
        // OCCT L483: Axe.XDirection().Dot(myBasisCurve->Line().Direction()).
        if x_dir.dot(self.my_basis_curve.line().direction) >= -CONFUSION {
            // OCCT: Axe.XReverse() (gp_Ax3.hxx L125 — vxdir.Reverse()).
            x_dir = -x_dir;
        }
        Plane {
            origin: p,
            normal: dir,
            u_dir: x_dir,
            v_dir: y_dir,
        }
    }

    /// OCCT GeomAdaptor_SurfaceOfRevolution::Cylinder() (cxx L493-501).
    pub fn cylinder(&self) -> CylindricalSurface {
        let p = self.value00(0.0);
        let axe_rev = self.my_axe_rev;
        let r = (p - axe_rev.0).dot(axe_rev.2);
        CylindricalSurface {
            origin: axe_rev.0,
            axis: axe_rev.1,
            radius: r,
            ref_dir: axe_rev.2,
            y_dir: None,
        }
    }

    /// OCCT GeomAdaptor_SurfaceOfRevolution::Cone() (cxx L505-527).
    pub fn cone(&self) -> ConicalSurface {
        let axe_rev = self.my_axe_rev;
        // OCCT L510: gp_Dir ldir = (myBasisCurve->Line()).Direction().
        let ldir = self.my_basis_curve.line().direction;
        let mut angle = dir_angle(axe_rev.1, ldir);
        let p0 = self.value00(0.0);
        let r = axe_rev.0.distance(p0);
        if r >= CONFUSION {
            let o = axe_rev.0;
            let mut op0 = p0 - o;
            let mut t = op0.dot(axe_rev.2);
            t /= ldir.dot(axe_rev.2);
            op0 = op0 - t * ldir;
            if op0.dot(axe_rev.1) > 0.0 {
                angle = -angle;
            }
        }
        // OCCT gp_Cone(Axe, Angle, Radius): the reference circle of Radius
        // sits in the frame origin plane (the rcad apex field is the point
        // on the axis where the surface radius equals `radius`).
        ConicalSurface {
            apex: axe_rev.0,
            axis: axe_rev.1,
            radius: r,
            half_angle_rad: angle,
            ref_dir: axe_rev.2,
        }
    }

    /// OCCT GeomAdaptor_SurfaceOfRevolution::Sphere() (cxx L531-540).
    pub fn sphere(&self) -> SphericalSurface {
        let c = self.my_basis_curve.circle();
        let (center, _normal, _x_dir, radius) = (c.center, c.normal, c.x_dir, c.radius);
        let axe_rev = self.my_axe_rev;
        SphericalSurface {
            center,
            axis: axe_rev.1,
            radius,
            ref_dir: axe_rev.2,
        }
    }

    /// OCCT GeomAdaptor_SurfaceOfRevolution::Torus() (cxx L544-552).
    pub fn torus(&self) -> ToroidalSurface {
        let c = self.my_basis_curve.circle();
        let (_center, _normal, _x_dir, radius) = (c.center, c.normal, c.x_dir, c.radius);
        let major_radius = line_point_distance(self.my_axis.0, self.my_axis.1, c.center);
        let axe_rev = self.my_axe_rev;
        ToroidalSurface {
            center: axe_rev.0,
            axis: axe_rev.1,
            ref_dir: axe_rev.2,
            major_radius,
            minor_radius: radius,
        }
    }

    /// OCCT `new Geom_SurfaceOfRevolution(C, myAxe)`
    /// (BRepSweep_Rotation.cxx L363-366) — the swept surface value over the
    /// trimmed generator curve C kept at the caller (NOT the unwrapped basis
    /// of the adaptor).
    pub fn surface_value(&self, the_curve: Curve3) -> Surface3 {
        Surface3::Revolution(RevolutionSurface {
            profile: Box::new(the_curve),
            axis_origin: self.my_axis.0,
            axis_dir: self.my_axis.1,
        })
    }
}

// ---------------------------------------------------------------------------
// Geom curve value re-hosts consumed by Translation/Rotation.
// ---------------------------------------------------------------------------

/// OCCT Geom2d_Line(gp_Lin2d) construction — the Curve2d line value
/// (direction normalized per the gp_Dir2d invariant).
pub fn geom2d_line(origin: glam::DVec2, direction: glam::DVec2) -> Curve2d {
    Curve2d::Line(Line2d::new(origin, direction))
}

/// OCCT Geom2d_Circle(gp_Circ2d) construction.
pub fn geom2d_circle(center: glam::DVec2, radius: f64) -> Curve2d {
    Curve2d::Circle(Circle2d {
        center,
        x_dir: glam::DVec2::X,
        y_dir: glam::DVec2::Y,
        radius,
    })
}

/// OCCT Geom_Line(gp_Lin) construction.
pub fn geom_line(origin: glam::DVec3, direction: glam::DVec3) -> Curve3 {
    Curve3::Line(Line3 {
        origin,
        direction: direction.normalize_or_zero(),
    })
}

/// OCCT Geom_Circle(gp_Ax2, Radius) construction (axis = Location,
/// Direction(N), XDirection).
pub fn geom_circle(axis_loc: glam::DVec3, axis_dir: glam::DVec3, x_dir: glam::DVec3, radius: f64) -> Curve3 {
    let n = axis_dir.normalize_or_zero();
    let y = n.cross(x_dir);
    Curve3::Circle(Circle3 {
        center: axis_loc,
        normal: n,
        x_dir,
        y_dir: y.normalize_or_zero(),
        radius,
    })
}

/// OCCT Geom_TrimmedCurve(BasisCurve, U1, U2) construction.
pub fn geom_trimmed_curve(basis: Curve3, first: f64, last: f64) -> Curve3 {
    Curve3::Trimmed(TrimmedCurve3 {
        curve: Box::new(basis),
        first,
        last,
    })
}

/// OCCT Geom_Curve::Copy() + Transform(Tr) — the copied curve with the
/// transform applied (the rcad transform_curve; rigid transforms preserve
/// the parameterization).
pub fn geom_curve_transformed(c: &Curve3, t: &glam::DAffine3) -> Curve3 {
    transform_curve(c, t)
}

// ---------------------------------------------------------------------------
// In-place TShape mutation helpers (the kernel edge_mut_inplace /
// wire_mut pattern): OCCT BRep_Builder edits TShapes in place — the shared
// Arc must not split (the sweep stores the same handle in myShapes slots and
// in parent containers).
// ---------------------------------------------------------------------------

/// The TShape payload of a shape as a mutable reference, in place.
macro_rules! tshape_mut {
    ($shape:expr, $variant:path, $name:literal) => {
        match unsafe { &mut *(Arc::as_ptr(&$shape.data) as *mut TShape) } {
            $variant(v) => v,
            _ => panic!(concat!($name, ": shape is not the expected type")),
        }
    };
}

pub(crate) fn vertex_data_mut(s: &Shape) -> &mut rcad_kernel::topods::TVertexData {
    tshape_mut!(s, TShape::Vertex, "vertex_mut")
}

pub(crate) fn edge_data_mut(s: &Shape) -> &mut rcad_kernel::topods::TEdgeData {
    tshape_mut!(s, TShape::Edge, "edge_mut")
}

pub(crate) fn wire_data_mut(s: &Shape) -> &mut rcad_kernel::topods::TWireData {
    tshape_mut!(s, TShape::Wire, "wire_mut")
}

pub(crate) fn face_data_mut(s: &Shape) -> &mut rcad_kernel::topods::TFaceData {
    tshape_mut!(s, TShape::Face, "face_mut")
}

pub(crate) fn shell_data_mut(s: &Shape) -> &mut rcad_kernel::topods::TShellData {
    tshape_mut!(s, TShape::Shell, "shell_mut")
}

pub(crate) fn solid_data_mut(s: &Shape) -> &mut rcad_kernel::topods::TSolidData {
    tshape_mut!(s, TShape::Solid, "solid_mut")
}
/// OCCT TopoDS_Shape::Closed(theIsClosed) — the CLOSED flag write on the
/// shape's own TShape.
pub(crate) fn set_closed_flag(s: &Shape, flag: bool) {
    let flags = match unsafe { &mut *(Arc::as_ptr(&s.data) as *mut TShape) } {
        TShape::Vertex(v) => &mut v.flags,
        TShape::Edge(e) => &mut e.flags,
        TShape::Wire(w) => &mut w.flags,
        TShape::Face(f) => &mut f.flags,
        TShape::Shell(sh) => &mut sh.flags,
        TShape::Solid(so) => &mut so.flags,
        TShape::CompSolid(_) | TShape::Compound(_) => return,
    };
    if flag {
        *flags |= tshape_flags::CLOSED;
    } else {
        *flags &= !tshape_flags::CLOSED;
    }
}

/// OCCT TopoDS_Iterator (cumOri=true): the stored sub-shapes with
/// TopAbs::Compose(parent.Orientation(), child.Orientation()) applied
/// (TopoDS_Iterator.cxx L72-80; the brep_algo::tool re-host).
pub(crate) fn topods_sub_shapes(sh: &Shape) -> Vec<Shape> {
    crate::brep_algo::tool::sub_shapes(sh)
}

// ---------------------------------------------------------------------------
// BRepAdaptor_Curve::Continuity equivalent + BRepLProp::Continuity re-host
// (TKBRep/BRepLProp, BRepLProp.cxx L28-131).
// ---------------------------------------------------------------------------

/// The GeomAbs_Shape continuity ranking (GeomAbs_Shape.hxx L47-56):
/// C0 < G1 < C1 < G2 < C2 < C3 < CN.
pub const GEOM_ABS_C0: i32 = 0;
pub const GEOM_ABS_G1: i32 = 1;
pub const GEOM_ABS_C1: i32 = 2;
pub const GEOM_ABS_C2: i32 = 4;
pub const GEOM_ABS_C3: i32 = 5;
pub const GEOM_ABS_CN: i32 = 6;

/// OCCT BRepAdaptor_Curve::Continuity() (BRepAdaptor_Curve.cxx Continuity):
/// the analytic curves and Bezier are CN; the B-spline continuity is the
/// minimum over the knot spans, derived from the interior-knot
/// multiplicities C(D - M) (Geom_BSplineCurve::Continuity).  The rcad
/// BSplineCurve3 stores the flat knot vector, so the multiplicities are
/// recomputed from the repeats.
pub fn brep_adaptor_curve_continuity(c: &Curve3) -> i32 {
    match c {
        Curve3::BSpline(b) => {
            // Geom_BSplineCurve::Continuity: for every interior knot of
            // multiplicity M the continuity degrades to C(D - M); the flat
            // knot vector carries the multiplicities as repeats.
            let deg = b.degree as i32;
            let mut cont = GEOM_ABS_CN;
            let knots = &b.knots;
            if knots.len() >= 2 {
                let mut i = 1usize;
                while i < knots.len() - 1 {
                    let mut m = 1usize;
                    while i + m < knots.len() - 1 && (knots[i + m] - knots[i]).abs() <= 1e-12 {
                        m += 1;
                    }
                    cont = cont.min(deg - m as i32);
                    i += m;
                }
                if cont < 0 {
                    cont = 0;
                }
            }
            cont
        }
        // Line / Circle / Ellipse / Hyperbola / Parabola / Bezier / Offset /
        // Trimmed(basis): GeomAbs_CN.
        _ => GEOM_ABS_CN,
    }
}

/// The BRepLProp_CLProps subset consumed by BRepLProp::Continuity: value, D1,
/// D2 and tangent-defined at one parameter (BRepLProp_CLProps.hxx/.cxx: the
/// tangent is defined when |D1| > Resolution).
struct ClProps {
    curve: Curve3,
    resolution: f64,
    u: f64,
    n: i32,
}

impl ClProps {
    fn new(curve: &Curve3, u: f64, n: i32, resolution: f64) -> Self {
        ClProps {
            curve: curve.clone(),
            resolution,
            u,
            n,
        }
    }
    fn value(&self) -> glam::DVec3 {
        self.curve.point_at(self.u)
    }
    fn d1(&self) -> glam::DVec3 {
        if self.n >= 1 {
            self.curve.derivative_at(self.u)
        } else {
            glam::DVec3::ZERO
        }
    }
    fn d2(&self) -> glam::DVec3 {
        if self.n >= 2 {
            self.curve.derivative2_at(self.u)
        } else {
            glam::DVec3::ZERO
        }
    }
    fn is_tangent_defined(&self) -> bool {
        let d1 = self.curve.derivative_at(self.u);
        d1.length_squared() > self.resolution * self.resolution
    }
    fn tangent(&self) -> glam::DVec3 {
        self.curve.derivative_at(self.u).normalize_or_zero()
    }
}

/// OCCT BRepLProp::Continuity(C1, C2, u1, u2, tl, ta) (BRepLProp.cxx L28-131)
/// — the joint continuity classification of two edge curves at their meeting
/// parameters.  The orientation flags carry the OCCT C1.Edge().Orientation()
/// REVERSED checks; `same_edge_periodic` is the OCCT E1.IsSame(E2) &&
/// C1.IsPeriodic() final condition.
#[allow(clippy::too_many_arguments)]
pub fn brep_lprop_continuity(
    c1: &Curve3,
    c2: &Curve3,
    u1: f64,
    u2: f64,
    tl: f64,
    ta: f64,
    c1_reversed: bool,
    c2_reversed: bool,
    same_edge_periodic: bool,
) -> i32 {
    let mut cont = GEOM_ABS_C0;
    let mut fini = false;
    let cont1 = brep_adaptor_curve_continuity(c1);
    let cont2 = brep_adaptor_curve_continuity(c2);
    let mut n1 = 0;
    let mut n2 = 0;
    if cont1 >= 5 {
        n1 = 3;
    } else if cont1 == 4 {
        n1 = 2;
    } else if cont1 == 2 {
        n1 = 1;
    }
    if cont2 >= 5 {
        n2 = 3;
    } else if cont2 == 4 {
        n2 = 2;
    } else if cont2 == 2 {
        n2 = 1;
    }
    let clp1 = ClProps::new(c1, u1, n1, tl);
    let clp2 = ClProps::new(c2, u2, n2, tl);
    // OCCT: if (!(clp1.Value().IsEqual(clp2.Value(), tl)))
    //   throw Standard_Failure("Courbes non jointives");
    if !gp_vec_is_equal(clp1.value(), clp2.value(), tl, ta) {
        panic!("Standard_Failure: Courbes non jointives");
    }
    let min = n1.min(n2);
    if min >= 1 {
        let mut d1 = clp1.d1();
        let mut d2 = clp2.d1();
        if c1_reversed {
            d1 = -d1;
        }
        if c2_reversed {
            d2 = -d2;
        }
        if gp_vec_is_equal(d1, d2, tl, ta) {
            cont = GEOM_ABS_C1;
        } else if clp1.is_tangent_defined() && clp2.is_tangent_defined() {
            let mut dir1 = clp1.tangent();
            let mut dir2 = clp2.tangent();
            if c1_reversed {
                dir1 = -dir1;
            }
            if c2_reversed {
                dir2 = -dir2;
            }
            if gp_dir_is_equal(dir1, dir2, ta) {
                cont = GEOM_ABS_G1;
            }
            fini = true;
        } else {
            fini = true;
        }
    }
    if min >= 2 && !fini {
        let d1 = clp1.d2();
        let d2 = clp2.d2();
        if gp_vec_is_equal(d1, d2, tl, ta) {
            cont = GEOM_ABS_C2;
        }
    }
    // OCCT: if (E1.IsSame(E2) && C1.IsPeriodic() && cont >= GeomAbs_G1)
    //   cont = GeomAbs_CN.
    if same_edge_periodic && cont >= GEOM_ABS_G1 {
        cont = GEOM_ABS_CN;
    }
    cont
}

/// An arbitrary unit vector perpendicular to the input (the gp construction
/// of a reference direction when none is carried).
fn any_perpendicular(n: glam::DVec3) -> glam::DVec3 {
    let a = if n.x.abs() < 0.9 {
        glam::DVec3::X
    } else {
        glam::DVec3::Y
    };
    a.cross(n).normalize_or_zero()
}
