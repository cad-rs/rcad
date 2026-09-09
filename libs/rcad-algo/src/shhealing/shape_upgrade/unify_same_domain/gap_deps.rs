//! GAP carriers and small re-hosts local to ShapeUpgrade_UnifySameDomain
//! (W1-6).
//!
//! Dependencies owned by *other* packages/batches are carried here with the
//! OCCT anchor and the owning batch; pure math re-hosts (gp/ElCLib statics)
//! are architecture bridges.  Every carrier keeps OCCT's failure path where
//! the real algorithm is missing.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::CurveEval;
use rcad_kernel::geom::{
    BSplineCurve2, BSplineCurve3, BezierCurve2, Circle2d, Circle3, Curve2d, Curve2dEval, Curve3,
    Line2d, Line3, Surface3, TrimmedCurve2, TrimmedCurve3,
};
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, BRepTool, GeomAbsShape, TShape};
use rcad_kernel::{base::convert as math_convert, geom};

/// OCCT GeomAbs_CurveType (GeomAbs_CurveType.hxx L26-66) — the members the
/// class switches on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomAbsCurveType {
    Line,
    Circle,
    Ellipse,
    BSplineCurve,
    BezierCurve,
    OtherCurve,
}

/// OCCT GeomAbs_SurfaceType — the Plane check of BRepAdaptor_Surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomAbsSurfaceType {
    Plane,
    Other,
}

// ---------------------------------------------------------------------------
// OCCT BRepAdaptor_Curve (TKTopAlgo/BRepAdaptor) — reduced re-host over the
// rcad edge payload (architecture bridge; the full adaptor is a kernel
// front-batch candidate).
// ---------------------------------------------------------------------------

/// OCCT `BRepAdaptor_Curve` reduced to the members UnifySameDomain consumes
/// (GetType, FirstParameter, LastParameter, Value, D1, Line, Circle,
/// NbPoles).  The adaptor keeps the edge 3D curve parametrization (the OCCT
/// callers of this class handle edge orientation manually — see
/// IsMergingPossible cxx L2590-2608).
pub struct BRepAdaptorCurve {
    c3d: Option<Curve3>,
    range: [f64; 2],
}

impl BRepAdaptorCurve {
    /// OCCT BRepAdaptor_Curve(E) — the location-applied curve over the edge
    /// range (BRepAdaptor_Curve.cxx L96-131; the location rescale follows
    /// BRep_Tool::Curve(E, f, l)).
    pub fn new(brep: &BRep, e: &Shape) -> Self {
        let (c, range) = match super::topexp::brep_tool_curve(brep, e) {
            Some((c, f, l)) => (Some(c), [f, l]),
            None => (None, [0.0, 0.0]),
        };
        BRepAdaptorCurve { c3d: c, range }
    }

    /// OCCT GetType.
    pub fn get_type(&self) -> GeomAbsCurveType {
        match self.c3d.as_ref() {
            Some(Curve3::Line(_)) => GeomAbsCurveType::Line,
            Some(Curve3::Circle(_)) => GeomAbsCurveType::Circle,
            Some(Curve3::Ellipse(_)) => GeomAbsCurveType::Ellipse,
            Some(Curve3::BSpline(_)) => GeomAbsCurveType::BSplineCurve,
            Some(Curve3::Bezier(_)) => GeomAbsCurveType::BezierCurve,
            _ => GeomAbsCurveType::OtherCurve,
        }
    }

    /// OCCT FirstParameter — the edge range first (the trimmed bound).
    pub fn first_parameter(&self) -> f64 {
        self.range[0]
    }

    /// OCCT LastParameter.
    pub fn last_parameter(&self) -> f64 {
        self.range[1]
    }

    /// OCCT Value(U).
    pub fn value(&self, u: f64) -> DVec3 {
        match self.c3d.as_ref() {
            Some(c) => c.point_at(u),
            None => DVec3::ZERO,
        }
    }

    /// OCCT D1(U, P, V).
    pub fn d1(&self, u: f64) -> (DVec3, DVec3) {
        match self.c3d.as_ref() {
            Some(c) => (c.point_at(u), c.derivative_at(u)),
            None => (DVec3::ZERO, DVec3::ZERO),
        }
    }

    /// OCCT Line() — the line geometry (gp_Lin over Line3).
    pub fn line(&self) -> Option<&Line3> {
        match self.c3d.as_ref() {
            Some(Curve3::Line(l)) => Some(l),
            _ => None,
        }
    }

    /// OCCT Circle() — the circle geometry.
    pub fn circle(&self) -> Option<&Circle3> {
        match self.c3d.as_ref() {
            Some(Curve3::Circle(c)) => Some(c),
            _ => None,
        }
    }

    /// OCCT NbPoles — the pole count of the BSpline/Bezier basis.
    pub fn nb_poles(&self) -> usize {
        match self.c3d.as_ref() {
            Some(Curve3::BSpline(b)) => b.control_points.len(),
            Some(Curve3::Bezier(b)) => b.control_points.len(),
            _ => 0,
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT BRepAdaptor_Curve2d (TKTopAlgo/BRepAdaptor) — reduced re-host over
// BRep_Tool::CurveOnSurface.
// ---------------------------------------------------------------------------

/// OCCT `BRepAdaptor_Curve2d` reduced to the members UnifySameDomain
/// consumes (Curve, FirstParameter, LastParameter, Value, D1).  The pcurve
/// comes from the kernel BRep_Tool::CurveOnSurface port (seam-side aware).
pub struct BRepAdaptorCurve2d {
    pcurve: Option<Curve2d>,
    range: [f64; 2],
}

impl BRepAdaptorCurve2d {
    /// OCCT BRepAdaptor_Curve2d(E, F) (BRepAdaptor_Curve2d.cxx L69-110).
    pub fn new(brep: &BRep, e: &Shape, f: &Shape) -> Self {
        let (pcurve, range) = match brep.curve_on_surface(e, f) {
            Some((pc, first, last)) => (Some(pc), [first, last]),
            None => (None, [0.0, 0.0]),
        };
        BRepAdaptorCurve2d { pcurve, range }
    }

    /// OCCT Curve() — null when the edge carries no pcurve on the face
    /// (the IsNull checks at cxx L358 and L3839).
    pub fn curve(&self) -> Option<&Curve2d> {
        self.pcurve.as_ref()
    }

    /// OCCT FirstParameter.
    pub fn first_parameter(&self) -> f64 {
        self.range[0]
    }

    /// OCCT LastParameter.
    pub fn last_parameter(&self) -> f64 {
        self.range[1]
    }

    /// OCCT Value(U).
    pub fn value(&self, u: f64) -> DVec2 {
        match self.pcurve.as_ref() {
            Some(pc) => pc.point_at(u),
            None => DVec2::ZERO,
        }
    }

    /// OCCT D1(U, P, V).
    pub fn d1(&self, u: f64) -> (DVec2, DVec2) {
        match self.pcurve.as_ref() {
            Some(pc) => (pc.point_at(u), pc.derivative_at(u)),
            None => (DVec2::ZERO, DVec2::ZERO),
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT GeomConvert / Geom2dConvert (TKMath) — GAP carriers.  The C1
// concatenation family (ConcatC1 / C0BSplineToC1BSplineCurve /
// CompCurveToBSplineCurve) has no rcad port; the carriers keep the
// effective contract of the three GlueEdgesWith3DCurves / UnionPCurves call
// shapes (convert -> C1 raise -> concatenate into one curve).
// ---------------------------------------------------------------------------

/// OCCT GeomConvert::CurveToBSplineCurve(T) (GeomConvert.cxx L42-92) — GAP
/// carrier: exact conversions for the analytic bases the chains carry
/// (Line/BSpline), sampling conversion for the rest.
pub fn geom_convert_curve_to_bspline(c: &Curve3, first: f64, last: f64) -> BSplineCurve3 {
    let basis = match c {
        Curve3::Trimmed(t) => Some((*t.curve).clone()),
        _ => None,
    };
    match c {
        Curve3::Line(l) => math_convert::line_to_bspline_range(l, first, last),
        Curve3::BSpline(b) => b.clone(),
        _ => {
            let trimmed = Curve3::Trimmed(TrimmedCurve3::new(
                basis.unwrap_or_else(|| c.clone()),
                first,
                last,
            ));
            math_convert::curve_to_bspline(&trimmed, 33)
        }
    }
}

/// OCCT Geom2dConvert::CurveToBSplineCurve(T) (Geom2dConvert.cxx) — GAP
/// carrier (see [`geom_convert_curve_to_bspline`]).
pub fn geom2d_convert_curve_to_bspline(c: &Curve2d, first: f64, last: f64) -> BSplineCurve2 {
    let basis = match c {
        Curve2d::Trimmed(t) => Some((*t.curve).clone()),
        _ => None,
    };
    match c {
        Curve2d::Line(l) => {
            // The linear BSpline over [first, last] (degree 1, two poles).
            BSplineCurve2 {
                degree: 1,
                knots: vec![first, first, last, last],
                control_points: vec![l.point_at(first), l.point_at(last)],
                weights: vec![1.0, 1.0],
            }
        }
        Curve2d::BSpline(b) => b.clone(),
        _ => {
            // Sampling approximation over the trimmed range (the
            // BSplineCurve2::approximate fit; documented reduction).
            let n = 33usize;
            let mut pts = Vec::with_capacity(n);
            for i in 0..n {
                let t = first + (last - first) * (i as f64) / ((n - 1) as f64);
                pts.push(c.point_at(t));
            }
            let _ = basis;
            BSplineCurve2::approximate(&pts)
        }
    }
}

/// OCCT GeomConvert::C0BSplineToC1BSplineCurve(C, Confusion())
/// (GeomConvert.cxx) — GAP carrier: no-op.  The converter outputs of the
/// carried chains carry no C0 junctions; the C1 raise is TKMath scope.
pub fn geom_convert_c0_to_c1(_c: &mut BSplineCurve3) {}

/// OCCT Geom2dConvert::C0BSplineToC1BSplineCurve — see
/// [`geom_convert_c0_to_c1`].
pub fn geom2d_convert_c0_to_c1(_c: &mut BSplineCurve2) {}

/// OCCT GeomConvert::ConcatC1(tab_c, tabtolvertex, ArrayOfInd, tab_c, flag,
/// prec) (GeomConvert.cxx) — GAP carrier: concatenates the chain into a
/// single BSpline by pole/knot merging (same degree), returning the
/// concatenated list (OCCT returns 1..n curves; the >1 tail is merged again
/// through CompCurveToBSplineCurve at GlueEdgesWith3DCurves cxx L1727-1737,
/// so the single-curve result reproduces the effective outcome).
pub fn geom_convert_concat_c1(tab_c: &[BSplineCurve3], _tol: f64) -> Vec<BSplineCurve3> {
    if tab_c.len() <= 1 {
        return tab_c.to_vec();
    }
    let mut acc = tab_c[0].clone();
    for nxt in &tab_c[1..] {
        match bspline3_concat(&acc, nxt) {
            Some(m) => acc = m,
            None => return tab_c.to_vec(),
        }
    }
    vec![acc]
}

/// OCCT Geom2dConvert::ConcatC1 — see [`geom_convert_concat_c1`].
pub fn geom2d_convert_concat_c1(tab_c: &[BSplineCurve2], _tol: f64) -> Vec<BSplineCurve2> {
    if tab_c.len() <= 1 {
        return tab_c.to_vec();
    }
    let mut acc = tab_c[0].clone();
    for nxt in &tab_c[1..] {
        match bspline2_concat(&acc, nxt) {
            Some(m) => acc = m,
            None => return tab_c.to_vec(),
        }
    }
    vec![acc]
}

/// OCCT GeomConvert_CompCurveToBSplineCurve::Add(B, Tol, After) chains —
/// GAP carrier folding the tail curves into the head (the
/// GlueEdgesWith3DCurves cxx L1729-1736 / UnionPCurves cxx L2033-2040
/// fallbacks).
pub fn geom_convert_comp_curve_add(head: &mut BSplineCurve3, next: &BSplineCurve3, _tol: f64) {
    if let Some(m) = bspline3_concat(head, next) {
        *head = m;
    }
}

/// OCCT Geom2dConvert_CompCurveToBSplineCurve::Add — see
/// [`geom_convert_comp_curve_add`].
pub fn geom2d_convert_comp_curve_add(head: &mut BSplineCurve2, next: &BSplineCurve2, _tol: f64) {
    if let Some(m) = bspline2_concat(head, next) {
        *head = m;
    }
}

/// Pole/knot merge of two clamped BSplines of equal degree (the carrier
/// concatenation primitive; parameters continue across the junction).
fn bspline3_concat(a: &BSplineCurve3, b: &BSplineCurve3) -> Option<BSplineCurve3> {
    if a.degree != b.degree || a.degree + 1 > a.knots.len() || b.degree + 1 > b.knots.len() {
        return None;
    }
    let a_last = a.knots[a.knots.len() - 1 - a.degree];
    let b_first = b.knots[b.degree];
    let mut knots = a.knots.clone();
    let shift = a_last - b_first;
    for k in b.knots.iter().skip(b.degree + 1) {
        knots.push(k + shift);
    }
    let mut control_points = a.control_points.clone();
    control_points.extend(b.control_points.iter().cloned());
    let mut weights = a.weights.clone();
    weights.extend(b.weights.iter().cloned());
    Some(BSplineCurve3 {
        degree: a.degree,
        knots,
        control_points,
        weights,
        is_periodic: false,
    })
}

/// 2D variant of [`bspline3_concat`].
fn bspline2_concat(a: &BSplineCurve2, b: &BSplineCurve2) -> Option<BSplineCurve2> {
    if a.degree != b.degree || a.degree + 1 > a.knots.len() || b.degree + 1 > b.knots.len() {
        return None;
    }
    let a_last = a.knots[a.knots.len() - 1 - a.degree];
    let b_first = b.knots[b.degree];
    let mut knots = a.knots.clone();
    let shift = a_last - b_first;
    for k in b.knots.iter().skip(b.degree + 1) {
        knots.push(k + shift);
    }
    let mut control_points = a.control_points.clone();
    control_points.extend(b.control_points.iter().cloned());
    let mut weights = a.weights.clone();
    weights.extend(b.weights.iter().cloned());
    Some(BSplineCurve2 {
        degree: a.degree,
        knots,
        control_points,
        weights,
    })
}

// ---------------------------------------------------------------------------
// OCCT BRepLib statics (BRepLib_1.cxx) — the docket section 4 gap 3 family.
// ---------------------------------------------------------------------------

/// OCCT BRepLib::BuildPCurveForEdgeOnPlane(E, F) (BRepLib_1.cxx L246-311):
/// for an edge on a planar face without a stored pcurve, computes the
/// on-the-fly plane projection (BRep_Tool::CurveOnPlane) and stores it.
pub fn brep_lib_build_pcurve_for_edge_on_plane(brep: &mut BRep, edge: &Shape, face: &Shape) {
    // OCCT L248-250: the pcurve is computed only when absent.
    let key_loc = brep.compose_pcurve_location(face.location, edge.location);
    let key = (face.ptr_id(), key_loc);
    let already_stored = match &*brep.tshapes[edge.index] {
        TShape::Edge(ed) => ed.pcurves.contains_key(&key) || ed.representations.iter().any(|r| {
            matches!(
                r,
                rcad_kernel::topods::CurveRepresentation::CurveOnSurface { face: f, .. }
                    | rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface { face: f, .. }
                    if *f == key
            )
        }),
        _ => false,
    };
    if already_stored {
        return;
    }
    // OCCT L252-303: the plane projection; the kernel CurveOnSurface port
    // falls back to CurveOnPlane exactly when nothing is stored
    // (BRep_Tool.cxx L367-372).
    if let Some((pc, _first, _last)) = brep.curve_on_surface(edge, face) {
        let mut builder = BRepBuilder::new();
        builder.update_edge_pcurve(brep, edge.clone(), pc, face.clone(), brep.tolerance(edge));
    }
}

/// OCCT BRepLib::BuildPCurveForEdgesOnPlane(const TopTools_ListOfShape& LE,
/// const TopoDS_Face& F) (BRepLib_1.cxx L313-339): the list walk over
/// BuildPCurveForEdgeOnPlane.
pub fn brep_lib_build_pcurve_for_edges_on_plane(brep: &mut BRep, le: &[Shape], face: &Shape) {
    for le_it in le {
        brep_lib_build_pcurve_for_edge_on_plane(brep, le_it, face);
    }
}

/// OCCT BRepLib::ContinuityOfFaces(E, F1, F2, Tol) — GAP carrier: the
/// face-face continuity evaluation has no rcad port; the carrier keeps
/// OCCT's not-continuous path (GeomAbs_C0), which drives the
/// `anOrderOfCont >= GeomAbs_G1` smoothness test at cxx L3595 to false.
pub fn brep_lib_continuity_of_faces(
    _brep: &BRep,
    _edge: &Shape,
    _f1: &Shape,
    _f2: &Shape,
    _tol: f64,
) -> GeomAbsShape {
    GeomAbsShape::C0
}

/// OCCT BRepLib_MakeEdge(Curve, V1, V2, p1, p2) (BRepLib_MakeEdge.cxx
/// L287-330) — the architecture re-host over add_tedge + UpdateVertex: the
/// edge carries the curve over [p1, p2] with the given extremity vertices.
pub fn brep_lib_make_edge(
    brep: &mut BRep,
    curve: &Curve3,
    v1: &Shape,
    v2: &Shape,
    p1: f64,
    p2: f64,
) -> Shape {
    let mut builder = BRepBuilder::new();
    let e = builder.add_edge(brep, Some(curve.clone()), v1.clone(), v2.clone(), [p1, p2]);
    // OCCT: the vertices get their parameters on the edge.
    builder.update_vertex_on_edge(brep, v1.clone(), p1, e.clone(), CONFUSION);
    builder.update_vertex_on_edge(brep, v2.clone(), p2, e.clone(), CONFUSION);
    e
}

// ---------------------------------------------------------------------------
// OCCT GC_MakeCircle / GC_MakeLine2d (TKTopAlgo/GC) — pure-math re-hosts.
// ---------------------------------------------------------------------------

/// OCCT GC_MakeCircle(P1, P2, P3) (GC_MakeCircle.cxx L58-77): the circle
/// through three points (gce_MakeCirc).  None keeps OCCT's NotDone path.
pub fn gc_make_circle_3_points(p1: DVec3, p2: DVec3, p3: DVec3) -> Option<Circle3> {
    // Circumcenter: the standard perpendicular-bisector solution.
    let u = p2 - p1;
    let v = p3 - p1;
    let n = u.cross(v);
    if n.length_squared() < 1e-24 {
        return None;
    }
    let offset = (u.length_squared() * v.cross(n) + v.length_squared() * n.cross(u))
        / (2.0 * n.length_squared());
    let center = p1 + offset;
    let radius = (p1 - center).length();
    if radius <= 0.0 {
        return None;
    }
    Some(Circle3::new(center, n.normalize_or_zero(), radius))
}

/// OCCT GC_MakeLine2d(P1, P2) (GC_MakeLine2d.cxx): the 2D line.
pub fn gc_make_line2d(p1: DVec2, p2: DVec2) -> Line2d {
    Line2d::new(p1, p2 - p1)
}

// ---------------------------------------------------------------------------
// gp / ElCLib pure-math re-hosts.
// ---------------------------------------------------------------------------

/// OCCT gp_Dir::Angle(Other) (gp_Dir.cxx): the angle in [0, PI] — 2D.
pub fn dir_angle_2d(a: DVec2, b: DVec2) -> f64 {
    let mut ang = (a.dot(b) / (a.length() * b.length()).max(1e-300))
        .clamp(-1.0, 1.0)
        .acos();
    if ang.is_nan() {
        ang = 0.0;
    }
    ang
}

/// OCCT gp_Dir::Angle(Other) — the 3D angle in [0, PI].
pub fn dir_angle_3d(a: DVec3, b: DVec3) -> f64 {
    let cos = (a.dot(b) / (a.length() * b.length()).max(1e-300)).clamp(-1.0, 1.0);
    let mut ang = cos.acos();
    if ang.is_nan() {
        ang = 0.0;
    }
    ang
}

/// OCCT gp_Dir::AngleWithRef(Other, V) — the 2D signed angle around the
/// reference direction: atan2((A x B) . V, A . B).
pub fn dir_angle_with_ref_2d(a: DVec2, b: DVec2, ref_dir: DVec2) -> f64 {
    let sin = ref_dir.dot(DVec2::new(a.x * b.y - a.y * b.x, a.y * b.x - a.x * b.y));
    let cos = (a.dot(b) / (a.length() * b.length()).max(1e-300)).clamp(-1.0, 1.0);
    sin.atan2(cos)
}

/// OCCT gp_Dir::AngleWithRef(Other, V) — the 3D form:
/// atan2((A x B) . V, A . B).
pub fn dir_angle_with_ref_3d(a: DVec3, b: DVec3, ref_dir: DVec3) -> f64 {
    let sin = ref_dir.dot(a.cross(b));
    let cos = (a.dot(b) / (a.length() * b.length()).max(1e-300)).clamp(-1.0, 1.0);
    sin.atan2(cos)
}

/// OCCT gp_Dir::IsParallel(Other, AngularTolerance) (gp_Dir.cxx) — 2D.
pub fn dir_is_parallel_2d(a: DVec2, b: DVec2, tol: f64) -> bool {
    let ang = dir_angle_2d(a, b);
    ang <= tol || (std::f64::consts::PI - ang) <= tol
}

/// OCCT gp_Dir::IsParallel — the 3D form.
pub fn dir_is_parallel_3d(a: DVec3, b: DVec3, tol: f64) -> bool {
    let ang = dir_angle_3d(a, b);
    ang <= tol || (std::f64::consts::PI - ang) <= tol
}

/// OCCT ElCLib::LineParameter(Axis, P) (ElCLib.cxx): the projection
/// parameter onto the axis direction.
pub fn elclib_line_parameter_3d(dir: DVec3, loc: DVec3, p: DVec3) -> f64 {
    (p - loc).dot(dir)
}

/// OCCT ElCLib::Parameter(gp_Lin2d, P) (ElCLib.cxx).
pub fn elclib_parameter_lin2d(origin: DVec2, dir: DVec2, p: DVec2) -> f64 {
    (p - origin).dot(dir)
}

/// OCCT ElCLib::Parameter(gp_Circ2d, P) (ElCLib.cxx) — the angle parameter
/// from the X axis of the circle frame; rcad Circle2d carries no frame, so
/// the default X-direction frame is used (documented reduction).
pub fn elclib_parameter_circ2d(center: DVec2, p: DVec2) -> f64 {
    (p - center).y.atan2((p - center).x)
}

/// OCCT gp_Lin2d::Contains(P, Tolerance) (gp_Lin2d.hxx): the distance of the
/// point to the line within the tolerance.
pub fn gp_lin2d_contains(origin: DVec2, dir: DVec2, p: DVec2, tol: f64) -> bool {
    let v = p - origin;
    let t = v.dot(dir);
    let proj = origin + dir * t;
    p.distance(proj) <= tol
}

// ---------------------------------------------------------------------------
// Geom2d curve mutation re-hosts (Geom2d_Curve::Translate / Mirror).
// ---------------------------------------------------------------------------

/// OCCT Geom2d_Curve::Translate(V) — the value-type re-host over the rcad
/// Curve2d variants (the pole/origin/center shift; exotic types keep the
/// unmodified copy, documented GAP).
pub fn geom2d_translate(c: &Curve2d, t: DVec2) -> Curve2d {
    match c {
        Curve2d::Line(l) => Curve2d::Line(Line2d {
            origin: l.origin + t,
            direction: l.direction,
        }),
        Curve2d::Circle(ci) => Curve2d::Circle(Circle2d {
            center: ci.center + t,
            x_dir: ci.x_dir,
            y_dir: ci.y_dir,
            radius: ci.radius,
        }),
        Curve2d::BSpline(b) => Curve2d::BSpline(BSplineCurve2 {
            degree: b.degree,
            knots: b.knots.clone(),
            control_points: b.control_points.iter().map(|p| *p + t).collect(),
            weights: b.weights.clone(),
        }),
        Curve2d::Bezier(b) => Curve2d::Bezier(BezierCurve2 {
            control_points: b.control_points.iter().map(|p| *p + t).collect(),
            weights: b.weights.clone(),
        }),
        Curve2d::Trimmed(tr) => Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(geom2d_translate(&tr.curve, t)),
            t_min: tr.t_min,
            t_max: tr.t_max,
        }),
        other => other.clone(),
    }
}

/// OCCT Geom2d_Curve::Mirror(gp::OX2d()) — the mirroring about the X axis
/// ((x, y) -> (x, -y)); see [`geom2d_translate`] for the type coverage.
pub fn geom2d_mirror_ox2d(c: &Curve2d) -> Curve2d {
    let mir = |p: DVec2| DVec2::new(p.x, -p.y);
    match c {
        Curve2d::Line(l) => Curve2d::Line(Line2d {
            origin: mir(l.origin),
            direction: mir(l.direction),
        }),
        Curve2d::Circle(ci) => Curve2d::Circle(Circle2d {
            center: mir(ci.center),
            x_dir: mir(ci.x_dir),
            y_dir: mir(ci.y_dir),
            radius: ci.radius,
        }),
        Curve2d::BSpline(b) => Curve2d::BSpline(BSplineCurve2 {
            degree: b.degree,
            knots: b.knots.clone(),
            control_points: b.control_points.iter().map(|p| mir(*p)).collect(),
            weights: b.weights.clone(),
        }),
        Curve2d::Bezier(b) => Curve2d::Bezier(BezierCurve2 {
            control_points: b.control_points.iter().map(|p| mir(*p)).collect(),
            weights: b.weights.clone(),
        }),
        Curve2d::Trimmed(tr) => Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(geom2d_mirror_ox2d(&tr.curve)),
            t_min: tr.t_min,
            t_max: tr.t_max,
        }),
        other => other.clone(),
    }
}

/// OCCT Geom2d_Curve::Mirror(gp::OY2d()) — the mirroring about the Y axis
/// ((x, y) -> (-x, y)); see [`geom2d_translate`] for the type coverage.
pub fn geom2d_mirror_oy2d(c: &Curve2d) -> Curve2d {
    let mir = |p: DVec2| DVec2::new(-p.x, p.y);
    match c {
        Curve2d::Line(l) => Curve2d::Line(Line2d {
            origin: mir(l.origin),
            direction: mir(l.direction),
        }),
        Curve2d::Circle(ci) => Curve2d::Circle(Circle2d {
            center: mir(ci.center),
            x_dir: mir(ci.x_dir),
            y_dir: mir(ci.y_dir),
            radius: ci.radius,
        }),
        Curve2d::BSpline(b) => Curve2d::BSpline(BSplineCurve2 {
            degree: b.degree,
            knots: b.knots.clone(),
            control_points: b.control_points.iter().map(|p| mir(*p)).collect(),
            weights: b.weights.clone(),
        }),
        Curve2d::Bezier(b) => Curve2d::Bezier(BezierCurve2 {
            control_points: b.control_points.iter().map(|p| mir(*p)).collect(),
            weights: b.weights.clone(),
        }),
        Curve2d::Trimmed(tr) => Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(geom2d_mirror_oy2d(&tr.curve)),
            t_min: tr.t_min,
            t_max: tr.t_max,
        }),
        other => other.clone(),
    }
}

/// OCCT Geom2d_Curve::Copy() unwrapped to the Geom2d_TrimmedCurve basis —
/// the Copy + IsKind(Geom2d_TrimmedCurve) + BasisCurve chain (cxx
/// L1873-1877, L1958-1962).
pub fn geom2d_copy_untrim(c: &Curve2d) -> Curve2d {
    match c {
        Curve2d::Trimmed(t) => (*t.curve).clone(),
        other => other.clone(),
    }
}

/// OCCT GeomConvert_ApproxSurface(Surface, Tol, C1, C1, 14, 14, 16, 1)
/// (cxx L3629-3631) — GAP carrier: the kernel analytic-to-BSpline
/// conversion (`surface_to_bspline`); the tolerance/continuity controlled
/// approximation is TKMath scope.
pub fn geom_convert_approx_surface(s: &Surface3) -> geom::BSplineSurface {
    math_convert::surface_to_bspline(s, 16, 16)
}

/// OCCT Geom_BSplineSurface::SetUPeriodic/SetVPeriodic — GAP carrier: the
/// rcad BSplineSurface value model carries no periodic flag, so the
/// periodicity re-read falls back to the bounds span (the Uperiod/Vperiod
/// updates at cxx L3634-3643 read the bounds).
pub fn geom_bspline_surface_period_span_u(s: &geom::BSplineSurface) -> f64 {
    let n = s.knots_u.len();
    if n == 0 {
        return 0.0;
    }
    s.knots_u[n - 1] - s.knots_u[0]
}

/// The V span — see [`geom_bspline_surface_period_span_u`].
pub fn geom_bspline_surface_period_span_v(s: &geom::BSplineSurface) -> f64 {
    let n = s.knots_v.len();
    if n == 0 {
        return 0.0;
    }
    s.knots_v[n - 1] - s.knots_v[0]
}

/// OCCT `Geom_TrimmedCurve` construction (the new Geom_TrimmedCurve(C, f, l)
/// chain) over the rcad value model.
pub fn geom_trimmed_3d(c: &Curve3, first: f64, last: f64) -> Curve3 {
    match c {
        Curve3::Trimmed(t) => {
            let b = (*t.curve).clone();
            Curve3::Trimmed(TrimmedCurve3::new(b, first, last))
        }
        other => Curve3::Trimmed(TrimmedCurve3::new(other.clone(), first, last)),
    }
}

/// OCCT `Geom2d_TrimmedCurve` construction — see [`geom_trimmed_3d`].
pub fn geom2d_trimmed(c: &Curve2d, first: f64, last: f64) -> Curve2d {
    match c {
        Curve2d::Trimmed(t) => {
            let b = (*t.curve).clone();
            Curve2d::Trimmed(TrimmedCurve2 {
                curve: Box::new(b),
                t_min: first,
                t_max: last,
            })
        }
        other => Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(other.clone()),
            t_min: first,
            t_max: last,
        }),
    }
}

/// The pcurve handle-identity stand-in: OCCT compares Geom2d_Curve HANDLES
/// (cxx L1175, L1187, L1849, L1891); rcad curves travel as values, so the
/// structural Debug image is the identity proxy (the same fetched row
/// prints identically; distinct rows differ).
pub fn pcurve_handle_same(a: &Curve2d, b: &Curve2d) -> bool {
    format!("{:?}", a) == format!("{:?}", b)
}

/// The Geom_Surface handle-identity proxy (the `S1 == S2` form of cxx
/// L1531; the structural Debug image over the rcad value model).
pub fn surface_handle_same(a: &Surface3, b: &Surface3) -> bool {
    format!("{:?}", a) == format!("{:?}", b)
}

/// OCCT Precision::IsNegativeInfinite / IsPositiveInfinite re-hosts over the
/// kernel constants.
pub fn precision_is_neg_inf(v: f64) -> bool {
    v <= -rcad_kernel::precision::INFINITE_VALUE
}
pub fn precision_is_pos_inf(v: f64) -> bool {
    v >= rcad_kernel::precision::INFINITE_VALUE
}
