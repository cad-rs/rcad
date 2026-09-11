//! OCCT ChFi2d_Builder — 1:1 translation.
//!
//! Sources (TKFillet/ChFi2d):
//!   - ChFi2d_Builder.hxx L41-281 (class declaration + members)
//!   - ChFi2d_Builder.lxx L26-87 (inline methods)
//!   - ChFi2d_Builder.cxx L57-1289 (constructors, Init, AddFillet,
//!     ModifyFillet, RemoveFillet, ComputeFillet, BuildNewWire,
//!     BuildNewEdge, UpDateHistory, BasisEdge, BuildFilletEdge,
//!     IsAFillet, IsAChamfer, IsIssuedFrom, IsLineOrCircle)
//!
//! OCCT class hierarchy: ChFi2d_Builder is a standalone class (no base).
//! Architecture differences (rcad <-> OCCT):
//!   - OCCT carries a global TShape graph; rcad shapes live in a
//!     [`BRep`] table, so the builder owns `my_brep` in addition to the
//!     OCCT member set.
//!   - `history` (NCollection_DataMap keyed by TopTools_ShapeMapHasher,
//!     i.e. TShape+Location+Orientation identity) is keyed here by the
//!     Shape TShape pointer (`ptr_id`) — the rcad stand-in for shape
//!     identity (same convention as ChFi3d's myEdgeFirstFace map).
//!   - OCCT `Geom_Curve` handle identity (`c1 == c2`) is replaced by
//!     value comparison (`geom_curve_same`), mirroring the rcad
//!     `surface_same` convention for Geom_Surface handles.
//!   - OCCT const methods that build edges (BuildNewEdge) take
//!     `&mut self` in rcad: the BRep pool mutation requires it.

use std::collections::HashMap;

use glam::{DVec2, DVec3};
use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_ext_pc::ExtremaExtPC;
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor;
use rcad_kernel::core::precision::{ANGULAR, CONFUSION, PCONFUSION};
use rcad_kernel::geom::{Circle2d, Curve2d, Curve3, CurveEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topo::topods::{
    tshape_flags, BRep, BRepBuilder, BRepTool as _, Orientation, TShape,
};

use crate::geomalgo::geom2d_gcc::{Circ2d2TanRad, GccEntPosition, QualifiedCurve};
use crate::geomalgo::geom2d_int::elclib2d;

// ---------------------------------------------------------------------------
// Cross-file types translated by the parallel Stage 1a agents (pending
// integration — these imports compile once `fillet/chfi2d.rs` lands):
//   - ChFi2dConstructionError  (OCCT ChFi2d_ConstructionError.hxx L28-44)
//   - common_vertex            (OCCT ChFi2d.hxx L47-49, ChFi2d.cxx)
//   - find_connected_edges     (OCCT ChFi2d.hxx L51-54, ChFi2d.cxx)
// ---------------------------------------------------------------------------
use super::chfi2d::{
    chfi2d_common_vertex, chfi2d_find_connected_edges, ChFi2dConstructionError,
};
use super::chfi2d_builder_0::brep_tool_parameter;

/// OCCT ChFi2d_Builder (ChFi2d_Builder.hxx L41-281): builds fillets and
/// chamfers on the wires of a planar face.
#[derive(Debug, Clone)]
pub struct ChFi2dBuilder {
    /// OCCT: ChFi2d_ConstructionError status (hxx L275)
    pub status: ChFi2dConstructionError,
    /// OCCT: TopoDS_Face refFace (hxx L276)
    pub ref_face: Shape,
    /// OCCT: TopoDS_Face newFace (hxx L277)
    pub new_face: Shape,
    /// OCCT: NCollection_Sequence<TopoDS_Shape> fillets (hxx L278)
    pub fillets: Vec<Shape>,
    /// OCCT: NCollection_Sequence<TopoDS_Shape> chamfers (hxx L279)
    pub chamfers: Vec<Shape>,
    /// OCCT: NCollection_DataMap<TopoDS_Shape, TopoDS_Shape,
    /// TopTools_ShapeMapHasher> history (hxx L280). Keyed by the key
    /// shape's ptr_id (architecture note in the module docs); the key
    /// Shape is kept alongside the value because OCCT BasisEdge returns
    /// the key edge to the caller.
    pub history: HashMap<u64, (Shape, Shape)>,
    /// The BRep the shapes belong to (rcad architecture: TopoDS_Shape
    /// lives inside a BRep TShape table; OCCT has global handle graphs).
    pub my_brep: BRep,
}

// ===========================================================================
// File-level helpers (OCCT ChFi2d_Builder.cxx)
// ===========================================================================

/// OCCT BRep_Builder::MakeVertex(V, P, Tol)
/// (BRepLib_MakeVertex.cxx: `B.MakeVertex(V, Point(P), Tol)`) — a fresh
/// vertex TShape, never a position-cache lookup.
pub(crate) fn brep_builder_make_vertex(brep: &mut BRep, p: DVec3, tol: f64) -> Shape {
    let v = brep.add_tvertex_unique(p);
    brep.vertex_mut(v.clone()).tolerance = tol;
    v
}

/// OCCT BRepLib::BuildCurves3d(S) — computes the 3D curves of the edges
/// from their pcurves. rcad architecture gap: edges carry the pcurve-only
/// representation natively and the 3D-curve regeneration pass is not yet
/// translated (TKBRep BRepLib::BuildCurves3d); the call site is kept with
/// a no-op body so the OCCT form stays visible.
pub(crate) fn brep_lib_build_curves3d(_brep: &mut BRep, _s: &Shape) {
    // OCCT BRepLib::BuildCurves3d — not yet translated (rcad gap).
}

/// OCCT TopExp_Explorer(F, TopAbs_WIRE) — the wires of a face in storage
/// order. rcad architecture: TFaceData carries outer_wire + inner_wires;
/// the explorer walks the outer wire first, then the inner wires.
pub(crate) fn topexp_explore_face_wires(brep: &BRep, f: &Shape) -> Vec<Shape> {
    let fd = brep.face(f.clone());
    let mut wires = Vec::new();
    if !fd.outer_wire.is_null() {
        wires.push(fd.outer_wire.clone());
    }
    wires.extend(fd.inner_wires.iter().cloned());
    wires
}

/// OCCT TopExp_Explorer(F, TopAbs_EDGE) — the edges of the face through
/// its wires (FACE -> WIRE -> EDGE depth-first order, every occurrence).
pub(crate) fn topexp_explore_face_edges(brep: &BRep, f: &Shape) -> Vec<Shape> {
    let mut edges = Vec::new();
    for w in topexp_explore_face_wires(brep, f) {
        if let TShape::Wire(wd) = brep.tshapes[w.index].as_ref() {
            edges.extend(wd.edges.iter().cloned());
        }
    }
    edges
}

/// OCCT TopTools_ShapeMapHasher Contains stand-in: TShape identity via
/// ptr_id (the rcad surface_same convention; location/orientation are not
/// part of the key — architecture note in the module docs).
pub(crate) fn shape_map_contains(map: &[Shape], s: &Shape) -> bool {
    map.iter().any(|m| m.ptr_id() == s.ptr_id())
}

/// OCCT Geom_Curve handle identity (`c1 == c2`) stand-in: value
/// comparison with Precision::Confusion (mirrors rcad `surface_same` for
/// Geom_Surface handles). Only the kinds ChFi2d compares are covered
/// (Geom_Line / Geom_Circle).
fn geom_curve_same(a: &Curve3, b: &Curve3) -> bool {
    const T: f64 = CONFUSION;
    let v = |x: DVec3, y: DVec3| (x - y).length() < T;
    match (a, b) {
        (Curve3::Line(a), Curve3::Line(b)) => {
            v(a.origin, b.origin) && v(a.direction, b.direction)
        }
        (Curve3::Circle(a), Curve3::Circle(b)) => {
            v(a.center, b.center)
                && v(a.normal, b.normal)
                && v(a.x_dir, b.x_dir)
                && v(a.y_dir, b.y_dir)
                && (a.radius - b.radius).abs() < T
        }
        _ => false,
    }
}

/// OCCT gp_Dir2d::Angle between two 2D directions
/// (gp_Dir2d.hxx L400-407): atan2(cross, dot).
fn gp_dir2d_angle(a: DVec2, b: DVec2) -> f64 {
    (a.x * b.y - a.y * b.x).atan2(a.dot(b))
}

/// OCCT gp_Vec2d::IsParallel(V, AngularTolerance) via gp_Dir2d::IsParallel
/// (gp_Dir2d.hxx L422-430): |angle| <= tol or PI - |angle| <= tol.
fn gp_vec2d_is_parallel(a: DVec2, b: DVec2, angular_tolerance: f64) -> bool {
    let mut an_ang = gp_dir2d_angle(a, b);
    if an_ang < 0.0 {
        an_ang = -an_ang;
    }
    an_ang <= angular_tolerance || std::f64::consts::PI - an_ang <= angular_tolerance
}

/// OCCT gp_Vec2d::IsOpposite(V, AngularTolerance) via gp_Dir2d::IsOpposite
/// (gp_Dir2d.hxx L410-418): PI - |angle| <= tol.
fn gp_vec2d_is_opposite(a: DVec2, b: DVec2, angular_tolerance: f64) -> bool {
    let mut an_ang = gp_dir2d_angle(a, b);
    if an_ang < 0.0 {
        an_ang = -an_ang;
    }
    std::f64::consts::PI - an_ang <= angular_tolerance
}

/// OCCT ElCLib::AdjustPeriodic(UFirst, ULast, Preci, U1, U2)
/// (ElCLib.cxx L115-151).
pub(crate) fn elclib_adjust_periodic(
    u_first: f64,
    u_last: f64,
    preci: f64,
    u1: &mut f64,
    u2: &mut f64,
) {
    // OCCT ElCLib.cxx L121: Precision::IsInfinite(UFirst) || Precision::IsInfinite(ULast)
    // (Precision.hxx L350-353).
    if rcad_kernel::precision::is_infinite_value(u_first)
        || rcad_kernel::precision::is_infinite_value(u_last)
    {
        *u1 = u_first;
        *u2 = u_last;
        return;
    }
    let a_period = u_last - u_first;
    if a_period < u_last.abs() * f64::EPSILON {
        // In order to avoid FLT_Overflow exception
        // (test bugs moddata_1 bug22757)
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

/// OCCT BRepLib_MakeEdge::Init(C, V1, V2) + Edge() (BRepLib_MakeEdge.cxx
/// L543-567): projects the vertex points onto the curve (Project
/// L73-117 via Extrema_ExtPC), then Init(C, V1, V2, p1, p2). On error
/// (PointProjectionFailed / ParameterOutOfRange / LineThroughIdenticPoints)
/// Edge() is a null edge.
/// rcad architecture: expressed over BRepBuilder::add_edge (the
/// BRep_Builder equivalent).
pub(crate) fn brep_lib_make_edge_init_vertices(
    brep: &mut BRep,
    c: &Curve3,
    vv1: &Shape,
    vv2: &Shape,
) -> Shape {
    // OCCT L546-566: project the vertices on the curve.
    let domain = c.default_domain();
    let cf = domain[0];
    let cl = domain[1];
    let project = |brep: &BRep, v: &Shape| -> Option<f64> {
        // OCCT Project(C, V, p) (BRepLib_MakeEdge.cxx L73-117).
        let p = brep.vertex_position(v);
        let eps2 = brep.vertex_tolerance(v);
        let eps2 = eps2 * eps2;
        // OCCT L78: GeomAdaptor_Curve GAC(C); L80-93: the endpoint
        // shortcut — the closer endpoint wins when within Eps2.
        let p1 = c.point_at(cf);
        let p2 = c.point_at(cl);
        let d1 = p1.distance_squared(p);
        let d2 = p2.distance_squared(p);
        if d1 < d2 && d1 <= eps2 {
            return Some(cf);
        } else if d2 < d1 && d2 <= eps2 {
            return Some(cl);
        }
        // OCCT L95: Extrema_ExtPC extrema(P, GAC) — the two-arg ctor over
        // the full domain, the default theTolF is 1.0e-10.
        let a_adaptor = GeomCurveAdaptor::new(c.clone());
        let a_tool = CurveToolHandle::for_curve3(c, &a_adaptor, &a_adaptor);
        let extrema = ExtremaExtPC::new_point_curve(p, &a_tool, 1.0e-10);
        // OCCT L96-113: the nearest-extremum pick, accepted only when its
        // square distance is within Eps2.
        if extrema.is_done() {
            let mut index: usize = 0;
            let mut dist2 = f64::MAX; // RealLast()
            for i in 1..=extrema.nb_ext() {
                let dist2min = extrema.square_distance(i);
                if dist2min < dist2 {
                    index = i;
                    dist2 = dist2min;
                }
            }
            if index != 0 && dist2 <= eps2 {
                return Some(extrema.point(index).param);
            }
        }
        None
    };
    let p1 = match project(brep, vv1) {
        Some(p) => p,
        None => return Shape::null(), // BRepLib_PointProjectionFailed
    };
    let p2 = match project(brep, vv2) {
        Some(p) => p,
        None => return Shape::null(), // BRepLib_PointProjectionFailed
    };
    brep_lib_make_edge_init_range(brep, c, vv1, vv2, p1, p2)
}

/// OCCT BRepLib_MakeEdge::Init(CC, VV1, VV2, pp1, pp2) + Edge()
/// (BRepLib_MakeEdge.cxx L603-700): the "really makes the job" overload —
/// trims trimmed curves, orders/adjusts the parameters, then builds the
/// edge with the computed range. A null edge is returned on
/// BRepLib_ParameterOutOfRange / BRepLib_LineThroughIdenticPoints.
pub(crate) fn brep_lib_make_edge_init_range(
    brep: &mut BRep,
    cc: &Curve3,
    vv1: &Shape,
    vv2: &Shape,
    pp1: f64,
    pp2: f64,
) -> Shape {
    // OCCT L609-616: kill trimmed curves.
    let mut c = cc.clone();
    loop {
        match &c {
            Curve3::Trimmed(ct) => {
                let basis = (*ct.curve).clone();
                c = basis;
            }
            _ => break,
        }
    }
    // OCCT L619-623: check parameters.
    let mut p1 = pp1;
    let mut p2 = pp2;
    let domain = c.default_domain();
    let cf = domain[0];
    let cl = domain[1];
    let epsilon = PCONFUSION;
    let periodic = c.is_periodic();
    let (mut v1, mut v2) = (vv1.clone(), vv2.clone());
    if periodic {
        // OCCT L629-635: adjust in period.
        elclib_adjust_periodic(cf, cl, epsilon, &mut p1, &mut p2);
    } else {
        // OCCT L640-651: reordonate.
        if p1 < p2 {
            // V1 = VV1; V2 = VV2 (already set)
        } else {
            v2 = vv1.clone();
            v1 = vv2.clone();
            let x = p1;
            p1 = p2;
            p2 = x;
        }
        // OCCT L654-658: check range.
        if (cf - p1 > epsilon) || (p2 - cl > epsilon) {
            return Shape::null(); // BRepLib_ParameterOutOfRange
        }
        // OCCT L661-666: check punctuality.
        if (p2 - p1) <= f64::EPSILON {
            // gp::Resolution()
            return Shape::null(); // BRepLib_LineThroughIdenticPoints
        }
    }
    // OCCT L668-760 (tail of Init): compute the points on the curve, check
    // closedness, and make the edge over [p1, p2] (BRep_Builder::UpdateEdge
    // + Add + UpdateVertex). rcad architecture: BRepBuilder::add_edge
    // carries curve + vertices + range in one step.
    let mut b = BRepBuilder::new();
    b.add_edge(brep, Some(c), v1, v2, [p1, p2])
}

/// OCCT BRepLib_MakeEdge(C, S, V1, V2, p1, p2) (BRepLib_MakeEdge.cxx
/// L907+): an edge carrying the pcurve C on surface S between V1/V2 with
/// the parameter range [p1, p2]. rcad architecture: the pcurve is stored
/// keyed by the owning face (`face`) via BRepBuilder::add_pcurve;
/// BRepTool::curve_on_surface resolves it back by surface value (the
/// surface_same stand-in for handle identity).
pub(crate) fn brep_lib_make_edge_pcurve(
    brep: &mut BRep,
    c: &Curve2d,
    face: &Shape,
    v1: &Shape,
    v2: &Shape,
    p1: f64,
    p2: f64,
) -> Shape {
    let mut b = BRepBuilder::new();
    let e = b.add_edge(brep, None, v1.clone(), v2.clone(), [p1, p2]);
    b.add_pcurve(brep, e.clone(), face.clone(), c.clone(), p1, p2);
    b.set_vertex_param(brep, e.clone(), v1.clone(), p1);
    b.set_vertex_param(brep, e.clone(), v2.clone(), p2);
    e
}

// ===========================================================================
// OCCT ChFi2d_Builder.cxx L57-59 / L210-239 — file-static IsIssuedFrom
// ===========================================================================

/// OCCT ChFi2d_Builder.cxx L210-239 — Search in <map> if <e> has a parent
/// edge. If a parent has been found, this edge is returned in
/// <basis_edge>, else <e> is returned in <basis_edge>.
fn is_issued_from(brep: &BRep, e: &Shape, map: &[Shape], basis_edge: &mut Shape) -> bool {
    // OCCT L219: occ::handle<Geom_Curve> c1 = BRep_Tool::Curve(E, loc1, f1, L1);
    let (c1, f1l1) = match brep.edge_curve_world(e) {
        Some((c, r)) => (c, r),
        None => return false,
    };
    let (f1, l1) = (f1l1[0], f1l1[1]);

    for item in map.iter() {
        let current_edge = item; // OCCT L223: TopoDS::Edge(Map.FindKey(i))
        // OCCT L227: occ::handle<Geom_Curve> c2 = BRep_Tool::Curve(currentEdge, ...);
        let (c2, f2l2) = match brep.edge_curve_world(current_edge) {
            Some((c, r)) => (c, r),
            None => continue,
        };
        let (f2, l2) = (f2l2[0], f2l2[1]);
        // OCCT L228-230:
        // if (c1 == c2 && (((f1 > f2 && f1 < L2) || (L1 > f2 && L1 < L2))
        //                || ((f1 > L2 && f1 < f2) || (L1 > L2 && L1 < f2))))
        if geom_curve_same(&c1, &c2)
            && (((f1 > f2 && f1 < l2) || (l1 > f2 && l1 < l2))
                || ((f1 > l2 && f1 < f2) || (l1 > l2 && l1 < f2)))
        {
            *basis_edge = current_edge.clone();
            basis_edge.orientation = e.orientation;
            return true;
        } // if (c1 == c2
    } // for (int i ...
    false
} // IsIssuedFrom

// ===========================================================================
// OCCT ChFi2d_Builder.cxx L61 / L1268-1289 — file-static IsLineOrCircle
// ===========================================================================

/// OCCT ChFi2d_Builder.cxx L1268-1289 — the pcurve basis of <e> on <f>
/// must be a Geom2d_Circle or a Geom2d_Line.
fn is_line_or_circle(brep: &BRep, e: &Shape, f: &Shape) -> bool {
    // OCCT L1275: occ::handle<Geom2d_Curve> C = BRep_Tool::CurveOnSurface(E, F, first, last);
    let c = match brep.curve_on_surface(e, f) {
        Some((c, _first, _last)) => c,
        None => return false,
    };
    // OCCT L1277-1285: down_cast<Geom2d_TrimmedCurve> -> BasisCurve.
    let basis_c = match &c {
        Curve2d::Trimmed(tc) => (*tc.curve).clone(),
        _ => c,
    };
    matches!(basis_c, Curve2d::Circle(_) | Curve2d::Line(_))
}

impl ChFi2dBuilder {
    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L65-68 — default constructor
    // =======================================================================

    /// OCCT ChFi2d_Builder() : status(ChFi2d_NotPlanar) (L65-68).
    pub fn new() -> Self {
        ChFi2dBuilder {
            status: ChFi2dConstructionError::NotPlanar,
            ref_face: Shape::null(),
            new_face: Shape::null(),
            fillets: Vec::new(),
            chamfers: Vec::new(),
            history: HashMap::new(),
            my_brep: BRep::new(),
        }
    }

    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L72-94 — constructor from a face
    // =======================================================================

    /// OCCT ChFi2d_Builder(const TopoDS_Face& F) (L72-94). The face <f>
    /// can be built on a closed or an open wire.
    /// rcad architecture: the owning BRep is passed alongside the face.
    pub fn new_with_face(brep: &BRep, f: Shape) -> Self {
        let mut b = ChFi2dBuilder::new();
        b.my_brep = brep.clone();
        if f.is_null() {
            b.status = ChFi2dConstructionError::NoFace;
            return b;
        }
        // OCCT L83: if (BRep_Tool::Surface(F, Loc)->IsKind(STANDARD_TYPE(Geom_Plane)))
        // rcad: no RTTI IsKind — planarity is the Surface3::Plane variant test.
        if matches!(brep.face_surface(&f), Some(Surface3::Plane(_))) {
            b.ref_face = f.clone();
            b.new_face = f;
            b.new_face.orientation = Orientation::Forward;
            brep_lib_build_curves3d(&mut b.my_brep, &b.new_face); // OCCT L87
            b.status = ChFi2dConstructionError::Ready;
        } else {
            b.status = ChFi2dConstructionError::NotPlanar;
        }
        b
    } // ChFi2d_Builder

    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L98-122 — Init(F)
    // =======================================================================

    /// OCCT ChFi2d_Builder::Init(const TopoDS_Face& F) (L98-122).
    pub fn init(&mut self, brep: &BRep, f: Shape) {
        if f.is_null() {
            self.status = ChFi2dConstructionError::NoFace;
            return;
        }
        self.fillets.clear(); // OCCT L105: fillets.Clear();
        self.chamfers.clear(); // OCCT L106: chamfers.Clear();
        self.history.clear(); // OCCT L107: history.Clear();
        // OCCT L112: if (BRep_Tool::Surface(F, Loc)->IsKind(STANDARD_TYPE(Geom_Plane)))
        if matches!(brep.face_surface(&f), Some(Surface3::Plane(_))) {
            self.ref_face = f.clone();
            self.new_face = f;
            self.new_face.orientation = Orientation::Forward;
            self.status = ChFi2dConstructionError::Ready;
        } else {
            self.status = ChFi2dConstructionError::NotPlanar;
        }
        self.my_brep = brep.clone();
    } // Init

    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L126-202 — Init(RefFace, ModFace)
    // =======================================================================

    /// OCCT ChFi2d_Builder::Init(const TopoDS_Face& RefFace,
    /// const TopoDS_Face& ModFace) (L126-202). Rust overload name:
    /// `init_modified` (OCCT second Init overload).
    pub fn init_modified(&mut self, brep: &BRep, ref_face: Shape, mod_face: Shape) {
        if ref_face.is_null() || mod_face.is_null() {
            self.status = ChFi2dConstructionError::NoFace;
            return;
        }
        self.fillets.clear();
        self.chamfers.clear();
        self.history.clear();
        // OCCT L140: if (!BRep_Tool::Surface(RefFace, loc)->IsKind(STANDARD_TYPE(Geom_Plane)))
        if !matches!(brep.face_surface(&ref_face), Some(Surface3::Plane(_))) {
            self.status = ChFi2dConstructionError::NotPlanar;
            return;
        }

        self.ref_face = ref_face.clone();
        self.new_face = mod_face;
        self.new_face.orientation = Orientation::Forward;
        self.status = ChFi2dConstructionError::Ready;
        self.my_brep = brep.clone();

        // OCCT L151-155: Research in newFace all the edges which not appear
        // in RefFace. The sequence newEdges will contain these edges.
        let mut new_edges: Vec<Shape> = Vec::new(); // NCollection_Sequence<TopoDS_Shape>
        let mut ref_edges_map: Vec<Shape> = Vec::new(); // NCollection_IndexedMap
        // OCCT L156: TopExp::MapShapes(refFace, TopAbs_EDGE, refEdgesMap);
        for e in topexp_explore_face_edges(&self.my_brep, &self.ref_face) {
            if !shape_map_contains(&ref_edges_map, &e) {
                ref_edges_map.push(e);
            }
        }
        // OCCT L157-166: TopExp_Explorer ex(newFace, TopAbs_EDGE);
        // while (ex.More()) { ... ex.Next(); }
        for e in topexp_explore_face_edges(&self.my_brep, &self.new_face) {
            let current_edge = e; // TopoDS::Edge(ex.Current())
            if !shape_map_contains(&ref_edges_map, &current_edge) {
                new_edges.push(current_edge);
            }
        } // while (ex ...

        // OCCT L168-201: update of history, fillets and chamfers fields
        let mut i = 1usize; // OCCT: int i = 1;
        let mut basis_edge = Shape::null();
        while i <= new_edges.len() {
            let current_edge = new_edges[i - 1].clone(); // TopoDS::Edge(newEdges.Value(i))
            if is_issued_from(&self.my_brep, &current_edge, &ref_edges_map, &mut basis_edge) {
                // OCCT L177: history.Bind(basisEdge, currentEdge);
                self.history
                    .insert(basis_edge.ptr_id(), (basis_edge.clone(), current_edge.clone()));
            } else {
                // this edge is a chamfer or a fillet
                // OCCT L185: occ::handle<Geom_Curve> curve =
                //   BRep_Tool::Curve(currentEdge, loc, first, last);
                let (curve, range) = match self.my_brep.edge_curve_world(&current_edge) {
                    Some(x) => x,
                    None => {
                        self.status = ChFi2dConstructionError::InitialisationError;
                        return;
                    }
                };
                let _first = range[0];
                let _last = range[1];
                if matches!(curve, Curve3::Circle(_)) {
                    self.fillets.push(current_edge);
                } else if matches!(curve, Curve3::Line(_)) {
                    self.chamfers.push(current_edge);
                } else {
                    self.status = ChFi2dConstructionError::InitialisationError;
                    return;
                } // else ...
            } // this edge is ...
            i += 1;
        } // while ...
    } // Init

    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L243-277 — AddFillet
    // =======================================================================

    /// OCCT ChFi2d_Builder::AddFillet(const TopoDS_Vertex& V,
    /// const double Radius) (L243-277).
    pub fn add_fillet(&mut self, v: &Shape, radius: f64) -> Shape {
        let mut adj_edge1 = Shape::null();
        let mut adj_edge2 = Shape::null();
        let mut adj_edge1_mod = Shape::null();
        let mut adj_edge2_mod = Shape::null();
        let mut fillet = Shape::null();
        self.status =
            chfi2d_find_connected_edges(&self.new_face, v, &mut adj_edge1, &mut adj_edge2);
        if self.status == ChFi2dConstructionError::ConnexionError {
            return fillet;
        }

        if self.is_a_fillet(&adj_edge1)
            || self.is_a_chamfer(&adj_edge1)
            || self.is_a_fillet(&adj_edge2)
            || self.is_a_chamfer(&adj_edge2)
        {
            self.status = ChFi2dConstructionError::NotAuthorized;
            return fillet;
        } // if (IsAFillet ...

        if !is_line_or_circle(&self.my_brep, &adj_edge1, &self.new_face)
            || !is_line_or_circle(&self.my_brep, &adj_edge2, &self.new_face)
        {
            self.status = ChFi2dConstructionError::NotAuthorized;
            return fillet;
        } // if (!IsLineOrCircle ...

        self.compute_fillet(
            v,
            &adj_edge1,
            &adj_edge2,
            radius,
            &mut adj_edge1_mod,
            &mut adj_edge2_mod,
            &mut fillet,
        );
        if self.status == ChFi2dConstructionError::IsDone
            || self.status == ChFi2dConstructionError::FirstEdgeDegenerated
            || self.status == ChFi2dConstructionError::LastEdgeDegenerated
            || self.status == ChFi2dConstructionError::BothEdgesDegenerated
        {
            self.build_new_wire(&adj_edge1, &adj_edge2, &adj_edge1_mod, &fillet, &adj_edge2_mod);
            let basis_edge1 = self.basis_edge(&adj_edge1);
            let basis_edge2 = self.basis_edge(&adj_edge2);
            self.up_date_history_new_edge(
                &basis_edge1,
                &basis_edge2,
                &adj_edge1_mod,
                &adj_edge2_mod,
                &fillet,
                1,
            );
            self.status = ChFi2dConstructionError::IsDone;
            return self.fillets[self.fillets.len() - 1].clone();
        }
        fillet
    } // AddFillet

    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L281-286 — ModifyFillet
    // =======================================================================

    /// OCCT ChFi2d_Builder::ModifyFillet(const TopoDS_Edge& Fillet,
    /// const double Radius) (L281-286).
    pub fn modify_fillet(&mut self, fillet: &Shape, radius: f64) -> Shape {
        let a_vertex = self.remove_fillet(fillet);
        let a_fillet = self.add_fillet(&a_vertex, radius);
        a_fillet
    } // ModifyFillet

    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L290-497 — RemoveFillet
    // =======================================================================

    /// OCCT ChFi2d_Builder::RemoveFillet(const TopoDS_Edge& Fillet)
    /// (L290-497).
    pub fn remove_fillet(&mut self, fillet: &Shape) -> Shape {
        let mut common_vertex_shape = Shape::null();
        let mut i = 1usize; // OCCT: int i = 1;
        let mut is_find = false; // OCCT: int IsFind = false;
        while i <= self.fillets.len() {
            let a_fillet = self.fillets[i - 1].clone(); // TopoDS::Edge(fillets.Value(i))
            if a_fillet.is_same(fillet) {
                self.fillets.remove(i - 1);
                is_find = true;
                break;
            }
            i += 1;
        }
        if !is_find {
            return common_vertex_shape;
        }

        let mut first_vertex;
        let mut last_vertex;
        // OCCT L312: TopExp::Vertices(Fillet, firstVertex, lastVertex);
        first_vertex = self.my_brep.first_vertex(fillet);
        last_vertex = self.my_brep.last_vertex(fillet);

        let mut adj_edge1 = Shape::null();
        let mut adj_edge2 = Shape::null();
        self.status = chfi2d_find_connected_edges(
            &self.new_face,
            &first_vertex,
            &mut adj_edge1,
            &mut adj_edge2,
        );
        if self.status == ChFi2dConstructionError::ConnexionError {
            return common_vertex_shape;
        }

        let basis_edge1;
        let basis_edge2;
        let e1;
        let mut e2;
        // E1 and E2 are the adjacent edges to Fillet

        if adj_edge1.is_same(fillet) {
            e1 = adj_edge2.clone();
        } else {
            e1 = adj_edge1.clone();
        }
        basis_edge1 = self.basis_edge(&e1);
        self.status = chfi2d_find_connected_edges(
            &self.new_face,
            &last_vertex,
            &mut adj_edge1,
            &mut adj_edge2,
        );
        if self.status == ChFi2dConstructionError::ConnexionError {
            return common_vertex_shape;
        }
        if adj_edge1.is_same(fillet) {
            e2 = adj_edge2.clone();
        } else {
            e2 = adj_edge1.clone();
        }
        basis_edge2 = self.basis_edge(&e2);
        let mut connection_e1_fillet = Shape::null();
        let mut connection_e2_fillet = Shape::null();
        // OCCT: bool hasConnection = ChFi2d::CommonVertex(basisEdge1, basisEdge2, commonVertex);
        let (cv, mut has_connection) = chfi2d_common_vertex(&basis_edge1, &basis_edge2);
        common_vertex_shape = cv;
        if !has_connection {
            self.status = ChFi2dConstructionError::ConnexionError;
            return common_vertex_shape;
        }
        // OCCT: hasConnection = ChFi2d::CommonVertex(E1, Fillet, connectionE1Fillet);
        let (cv1, hc1) = chfi2d_common_vertex(&e1, fillet);
        connection_e1_fillet = cv1;
        has_connection = hc1;
        if !has_connection {
            self.status = ChFi2dConstructionError::ConnexionError;
            return common_vertex_shape;
        }
        // OCCT: hasConnection = ChFi2d::CommonVertex(E2, Fillet, connectionE2Fillet);
        let (cv2, hc2) = chfi2d_common_vertex(&e2, fillet);
        connection_e2_fillet = cv2;
        has_connection = hc2;
        if !has_connection {
            self.status = ChFi2dConstructionError::ConnexionError;
            return common_vertex_shape;
        }

        // rebuild edges on wire
        let mut new_edge1 = Shape::null();
        let mut new_edge2 = Shape::null();
        let mut v;
        let mut v1;
        let mut v2;

        // OCCT L374: TopExp::Vertices(E1, firstVertex, lastVertex);
        first_vertex = self.my_brep.first_vertex(&e1);
        last_vertex = self.my_brep.last_vertex(&e1);
        // OCCT L375: TopExp::Vertices(basisEdge1, v1, v2);
        v1 = self.my_brep.first_vertex(&basis_edge1);
        v2 = self.my_brep.last_vertex(&basis_edge1);
        if v1.is_same(&common_vertex_shape) {
            v = v2.clone();
        } else {
            v = v1.clone();
        }

        if first_vertex.is_same(&v) || last_vertex.is_same(&v) {
            // It means the edge support only one fillet. In this case
            // the new edge must be the basis edge.
            new_edge1 = basis_edge1.clone();
        } else {
            // It means the edge support one fillet on each end.
            if first_vertex.is_same(&connection_e1_fillet) {
                // OCCT L399: occ::handle<Geom_Curve> curve = BRep_Tool::Curve(E1, loc, first, last);
                let (curve, _range) = match self.my_brep.edge_curve_world(&e1) {
                    Some(x) => x,
                    None => return common_vertex_shape,
                };
                new_edge1 = brep_lib_make_edge_init_vertices(
                    &mut self.my_brep,
                    &curve,
                    &common_vertex_shape,
                    &last_vertex,
                );
                new_edge1.orientation = e1.orientation;
                new_edge1.location = e1.location;
            } // if (firstVertex ...
            else if last_vertex.is_same(&connection_e1_fillet) {
                // OCCT L410: occ::handle<Geom_Curve> curve = BRep_Tool::Curve(E1, loc, first, last);
                let (curve, _range) = match self.my_brep.edge_curve_world(&e1) {
                    Some(x) => x,
                    None => return common_vertex_shape,
                };
                new_edge1 = brep_lib_make_edge_init_vertices(
                    &mut self.my_brep,
                    &curve,
                    &first_vertex,
                    &common_vertex_shape,
                );
                new_edge1.orientation = e1.orientation;
                new_edge1.location = e1.location;
            } // else if (lastVertex ...
        } // else ...

        // OCCT L418: TopExp::Vertices(basisEdge2, v1, v2);
        v1 = self.my_brep.first_vertex(&basis_edge2);
        v2 = self.my_brep.last_vertex(&basis_edge2);
        if v1.is_same(&common_vertex_shape) {
            v = v2.clone();
        } else {
            v = v1.clone();
        }

        // OCCT L428: TopExp::Vertices(E2, firstVertex, lastVertex);
        first_vertex = self.my_brep.first_vertex(&e2);
        last_vertex = self.my_brep.last_vertex(&e2);
        if first_vertex.is_same(&v) || last_vertex.is_same(&v) {
            // It means the edge support only one fillet. In this case
            // the new edge must be the basis edge.
            new_edge2 = basis_edge2.clone();
        } else {
            // It means the edge support one fillet on each end.
            if first_vertex.is_same(&connection_e2_fillet) {
                // OCCT L443: occ::handle<Geom_Curve> curve = BRep_Tool::Curve(E2, loc, first, last);
                let (curve, _range) = match self.my_brep.edge_curve_world(&e2) {
                    Some(x) => x,
                    None => return common_vertex_shape,
                };
                new_edge2 = brep_lib_make_edge_init_vertices(
                    &mut self.my_brep,
                    &curve,
                    &common_vertex_shape,
                    &last_vertex,
                );
                new_edge2.orientation = e2.orientation;
                new_edge2.location = e2.location;
            } // if (firstVertex ...
            else if last_vertex.is_same(&connection_e2_fillet) {
                // OCCT L454: occ::handle<Geom_Curve> curve = BRep_Tool::Curve(E2, loc, first, last);
                let (curve, _range) = match self.my_brep.edge_curve_world(&e2) {
                    Some(x) => x,
                    None => return common_vertex_shape,
                };
                new_edge2 = brep_lib_make_edge_init_vertices(
                    &mut self.my_brep,
                    &curve,
                    &first_vertex,
                    &common_vertex_shape,
                );
                new_edge2.orientation = e2.orientation;
                new_edge2.location = e2.location;
            } // else if (lastVertex ...
        } // else ...

        // rebuild the newFace
        // OCCT L463: TopExp_Explorer Ex(newFace, TopAbs_EDGE);
        let ex = topexp_explore_face_edges(&self.my_brep, &self.new_face);
        // OCCT L466-467: BRep_Builder B; B.MakeWire(newWire);
        let mut b = BRepBuilder::new();
        let new_wire = b.make_wire(&mut self.my_brep);

        for the_edge in ex {
            // OCCT L471: const TopoDS_Edge& theEdge = TopoDS::Edge(Ex.Current());
            if !the_edge.is_same(&e1) && !the_edge.is_same(&e2) && !the_edge.is_same(fillet) {
                b.add_to_wire(&mut self.my_brep, new_wire.clone(), the_edge);
            } else {
                // OCCT L478: if (theEdge == E1) — TopoDS_Shape operator== is
                // IsEqual (TShape+Location+Orientation) -> Shape::is_equal.
                if the_edge.is_equal(&e1) {
                    b.add_to_wire(&mut self.my_brep, new_wire.clone(), new_edge1.clone());
                } else if the_edge.is_equal(&e2) {
                    b.add_to_wire(&mut self.my_brep, new_wire.clone(), new_edge2.clone());
                }
            } // else
        } // while ...
        // OCCT L489-492:
        // BRepAdaptor_Surface Adaptor3dSurface(refFace);
        // BRepLib_MakeFace mFace(Adaptor3dSurface.Plane(), newWire);
        // newFace.Nullify(); newFace = mFace;
        let ref_plane = match self.my_brep.face_surface_world(&self.ref_face) {
            Some(Surface3::Plane(pl)) => pl,
            _ => panic!("Standard_TypeMismatch: refFace is not a plane"),
        };
        let m_face = b.make_face(
            &mut self.my_brep,
            Some(Surface3::Plane(ref_plane)),
            new_wire.clone(),
        );
        self.new_face = m_face;

        self.up_date_history(&basis_edge1, &basis_edge2, &new_edge1, &new_edge2);

        common_vertex_shape
    } // RemoveFillet

    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L501-530 — ComputeFillet
    // =======================================================================

    /// OCCT ChFi2d_Builder::ComputeFillet (L501-530). Is internally used
    /// by <add_fillet>: <trim_e1>, <trim_e2>, <fillet> have sense only if
    /// the status <status> is equal to <IsDone>.
    #[allow(clippy::too_many_arguments)]
    pub fn compute_fillet(
        &mut self,
        v: &Shape,
        e1: &Shape,
        e2: &Shape,
        radius: f64,
        trim_e1: &mut Shape,
        trim_e2: &mut Shape,
        fillet: &mut Shape,
    ) {
        let mut new_extr1 = Shape::null();
        let mut new_extr2 = Shape::null();
        let mut degen1 = false;
        let mut degen2 = false;
        *fillet = self.build_fillet_edge(v, e1, e2, radius, &mut new_extr1, &mut new_extr2);
        if self.status != ChFi2dConstructionError::IsDone {
            return;
        }
        *trim_e1 = self.build_new_edge_degenerated(e1, v, &new_extr1, &mut degen1);
        *trim_e2 = self.build_new_edge_degenerated(e2, v, &new_extr2, &mut degen2);
        if degen1 && degen2 {
            self.status = ChFi2dConstructionError::BothEdgesDegenerated;
        }
        if degen1 && !degen2 {
            self.status = ChFi2dConstructionError::FirstEdgeDegenerated;
        }
        if !degen1 && degen2 {
            self.status = ChFi2dConstructionError::LastEdgeDegenerated;
        }
    } // ComputeFillet

    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L534-600 — BuildNewWire
    // =======================================================================

    /// OCCT ChFi2d_Builder::BuildNewWire (L534-600): replaces in the new
    /// face <new_face> <old_e1> and <old_e2> by <e1>, <fillet> and <e2>
    /// (degenerate cases drop the corresponding trimmed edge).
    pub fn build_new_wire(
        &mut self,
        old_e1: &Shape,
        old_e2: &Shape,
        e1: &Shape,
        fillet: &Shape,
        e2: &Shape,
    ) {
        let mut a_closed_status = true;

        // OCCT L543-549: TopExp_Explorer Ex(refFace, TopAbs_WIRE);
        for a_wire in topexp_explore_face_wires(&self.my_brep, &self.ref_face) {
            a_closed_status = self.my_brep.has_flag(a_wire.clone(), tshape_flags::CLOSED);
            break;
        }

        let mut fillet_is_added = false;

        // OCCT L553: Ex.Init(newFace, TopAbs_EDGE);
        let ex = topexp_explore_face_edges(&self.my_brep, &self.new_face);
        // OCCT L555-556: BRep_Builder B; B.MakeWire(newWire);
        let mut b = BRepBuilder::new();
        let new_wire = b.make_wire(&mut self.my_brep);

        for the_edge in ex {
            // OCCT L560: const TopoDS_Edge& theEdge = TopoDS::Edge(Ex.Current());
            if !the_edge.is_same(old_e1) && !the_edge.is_same(old_e2) {
                b.add_to_wire(&mut self.my_brep, new_wire.clone(), the_edge);
            } else {
                // OCCT L567: if (theEdge == OldE1) — operator== is IsEqual.
                if the_edge.is_equal(old_e1) {
                    if self.status != ChFi2dConstructionError::FirstEdgeDegenerated
                        && self.status != ChFi2dConstructionError::BothEdgesDegenerated
                    {
                        b.add_to_wire(&mut self.my_brep, new_wire.clone(), e1.clone());
                    }
                    if !fillet_is_added {
                        b.add_to_wire(&mut self.my_brep, new_wire.clone(), fillet.clone());
                        fillet_is_added = true;
                    } // if ( !filletIsAdded ...
                } // if (theEdge == ...
                else {
                    if self.status != ChFi2dConstructionError::LastEdgeDegenerated
                        && self.status != ChFi2dConstructionError::BothEdgesDegenerated
                    {
                        b.add_to_wire(&mut self.my_brep, new_wire.clone(), e2.clone());
                    }
                    if !fillet_is_added {
                        b.add_to_wire(&mut self.my_brep, new_wire.clone(), fillet.clone());
                        fillet_is_added = true;
                    } // if ( !filletIsAdded ...
                } // else ...
            } // else ...
        } // while ...

        // OCCT L595: newWire.Closed(aClosedStatus);
        {
            let wd = self.my_brep.wire_mut(new_wire.clone());
            if a_closed_status {
                wd.flags |= tshape_flags::CLOSED;
            } else {
                wd.flags &= !tshape_flags::CLOSED;
            }
        }
        // OCCT L596-598:
        // BRepAdaptor_Surface Adaptor3dSurface(refFace);
        // BRepLib_MakeFace mFace(Adaptor3dSurface.Plane(), newWire);
        // newFace = mFace;
        let ref_plane = match self.my_brep.face_surface_world(&self.ref_face) {
            Some(Surface3::Plane(pl)) => pl,
            _ => panic!("Standard_TypeMismatch: refFace is not a plane"),
        };
        let m_face =
            b.make_face(&mut self.my_brep, Some(Surface3::Plane(ref_plane)), new_wire.clone());
        self.new_face = m_face;
    } // BuildNewWire

    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L604-629 — BuildNewEdge (3-arg)
    // =======================================================================

    /// OCCT ChFi2d_Builder::BuildNewEdge(const TopoDS_Edge& E1,
    /// const TopoDS_Vertex& OldExtr, const TopoDS_Vertex& NewExtr) const
    /// (L604-629): changes <old_extr> of <e1> by <new_extr>.
    /// rcad architecture: OCCT declares the method const; the rcad BRep
    /// pool mutation (edge registration inside BRepLib_MakeEdge) requires
    /// &mut self.
    pub fn build_new_edge(&mut self, e1: &Shape, old_extr: &Shape, new_extr: &Shape) -> Shape {
        // OCCT L616: occ::handle<Geom_Curve> curve = BRep_Tool::Curve(E1, first, last);
        let (curve, _range) = match self.my_brep.edge_curve_world(e1) {
            Some(x) => x,
            None => return Shape::null(),
        };
        let first_vertex = self.my_brep.first_vertex(e1);
        let last_vertex = self.my_brep.last_vertex(e1);
        // OCCT L617-624: BRepLib_MakeEdge makeEdge; makeEdge.Init(curve, ...);
        let made = if first_vertex.is_same(old_extr) {
            brep_lib_make_edge_init_vertices(&mut self.my_brep, &curve, new_extr, &last_vertex)
        } else {
            brep_lib_make_edge_init_vertices(&mut self.my_brep, &curve, &first_vertex, new_extr)
        };
        // OCCT L625: TopoDS_Edge anEdge = makeEdge;  (implicit Edge() conversion)
        let mut an_edge = made;
        an_edge.orientation = e1.orientation;
        //  anEdge.Location(E1.Location());
        an_edge
    }

    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L636-680 — BuildNewEdge (4-arg)
    // =======================================================================

    /// OCCT ChFi2d_Builder::BuildNewEdge(E1, OldExtr, NewExtr,
    /// bool& IsDegenerated) const (L636-680) — special flag if the new
    /// edge is degenerated. Rust overload name:
    /// `build_new_edge_degenerated` (OCCT second BuildNewEdge overload).
    /// rcad architecture: &mut self per the note on `build_new_edge`.
    pub fn build_new_edge_degenerated(
        &mut self,
        e1: &Shape,
        old_extr: &Shape,
        new_extr: &Shape,
        is_degenerated: &mut bool,
    ) -> Shape {
        *is_degenerated = false; // OCCT L644
        let first_vertex = self.my_brep.first_vertex(e1);
        let last_vertex = self.my_brep.last_vertex(e1);
        let pnew = self.my_brep.vertex_position(new_extr); // OCCT L647: gp_Pnt Pnew = BRep_Tool::Pnt(NewExtr);
        let mut ponctual_edge = false; // OCCT L648: bool PonctualEdge = false;
        let tol = CONFUSION; // OCCT L649: double Tol = Precision::Confusion();
        // OCCT L653: occ::handle<Geom_Curve> curve = BRep_Tool::Curve(E1, first, last);
        let (curve, _range) = match self.my_brep.edge_curve_world(e1) {
            Some(x) => x,
            None => return Shape::null(),
        };
        let mut an_edge;
        if first_vertex.is_same(old_extr) {
            // OCCT L656: makeEdge.Init(curve, NewExtr, lastVertex);
            let made =
                brep_lib_make_edge_init_vertices(&mut self.my_brep, &curve, new_extr, &last_vertex);
            // OCCT L657-658: gp_Pnt PV = BRep_Tool::Pnt(lastVertex);
            //                 PonctualEdge = (Pnew.Distance(PV) < Tol);
            let pv = self.my_brep.vertex_position(&last_vertex);
            ponctual_edge = pnew.distance(pv) < tol;
            // OCCT L667: BRepLib_EdgeError error = makeEdge.Error();
            // rcad architecture: BRepLib_MakeEdge carries no error state in
            // the rcad make helper — a failed make returns a null edge; the
            // BRepLib_LineThroughIdenticPoints condition (two construction
            // vertices at the same point within tolerance) coincides with
            // the PonctualEdge test.
            let error_line_through_identic_points = made.is_null();
            // OCCT L668-672:
            // if (error == BRepLib_LineThroughIdenticPoints || PonctualEdge)
            // { IsDegenerated = true; anEdge = E1; }
            if error_line_through_identic_points || ponctual_edge {
                *is_degenerated = true;
                return e1.clone();
            }
            an_edge = made;
            an_edge.orientation = e1.orientation;
            //    anEdge.Location(E1.Location());
        } else {
            // OCCT L662: makeEdge.Init(curve, firstVertex, NewExtr);
            let made = brep_lib_make_edge_init_vertices(
                &mut self.my_brep,
                &curve,
                &first_vertex,
                new_extr,
            );
            // OCCT L663-664: gp_Pnt PV = BRep_Tool::Pnt(firstVertex);
            //                 PonctualEdge = (Pnew.Distance(PV) < Tol);
            let pv = self.my_brep.vertex_position(&first_vertex);
            ponctual_edge = pnew.distance(pv) < tol;
            let error_line_through_identic_points = made.is_null();
            if error_line_through_identic_points || ponctual_edge {
                *is_degenerated = true;
                return e1.clone();
            }
            an_edge = made;
            an_edge.orientation = e1.orientation;
            //    anEdge.Location(E1.Location());
        }
        let _ = ponctual_edge;
        an_edge
    }

    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L684-716 — UpDateHistory (6-arg)
    // =======================================================================

    /// OCCT ChFi2d_Builder::UpDateHistory(E1, E2, TrimE1, TrimE2, NewEdge,
    /// Id) (L684-716): writes <new_edge> in <fillets> if <id> is equal to
    /// 1, or in <chamfers> if <id> is equal to 2; writes the modifications
    /// in <history>: <trim_e1> is given by <e1>, <trim_e2> by <e2> if
    /// <trim_e1> and <trim_e2> are not degenerated. Rust overload name:
    /// `up_date_history_new_edge` (OCCT 6-arg UpDateHistory overload).
    #[allow(clippy::too_many_arguments)]
    pub fn up_date_history_new_edge(
        &mut self,
        e1: &Shape,
        e2: &Shape,
        trim_e1: &Shape,
        trim_e2: &Shape,
        new_edge: &Shape,
        id: i32,
    ) {
        if id == 1 {
            // the new edge is a fillet
            self.fillets.push(new_edge.clone());
        } else {
            // the new edge is a chamfer
            self.chamfers.push(new_edge.clone());
        }

        self.history.remove(&e1.ptr_id()); // OCCT L700: history.UnBind(E1);
        if self.status != ChFi2dConstructionError::FirstEdgeDegenerated
            && self.status != ChFi2dConstructionError::BothEdgesDegenerated
        {
            if !e1.is_same(trim_e1) {
                // history.Bind(E1, TrimE1);
                self.history.insert(e1.ptr_id(), (e1.clone(), trim_e1.clone()));
            }
        }
        self.history.remove(&e2.ptr_id()); // OCCT L708: history.UnBind(E2);
        if self.status != ChFi2dConstructionError::LastEdgeDegenerated
            && self.status != ChFi2dConstructionError::BothEdgesDegenerated
        {
            if !e2.is_same(trim_e2) {
                // history.Bind(E2, TrimE2);
                self.history.insert(e2.ptr_id(), (e2.clone(), trim_e2.clone()));
            }
        }
    } // UpDateHistory

    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L720-742 — UpDateHistory (4-arg)
    // =======================================================================

    /// OCCT ChFi2d_Builder::UpDateHistory(E1, E2, TrimE1, TrimE2)
    /// (L720-742): writes the modifications in <history>. <trim_e1> is
    /// given by <e1>, <trim_e2> by <e2>.
    pub fn up_date_history(&mut self, e1: &Shape, e2: &Shape, trim_e1: &Shape, trim_e2: &Shape) {
        if self.history.contains_key(&e1.ptr_id()) {
            self.history.remove(&e1.ptr_id()); // history.UnBind(E1);
        }
        if !e1.is_same(trim_e1) {
            // history.Bind(E1, TrimE1);
            self.history.insert(e1.ptr_id(), (e1.clone(), trim_e1.clone()));
        }
        if self.history.contains_key(&e2.ptr_id()) {
            self.history.remove(&e2.ptr_id()); // history.UnBind(E2);
        }
        if !e2.is_same(trim_e2) {
            // history.Bind(E2, TrimE2);
            self.history.insert(e2.ptr_id(), (e2.clone(), trim_e2.clone()));
        }
    } // UpDateHistory

    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L746-762 — BasisEdge
    // =======================================================================

    /// OCCT ChFi2d_Builder::BasisEdge(const TopoDS_Edge& E) const
    /// (L746-762): returns the modified edge if <e> has descendant or
    /// <e> in the other case. Architecture note: the OCCT reference
    /// return becomes a clone (Rust cannot return a reference into the
    /// map).
    pub fn basis_edge(&self, e: &Shape) -> Shape {
        for (_key, (key_shape, value_shape)) in &self.history {
            let an_edge = value_shape; // TopoDS::Edge(iterator.Value());
            if an_edge.is_same(e) {
                let another_edge = key_shape; // TopoDS::Edge(iterator.Key());
                return another_edge.clone();
            } // if (anEdge.IsSame ...
        } // while (Iterator.More ...
        e.clone()
    } // BasisEdge

    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L766-1230 — BuildFilletEdge
    // =======================================================================

    /// OCCT ChFi2d_Builder::BuildFilletEdge (L766-1230). Is internally
    /// used by <compute_fillet>: <new_extr1> and <new_extr2> will contain
    /// the new extremities of <adj_edge1> and <adj_edge2>.
    #[allow(clippy::too_many_arguments)]
    pub fn build_fillet_edge(
        &mut self,
        v: &Shape,
        adj_edge1: &Shape,
        adj_edge2: &Shape,
        radius: f64,
        new_extr1: &mut Shape,
        new_extr2: &mut Shape,
    ) -> Shape {
        let mut e1 = adj_edge1.clone(); // OCCT L774: E1 = AdjEdge1;
        let mut e2 = adj_edge2.clone(); // OCCT L775: E2 = AdjEdge2;
        // OCCT L776-779: TopExp::FirstVertex/LastVertex on E1/E2.
        let v1 = self.my_brep.first_vertex(&e1);
        let v2 = self.my_brep.last_vertex(&e1);
        let v3 = self.my_brep.first_vertex(&e2);
        let v4 = self.my_brep.last_vertex(&e2);

        // ====================================================================
        //    The first arc is found.                                        +
        // ====================================================================

        let o1: Orientation; // OCCT L785: TopAbs_Orientation O1;
        let oe1 = e1.orientation; // OCCT L787: OE1 = E1.Orientation();
        e1.orientation = Orientation::Forward; // OCCT L788: E1.Orientation(TopAbs_FORWARD);
        e2.orientation = Orientation::Forward; // OCCT L789: E2.Orientation(TopAbs_FORWARD);
        // OCCT L790-793: Ebid1 = TopoDS::Edge(E1.EmptyCopied());
        //                Ebid2 = TopoDS::Edge(E2.EmptyCopied());
        // (created and left unused in OCCT as well)
        let _ebid1 = self.my_brep.empty_copied(&e1);
        let _ebid2 = self.my_brep.empty_copied(&e2);

        // ====================================================================
        //    Save non-modified parts of edges concerned.      +
        // ====================================================================

        let mut param1;
        let mut param2;
        let mut param3;
        let mut param4;

        if v1.is_same(v) {
            param1 = brep_tool_parameter(&self.my_brep,&v1, &e1); // OCCT L804: BRep_Tool::Parameter(V1, E1)
            param2 = brep_tool_parameter(&self.my_brep,&v2, &e1); // OCCT L805: BRep_Tool::Parameter(V2, E1)
            o1 = v2.orientation; // OCCT L806: O1 = V2.Orientation();
        } else {
            param1 = brep_tool_parameter(&self.my_brep,&v2, &e1); // OCCT L810
            param2 = brep_tool_parameter(&self.my_brep,&v1, &e1); // OCCT L811
            o1 = v1.orientation; // OCCT L812
        }
        if v3.is_same(v) {
            param3 = brep_tool_parameter(&self.my_brep,&v3, &e2); // OCCT L816
            param4 = brep_tool_parameter(&self.my_brep,&v4, &e2); // OCCT L817
        } else {
            param3 = brep_tool_parameter(&self.my_brep,&v4, &e2); // OCCT L821
            param4 = brep_tool_parameter(&self.my_brep,&v3, &e2); // OCCT L822
        }

        // ====================================================================
        //    Restore geometric supports.                            +
        // ====================================================================

        // OCCT L832-833:
        // C1 = BRep_Tool::CurveOnSurface(E1, newFace, ufirst1, ulast1);
        // C2 = BRep_Tool::CurveOnSurface(E2, newFace, ufirst2, ulast2);
        let (c1, ufirst1, ulast1) = match self.my_brep.curve_on_surface(&e1, &self.new_face) {
            Some((c, f, l)) => (c, f, l),
            None => return Shape::null(),
        };
        let (c2, ufirst2, ulast2) = match self.my_brep.curve_on_surface(&e2, &self.new_face) {
            Some((c, f, l)) => (c, f, l),
            None => return Shape::null(),
        };
        let mut u1;
        let mut u2;
        let mut pu1;
        let mut pu2;
        let mut vv1;
        let mut vv2;
        let mut ppu1;
        let mut ppu2;

        // ====================================================================
        //   Determination of the face for fillet.                                +
        // ====================================================================

        let mut p = DVec2::ZERO; // OCCT L839: gp_Pnt2d p;
        let ve1;
        let ve2;
        let mut ve3;
        let mut ve4;
        let sens1: bool;
        let sens2: bool;

        // OCCT L844-862: basisC1/basisC2 = the un-trimmed basis curves.
        let basis_c1 = match &c1 {
            // OCCT L845-849: T1 = down_cast<Geom2d_TrimmedCurve>(C1);
            // basisC1 = T1->BasisCurve();
            Curve2d::Trimmed(t1) => (*t1.curve).clone(),
            _ => c1.clone(),
        };
        let basis_c2 = match &c2 {
            Curve2d::Trimmed(t2) => (*t2.curve).clone(),
            _ => c2.clone(),
        };

        if matches!(basis_c1, Curve2d::Circle(_)) {
            // OCCT L866-869:
            // occ::handle<Geom2d_Circle> CC1 = down_cast<Geom2d_Circle>(basisC1);
            // ElCLib::D1(param1, CC1->Circ2d(), p, Ve1);
            // Sens1 = (CC1->Circ2d()).IsDirect();
            if let Curve2d::Circle(cc1) = &basis_c1 {
                let (p_val, ve1_val) =
                    elclib2d::circle_d1(cc1.center, cc1.x_dir, cc1.y_dir, cc1.radius, param1);
                p = p_val;
                ve1 = ve1_val;
                // OCCT gp_Circ2d::IsDirect — the frame is direct when
                // x_dir x y_dir is positive (gp_Ax22d sense flag).
                sens1 = cc1.x_dir.x * cc1.y_dir.y - cc1.x_dir.y * cc1.y_dir.x > 0.0;
            } else {
                unreachable!()
            }
        } // if (C1->DynamicType() ...
        else {
            // OCCT L872-874:
            // occ::handle<Geom2d_Line> CC1 = down_cast<Geom2d_Line>(basisC1);
            // ElCLib::D1(param1, CC1->Lin2d(), p, Ve1);
            // Sens1 = true;
            if let Curve2d::Line(cc1) = &basis_c1 {
                let (p_val, ve1_val) = elclib2d::line_d1(cc1.origin, cc1.direction, param1);
                p = p_val;
                ve1 = ve1_val;
                sens1 = true;
            } else {
                unreachable!()
            }
        } // else ...
        if matches!(basis_c2, Curve2d::Circle(_)) {
            // OCCT L878-880: same for C2 with param3.
            if let Curve2d::Circle(cc2) = &basis_c2 {
                let (p_val, ve2_val) =
                    elclib2d::circle_d1(cc2.center, cc2.x_dir, cc2.y_dir, cc2.radius, param3);
                p = p_val;
                ve2 = ve2_val;
                sens2 = cc2.x_dir.x * cc2.y_dir.y - cc2.x_dir.y * cc2.y_dir.x > 0.0;
            } else {
                unreachable!()
            }
        } // if (C2->DynamicType() ...
        else {
            // OCCT L884-886: same for C2 line branch.
            if let Curve2d::Line(cc2) = &basis_c2 {
                let (p_val, ve2_val) = elclib2d::line_d1(cc2.origin, cc2.direction, param3);
                p = p_val;
                ve2 = ve2_val;
                sens2 = true;
            } else {
                unreachable!()
            }
        } // else ...

        let mut fillet_edge = Shape::null(); // OCCT L889: TopoDS_Edge filletEdge;

        // OCCT L891-893: cross = Ve1.Crossed(Ve2); Ve3 = Ve1; Ve4 = Ve2;
        let mut cross = ve1.x * ve2.y - ve1.y * ve2.x;
        ve3 = ve1;
        ve4 = ve2;

        // processing of tangency or downcast point
        // OCCT L896: if (Ve1.IsParallel(Ve2, Precision::Angular()))
        if gp_vec2d_is_parallel(ve1, ve2, ANGULAR) {
            // Ve1 and Ve2 are parallel : cross at 0
            cross = 0.0; // OCCT L899: cross = 0.;
            if param1 < param2 {
                ve3 = -ve1; // OCCT L902: Ve3 = -Ve1;
            }
            if param3 > param4 {
                ve4 = -ve2; // OCCT L905: Ve4 = -Ve2;
            }

            // OCCT L909: if (!Ve4.IsOpposite(Ve3, Precision::Angular()))
            if !gp_vec2d_is_opposite(ve4, ve3, ANGULAR) {
                // There is a true tangency point and the calculation is stopped
                self.status = ChFi2dConstructionError::TangencyError;
                return fillet_edge;
            }
            // Otherwise this is a downcast point, and the calculation is continued
        }

        // OCCT L918: GccEnt_Position Qual1, Qual2;
        let qual1: GccEntPosition;
        let qual2: GccEntPosition;
        if cross < 0.0 {
            if param3 > param4 {
                if sens1 {
                    qual1 = GccEntPosition::Enclosed; // GccEnt_enclosed
                } else {
                    qual1 = GccEntPosition::Outside; // GccEnt_outside
                }
            } else {
                if sens1 {
                    qual1 = GccEntPosition::Outside;
                } else {
                    qual1 = GccEntPosition::Enclosed;
                }
            }
            if param1 > param2 {
                if sens2 {
                    qual2 = GccEntPosition::Outside;
                } else {
                    qual2 = GccEntPosition::Enclosed;
                }
            } else {
                if sens2 {
                    qual2 = GccEntPosition::Enclosed;
                } else {
                    qual2 = GccEntPosition::Outside;
                }
            }
        } // if (cross < 0 ...
        else {
            if param3 > param4 {
                if sens1 {
                    qual1 = GccEntPosition::Outside;
                } else {
                    qual1 = GccEntPosition::Enclosed;
                }
            } else {
                if sens1 {
                    qual1 = GccEntPosition::Enclosed;
                } else {
                    qual1 = GccEntPosition::Outside;
                }
            }
            if param1 > param2 {
                if sens2 {
                    qual2 = GccEntPosition::Enclosed;
                } else {
                    qual2 = GccEntPosition::Outside;
                }
            } else {
                if sens2 {
                    qual2 = GccEntPosition::Outside;
                } else {
                    qual2 = GccEntPosition::Enclosed;
                }
            }
        } // else ...

        let tol = CONFUSION; // OCCT L1014: double Tol = Precision::Confusion();
        // OCCT L1015-1018:
        // Geom2dGcc_Circ2d2TanRad Fillet(Geom2dGcc_QualifiedCurve(basisC1, Qual1),
        //                                Geom2dGcc_QualifiedCurve(basisC2, Qual2),
        //                                Radius, Tol);
        // (the local OCCT variable name `Fillet` collides with the output
        //  parameter in Rust; renamed a_fillet)
        let a_fillet = Circ2d2TanRad::new_curve_curve(
            &QualifiedCurve::new(basis_c1.clone(), qual1),
            &QualifiedCurve::new(basis_c2.clone(), qual2),
            radius,
            tol,
        );
        // OCCT L1019-1023: if (!Fillet.IsDone() || Fillet.NbSolutions() == 0)
        if !a_fillet.is_done() || a_fillet.nb_solutions() == 0 {
            self.status = ChFi2dConstructionError::ComputationError;
            return fillet_edge;
        }
        // OCCT L1024: else if (Fillet.NbSolutions() >= 1)
        self.status = ChFi2dConstructionError::IsDone;
        let mut numsol = 1usize; // OCCT L1027: int numsol = 1;
        let mut nsol = 1usize; // OCCT L1028: int nsol = 1;
        let mut ptg1 = DVec2::ZERO; // OCCT L1030: gp_Pnt2d Ptg1, Ptg2;
        let mut ptg2 = DVec2::ZERO;
        let mut dist;
        let mut dist1 = 1.0e40f64; // OCCT L1032: double dist1 = 1.e40;
        let mut inside = false; // OCCT L1033: bool inside = false;
        let two_pi = std::f64::consts::PI * 2.0;
        while nsol <= a_fillet.nb_solutions() {
            // OCCT L1036: Fillet.Tangency1(nsol, PU1, PU2, Ptg1);
            let (pu1_v, pu2_v, ptg1_v) = a_fillet.tangency1(nsol);
            pu1 = pu1_v;
            pu2 = pu2_v;
            ptg1 = ptg1_v;
            let _ = pu1;
            dist = ptg1.distance(p);
            if matches!(basis_c1, Curve2d::Line(_)) {
                // OCCT L1040: inside = (PU2 < param1 && PU2 > param2) || (PU2 < param2 && PU2 > param1);
                inside = (pu2 < param1 && pu2 > param2) || (pu2 < param2 && pu2 > param1);
                if inside && dist < dist1 {
                    numsol = nsol;
                    dist1 = dist;
                } // if ((((inside && ...
            } // if (C1->DynamicType( ...
            else {
                // OCCT L1049: Fillet.Tangency2(nsol, PPU1, PPU2, Ptg2);
                let (ppu1_v, ppu2_v, ptg2_v) = a_fillet.tangency2(nsol);
                ppu1 = ppu1_v;
                ppu2 = ppu2_v;
                ptg2 = ptg2_v;
                let _ = ppu1;
                dist = ptg2.distance(p);
                // OCCT L1051: inside = (PPU2 < param3 && PPU2 > param4) || (PPU2 < param4 && PPU2 > param3);
                inside = (ppu2 < param3 && ppu2 > param4) || (ppu2 < param4 && ppu2 > param3);
                // OCCT L1053-1055: case of arc of circle passing on the sewing
                if matches!(basis_c2, Curve2d::Circle(_))
                    && ((two_pi < param3 && two_pi > param4)
                        || (two_pi < param4 && two_pi > param3))
                {
                    // cas param3<param4
                    // OCCT L1058: inside = (param3 < PPU2 && PPU2 < 2*M_PI) || (0 <= PPU2 && PPU2 < param4 - 2*M_PI);
                    inside =
                        (param3 < ppu2 && ppu2 < two_pi) || (0.0 <= ppu2 && ppu2 < param4 - two_pi);
                    // cas param4<param3
                    // OCCT L1060-1061: inside = inside || (param4 < PPU2 && PPU2 < 2*M_PI) || (0 <= PPU2 && PPU2 < param3 - 2*M_PI);
                    inside = inside
                        || (param4 < ppu2 && ppu2 < two_pi)
                        || (0.0 <= ppu2 && ppu2 < param3 - two_pi);
                }
                if inside && dist < dist1 {
                    numsol = nsol;
                    dist1 = dist;
                } // if ((((param3 ...
            } // else ...
            nsol += 1;
        } // while (nsol ...
        // OCCT L1071-1072:
        // gp_Circ2d cir(Fillet.ThisSolution(numsol));
        // occ::handle<Geom2d_Circle> circle = new Geom2d_Circle(cir);
        let cir: Circle2d = a_fillet.this_solution(numsol);
        let circle = Curve2d::Circle(cir);

        // OCCT L1074-1076:
        // BRep_Builder B;
        // BRepAdaptor_Surface Adaptor3dSurface(refFace);
        // occ::handle<Geom_Plane> refSurf = new Geom_Plane(Adaptor3dSurface.Plane());
        let ref_plane = match self.my_brep.face_surface_world(&self.ref_face) {
            Some(Surface3::Plane(pl)) => pl,
            _ => panic!("Standard_TypeMismatch: refFace is not a plane"),
        };
        let _ = &ref_plane;
        // OCCT L1077-1078: Fillet.Tangency1(numsol, U1, U2, Ptg1);
        //                  Fillet.Tangency2(numsol, Vv1, Vv2, Ptg2);
        let (u1_v, u2_v, ptg1_v) = a_fillet.tangency1(numsol);
        u1 = u1_v;
        u2 = u2_v;
        ptg1 = ptg1_v;
        let (vv1_v, vv2_v, ptg2_v) = a_fillet.tangency2(numsol);
        vv1 = vv1_v;
        vv2 = vv2_v;
        ptg2 = ptg2_v;

        // check the validity of parameters
        //// modified by jgv, 08.08.2011 for bug 0022695 ////
        // OCCT L1083: inside = (U2 < param1 && U2 >= param2) || (U2 <= param2 && U2 > param1);
        inside = (u2 < param1 && u2 >= param2) || (u2 <= param2 && u2 > param1);
        /////////////////////////////////////////////////////
        // OCCT L1085-1086:
        if matches!(basis_c1, Curve2d::Circle(_))
            && ((two_pi < param1 && two_pi > param2) || (two_pi < param2 && two_pi > param1))
        {
            // arc of circle containing the circle origin
            // case param1<param2
            // OCCT L1090: inside = (param1 < U2 && U2 < 2*M_PI) || (0 <= U2 && U2 < param2 - 2*M_PI);
            inside = (param1 < u2 && u2 < two_pi) || (0.0 <= u2 && u2 < param2 - two_pi);
            // case param2<param1
            // OCCT L1092: inside = inside || (param2 < U2 && U2 < 2*M_PI) || (0 <= U2 && U2 < param1 - 2*M_PI);
            inside = inside || (param2 < u2 && u2 < two_pi) || (0.0 <= u2 && u2 < param1 - two_pi);
        }
        if !inside {
            self.status = ChFi2dConstructionError::ComputationError;
            return fillet_edge;
        }

        //// modified by jgv, 08.08.2011 for bug 0022695 ////
        // OCCT L1102: inside = (Vv2 < param3 && Vv2 >= param4) || (Vv2 <= param4 && Vv2 > param3);
        inside = (vv2 < param3 && vv2 >= param4) || (vv2 <= param4 && vv2 > param3);
        /////////////////////////////////////////////////////
        // OCCT L1104-1105:
        if matches!(basis_c2, Curve2d::Circle(_))
            && ((two_pi < param3 && two_pi > param4) || (two_pi < param4 && two_pi > param3))
        {
            // arc of circle containing the circle origin
            // cas param3<param4
            // OCCT L1109: inside = (param3 < Vv2 && Vv2 < 2*M_PI) || (0 <= Vv2 && Vv2 < param4 - 2*M_PI);
            inside = (param3 < vv2 && vv2 < two_pi) || (0.0 <= vv2 && vv2 < param4 - two_pi);
            // cas param4<param3
            // OCCT L1110-1111: inside = inside || (param4 < Vv2 && Vv2 < 2*M_PI) || (0 <= Vv2 && Vv2 < param3 - 2*M_PI);
            inside =
                inside || (param4 < vv2 && vv2 < two_pi) || (0.0 <= vv2 && vv2 < param3 - two_pi);
        }
        if !inside {
            self.status = ChFi2dConstructionError::ComputationError;
            return fillet_edge;
        }

        // OCCT L1119-1120:
        // gp_Pnt p1 = Adaptor3dSurface.Value(Ptg1.X(), Ptg1.Y());
        // gp_Pnt p2 = Adaptor3dSurface.Value(Ptg2.X(), Ptg2.Y());
        let p1 = ref_plane.point_at(ptg1.x, ptg1.y);
        let p2 = ref_plane.point_at(ptg2.x, ptg2.y);
        // OCCT L1121-1122: B.MakeVertex(Vertex1, p1, Tol); NewExtr1 = Vertex1;
        let vertex1 = brep_builder_make_vertex(&mut self.my_brep, p1, tol);
        *new_extr1 = vertex1.clone();
        // OCCT L1123-1126:
        if (u2 - ufirst1).abs() <= PCONFUSION {
            *new_extr1 = v1.clone();
        }
        if (u2 - ulast1).abs() <= PCONFUSION {
            *new_extr1 = v2.clone();
        }

        // OCCT L1132-1133: B.MakeVertex(Vertex2, p2, Tol); NewExtr2 = Vertex2;
        let vertex2 = brep_builder_make_vertex(&mut self.my_brep, p2, tol);
        *new_extr2 = vertex2.clone();
        // OCCT L1134-1137:
        if (vv2 - ufirst2).abs() <= PCONFUSION {
            *new_extr2 = v3.clone();
        }
        if (vv2 - ulast2).abs() <= PCONFUSION {
            *new_extr2 = v4.clone();
        }

        // ====================================================================
        //   Update tops of the fillet.                                  +
        // ====================================================================
        // OCCT L1146-1149:
        // gp_Pnt Pntbid; gp_Pnt2d sommet;
        // Pntbid = BRep_Tool::Pnt(V); sommet = gp_Pnt2d(Pntbid.X(), Pntbid.Y());
        let pntbid = self.my_brep.vertex_position(v);
        let sommet = DVec2::new(pntbid.x, pntbid.y);

        // OCCT L1151-1162: gp_Pnt pntBid; gp_Pnt2d somBid;
        let pnt_bid;
        let som_bid;
        if v1.is_same(v) {
            pnt_bid = self.my_brep.vertex_position(&v2);
            som_bid = DVec2::new(pnt_bid.x, pnt_bid.y);
        } else {
            pnt_bid = self.my_brep.vertex_position(&v1);
            som_bid = DVec2::new(pnt_bid.x, pnt_bid.y);
        }

        // OCCT L1164-1165: gp_Vec2d vec; ElCLib::D1(U1, cir, Ptg1, vec);
        let (_ptg1_d1, vec) = elclib2d::circle_d1(cir.center, cir.x_dir, cir.y_dir, cir.radius, u1);
        let _ = _ptg1_d1;

        let mut vec1;
        // OCCT L1168-1175:
        if matches!(basis_c1, Curve2d::Circle(_)) {
            if let Curve2d::Circle(cc1) = &basis_c1 {
                // OCCT L1171: gp_Circ2d cir2d(CC1->Circ2d());
                let cir2d = *cc1;
                // OCCT L1172: double par = ElCLib::Parameter(cir2d, Ptg1);
                let par = elclib2d::circle_parameter(cir2d.center, cir2d.x_dir, cir2d.y_dir, ptg1);
                // OCCT L1174: ElCLib::D1(par, cir2d, Pd, vec1);
                let (_pd, vec1_v) =
                    elclib2d::circle_d1(cir2d.center, cir2d.x_dir, cir2d.y_dir, cir2d.radius, par);
                vec1 = vec1_v;
            } else {
                unreachable!()
            }
        } // if (C1->DynamicType() ...
        else if matches!(basis_c1, Curve2d::Line(_)) {
            if let Curve2d::Line(cc1) = &basis_c1 {
                // OCCT L1179: gp_Lin2d lin2d(CC1->Lin2d());
                let lin2d = *cc1;
                // OCCT L1180: double par = ElCLib::Parameter(lin2d, sommet);
                let par = elclib2d::line_parameter(lin2d.origin, lin2d.direction, sommet);
                // OCCT L1181: vec1 = gp_Vec2d(sommet.X() - somBid.X(), sommet.Y() - somBid.Y());
                vec1 = DVec2::new(sommet.x - som_bid.x, sommet.y - som_bid.y);
                // OCCT L1183: ElCLib::D1(par, lin2d, Pd, vec1);
                let (_pd, vec1_v) = elclib2d::line_d1(lin2d.origin, lin2d.direction, par);
                vec1 = vec1_v;
            } else {
                unreachable!()
            }
        } // else if ...
        else {
            vec1 = DVec2::ZERO;
        }

        // OCCT L1186-1189: if (OE1 == TopAbs_REVERSED) vec1.Reverse();
        if oe1 == Orientation::Reversed {
            vec1 = -vec1;
        } // if (OE1 ...
        // OCCT L1190: bool Sense = (vec1 * vec) > 0.;
        // (gp_Vec2d operator*(gp_Vec2d) is the dot product)
        let sense = vec1.dot(vec) > 0.0;
        // OCCT L1191-1194: if (U1 > Vv1 && U1 > 2. * M_PI)
        //   ElCLib::AdjustPeriodic(0., 2. * M_PI, Precision::Confusion(), U1, Vv1);
        if u1 > vv1 && u1 > 2.0 * std::f64::consts::PI {
            elclib_adjust_periodic(0.0, 2.0 * std::f64::consts::PI, CONFUSION, &mut u1, &mut vv1);
        } // if (U1 ...
        // OCCT L1195-1203:
        if (o1 == Orientation::Forward && oe1 == Orientation::Forward)
            || (o1 == Orientation::Reversed && oe1 == Orientation::Reversed)
        {
            // OCCT L1198: filletEdge = BRepLib_MakeEdge(circle, refSurf, NewExtr1, NewExtr2, U1, Vv1);
            fillet_edge = brep_lib_make_edge_pcurve(
                &mut self.my_brep,
                &circle,
                &self.ref_face,
                new_extr1,
                new_extr2,
                u1,
                vv1,
            );
        } // if (O1 == ...
        else {
            // OCCT L1202: filletEdge = BRepLib_MakeEdge(circle, refSurf, NewExtr2, NewExtr1, Vv1, U1);
            fillet_edge = brep_lib_make_edge_pcurve(
                &mut self.my_brep,
                &circle,
                &self.ref_face,
                new_extr2,
                new_extr1,
                vv1,
                u1,
            );
        } // else ...
        // OCCT L1204-1224: if (!Sense)
        if !sense {
            let s1 = fillet_edge.orientation; // OCCT L1206: TopAbs_Orientation S1 = filletEdge.Orientation();
            if (o1 == Orientation::Forward && oe1 == Orientation::Forward)
                || (o1 == Orientation::Reversed && oe1 == Orientation::Reversed)
            {
                // OCCT L1210: filletEdge = BRepLib_MakeEdge(circle, refSurf, NewExtr2, NewExtr1, Vv1, U1);
                fillet_edge = brep_lib_make_edge_pcurve(
                    &mut self.my_brep,
                    &circle,
                    &self.ref_face,
                    new_extr2,
                    new_extr1,
                    vv1,
                    u1,
                );
            } else {
                // OCCT L1214: filletEdge = BRepLib_MakeEdge(circle, refSurf, NewExtr1, NewExtr2, U1, Vv1);
                fillet_edge = brep_lib_make_edge_pcurve(
                    &mut self.my_brep,
                    &circle,
                    &self.ref_face,
                    new_extr1,
                    new_extr2,
                    u1,
                    vv1,
                );
            }
            if s1 == Orientation::Forward {
                fillet_edge.orientation = Orientation::Reversed;
            } else {
                fillet_edge.orientation = Orientation::Forward;
            }
        } // if (!Sense

        brep_lib_build_curves3d(&mut self.my_brep, &fillet_edge); // OCCT L1228
        fillet_edge
    } // BuildFilletEdge

    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L1234-1247 — IsAFillet
    // =======================================================================

    /// OCCT ChFi2d_Builder::IsAFillet(const TopoDS_Edge& E) const
    /// (L1234-1247).
    pub fn is_a_fillet(&self, e: &Shape) -> bool {
        let mut i = 1usize;
        while i <= self.fillets.len() {
            let current_edge = &self.fillets[i - 1]; // TopoDS::Edge(fillets.Value(i))
            if current_edge.is_same(e) {
                return true;
            }
            i += 1;
        }
        false
    } // IsAFillet

    // =======================================================================
    // OCCT ChFi2d_Builder.cxx L1251-1264 — IsAChamfer
    // =======================================================================

    /// OCCT ChFi2d_Builder::IsAChamfer(const TopoDS_Edge& E) const
    /// (L1251-1264).
    pub fn is_a_chamfer(&self, e: &Shape) -> bool {
        let mut i = 1usize;
        while i <= self.chamfers.len() {
            let current_edge = &self.chamfers[i - 1]; // TopoDS::Edge(chamfers.Value(i))
            if current_edge.is_same(e) {
                return true;
            }
            i += 1;
        }
        false
    } // IsAChamfer

    // =======================================================================
    // OCCT ChFi2d_Builder.lxx L26-31 — Result()
    // =======================================================================

    /// OCCT ChFi2d_Builder::Result() const (lxx L26-31): returns the
    /// modified face with refFace's orientation.
    pub fn result(&self) -> Shape {
        let mut a_face = self.new_face.clone();
        a_face.orientation = self.ref_face.orientation;
        a_face
    }

    // =======================================================================
    // OCCT ChFi2d_Builder.lxx L35-38 — IsModified()
    // =======================================================================

    /// OCCT ChFi2d_Builder::IsModified(const TopoDS_Edge& E) const
    /// (lxx L35-38): history.IsBound(E).
    pub fn is_modified(&self, e: &Shape) -> bool {
        self.history.contains_key(&e.ptr_id())
    }

    // =======================================================================
    // OCCT ChFi2d_Builder.lxx L42-45 — FilletEdges()
    // =======================================================================

    /// OCCT ChFi2d_Builder::FilletEdges() const (lxx L42-45): returns the
    /// list of new edges.
    pub fn fillet_edges(&self) -> &[Shape] {
        &self.fillets
    }

    // =======================================================================
    // OCCT ChFi2d_Builder.lxx L49-52 — ChamferEdges()
    // =======================================================================

    /// OCCT ChFi2d_Builder::ChamferEdges() const (lxx L49-52): returns the
    /// list of new edges.
    pub fn chamfer_edges(&self) -> &[Shape] {
        &self.chamfers
    }

    // =======================================================================
    // OCCT ChFi2d_Builder.lxx L56-59 — NbFillet()
    // =======================================================================

    /// OCCT ChFi2d_Builder::NbFillet() const (lxx L56-59).
    pub fn nb_fillet(&self) -> usize {
        self.fillets.len() // fillets.Length()
    }

    // =======================================================================
    // OCCT ChFi2d_Builder.lxx L63-66 — NbChamfer()
    // =======================================================================

    /// OCCT ChFi2d_Builder::NbChamfer() const (lxx L63-66).
    pub fn nb_chamfer(&self) -> usize {
        self.chamfers.len() // chamfers.Length()
    }

    // =======================================================================
    // OCCT ChFi2d_Builder.lxx L70-73 — HasDescendant()
    // =======================================================================

    /// OCCT ChFi2d_Builder::HasDescendant(const TopoDS_Edge& E) const
    /// (lxx L70-73): history.IsBound(E).
    pub fn has_descendant(&self, e: &Shape) -> bool {
        self.history.contains_key(&e.ptr_id())
    }

    // =======================================================================
    // OCCT ChFi2d_Builder.lxx L77-80 — DescendantEdge()
    // =======================================================================

    /// OCCT ChFi2d_Builder::DescendantEdge(const TopoDS_Edge& E) const
    /// (lxx L77-80): TopoDS::Edge(history.Find(E)).
    pub fn descendant_edge(&self, e: &Shape) -> Shape {
        self.history[&e.ptr_id()].1.clone()
    }

    // =======================================================================
    // OCCT ChFi2d_Builder.lxx L84-87 — Status()
    // =======================================================================

    /// OCCT ChFi2d_Builder::Status() const (lxx L84-87).
    pub fn status(&self) -> ChFi2dConstructionError {
        self.status
    }
}
