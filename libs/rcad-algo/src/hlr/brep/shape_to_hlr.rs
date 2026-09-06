// OCCT HLRBRep_ShapeToHLR (TKHLR/HLRBRep/HLRBRep_ShapeToHLR.hxx L1-68
// + .cxx L1-313) — the OutLiner+Projector to Data loader.
//
// Translation notes:
// - OCCT NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> is
//   the insertion-ordered Vec<Shape> keyed by Shape::ptr_id (the
//   HLRTopoBRep_Data / HLRBRep_Data precedent); FindIndex returns the
//   1-based position, 0 when absent.
// - OCCT TopExp_Explorer / TopExp::MapShapes / MapShapesAndAncestors /
//   TopExp::Vertices are translated locally (TopExp_Explorer.cxx L34-244,
//   TopExp.cxx L39-51, L84-120, L214-251) — the out_liner.rs precedent:
//   rcad has no shared TopExp module and the landed helper set there is
//   private to that file.
// - The OCCT handle<HLRBRep_Data> return maps to Box<Data<'static>> (the
//   InternalAlgo storage shape — the myDS unique-ownership Box).
// - The kernel-context deviations (no OCCT statement): the BRep arena
//   travels as an explicit parameter (the out_liner::fill precedent) and
//   the per-face Set / per-edge curve loads need a 'static BRep view —
//   one clone is leaked per Load call (the Data<'static> leak route
//   sanctioned for this stage), and every edge curve adaptor is leaked
//   next to it (the &'static CurveView of the Data object graph).
// - OCCT BRep_Tool::Continuity reads the regularity stored on the edge
//   curve representations (the BRep_CurveOn2Surfaces records); the rcad
//   equivalent is the CurveOn2Surfaces representation written through the
//   BRepBuilder::continuity (the BRep_Builder::Continuity form).

use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::geom::Curve3;
use rcad_kernel::topods::{
    compose_pcurve_location, face_surface_value, tshape_flags, CurveRepresentation, GeomAbsShape,
    BRep, BRepTool, Orientation, Shape, ShapeType,
};

use crate::hlr::algo::projector::Projector;
use crate::hlr::brep::b_curve_tool::CurveView;
use crate::hlr::brep::curve::Curve;
use crate::hlr::brep::data::Data;
use crate::hlr::topo_brep::out_liner::OutLiner;
use crate::shhealing::shape_build::brep_tool::raw_subshapes;
use crate::topalgo::brep_top_adaptor::tool::BRepTopAdaptorTool;

/// OCCT TopAbs_ShapeEnum ranking (TopAbs_ShapeEnum.hxx L48-58): COMPOUND = 0
/// < ... < VERTEX = 7 < SHAPE = 8 — "more complex" is the LOWER value.  The
/// rcad ShapeType enum is declared in the opposite order, so the OCCT value
/// is mapped explicitly.
fn occt_shape_rank(t: ShapeType) -> i32 {
    match t {
        ShapeType::Compound => 0,
        ShapeType::CompSolid => 1,
        ShapeType::Solid => 2,
        ShapeType::Shell => 3,
        ShapeType::Face => 4,
        ShapeType::Wire => 5,
        ShapeType::Edge => 6,
        ShapeType::Vertex => 7,
        ShapeType::Shape => 8,
    }
}

/// OCCT NCollection_IndexedMap::Contains — the TShape-pointer identity key.
fn indexed_map_contains(map: &[Shape], s: &Shape) -> bool {
    let id = s.ptr_id();
    map.iter().any(|k| k.ptr_id() == id)
}

/// OCCT NCollection_IndexedMap::Add — append when not yet present.
fn indexed_map_add(map: &mut Vec<Shape>, s: &Shape) {
    if !indexed_map_contains(map, s) {
        map.push(s.clone());
    }
}

/// OCCT NCollection_IndexedMap::FindIndex — the 1-based position of the
/// key, 0 when absent.
fn indexed_map_find_index(map: &[Shape], s: &Shape) -> i32 {
    let id = s.ptr_id();
    map.iter()
        .position(|k| k.ptr_id() == id)
        .map(|p| p as i32 + 1)
        .unwrap_or(0)
}

/// OCCT NCollection_IndexedDataMap::FindIndex + Add (TopExp.cxx L98-102,
/// L113-117): the 1-based position of the key, appending the key with an
/// empty list when not yet present.
fn ve_map_pos_or_add(map: &mut Vec<(Shape, Vec<Shape>)>, s: &Shape) -> usize {
    let id = s.ptr_id();
    if let Some(pos) = map.iter().position(|(k, _)| k.ptr_id() == id) {
        pos
    } else {
        map.push((s.clone(), Vec::new()));
        map.len() - 1
    }
}

/// OCCT NCollection_IndexedDataMap::FindIndex — the 1-based position, 0
/// when absent.
fn ve_map_find_index(map: &[(Shape, Vec<Shape>)], s: &Shape) -> i32 {
    let id = s.ptr_id();
    map.iter()
        .position(|(k, _)| k.ptr_id() == id)
        .map(|p| p as i32 + 1)
        .unwrap_or(0)
}

/// OCCT NCollection_Map::Add — true when the key was not yet present
/// (ExploreShape cxx L299 / L308).
fn shape_map_add(shape_map: &mut Vec<Shape>, s: &Shape) -> bool {
    let id = s.ptr_id();
    if shape_map.iter().any(|x| x.ptr_id() == id) {
        false
    } else {
        shape_map.push(s.clone());
        true
    }
}

/// OCCT Standard_Real.hxx Epsilon(V) (L242-246): the distance to the next
/// representable value away from zero (nextafter toward RealLast for
/// V >= 0, toward RealFirst for V < 0).  For V = RealFirst / RealLast the
/// result overflows to +infinity (the (float) cast keeps it).
fn epsilon(v: f64) -> f64 {
    // the IEEE bit-increment walks away from zero on both signs.
    let na = f64::from_bits(v.to_bits() + 1);
    if v >= 0.0 {
        na - v
    } else {
        v - na
    }
}

/// OCCT BRep_Tool::Continuity(E, F1, F2) (BRep_Tool.cxx L1180-1188) ->
/// Continuity(E, S1, S2, L1, L2) (BRep_Tool.cxx L1223-1246) — the regularity
/// stored on the (E, F1/F2) curve representations.  The OCCT default is
/// GeomAbs_C0 when no representation matches (L1245).
fn brep_tool_continuity(brep: &BRep, e: &Shape, f1: &Shape, f2: &Shape) -> GeomAbsShape {
    // OCCT L1184-1186: S1 = Surface(F1, l1); S2 = Surface(F2, l2).
    let (Some(s1), Some(s2)) = (
        face_surface_value(brep, f1),
        face_surface_value(brep, f2),
    ) else {
        // the OCCT null-surface case has no matching representation.
        return GeomAbsShape::C0;
    };
    // OCCT L1229-1230: l1 = L1.Predivided(E.Location());
    // l2 = L2.Predivided(E.Location());
    let l1 = compose_pcurve_location(f1.location, e.location, &brep.locations);
    let l2 = compose_pcurve_location(f2.location, e.location, &brep.locations);

    // the representation scan (OCCT L1236-1244).
    let ed = match e.as_edge() {
        Some(ed) => ed,
        None => return GeomAbsShape::C0,
    };
    for cr in &ed.representations {
        // OCCT L1239: cr->IsRegularity(S1, S2, l1, l2).
        if cr.is_regularity_on(&s1, &s2, l1, l2) {
            // OCCT L1241: return cr->Continuity();
            if let CurveRepresentation::CurveOn2Surfaces { continuity, .. } = cr {
                return *continuity;
            }
        }
    }
    // OCCT L1245: return GeomAbs_C0;
    GeomAbsShape::C0
}

/// OCCT BRep_Tool::IsClosed(E, F) (BRep_Tool.cxx L795-806 + L812-860) —
/// true when the edge has two pcurves on the face: the plane shortcut
/// (IsPlane(S) -> false) then the CurveOnClosedSurface representation
/// scan.  The rcad TEdgeData records the closed-surface pcurve pair as the
/// CurveOnClosedSurface representation keyed by the owning face; the
/// triangulation branch has no rcad reader (no polygon-on-triangulation).
fn brep_tool_is_closed_edge_face(e: &Shape, f: &Shape) -> bool {
    let ed = match e.as_edge() {
        Some(ed) => ed,
        None => return false,
    };
    let fd = match f.as_face() {
        Some(fd) => fd,
        None => return false,
    };
    // IsClosed(E, S, L): if (IsPlane(S)) return false;
    if let Some(rcad_kernel::geom::Surface3::Plane(_)) = fd.surface.as_ref() {
        return false;
    }
    // the representation scan: a CurveOnClosedSurface for this face.
    ed.representations.iter().any(|cr| match cr {
        rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface { face, .. } => {
            face.0 == f.ptr_id()
        }
        _ => false,
    })
}

/// OCCT BRepTools::IsReallyClosed(E, F) (BRepTools.cxx L1204-1221).
fn brep_tools_is_really_closed(brep: &mut BRep, e: &Shape, f: &Shape) -> bool {
    // if (!BRep_Tool::IsClosed(E, F)) return false;
    if !brep_tool_is_closed_edge_face(e, f) {
        return false;
    }
    // int nbocc = 0; for (exp.Init(F, TopAbs_EDGE); ...) if
    // (exp.Current().IsSame(E)) nbocc++;
    let mut nbocc = 0;
    let mut exp = TopExpExplorer::new();
    exp.init(brep, f, ShapeType::Edge, ShapeType::Shape);
    while exp.more() {
        if exp.current().is_same(e) {
            nbocc += 1;
        }
        exp.next(brep);
    }
    // return nbocc == 2;
    nbocc == 2
}

/// OCCT TopExp_Explorer (TopExp_Explorer.cxx L34-244) — the stack explorer
/// with ToFind/ToAvoid.  Each stack level is a TopoDS_Iterator
/// (TopoDS_Iterator.cxx L28-70, the defaults cumOri = cumLoc = true): the
/// composed children snapshot and the current position.  Local copy of the
/// out_liner.rs translation (private there).
struct TopExpExplorer {
    /// OCCT toFind.
    to_find: ShapeType,
    /// OCCT toAvoid.
    to_avoid: ShapeType,
    /// OCCT hasMore.
    has_more: bool,
    /// OCCT myShape — Current() while the stack is empty.
    my_shape: Shape,
    /// OCCT myStack[0..=myStackTop]: (the composed children, the position).
    my_stack: Vec<(Vec<Shape>, usize)>,
}

impl TopExpExplorer {
    /// OCCT TopExp_Explorer() (cxx L39-45): toFind = toAvoid = TopAbs_SHAPE,
    /// hasMore = false.
    fn new() -> Self {
        TopExpExplorer {
            to_find: ShapeType::Shape,
            to_avoid: ShapeType::Shape,
            has_more: false,
            my_shape: Shape::null(),
            my_stack: Vec::new(),
        }
    }

    /// OCCT shouldAvoid (cxx L25-28).
    fn should_avoid(&self, ty: ShapeType) -> bool {
        self.to_avoid != ShapeType::Shape && self.to_avoid == ty
    }

    /// OCCT pushIterator (cxx L182-193) with the TopoDS_Iterator(S)
    /// construction (TopoDS_Iterator.cxx L28-70): the children carry the
    /// composed orientation (TopAbs::Compose) and location (Move).
    fn push_iterator(&mut self, brep: &mut BRep, s: &Shape) {
        let mut children = raw_subshapes(brep, s);
        for c in children.iter_mut() {
            c.orientation = s.orientation.compose(c.orientation);
            if s.location != 0 {
                // OCCT updateCurrentShape: myShape.Move(myLocation, false).
                let composed = brep.get_location(s.location) * brep.get_location(c.location);
                c.location = brep.add_location(composed);
            }
        }
        self.my_stack.push((children, 0));
    }

    /// OCCT Init (cxx L66-108).
    fn init(&mut self, brep: &mut BRep, s: &Shape, to_find: ShapeType, to_avoid: ShapeType) {
        // OCCT Clear() (cxx L195-204).
        self.my_stack.clear();
        self.has_more = false;
        self.my_shape = s.clone();
        self.to_find = to_find;
        self.to_avoid = to_avoid;

        if s.is_null() {
            self.has_more = false;
            return;
        }

        if to_find == ShapeType::Shape {
            self.has_more = false;
        } else {
            let ty = s.shape_type();
            if occt_shape_rank(ty) > occt_shape_rank(to_find) {
                self.has_more = false;
            } else if ty != to_find {
                self.has_more = true;
                self.next(brep);
            } else {
                self.has_more = true;
            }
        }
    }

    /// OCCT More (TopExp_Explorer.hxx L155-158).
    fn more(&self) -> bool {
        self.has_more
    }

    /// OCCT Current (cxx L156-161): myStackTop < 0 ? myShape :
    /// myStack[myStackTop].Value().
    fn current(&self) -> &Shape {
        if self.my_stack.is_empty() {
            &self.my_shape
        } else {
            let top = self.my_stack.last().expect("TopExp_Explorer: stack");
            &top.0[top.1]
        }
    }

    /// OCCT Next (cxx L110-154).
    fn next(&mut self, brep: &mut BRep) {
        if self.my_stack.is_empty() {
            let ty = self.my_shape.shape_type();
            if ty == self.to_find {
                self.has_more = false;
                return;
            } else if self.should_avoid(ty) {
                self.has_more = false;
                return;
            } else {
                let s = self.my_shape.clone();
                self.push_iterator(brep, &s);
            }
        } else {
            // myStack[myStackTop].Next().
            let top = self.my_stack.len() - 1;
            self.my_stack[top].1 += 1;
        }

        loop {
            let top = self.my_stack.len() - 1;
            if self.my_stack[top].1 < self.my_stack[top].0.len() {
                let shap_top = self.my_stack[top].0[self.my_stack[top].1].clone();
                let ty = shap_top.shape_type();
                if ty == self.to_find {
                    self.has_more = true;
                    return;
                } else if occt_shape_rank(self.to_find) > occt_shape_rank(ty)
                    && !self.should_avoid(ty)
                {
                    // pushIterator(TopoDS_Iterator(aShapTop)) — aTopIter
                    // reference invalid after the push (cxx L152-155).
                    self.push_iterator(brep, &shap_top);
                } else {
                    self.my_stack[top].1 += 1;
                }
            } else {
                // popIterator (cxx L196-201).
                self.my_stack.pop();
                if self.my_stack.is_empty() {
                    break;
                }
                let top = self.my_stack.len() - 1;
                self.my_stack[top].1 += 1;
            }
        }
        self.has_more = false;
    }
}

/// OCCT TopExp::MapShapes(S, T, M) (TopExp.cxx L39-51) — the indexed map of
/// the sub-shapes of type T.
fn map_shapes(brep: &mut BRep, s: &Shape, t: ShapeType, m: &mut Vec<Shape>) {
    let mut ex = TopExpExplorer::new();
    ex.init(brep, s, t, ShapeType::Shape);
    while ex.more() {
        let cur = ex.current().clone();
        indexed_map_add(m, &cur);
        ex.next(brep);
    }
}

/// OCCT TopExp::MapShapesAndAncestors (TopExp.cxx L84-120) — the indexed
/// (sub-shape -> ancestors) map; rcad encodes the
/// NCollection_IndexedDataMap as Vec<(Shape, Vec<Shape>)> (insertion order,
/// ptr_id keys).
fn map_shapes_and_ancestors(
    brep: &mut BRep,
    s: &Shape,
    ts: ShapeType,
    ta: ShapeType,
    m: &mut Vec<(Shape, Vec<Shape>)>,
) {
    // visit ancestors (cxx L89-107).
    let mut exa = TopExpExplorer::new();
    exa.init(brep, s, ta, ShapeType::Shape);
    while exa.more() {
        // const TopoDS_Shape& anc = exa.Current();
        let anc = exa.current().clone();
        let mut exs = TopExpExplorer::new();
        exs.init(brep, &anc, ts, ShapeType::Shape);
        while exs.more() {
            // int index = M.FindIndex(exs.Current()); if (index == 0)
            // index = M.Add(exs.Current(), empty); M(index).Append(anc);
            let cur = exs.current().clone();
            let pos = ve_map_pos_or_add(m, &cur);
            m[pos].1.push(anc.clone());
            exs.next(brep);
        }
        exa.next(brep);
    }

    // visit shapes not under ancestors (cxx L109-119).
    let mut ex = TopExpExplorer::new();
    ex.init(brep, s, ts, ta);
    while ex.more() {
        // int index = M.FindIndex(ex.Current()); if (index == 0)
        // index = M.Add(ex.Current(), empty);
        let cur = ex.current().clone();
        ve_map_pos_or_add(m, &cur);
        ex.next(brep);
    }
}

/// OCCT TopExp::Vertices(E, Vfirst, Vlast) (TopExp.cxx L214-251, the default
/// CumOri = false): the FORWARD child is Vfirst, the REVERSED child is Vlast,
/// both nullified when undefined.
fn top_exp_vertices(e: &Shape) -> (Shape, Shape) {
    let mut v_first = Shape::null();
    let mut v_last = Shape::null();
    if let Some(ed) = e.as_edge() {
        for v in [&ed.first, &ed.last] {
            if v.orientation == Orientation::Forward {
                v_first = v.clone();
            } else if v.orientation == Orientation::Reversed {
                v_last = v.clone();
            }
        }
    }
    (v_first, v_last)
}

/// OCCT BRepAdaptor_Curve(E) (BRepAdaptor_Curve over the TopoDS edge) seen
/// through the [`CurveView`] interface of HLRBRep_BCurveTool — the
/// ChangeGeometry().Curve(EG) load of HLRBRep_EdgeData::Set (EdgeData.cxx
/// L41).  The adaptor reads the rcad TEdgeData 3D curve + range; the
/// per-type resolution follows the bop BRepAdaptorCurve
/// (curve_range::curve_resolution); IsClosed / IsPeriodic read the conic
/// and B-spline periodicity.
struct EdgeCurveAdaptor {
    curve: Curve3,
    first: f64,
    last: f64,
}

impl EdgeCurveAdaptor {
    /// OCCT BRepAdaptor_Curve::FirstParameter/LastParameter — the edge
    /// trimmed range.
    fn edge_view(brep: &BRep, e: &Shape) -> &'static dyn CurveView {
        let ed = e
            .as_edge()
            .expect("Standard_NoSuchObject: BRepAdaptor_Curve: not an edge");
        let curve = ed
            .curve
            .clone()
            .expect("StdFail_NotDone: BRepAdaptor_Curve: edge without 3D curve");
        let range = brep.edge_range(e);
        Box::leak(Box::new(EdgeCurveAdaptor {
            curve,
            first: range[0],
            last: range[1],
        }))
    }
}

impl CurveView for EdgeCurveAdaptor {
    fn first_parameter(&self) -> f64 {
        self.first
    }
    fn last_parameter(&self) -> f64 {
        self.last
    }
    fn d0(&self, u: f64) -> rcad_kernel::geom::Point3 {
        use rcad_kernel::geom::CurveEval;
        self.curve.point_at(u)
    }
    fn d1(&self, u: f64) -> (rcad_kernel::geom::Point3, rcad_kernel::geom::Vec3) {
        use rcad_kernel::geom::CurveEval;
        (self.curve.point_at(u), self.curve.derivative_at(u))
    }
    fn d2(
        &self,
        u: f64,
    ) -> (
        rcad_kernel::geom::Point3,
        rcad_kernel::geom::Vec3,
        rcad_kernel::geom::Vec3,
    ) {
        use rcad_kernel::geom::CurveEval;
        (
            self.curve.point_at(u),
            self.curve.derivative_at(u),
            self.curve.derivative2_at(u),
        )
    }
    fn get_type(&self) -> CurveType {
        match &self.curve {
            Curve3::Line(_) => CurveType::Line,
            Curve3::Circle(_) => CurveType::Circle,
            Curve3::Ellipse(_) => CurveType::Ellipse,
            Curve3::Hyperbola(_) => CurveType::Hyperbola,
            Curve3::Parabola(_) => CurveType::Parabola,
            Curve3::Bezier(_) => CurveType::Bezier,
            Curve3::BSpline(_) => CurveType::BSpline,
            _ => CurveType::Other,
        }
    }
    fn line(&self) -> rcad_kernel::geom::Line3 {
        match &self.curve {
            Curve3::Line(l) => *l,
            _ => panic!("Standard_NoSuchObject: BRepAdaptor_Curve::Line"),
        }
    }
    fn circle(&self) -> rcad_kernel::geom::Circle3 {
        match &self.curve {
            Curve3::Circle(c) => *c,
            _ => panic!("Standard_NoSuchObject: BRepAdaptor_Curve::Circle"),
        }
    }
    fn ellipse(&self) -> rcad_kernel::geom::Ellipse3 {
        match &self.curve {
            Curve3::Ellipse(e) => *e,
            _ => panic!("Standard_NoSuchObject: BRepAdaptor_Curve::Ellipse"),
        }
    }
    fn degree(&self) -> i32 {
        match &self.curve {
            Curve3::BSpline(b) => b.degree as i32,
            Curve3::Bezier(b) => b.control_points.len() as i32 - 1,
            _ => 0,
        }
    }
    fn nb_poles(&self) -> i32 {
        match &self.curve {
            Curve3::BSpline(b) => b.control_points.len() as i32,
            Curve3::Bezier(b) => b.control_points.len() as i32,
            _ => 0,
        }
    }
    fn nb_knots(&self) -> i32 {
        match &self.curve {
            Curve3::BSpline(b) => run_length_knots(&b.knots).len() as i32,
            _ => 0,
        }
    }
    fn is_closed(&self) -> bool {
        // OCCT BRepAdaptor_Curve::IsClosed — the periodic-curve seam test;
        // the conic/BSpline periodicity is the rcad carrier.
        self.is_periodic_curve()
            && (self.last - self.first) >= self.period() - 1.0e-10
    }
    fn is_periodic(&self) -> bool {
        self.is_periodic_curve()
    }
    fn period(&self) -> f64 {
        match &self.curve {
            Curve3::Circle(_) | Curve3::Ellipse(_) => 2.0 * std::f64::consts::PI,
            Curve3::BSpline(b) => {
                if b.is_periodic {
                    self.last - self.first
                } else {
                    0.0
                }
            }
            _ => 0.0,
        }
    }
    fn resolution(&self, r3d: f64) -> f64 {
        crate::bop::int_tools::curve_range::curve_resolution(
            &self.curve,
            0.5 * (self.first + self.last),
            r3d,
        )
    }
    fn parameter_3d(&self, p2d: f64) -> f64 {
        p2d
    }
    fn poles(&self) -> Vec<rcad_kernel::geom::Point3> {
        match &self.curve {
            Curve3::BSpline(b) => b.control_points.clone(),
            Curve3::Bezier(b) => b.control_points.clone(),
            _ => Vec::new(),
        }
    }
    fn knots(&self) -> Vec<f64> {
        match &self.curve {
            // OCCT Geom_BSplineCurve::Knots() — the distinct knots; the rcad
            // carrier stores the flat vector.
            Curve3::BSpline(b) => run_length_knots(&b.knots).into_iter().map(|(k, _)| k).collect(),
            _ => Vec::new(),
        }
    }
    fn multiplicities(&self) -> Vec<i32> {
        match &self.curve {
            // OCCT Geom_BSplineCurve::Multiplicities() — per distinct knot.
            Curve3::BSpline(b) => {
                run_length_knots(&b.knots).into_iter().map(|(_, m)| m).collect()
            }
            _ => Vec::new(),
        }
    }
    fn weights(&self) -> Vec<f64> {
        match &self.curve {
            Curve3::BSpline(b) => b.weights.clone(),
            Curve3::Bezier(b) => b.weights.clone(),
            _ => Vec::new(),
        }
    }
}

/// The OCCT-style (knot, multiplicity) pairs from the rcad flat knot vector
/// (each knot repeated its multiplicity — the expand_knots precedent).
fn run_length_knots(flat: &[f64]) -> Vec<(f64, i32)> {
    let mut out: Vec<(f64, i32)> = Vec::new();
    for &k in flat {
        if let Some(last) = out.last_mut() {
            if last.0 == k {
                last.1 += 1;
                continue;
            }
        }
        out.push((k, 1));
    }
    out
}

impl EdgeCurveAdaptor {
    fn is_periodic_curve(&self) -> bool {
        match &self.curve {
            Curve3::Circle(_) | Curve3::Ellipse(_) => true,
            Curve3::BSpline(b) => b.is_periodic,
            _ => false,
        }
    }
}

/// OCCT HLRBRep_ShapeToHLR::Load (hxx L45-49 + cxx L40-174) — creates a
/// DataStructure containing the OutLiner <s> depending on the projector <p>
/// and nbIso.
///
/// rcad deviations (reported): the BRep builder arena travels as the
/// leading parameter (the out_liner::fill precedent) and S is taken by
/// &mut — OCCT mutates the OutLiner through the handle (S->Fill, cxx L46).
/// The return is the OCCT handle as the unique-ownership Box<Data<'static>>
/// (the InternalAlgo storage shape).
pub fn load(
    brep: &mut BRep,
    s: &mut OutLiner,
    p: &Projector,
    mst: &mut Vec<(Shape, BRepTopAdaptorTool)>,
    nb_iso: usize,
) -> Box<Data<'static>> {
    // OCCT L46: S->Fill(P, MST, nbIso);
    s.fill(brep, p, mst, nb_iso);

    // The 'static BRep view for the HLRBRep_Data object graph (the
    // FaceData::Set / HLRBRep_Surface loads) — one leaked clone per Load
    // (see the module notes).
    let brep_st: &'static BRep = Box::leak(Box::new(brep.clone()));

    // OCCT L48: HLRTopoBRep_Data& TopDS = S->DataStructure();
    // (borrowed from <s> at each use site below.)
    // OCCT L49-50: NCollection_IndexedMap FM; NCollection_IndexedMap EM;
    let mut fm: Vec<Shape> = Vec::new();
    let mut em: Vec<Shape> = Vec::new();
    // OCCT L51-54: NCollection_IndexedDataMap VerticesToEdges;
    // NCollection_IndexedDataMap EdgesToFaces;
    let mut vertices_to_edges: Vec<(Shape, Vec<Shape>)> = Vec::new();
    let mut edges_to_faces: Vec<(Shape, Vec<Shape>)> = Vec::new();

    // OCCT L56: TopExp_Explorer exshell, exface;
    let mut exshell = TopExpExplorer::new();
    let mut exface = TopExpExplorer::new();

    // OCCT L58: for (exshell.Init(S->OutLinedShape(), TopAbs_SHELL); ...)
    let out_lined = s.out_lined_shape().clone();
    exshell.init(brep, &out_lined, ShapeType::Shell, ShapeType::Shape);
    while exshell.more() {
        // faces in a shell
        let shell = exshell.current().clone();

        // OCCT L61: for (exface.Init(exshell.Current(), TopAbs_FACE); ...)
        exface.init(brep, &shell, ShapeType::Face, ShapeType::Shape);
        while exface.more() {
            let cur = exface.current().clone();
            // OCCT L63-66: if (!FM.Contains(exface.Current())) FM.Add(...);
            indexed_map_add(&mut fm, &cur);
            exface.next(brep);
        }
        exshell.next(brep);
    }

    // OCCT L70: for (exface.Init(S->OutLinedShape(), TopAbs_FACE,
    // TopAbs_SHELL); ...) — faces not in a shell
    exface.init(brep, &out_lined, ShapeType::Face, ShapeType::Shell);
    while exface.more() {
        let cur = exface.current().clone();
        if !indexed_map_contains(&fm, &cur) {
            indexed_map_add(&mut fm, &cur);
        }
        exface.next(brep);
    }

    // OCCT L78: TopExp::MapShapes(S->OutLinedShape(), TopAbs_EDGE, EM);
    map_shapes(brep, &out_lined, ShapeType::Edge, &mut em);

    // OCCT L80-86: for (i = 1; i <= nbEdge; i++) — vertices back to edges
    let nb_edge = em.len();
    let mut i: usize = 1;
    while i <= nb_edge {
        let e = em[i - 1].clone();
        map_shapes_and_ancestors(
            brep,
            &e,
            ShapeType::Vertex,
            ShapeType::Edge,
            &mut vertices_to_edges,
        );
        i += 1;
    }

    // OCCT L88-89: int nbVert = VerticesToEdges.Extent(); int nbFace = FM.Extent();
    let nb_vert = vertices_to_edges.len();
    let nb_face = fm.len();

    // OCCT L100: occ::handle<HLRBRep_Data> DS = new HLRBRep_Data(nbVert, nbEdge, nbFace);
    let mut ds: Data<'static> = Data::new(nb_vert, nb_edge, nb_face);
    // rcad: the Data kernel context rides with the Load (the [Data::set_brep]
    // deviation accessor) — the leaked BRep view above is the per-Load image
    // of the session kernel context the OCCT Data would read through the
    // global TShape graph.
    ds.set_brep(brep_st);

    // OCCT L101-105: HLRBRep_EdgeData* ed = nullptr;
    // if (nbEdge != 0) { ed = &(DS->EDataArray().ChangeValue(1)); }
    // (the ChangeValue(1) slot with the ed++ walk maps onto the 0-based
    // iter_mut below — the first slot is 1-based edge 1.)
    // OCCT L106: //  ed++;

    // OCCT L108-111: for (i = 1; i <= nbFace; i++) — test of Double edges
    let mut i: usize = 1;
    while i <= nb_face {
        let f = fm[i - 1].clone();
        map_shapes_and_ancestors(
            brep,
            &f,
            ShapeType::Edge,
            ShapeType::Face,
            &mut edges_to_faces,
        );
        i += 1;
    }

    // OCCT L113-170: for (i = 1; i <= nbEdge; i++) — load the Edges
    let mut i: usize = 1;
    while i <= nb_edge {
        // OCCT L115: const TopoDS_Edge& Edg = TopoDS::Edge(EM(i));
        let edg = em[i - 1].clone();
        // OCCT L116: TopExp::Vertices(Edg, VF, VL);
        let (vf, vl) = top_exp_vertices(&edg);
        // OCCT L117: BRep_Tool::Range(Edg, pf, pl);
        let range = brep.edge_range(&edg);
        let mut pf = range[0];
        let mut pl = range[1];
        // OCCT L118-119: bool reg1 = false; bool regn = false;
        let mut reg1 = false;
        let mut regn = false;
        // OCCT L120: int inde = EdgesToFaces.FindIndex(Edg);
        let inde = ve_map_find_index(&edges_to_faces, &edg);
        if inde > 0 {
            // OCCT L123: if (EdgesToFaces(inde).Extent() == 2)
            let list = &edges_to_faces[(inde - 1) as usize].1;
            if list.len() == 2 {
                // OCCT L125-128: itn = EdgesToFaces(inde); F1 = itn.Value();
                // itn.Next(); F2 = itn.Value();
                let f1 = list[0].clone();
                let f2 = list[1].clone();
                // OCCT L129-131: GeomAbs_Shape rg = BRep_Tool::Continuity(...);
                // reg1 = rg >= GeomAbs_G1; regn = rg >= GeomAbs_G2;
                let rg = brep_tool_continuity(brep, &edg, &f1, &f2);
                reg1 = rg >= GeomAbsShape::G1;
                regn = rg >= GeomAbsShape::G2;
            }
        }

        // OCCT L135-149: if (VF.IsNull()) { i1 = 0; o1 = false; c1 = false;
        // pf = RealFirst(); tf = (float)Epsilon(pf); } else { ... }
        let i1: i32;
        let o1: bool;
        let c1: bool;
        let tf: f32;
        if vf.is_null() {
            i1 = 0;
            o1 = false;
            c1 = false;
            pf = f64::MIN; // RealFirst()
            tf = epsilon(pf) as f32;
        } else {
            i1 = ve_map_find_index(&vertices_to_edges, &vf);
            o1 = s.data_structure().is_out_v(&vf);
            c1 = s.data_structure().is_int_v(&vf);
            tf = brep.tolerance(&vf) as f32;
        }

        // OCCT L151-165: if (VL.IsNull()) { ... pl = RealLast();
        // tl = (float)Epsilon(pl); } else { ... }
        let i2: i32;
        let o2: bool;
        let c2: bool;
        let tl: f32;
        if vl.is_null() {
            i2 = 0;
            o2 = false;
            c2 = false;
            pl = f64::MAX; // RealLast()
            tl = epsilon(pl) as f32;
        } else {
            i2 = ve_map_find_index(&vertices_to_edges, &vl);
            o2 = s.data_structure().is_out_v(&vl);
            c2 = s.data_structure().is_int_v(&vl);
            tl = brep.tolerance(&vl) as f32;
        }

        // OCCT L167: ed->Set(reg1, regn, Edg, i1, i2, o1, o2, c1, c2, pf, tf, pl, tl);
        // — the rcad Set takes the loaded adaptor + the BRep_Tool tolerance
        // read at the kernel boundary (the EdgeData documented deferral):
        // ChangeGeometry().Curve(EG) lands first.
        let mut geometry = Curve::new();
        geometry.load(EdgeCurveAdaptor::edge_view(brep, &edg));
        let tolerance = brep.tolerance(&edg) as f32;
        ds.e_data_array_mut()[i - 1].set(
            reg1, regn, geometry, tolerance, i1, i2, o1, o2, c1, c2, pf, tf, pl, tl,
        );
        // OCCT L168: DS->EdgeMap().Add(Edg);
        indexed_map_add(ds.edge_map_mut(), &edg);
        // OCCT L169: ed++;
        i += 1;
    }

    // OCCT L172: ExploreShape(S, DS, FM, EM);
    explore_shape(brep, brep_st, s, &mut ds, &fm, &em);
    // OCCT L173: return DS;
    Box::new(ds)
}

/// OCCT HLRBRep_ShapeToHLR::ExploreFace (hxx L52-59 + cxx L178-240).
fn explore_face(
    brep: &mut BRep,
    brep_st: &'static BRep,
    s: &mut OutLiner,
    ds: &mut Data<'static>,
    fm: &[Shape],
    em: &[Shape],
    i: &mut usize,
    f: &Shape,
    closed: bool,
) {
    // OCCT L187: i++;
    *i += 1;
    let idx = *i; // the 1-based FM position (L190/L193 read FM(i)).
    // OCCT L188: TopExp_Explorer Ex1, Ex2;
    let mut ex1 = TopExpExplorer::new();
    let mut ex2 = TopExpExplorer::new();
    // OCCT L189: HLRTopoBRep_Data& TopDS = S->DataStructure();
    // (borrowed from <s> at the Int/Iso/Out reads below.)
    // OCCT L190: TopAbs_Orientation orient = FM(i).Orientation();
    let orient = fm[idx - 1].orientation;
    // OCCT L191-192: TopoDS_Face theFace = TopoDS::Face(FM(i));
    // theFace.Orientation(TopAbs_FORWARD);
    let mut the_face = fm[idx - 1].clone();
    the_face.orientation = Orientation::Forward;
    // OCCT L193: HLRBRep_FaceData& fd = DS->FDataArray().ChangeValue(i);
    let fd = &mut ds.f_data_array_mut()[idx - 1];

    // OCCT L195-200: int nw = 0; for (Ex1.Init(theFace, TopAbs_WIRE); ...) nw++;
    let mut nw: usize = 0;
    ex1.init(brep, &the_face, ShapeType::Wire, ShapeType::Shape);
    while ex1.more() {
        nw += 1;
        ex1.next(brep);
    }

    // OCCT L202: fd.Set(theFace, orient, closed, nw);
    fd.set(brep_st, &the_face, orient, closed, nw);

    // OCCT L203: nw = 0;
    let mut nw: usize = 0;

    // OCCT L205: for (Ex1.Init(theFace, TopAbs_WIRE); Ex1.More(); Ex1.Next())
    ex1.init(brep, &the_face, ShapeType::Wire, ShapeType::Shape);
    while ex1.more() {
        nw += 1;
        // OCCT L208: int ne = 0;
        let mut ne: usize = 0;

        // OCCT L210-217: count the non-degenerated edges of the wire.
        let wire = ex1.current().clone();
        ex2.init(brep, &wire, ShapeType::Edge, ShapeType::Shape);
        while ex2.more() {
            let an_edge = ex2.current().clone();
            if !an_edge
                .as_edge()
                .map(|ed| ed.degenerated)
                .unwrap_or(false)
            {
                ne += 1;
            }
            ex2.next(brep);
        }

        // OCCT L219: fd.SetWire(nw, ne);
        fd.set_wire(nw, ne);

        // OCCT L220: ne = 0;
        let mut ne: usize = 0;

        // OCCT L222-237: for (Ex2.Init(Ex1.Current(), TopAbs_EDGE); ...)
        ex2.init(brep, &wire, ShapeType::Edge, ShapeType::Shape);
        while ex2.more() {
            let e = ex2.current().clone();
            // OCCT L225-228: if (BRep_Tool::Degenerated(E)) continue;
            if e.as_edge().map(|ed| ed.degenerated).unwrap_or(false) {
                ex2.next(brep);
                continue;
            }
            ne += 1;
            // OCCT L230: int ie = EM.FindIndex(E);
            let ie = indexed_map_find_index(em, &e);
            // OCCT L231: TopAbs_Orientation anOrientE = E.Orientation();
            let an_orient_e = e.orientation;
            // OCCT L232-234: the TopDS reads take the ORIGINAL face F.
            let int = s.data_structure().is_int_l_face_edge(f, &e);
            let iso = s.data_structure().is_iso_l_face_edge(f, &e);
            let out = s.data_structure().is_out_l_face_edge(f, &e);
            // OCCT L235: bool Dbl = BRepTools::IsReallyClosed(E, theFace);
            let dbl = brep_tools_is_really_closed(brep, &e, &the_face);
            // OCCT L236: fd.SetWEdge(nw, ne, ie, anOrientE, Out, Int, Dbl, Iso);
            fd.set_w_edge(nw, ne, ie, an_orient_e, out, int, dbl, iso);
            ex2.next(brep);
        }
        ex1.next(brep);
    }
    // OCCT L239: DS->FaceMap().Add(theFace);
    indexed_map_add(ds.face_map_mut(), &the_face);
}

/// OCCT HLRBRep_ShapeToHLR::ExploreShape (hxx L61-65 + cxx L244-313).
fn explore_shape(
    brep: &mut BRep,
    brep_st: &'static BRep,
    s: &mut OutLiner,
    ds: &mut Data<'static>,
    fm: &[Shape],
    em: &[Shape],
) {
    // OCCT L250: NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> ShapeMap;
    let mut shape_map: Vec<Shape> = Vec::new();
    // OCCT L251: TopExp_Explorer exshell, exface, exedge;
    let mut exshell = TopExpExplorer::new();
    let mut exface = TopExpExplorer::new();
    let mut exedge = TopExpExplorer::new();
    // OCCT L252: int i = 0;
    let mut i: usize = 0;

    // OCCT L254: for (exshell.Init(S->OriginalShape(), TopAbs_SHELL); ...)
    let original = s.original_shape().clone();
    exshell.init(brep, &original, ShapeType::Shell, ShapeType::Shape);
    while exshell.more() {
        // faces in a shell (open or close)
        let cur = exshell.current().clone();

        // OCCT L257: bool closed = exshell.Current().Closed();
        let mut closed = cur
            .as_shell()
            .map(|sd| sd.flags & tshape_flags::CLOSED != 0)
            .unwrap_or(false);

        // OCCT L259-295: if (!closed) — recompute the closedness from the
        // FORWARD/REVERSED edge balance of the shell.
        if !closed {
            // OCCT L261-263: int ie; int nbEdge = EM.Extent();
            // int* flag = new int[nbEdge + 1];
            let nb_edge = em.len();
            let mut flag = vec![0i32; nb_edge + 1];

            // OCCT L265-268: for (ie = 1; ie <= nbEdge; ie++) flag[ie] = 0;

            // OCCT L270-286: for (exedge.Init(exshell.Current(), TopAbs_EDGE); ...)
            exedge.init(brep, &cur, ShapeType::Edge, ShapeType::Shape);
            while exedge.more() {
                let e = exedge.current().clone();
                // OCCT L273: ie = EM.FindIndex(E); — 0 lands in the unused
                // slot flag[0] for edges outside the OutLinedShape map.
                let ie = indexed_map_find_index(em, &e) as usize;
                // OCCT L274: TopAbs_Orientation orient = E.Orientation();
                let orient = e.orientation;
                // OCCT L275: if (!BRep_Tool::Degenerated(E))
                if !e.as_edge().map(|ed| ed.degenerated).unwrap_or(false) {
                    if orient == Orientation::Forward {
                        flag[ie] += 1;
                    } else if orient == Orientation::Reversed {
                        flag[ie] -= 1;
                    }
                }
                exedge.next(brep);
            }
            // OCCT L287: closed = true;
            closed = true;

            // OCCT L289-292: for (ie = 1; ie <= nbEdge && closed; ie++)
            // closed = (flag[ie] == 0);
            let mut ie: usize = 1;
            while ie <= nb_edge && closed {
                closed = flag[ie] == 0;
                ie += 1;
            }
            // OCCT L293-294: delete[] flag; flag = nullptr; (the Vec drops)
        }

        // OCCT L297-303: for (exface.Init(exshell.Current(), TopAbs_FACE); ...)
        exface.init(brep, &cur, ShapeType::Face, ShapeType::Shape);
        while exface.more() {
            let f = exface.current().clone();
            // OCCT L299-301: if (ShapeMap.Add(exface.Current()))
            // ExploreFace(S, DS, FM, EM, i, TopoDS::Face(exface.Current()), closed);
            if shape_map_add(&mut shape_map, &f) {
                explore_face(brep, brep_st, s, ds, fm, em, &mut i, &f, closed);
            }
            exface.next(brep);
        }
        exshell.next(brep);
    }

    // OCCT L306-312: for (exface.Init(S->OriginalShape(), TopAbs_FACE,
    // TopAbs_SHELL); ...) — faces not in a shell
    exface.init(brep, &original, ShapeType::Face, ShapeType::Shell);
    while exface.more() {
        let f = exface.current().clone();
        // OCCT L308-310: if (ShapeMap.Add(exface.Current()))
        // ExploreFace(S, DS, FM, EM, i, TopoDS::Face(exface.Current()), false);
        if shape_map_add(&mut shape_map, &f) {
            explore_face(brep, brep_st, s, ds, fm, em, &mut i, &f, false);
        }
        exface.next(brep);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
    use rcad_kernel::math::gp::Ax2;
    use rcad_kernel::topo::topods::BRepBuilder;
    use glam::DVec3;

    /// A shell of two unit squares (the z = 0 and z = 1 planes) — the
    /// "square-face box" fixture (the data/update.rs plane_brep pattern).
    /// Returns (brep, shell, face1, face2).
    fn square_shell_fixture() -> (BRep, Shape, Shape, Shape) {
        let mut brep = BRep::new();
        let mut b = BRepBuilder::new();

        // The z = 0 square: vertices (0,0), (1,0), (1,1), (0,1).
        let pts = [
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let vs: Vec<Shape> = pts.iter().map(|p| b.add_vertex(&mut brep, *p, 1e-7)).collect();
        let mut edges = Vec::new();
        for i in 0..4 {
            let a = pts[i];
            let bb = pts[(i + 1) % 4];
            let curve = Curve3::Line(Line3::new(a, (bb - a).normalize()));
            // OCCT BRep_Builder::Add(E, V) stores the start vertex FORWARD
            // and the end vertex REVERSED on the edge.
            let v_first = vs[i].clone();
            let mut v_last = vs[(i + 1) % 4].clone();
            v_last.orientation = Orientation::Reversed;
            edges.push(b.add_edge(&mut brep, Some(curve), v_first, v_last, [0.0, 1.0]));
        }
        let wire0 = brep.add_twire(edges.clone());
        let face1 = brep.add_tface(
            Some(Surface3::Plane(Plane {
                origin: DVec3::ZERO,
                normal: DVec3::new(0.0, 0.0, 1.0),
                u_dir: DVec3::new(1.0, 0.0, 0.0),
                v_dir: DVec3::new(0.0, 1.0, 0.0),
            })),
            wire0,
            Vec::new(),
            None,
            Some([0.0, 1.0, 0.0, 1.0]),
            Vec::new(),
            true,
        );

        // The z = 1 square.
        let pts1: Vec<DVec3> = pts.iter().map(|p| *p + DVec3::new(0.0, 0.0, 1.0)).collect();
        let vs1: Vec<Shape> = pts1.iter().map(|p| b.add_vertex(&mut brep, *p, 1e-7)).collect();
        let mut edges1 = Vec::new();
        for i in 0..4 {
            let a = pts1[i];
            let bb = pts1[(i + 1) % 4];
            let curve = Curve3::Line(Line3::new(a, (bb - a).normalize()));
            let v_first = vs1[i].clone();
            let mut v_last = vs1[(i + 1) % 4].clone();
            v_last.orientation = Orientation::Reversed;
            edges1.push(b.add_edge(&mut brep, Some(curve), v_first, v_last, [0.0, 1.0]));
        }
        let wire1 = brep.add_twire(edges1);
        let face2 = brep.add_tface(
            Some(Surface3::Plane(Plane {
                origin: DVec3::new(0.0, 0.0, 1.0),
                normal: DVec3::new(0.0, 0.0, 1.0),
                u_dir: DVec3::new(1.0, 0.0, 0.0),
                v_dir: DVec3::new(0.0, 1.0, 0.0),
            })),
            wire1,
            Vec::new(),
            None,
            Some([0.0, 1.0, 0.0, 1.0]),
            Vec::new(),
            true,
        );

        let shell = brep.add_tshell(vec![face1.clone(), face2.clone()]);
        (brep, shell, face1, face2)
    }

    /// The MST map pre-bound for the faces (what DSFiller::Insert leaves
    /// behind for already-visited faces).
    fn bound_mst(brep: &BRep, faces: [&Shape; 2]) -> Vec<(Shape, BRepTopAdaptorTool)> {
        let arc_brep = Arc::new(brep.clone());
        faces
            .into_iter()
            .map(|f| {
                (
                    f.clone(),
                    BRepTopAdaptorTool::new_face(arc_brep.clone(), f, 1.0e-7),
                )
            })
            .collect()
    }

    /// The identity top-view projector.
    fn top_view_projector() -> Projector {
        Projector::from_ax2(&Ax2::new(
            DVec3::ZERO,
            DVec3::new(0.0, 0.0, 1.0),
            DVec3::new(1.0, 0.0, 0.0),
        ))
    }

    /// OCCT anchor: Load (cxx L40-174) over the two-square shell — the
    /// counts (8 vertices / 8 edges / 2 faces), the EMap / FMap contents
    /// (the OutLinedShape copies), the EdgeData Set fills (VSta/VEnd, Rg1,
    /// Used, Selected) and the FaceData Wires structure
    /// (SetWire / SetWEdge of one wire of four edges per face).
    #[test]
    fn shape_to_hlr_load_square_box() {
        let (mut brep, shell, face1, face2) = square_shell_fixture();
        let mut mst = bound_mst(&brep, [&face1, &face2]);
        let proj = top_view_projector();

        let mut ol = OutLiner::with_original_shape(&shell);
        let ds = load(&mut brep, &mut ol, &proj, &mut mst, 0);

        // OCCT hand-calc: nbVert = 8 (4 per square), nbEdge = 8, nbFace = 2.
        assert_eq!(ds.nb_vertices(), 8);
        assert_eq!(ds.nb_edges(), 8);
        assert_eq!(ds.nb_faces(), 2);

        // The OutLinedShape: a compound of one shell of the two face copies
        // (the OutLiner::BuildShape precedent) — FM/FMap hold the copies.
        let result = ol.out_lined_shape().clone();
        assert_eq!(result.shape_type(), ShapeType::Compound);
        let children = raw_subshapes(&brep, &result);
        assert_eq!(children.len(), 1);
        let faces = raw_subshapes(&brep, &children[0]);
        assert_eq!(faces.len(), 2);

        // DS->FaceMap().Add(theFace) (ExploreFace cxx L239): the two
        // FORWARD copies, in shell order.
        assert_eq!(ds.face_map().len(), 2);
        for (k, f) in faces.iter().enumerate() {
            assert_eq!(ds.face_map()[k].ptr_id(), f.ptr_id());
            assert_eq!(ds.face_map()[k].orientation, Orientation::Forward);
        }

        // DS->EdgeMap().Add(Edg) (cxx L168): the 8 edges of the OutLined
        // faces — 4 per face in wire order.
        assert_eq!(ds.edge_map().len(), 8);
        for k in 0..4 {
            let e1 = raw_subshapes(&brep, &raw_subshapes(&brep, &faces[0])[0])[k].clone();
            assert_eq!(ds.edge_map()[k].ptr_id(), e1.ptr_id());
            let e2 = raw_subshapes(&brep, &raw_subshapes(&brep, &faces[1])[0])[k].clone();
            assert_eq!(ds.edge_map()[4 + k].ptr_id(), e2.ptr_id());
        }

        // The EdgeData Set fills (cxx L167): the first edge of face1 runs
        // (0,0,0) -> (1,0,0); its vertices are VerticesToEdges 1 and 2.
        let ed0 = &ds.e_data_array()[0];
        assert_eq!(ed0.v_sta(), 1);
        assert_eq!(ed0.v_end(), 2);
        assert!(!ed0.rg1_line());
        assert!(!ed0.rg_n_line());
        assert!(!ed0.used());
        assert!(ed0.selected());
        // the face2 edges keep the shifted vertex indices (5..8 -> 6..9?
        // no: 8 vertices total — face1 edges 1-4 use vertices 1-4, face2
        // edges use 5-8).
        let ed4 = &ds.e_data_array()[4];
        assert_eq!(ed4.v_sta(), 5);
        assert_eq!(ed4.v_end(), 6);

        // The FaceData Wires structure (cxx L202/L219/L236): one wire of
        // four non-degenerated edges per face; the wEdge indices are the
        // EM positions 1..4 / 5..8, all FORWARD, Out/Int/Dbl/Iso false
        // (no IntL/IsoL records, plane faces never report a closed edge).
        let fdata = ds.f_data_array();
        assert_eq!(fdata.len(), 2);
        for (fi, fd) in fdata.iter().enumerate() {
            let wb = fd.wires();
            assert_eq!(unsafe { (*wb).nb_wires() }, 1);
            let w = unsafe { (*wb).wire(1) };
            assert_eq!(w.nb_edges(), 4);
            for ie in 1..=4 {
                assert_eq!(w.edge(ie), (fi * 4 + ie) as i32);
                assert_eq!(w.orientation(ie), Orientation::Forward);
                assert!(!w.out_line(ie));
                assert!(!w.internal(ie));
                assert!(!w.double(ie));
                assert!(!w.iso_line(ie));
            }
            assert!(fd.selected());
            assert!(!fd.closed());
        }

        // The MST map after Fill: the two pre-bound faces (cache) plus the
        // two ProcessFace NF bindings.
        assert_eq!(mst.len(), 4);
    }

    /// OCCT anchor: Data::Update after Load — the BigSize over the 16
    /// directions of the unit box (the OCCT formulas, HLRBRep_Data.cxx
    /// L619-655): the 14 trig directions of the x/y spans plus the two z
    /// slots, precad = 0.0005 -> mySurD -> myDeca.
    #[test]
    fn shape_to_hlr_load_update_big_size() {
        let (mut brep, shell, face1, face2) = square_shell_fixture();
        let mut mst = bound_mst(&brep, [&face1, &face2]);
        let proj = top_view_projector();

        let mut ol = OutLiner::with_original_shape(&shell);
        let mut ds = load(&mut brep, &mut ol, &proj, &mut mst, 0);

        // the 'static BRep view for the Update kernel context.
        let brep_st: &'static BRep = Box::leak(Box::new(brep.clone()));
        ds.update(brep_st, &proj);

        // the unit box of the edge cloud: x[0,1], y[0,1], z[0,1] —
        // d[2k] = d[2k+1] = sin(a) + cos(a) over a = k*PI/14, d[14/15] = 1.
        let mut span = [0.0f64; 16];
        let mut big = 1.0f64; // the z span
        for k in 0..7usize {
            let a = (k as f64) * std::f64::consts::PI / 14.0;
            let (s, c) = (a.sin(), a.cos());
            span[2 * k] = s + c;
            span[2 * k + 1] = s + c;
            big = big.max(span[2 * k]).max(span[2 * k + 1]);
        }
        span[14] = 1.0;
        span[15] = 1.0;
        // OCCT hand-calc: max over sin+cos lands at k = 3, 4 (~1.405321).
        assert!((big - 1.405_321).abs() < 1e-5, "big={}", big);
        assert!(
            (ds.big_size() - big).abs() <= 1e-9 * big,
            "bigSize {} vs {}",
            ds.big_size(),
            big
        );
        // mySurD[i] = 0x00007fff / (d[i] + precad), precad = bigSize*0.0005.
        let p1 = big * 0.0005;
        for (i, si) in span.iter().enumerate() {
            let expected = 32767.0 / (si + p1);
            assert!(
                (ds.sur_d()[i] - expected).abs() <= 1e-6 * expected.abs().max(1.0),
                "surD[{}]",
                i
            );
        }
    }

    /// OCCT anchor: the second Load with the same MST map — the pre-bound
    /// face tools are reused by the Fill (the DSFiller MST cache) and the
    /// rebuilt Data carries the same counts.
    #[test]
    fn shape_to_hlr_load_mst_cache_second_load() {
        let (mut brep, shell, face1, face2) = square_shell_fixture();
        let mut mst = bound_mst(&brep, [&face1, &face2]);
        let proj = top_view_projector();

        let mut ol = OutLiner::with_original_shape(&shell);
        let ds1 = load(&mut brep, &mut ol, &proj, &mut mst, 0);
        assert_eq!(mst.len(), 4);

        // the second Load: a fresh OutLiner of the same shape — the Fill
        // hits the MST bindings of the two faces (no re-Bind of the
        // pre-bound tools) and the new NF copies are appended.
        let mut ol2 = OutLiner::with_original_shape(&shell);
        let ds2 = load(&mut brep, &mut ol2, &proj, &mut mst, 0);

        // the pre-bound keys are still the first two entries (cache hit,
        // insertion order preserved).
        assert_eq!(mst[0].0.ptr_id(), face1.ptr_id());
        assert_eq!(mst[1].0.ptr_id(), face2.ptr_id());
        // only the two new NF bindings were appended.
        assert_eq!(mst.len(), 6);

        // the same shape re-loads with the same structure: the edge
        // TShapes are shared with the first OutLinedShape (ProcessFace
        // reuses the original edges), while the face copies are fresh
        // EmptyCopies of the second Fill.
        assert_eq!(ds2.nb_vertices(), ds1.nb_vertices());
        assert_eq!(ds2.nb_edges(), ds1.nb_edges());
        assert_eq!(ds2.nb_faces(), ds1.nb_faces());
        assert_eq!(ds2.edge_map().len(), ds1.edge_map().len());
        for (a, b) in ds1.edge_map().iter().zip(ds2.edge_map().iter()) {
            assert_eq!(a.ptr_id(), b.ptr_id());
        }
        // the fresh face copies of the second OutLinedShape.
        let result2 = ol2.out_lined_shape().clone();
        let faces2 = raw_subshapes(&brep, &raw_subshapes(&brep, &result2)[0]);
        assert_eq!(ds2.face_map().len(), 2);
        for (k, f) in faces2.iter().enumerate() {
            assert_eq!(ds2.face_map()[k].ptr_id(), f.ptr_id());
        }
    }
}
