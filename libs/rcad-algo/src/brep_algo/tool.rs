//! Shared BRep_Tool / TopExp / BRep_Builder re-hosts for the BRepAlgo
//! package (OCCT TKBRep vehicle used by AsDes / Image / Loop /
//! FaceRestrictor / NormalProjection).
//!
//! Re-host style follows the feat package precedent
//! (feat/loc_ope_wires_on_shape_b.rs):
//! - Shape identity key (TopTools_ShapeMapHasher): TShape pointer +
//!   Location, orientation ignored.
//! - BRep_Builder mutations: in-place Arc::make_mut edits on the Shape
//!   payload (architecture difference: rcad TShape data travels inside the
//!   Shape's Arc, no builder object).
//! - Locations are identity in this pipeline (feat arch. diff. #1): the
//!   OCCT "apply loc.Transformation()" branches are no-ops.

use rcad_kernel::geom::{Curve2d, Curve3, Surface3};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{tshape_flags, Orientation, ShapeType, TShape};
use std::collections::HashMap;
use std::sync::Arc;

pub(crate) type ShapeKey = (u64, u32);

/// OCCT TopTools_ShapeMapHasher: TShape + Location, orientation ignored.
pub(crate) fn shape_key(s: &Shape) -> ShapeKey {
    (s.ptr_id(), s.location)
}

/// OCCT TopoDS_Shape::Oriented(theOr) — a copy carrying the orientation.
pub(crate) fn oriented(s: &Shape, the_or: Orientation) -> Shape {
    let mut c = s.clone();
    c.orientation = the_or;
    c
}

/// OCCT TopAbs::Reverse (TopAbs.hxx L80-91): FORWARD<->REVERSED,
/// INTERNAL/EXTERNAL unchanged.
pub(crate) fn top_abs_reverse(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        o => o,
    }
}

/// OCCT TopoDS_Shape::Reversed() — Oriented(TopAbs::Reverse(Orientation)).
pub(crate) fn reversed(s: &Shape) -> Shape {
    oriented(s, top_abs_reverse(s.orientation))
}

/// OCCT TopoDS_Shape::IsSame over two possibly-null shapes (TopoDS_Shape.hxx
/// L268-271): two null shapes share the same null TShape handle and are the
/// same; rcad Shape::null() allocates distinct Arcs so the null-null case is
/// handled explicitly.
pub(crate) fn is_same_opt(a: &Option<Shape>, b: &Option<Shape>) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => x.is_same(y),
        (None, None) => true,
        _ => false,
    }
}

/// OCCT theShape.IsSame(theOther) with theOther possibly null.
pub(crate) fn shape_same_opt(s: &Shape, o: &Option<Shape>) -> bool {
    match o {
        Some(other) => s.is_same(other),
        None => s.is_null(),
    }
}

/// OCCT theShape.IsEqual(theOther) with theOther possibly null
/// (TopoDS_Shape.hxx L276-280: TShape + Location + Orientation).
pub(crate) fn shape_is_equal_opt(s: &Shape, o: &Option<Shape>) -> bool {
    match o {
        Some(other) => s.is_equal(other),
        None => s.is_null(),
    }
}

// ---------------------------------------------------------------------------
// BRep_Tool re-hosts.
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Tolerance — the three shape-kind overloads, each of
/// which floors the stored tolerance at Precision::Confusion():
///   - BRep_Tool::Tolerance(const TopoDS_Vertex& V)  BRep_Tool.cxx L1314-1333
///   - BRep_Tool::Tolerance(const TopoDS_Edge& E)    BRep_Tool.cxx L881-895
///   - BRep_Tool::Tolerance(const TopoDS_Face& F)    BRep_Tool.cxx L137-149
/// All three bodies are identical: p = TE->Tolerance(); pMin =
/// Precision::Confusion(); if (p > pMin) return p; else return pMin.
/// OCCT has no generic BRep_Tool::Tolerance(const TopoDS_Shape&)
/// dispatcher (BRep_Tool.hxx declares only the three overloads at
/// L87/L256/L345), so the match over the shape kind is the rcad stand-in
/// for the C++ overload resolution.  The OCCT vertex body throws
/// Standard_NullObject when the TVertex is null; the rcad TShape::Vertex
/// payload is non-nullable, so that branch has no counterpart.
pub(crate) fn brep_tool_tolerance(s: &Shape) -> f64 {
    // OCCT `constexpr double pMin = Precision::Confusion();`
    const P_MIN: f64 = rcad_kernel::precision::CONFUSION;
    match s.data.as_ref() {
        // OCCT BRep_Tool.cxx L1314-1333.
        TShape::Vertex(vd) => {
            let p = vd.tolerance;
            if p > P_MIN {
                p
            } else {
                P_MIN
            }
        }
        // OCCT BRep_Tool.cxx L881-895.
        TShape::Edge(ed) => {
            let p = ed.tolerance;
            if p > P_MIN {
                p
            } else {
                P_MIN
            }
        }
        // OCCT BRep_Tool.cxx L137-149.
        TShape::Face(fd) => {
            let p = fd.tolerance;
            if p > P_MIN {
                p
            } else {
                P_MIN
            }
        }
        // No OCCT overload exists for wire/shell/solid/compound.
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::Pnt(V).
pub(crate) fn brep_tool_pnt(vtx: &Shape) -> Option<glam::DVec3> {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => Some(vd.point),
        _ => None,
    }
}

/// OCCT BRep_Tool::Curve(E, f, l) — the 3D curve and range (identity
/// location; feat loc_ope_find_edges.rs arch. diff. #1).
pub(crate) fn brep_tool_curve(edg: &Shape) -> Option<(Curve3, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => {
            let c = ed.curve.as_ref()?;
            Some((c.clone(), ed.range[0], ed.range[1]))
        }
        _ => None,
    }
}

/// OCCT BRep_Tool::Range(E, f, l) — the 3D range.
pub(crate) fn brep_tool_range(edg: &Shape) -> (f64, f64) {
    match edg.data.as_ref() {
        TShape::Edge(ed) => (ed.range[0], ed.range[1]),
        _ => (0.0, 0.0),
    }
}

/// OCCT BRep_Tool::Surface(F) (feat loc_ope_gluer.rs arch. diff. #1).
pub(crate) fn brep_tool_surface(face: &Shape) -> Option<Surface3> {
    match face.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT BRep_Tool::CurveOnSurface(E, F, f, l) — the pcurve of the edge on
/// the face with its range.
pub(crate) fn brep_tool_curve_on_surface(edg: &Shape, face: &Shape) -> Option<(Curve2d, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed
            .pcurves
            .get(&shape_key(face))
            .map(|(c, f, l)| (c.clone(), *f, *l)),
        _ => None,
    }
}

/// OCCT BRep_Tool::CurveOnSurface(E, C2, S, L, f, l, Index) with Index == 1 —
/// the FIRST pcurve whatever its surface (BRep_Tool.cxx L354-401 reduced to
/// the rcad pcurve map: the first entry).
pub(crate) fn brep_tool_first_curve_on_surface(edg: &Shape) -> Option<(Curve2d, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed.pcurves.values().next().map(|(c, f, l)| (c.clone(), *f, *l)),
        _ => None,
    }
}

/// OCCT BRep_Tool::Parameter(V, E) — the stored vertex parameter on the edge.
pub(crate) fn brep_tool_parameter(vtx: &Shape, edg: &Shape) -> f64 {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed
            .vertex_params
            .get(&vtx.ptr_id())
            .copied()
            .unwrap_or(0.0),
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::UVPoints(E, F, PFirst, PLast) — the pcurve values at the
/// edge range bounds (feat loc_ope_generator_b.rs precedent).
pub(crate) fn brep_tool_uv_points(edg: &Shape, face: &Shape) -> (glam::DVec2, glam::DVec2) {
    use rcad_kernel::geom::Curve2dEval;
    let (f, l) = brep_tool_range(edg);
    if let Some((c, _, _)) = brep_tool_curve_on_surface(edg, face) {
        return (c.point_at(f), c.point_at(l));
    }
    (glam::DVec2::ZERO, glam::DVec2::ZERO)
}

/// OCCT BRep_Tool::IsClosed(theShape) for TopAbs_EDGE
/// (BRep_Tool.cxx L1757-1763): the two extremity vertices are the same.
pub(crate) fn brep_tool_is_closed_edge(edg: &Shape) -> bool {
    let (v1, v2) = top_exp_vertices_raw(edg);
    match (&v1, &v2) {
        (Some(a), Some(b)) => a.is_same(b),
        _ => false,
    }
}

/// OCCT BRep_Tool::IsClosed(E, S, L) (BRep_Tool.cxx L814-841): a seam edge on
/// the closed surface — a CurveOnClosedSurface representation matching the
/// surface (rcad keys the pcurves by face; the caller passes the face whose
/// surface is S — architecture difference, the Geom_Surface handle identity
/// maps to the face shape key).
pub(crate) fn brep_tool_is_closed_on_surface(edg: &Shape, face: &Shape) -> bool {
    match edg.data.as_ref() {
        TShape::Edge(ed) => {
            let fkey = shape_key(face);
            ed.representations.iter().any(|cr| match cr {
                rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface { face, .. } => {
                    *face == fkey
                }
                _ => false,
            })
        }
        _ => false,
    }
}

// ---------------------------------------------------------------------------
// TopExp re-hosts.
// ---------------------------------------------------------------------------

/// OCCT TopExp::Vertices(E, Vfirst, Vlast) with the default CumOri=false
/// (TopExp.cxx L214-253): the RAW stored child orientation selects
/// first/last and the vertices are returned with that raw orientation.
pub(crate) fn top_exp_vertices_raw(edg: &Shape) -> (Option<Shape>, Option<Shape>) {
    let ed = match edg.data.as_ref() {
        TShape::Edge(ed) => ed,
        _ => return (None, None),
    };
    let mut vfirst: Option<Shape> = None;
    let mut vlast: Option<Shape> = None;
    // The rcad edge children (TopoDS_Iterator over TShape::Shapes order:
    // first, last).
    for child in [&ed.first, &ed.last] {
        if child.is_null() {
            continue;
        }
        match child.orientation {
            Orientation::Forward => {
                if vfirst.is_none() {
                    vfirst = Some(child.clone());
                }
            }
            Orientation::Reversed => {
                if vlast.is_none() {
                    vlast = Some(child.clone());
                }
            }
            _ => {}
        }
    }
    (vfirst, vlast)
}

/// OCCT TopExp::Vertices(const TopoDS_Wire& W, Vfirst, Vlast)
/// (TopExp.cxx L255-318): walk the wire edges (cumulated orientation),
/// add/remove the extremity vertices in a map; empty map = closed (both ends
/// on the last V2), extent 2 = open (the FORWARD-oriented key is Vfirst, the
/// REVERSED one is Vlast).  The NCollection_Map scan order maps to the rcad
/// insertion-ordered map (architecture difference: bucket order vs insertion
/// order; both hold exactly the same two vertices for the open case).
pub(crate) fn top_exp_vertices_wire(w: &Shape) -> (Option<Shape>, Option<Shape>) {
    let mut vfirst: Option<Shape> = None;
    let mut vlast: Option<Shape> = None;
    let mut v2: Option<Shape> = None; // OCCT: the V2 local, last written per edge
    let mut vmap: indexmap::IndexMap<ShapeKey, Shape> = indexmap::IndexMap::new();

    for e in sub_shapes(w) {
        let (mut v1, mut v2e) = if e.orientation == Orientation::Reversed {
            let (a, b) = top_exp_vertices_raw(&e);
            (b, a)
        } else {
            top_exp_vertices_raw(&e)
        };
        // add or remove in the vertex map.
        if let Some(vv) = &mut v1 {
            vv.orientation = Orientation::Forward;
        }
        if let Some(vv) = &mut v2e {
            vv.orientation = Orientation::Reversed;
        }
        v2 = v2e.clone();
        if let Some(vv) = v1 {
            // OCCT !Add -> Remove (the entry was already present).
            if vmap.insert(shape_key(&vv), vv.clone()).is_some() {
                vmap.swap_remove(&shape_key(&vv));
            }
        }
        if let Some(vv) = v2e {
            if vmap.insert(shape_key(&vv), vv.clone()).is_some() {
                vmap.swap_remove(&shape_key(&vv));
            }
        }
    }
    if vmap.is_empty() {
        // closed: both ends on the last processed V2.
        if let Some(v2v) = &v2 {
            vfirst = Some(oriented(v2v, Orientation::Forward));
            vlast = Some(oriented(v2v, Orientation::Reversed));
        }
    } else if vmap.len() == 2 {
        // open: the FORWARD-oriented key is Vfirst, the REVERSED one is Vlast.
        for (_, k) in vmap.iter() {
            if k.orientation == Orientation::Forward {
                vfirst = Some(k.clone());
                break;
            }
        }
        for (_, k) in vmap.iter() {
            if k.orientation == Orientation::Reversed {
                vlast = Some(k.clone());
                break;
            }
        }
    }
    (vfirst, vlast)
}

/// OCCT TopAbs::ShapeTypeValue order (TopAbs_ShapeEnum) — used by the
/// explorer for the "ty > toFind" pruning (TopExp_Explorer.cxx L90-95).
fn topabs_value(t: ShapeType) -> i32 {
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

/// OCCT TopoDS_Iterator (cumOri=true): the stored sub-shapes with
/// TopAbs::Compose(parent.Orientation(), child.Orientation()) applied
/// (TopoDS_Iterator.cxx L72-80; brep_feat_builder.rs sub_shapes model).
pub(crate) fn sub_shapes(sh: &Shape) -> Vec<Shape> {
    let composed = |c: &Shape| {
        let mut c = c.clone();
        c.orientation = sh.orientation.compose(c.orientation);
        c
    };
    match sh.data.as_ref() {
        TShape::Vertex(_) => Vec::new(),
        TShape::Edge(ed) => vec![composed(&ed.first), composed(&ed.last)],
        TShape::Wire(wd) => wd.edges.iter().map(composed).collect(),
        TShape::Face(fd) => {
            let mut out =
                Vec::with_capacity(1 + fd.inner_wires.len() + fd.internal_vertices.len());
            // OCCT BRep_TFace owns no children: TopoDS_Iterator walks
            // TopoDS_TShape::myShapes, which holds exactly what
            // BRep_Builder::Add appended (the wires and the internal
            // vertices).  rcad splits those into typed slots, so only the
            // slots that actually hold a child are yielded.  A face built
            // without a boundary wire (OCCT BRepLib_MakeFace(F, S) /
            // MakeFace(F, S, Tol) — a natural restriction) carries the NULL
            // placeholder in outer_wire, which is a Vertex, NOT a wire; and
            // because pool-free builder shapes also carry index == usize::MAX
            // (see Shape::is_null), the test has to be on the child type.
            if matches!(fd.outer_wire.data.as_ref(), TShape::Wire(_)) {
                out.push(composed(&fd.outer_wire));
            }
            out.extend(fd.inner_wires.iter().map(composed));
            out.extend(fd.internal_vertices.iter().map(composed));
            out
        }
        TShape::Shell(sd) => sd.faces.iter().map(composed).collect(),
        TShape::Solid(sd) => {
            let mut out = Vec::new();
            out.extend(sd.shells.iter().map(composed));
            out.extend(sd.internal_vertices.iter().map(composed));
            out.extend(sd.internal_edges.iter().map(composed));
            out
        }
        TShape::CompSolid(cs) => cs.iter().map(composed).collect(),
        TShape::Compound(cd) => cd.iter().map(composed).collect(),
    }
}

/// OCCT TopExp_Explorer (TopExp_Explorer.cxx L68-171): yields the sub-shapes
/// of `s` of type `to_find`, descending into coarser container types, never
/// entering types of `to_avoid` (ShapeType::Shape = "avoid nothing").
pub(crate) fn explorer(s: &Shape, to_find: ShapeType, to_avoid: ShapeType) -> Vec<Shape> {
    let mut out: Vec<Shape> = Vec::new();
    if to_find == ShapeType::Shape {
        return out;
    }
    if topabs_value(s.shape_type()) > topabs_value(to_find) {
        return out;
    }
    explore_rec(s, to_find, to_avoid, &mut out);
    out
}

fn explore_rec(sh: &Shape, to_find: ShapeType, to_avoid: ShapeType, out: &mut Vec<Shape>) {
    let ty = sh.shape_type();
    if ty == to_find {
        out.push(sh.clone());
        return;
    }
    if topabs_value(ty) > topabs_value(to_find) {
        return;
    }
    if to_avoid != ShapeType::Shape && ty == to_avoid {
        return;
    }
    for c in sub_shapes(sh) {
        explore_rec(&c, to_find, to_avoid, out);
    }
}

// ---------------------------------------------------------------------------
// TopoDS_Shape::EmptyCopied (no pool: the TShape travels in the Arc).
// ---------------------------------------------------------------------------

/// OCCT TopoDS_Shape::EmptyCopied (TopoDS_Shape.hxx L168-172): a new TShape
/// via TShape::EmptyCopy() carrying the ORIGINAL Location and Orientation.
/// The per-type field copy mirrors the OCCT BRep_Txx::EmptyCopy overrides
/// (rcad-kernel topods.rs BRep::empty_copy; architecture difference: rcad
/// copies from the Shape's own Arc, no BRep pool slot involved).
pub(crate) fn empty_copied(r: &Shape) -> Shape {
    let new_data = match r.data.as_ref() {
        TShape::Vertex(vd) => Arc::new(TShape::Vertex(rcad_kernel::topods::TVertexData {
            my_shapes: Vec::new(),
            flags: tshape_flags::DEFAULT,
            point: vd.point,
            tolerance: vd.tolerance,
            points: Vec::new(),
        })),
        TShape::Edge(ed) => Arc::new(TShape::Edge(rcad_kernel::topods::TEdgeData {
            my_shapes: Vec::new(),
            flags: tshape_flags::DEFAULT,
            curve: ed.curve.clone(),
            first: Shape::null(),
            last: Shape::null(),
            range: ed.range,
            degenerated: ed.degenerated,
            pcurves: ed.pcurves.clone(),
            representations: ed.representations.clone(),
            vertex_params: HashMap::new(),
            tolerance: ed.tolerance,
            same_parameter: ed.same_parameter,
            same_range: ed.same_range,
        })),
        TShape::Wire(_) => Arc::new(TShape::Wire(rcad_kernel::topods::TWireData {
            my_shapes: Vec::new(),
            flags: tshape_flags::DEFAULT,
            edges: Vec::new(),
        })),
        TShape::Face(fd) => Arc::new(TShape::Face(rcad_kernel::topods::TFaceData {
            my_shapes: Vec::new(),
            flags: tshape_flags::DEFAULT,
            surface: fd.surface.clone(),
            surface_location: fd.surface_location,
            outer_wire: Shape::null(),
            inner_wires: Vec::new(),
            sample_point: None,
            uv_domain: None,
            internal_vertices: Vec::new(),
            tolerance: fd.tolerance,
            natural_restriction: false,
        })),
        TShape::Shell(_) => Arc::new(TShape::Shell(rcad_kernel::topods::TShellData {
            my_shapes: Vec::new(),
            flags: tshape_flags::DEFAULT,
            faces: Vec::new(),
        })),
        TShape::Solid(_) => Arc::new(TShape::Solid(rcad_kernel::topods::TSolidData {
            my_shapes: Vec::new(),
            flags: tshape_flags::DEFAULT,
            shells: Vec::new(),
            internal_vertices: Vec::new(),
            internal_edges: Vec::new(),
        })),
        TShape::CompSolid(_) => Arc::new(TShape::CompSolid(Vec::new())),
        TShape::Compound(_) => Arc::new(TShape::Compound(Vec::new())),
    };
    Shape {
        data: new_data,
        index: usize::MAX,
        location: r.location,
        orientation: r.orientation,
    }
}

// ---------------------------------------------------------------------------
// BRep_Builder re-hosts (in-place Arc::make_mut edits).
// ---------------------------------------------------------------------------

/// OCCT BRep_Builder::MakeWire(W) — an empty wire TShape.
pub(crate) fn builder_make_wire() -> Shape {
    Shape {
        data: Arc::new(TShape::Wire(rcad_kernel::topods::TWireData {
            my_shapes: Vec::new(),
            flags: tshape_flags::DEFAULT,
            edges: Vec::new(),
        })),
        index: usize::MAX,
        location: 0,
        orientation: Orientation::Forward,
    }
}

/// OCCT BRep_Builder::MakeCompound(C) — an empty compound TShape.
pub(crate) fn builder_make_compound() -> Shape {
    Shape {
        data: Arc::new(TShape::Compound(Vec::new())),
        index: usize::MAX,
        location: 0,
        orientation: Orientation::Forward,
    }
}

/// OCCT BRep_Builder::MakeVertex(V) — an empty vertex TShape (no point yet).
pub(crate) fn builder_make_vertex() -> Shape {
    Shape {
        data: Arc::new(TShape::Vertex(rcad_kernel::topods::TVertexData {
            my_shapes: Vec::new(),
            flags: tshape_flags::DEFAULT,
            point: glam::DVec3::ZERO,
            tolerance: 0.0,
            points: Vec::new(),
        })),
        index: usize::MAX,
        location: 0,
        orientation: Orientation::Forward,
    }
}

/// OCCT BRep_Builder::MakeEdge(E) — an empty edge TShape.
pub(crate) fn builder_make_edge() -> Shape {
    Shape {
        data: Arc::new(TShape::Edge(rcad_kernel::topods::TEdgeData {
            my_shapes: Vec::new(),
            flags: tshape_flags::DEFAULT,
            curve: None,
            first: Shape::null(),
            last: Shape::null(),
            range: [0.0, 0.0],
            degenerated: false,
            pcurves: indexmap::IndexMap::new(),
            representations: Vec::new(),
            vertex_params: HashMap::new(),
            tolerance: 0.0,
            same_parameter: false,
            same_range: false,
        })),
        index: usize::MAX,
        location: 0,
        orientation: Orientation::Forward,
    }
}

/// OCCT BRep_Builder::Add(W, E) — append the edge to the wire.
pub(crate) fn builder_add_wire_edge(the_w: &mut Shape, the_e: &Shape) {
    if let TShape::Wire(wd) = Arc::make_mut(&mut the_w.data) {
        wd.edges.push(the_e.clone());
    }
}

/// OCCT BRep_Builder::Add(C, S) — append the shape to the compound.
pub(crate) fn builder_add_compound_shape(the_c: &mut Shape, the_s: &Shape) {
    if let TShape::Compound(cd) = Arc::make_mut(&mut the_c.data) {
        cd.push(the_s.clone());
    }
}

/// OCCT BRep_Builder::Add(F, W) — the first wire is the outer wire,
/// the following ones are inner wires (BRep_Builder.cxx Add Face branch).
pub(crate) fn builder_add_face_wire(the_f: &mut Shape, the_w: &Shape) {
    if let TShape::Face(fd) = Arc::make_mut(&mut the_f.data) {
        if fd.outer_wire.is_null() {
            fd.outer_wire = the_w.clone();
        } else {
            fd.inner_wires.push(the_w.clone());
        }
    }
}

/// OCCT BRep_Builder::Add(E, V) — attach the vertex by orientation
/// (FORWARD -> first, REVERSED -> last).
pub(crate) fn builder_add_edge_vertex(the_e: &mut Shape, the_v: &Shape) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        match the_v.orientation {
            Orientation::Reversed => ed.last = the_v.clone(),
            _ => ed.first = the_v.clone(),
        }
    }
}

/// OCCT BRep_Builder::Remove(E, V) — detach the vertex sub-shape from the
/// edge (BRep_Builder.cxx Remove: the matching extremity is nullified and
/// its stored parameter dropped).
pub(crate) fn builder_remove_edge_vertex(the_e: &mut Shape, the_v: &Shape) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        if the_v.is_same(&ed.first) {
            ed.first = Shape::null();
        }
        if the_v.is_same(&ed.last) {
            ed.last = Shape::null();
        }
        ed.vertex_params.remove(&the_v.ptr_id());
    }
}

/// OCCT BRep_Builder::UpdateVertex(V, Tol) — BRep_TVertex::UpdateTolerance
/// keeps the max.
pub(crate) fn builder_update_vertex_tol(the_v: &mut Shape, the_tol: f64) {
    if let TShape::Vertex(vd) = Arc::make_mut(&mut the_v.data) {
        vd.tolerance = vd.tolerance.max(the_tol);
    }
}

/// OCCT BRep_Builder::UpdateVertex(V, P, Tol) — move the vertex point and
/// update the tolerance (BRep_Builder.cxx L1069-1100).
pub(crate) fn builder_update_vertex_point_tol(
    the_v: &mut Shape,
    the_p: glam::DVec3,
    the_tol: f64,
) {
    if let TShape::Vertex(vd) = Arc::make_mut(&mut the_v.data) {
        vd.point = the_p;
        vd.tolerance = vd.tolerance.max(the_tol);
    }
}

/// OCCT BRep_Builder::UpdateEdge(E, C2d, F, Tol) — bind the pcurve on the
/// face (the pcurve range is the edge 3D range).
pub(crate) fn builder_update_edge_pcurve(
    the_e: &mut Shape,
    the_c2d: &Curve2d,
    the_f: &Shape,
    the_tol: f64,
) {
    let key = shape_key(the_f);
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        let (f0, l0) = (ed.range[0], ed.range[1]);
        ed.pcurves.insert(key, (the_c2d.clone(), f0, l0));
        ed.tolerance = ed.tolerance.max(the_tol);
    }
}

/// OCCT BRep_Builder::Range(E, First, Last).
pub(crate) fn builder_range_edge(the_e: &mut Shape, the_first: f64, the_last: f64) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        ed.range = [the_first, the_last];
    }
}

/// OCCT BRep_Builder::Degenerated(E, flag).
pub(crate) fn builder_set_degenerated(the_e: &mut Shape, flag: bool) {
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        ed.degenerated = flag;
    }
}

/// OCCT TopoDS_Shape::Closed(theIsClosed) — the CLOSED flag write on the
/// shape's own TShape (TopoDS_Shape.hxx L216-223).
pub(crate) fn builder_set_closed(the_s: &mut Shape, flag: bool) {
    let ts = Arc::make_mut(&mut the_s.data);
    let flags = match ts {
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

/// OCCT BRep_Builder::Free(S, theIsFree) — the FREE flag write
/// (BRep_Builder.cxx Free branch).
pub(crate) fn builder_set_free(the_s: &mut Shape, flag: bool) {
    let ts = Arc::make_mut(&mut the_s.data);
    let flags = match ts {
        TShape::Vertex(v) => &mut v.flags,
        TShape::Edge(e) => &mut e.flags,
        TShape::Wire(w) => &mut w.flags,
        TShape::Face(f) => &mut f.flags,
        TShape::Shell(sh) => &mut sh.flags,
        TShape::Solid(so) => &mut so.flags,
        TShape::CompSolid(_) | TShape::Compound(_) => return,
    };
    if flag {
        *flags |= tshape_flags::FREE;
    } else {
        *flags &= !tshape_flags::FREE;
    }
}

/// OCCT TopoDS_Shape::Closed() — the CLOSED flag read on the shape's own
/// TShape (TopoDS_Shape.hxx L216-223 getter).
pub(crate) fn shape_is_closed(the_s: &Shape) -> bool {
    let flags = match the_s.data.as_ref() {
        TShape::Vertex(v) => v.flags,
        TShape::Edge(e) => e.flags,
        TShape::Wire(w) => w.flags,
        TShape::Face(f) => f.flags,
        TShape::Shell(sh) => sh.flags,
        TShape::Solid(so) => so.flags,
        TShape::CompSolid(_) | TShape::Compound(_) => 0,
    };
    flags & tshape_flags::CLOSED != 0
}

/// OCCT BRep_Builder::Range(E, F, First, Last) — the pcurve range of the
/// edge on the face (BRep_Builder.cxx Range CurveOnSurface branch).
pub(crate) fn builder_range_edge_on_face(
    the_e: &mut Shape,
    the_f: &Shape,
    the_first: f64,
    the_last: f64,
) {
    let key = shape_key(the_f);
    if let TShape::Edge(ed) = Arc::make_mut(&mut the_e.data) {
        if let Some(entry) = ed.pcurves.get_mut(&key) {
            entry.1 = the_first;
            entry.2 = the_last;
        }
    }
}
