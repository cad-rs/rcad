//! TopExp / TopoDS / BRep_Tool re-hosts local to
//! ShapeUpgrade_UnifySameDomain (architecture support, TKBRep/TKTopAlgo).
//!
//! The shhealing convention keeps these small walk helpers next to their
//! consuming package (see shape_build/brep_tool.rs for the shared ones and
//! feat/loc_ope_glued_shape.rs for the same re-host pattern).  Only the
//! members UnifySameDomain actually consumes are carried.

use rcad_kernel::geom::{transform_curve, Curve3, CurveEval, Surface3};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepTool, Orientation, ShapeType, TShape};

use super::shape_key;
use super::DataMapOfShapeMapOfShape;
use super::IndexedDataMapOfShapeListOfShape;
use super::MapOfShape;
use crate::shhealing::shape_build::brep_tool::iter_subshapes;
pub(crate) use crate::shhealing::shape_build::brep_tool::topexp_explorer;

/// OCCT TopExp_Explorer(S, ToFind, ToAvoid) (TopExp_Explorer.cxx L107-136):
/// the depth-first walk that skips the ToAvoid sub-trees entirely — the
/// avoided occurrence is neither visited nor descended into.
pub fn topexp_explorer_avoid(
    brep: &mut BRep,
    shape: &Shape,
    to_find: ShapeType,
    to_avoid: ShapeType,
) -> Vec<Shape> {
    let mut out = Vec::new();
    explore_avoid(
        brep,
        shape,
        to_find,
        to_avoid,
        shape.orientation,
        shape.location,
        &mut out,
    );
    out
}

fn explore_avoid(
    brep: &mut BRep,
    current: &Shape,
    to_find: ShapeType,
    to_avoid: ShapeType,
    cum_ori: Orientation,
    cum_loc: u32,
    out: &mut Vec<Shape>,
) {
    for c in crate::shhealing::shape_build::brep_tool::raw_subshapes(brep, current) {
        if c.shape_type() == to_avoid {
            continue;
        }
        let ori = cum_ori.compose(c.orientation);
        let loc_val = brep.get_location(cum_loc) * brep.get_location(c.location);
        let loc = brep.add_location(loc_val);
        let child = Shape {
            orientation: ori,
            location: loc,
            ..c
        };
        if child.shape_type() == to_find {
            out.push(child);
        } else {
            explore_avoid(brep, &child, to_find, to_avoid, ori, loc, out);
        }
    }
}

/// OCCT TopExp::MapShapes(S, T, M) (TopExp.cxx L35-46) — the indexed-map
/// form over the explorer.
pub fn map_shapes(brep: &mut BRep, s: &Shape, t: ShapeType, m: &mut MapOfShape) {
    for e in topexp_explorer(brep, s, t) {
        super::map_add(m, &e);
    }
}

/// OCCT TopExp::MapShapes(S, T, M) keyed variant producing the ordered list
/// (the aFaceMap / aSeqEdges walks).
pub fn map_shapes_list(brep: &mut BRep, s: &Shape, t: ShapeType) -> Vec<Shape> {
    topexp_explorer(brep, s, t)
}

/// OCCT TopExp::MapShapes(S, M, cumOri, cumLoc) (TopExp.cxx L49-63): the
/// all-shapes form — S itself is added, then every sub-shape recursively.
pub fn map_shapes_all(brep: &mut BRep, s: &Shape, m: &mut MapOfShape) {
    super::map_add(m, s);
    for it in iter_subshapes(brep, s, true, true) {
        map_shapes_all(brep, &it, m);
    }
}

/// OCCT TopExp::MapShapesAndAncestors (TopExp.cxx L80-118).
pub fn map_shapes_and_ancestors(
    brep: &mut BRep,
    the_s: &Shape,
    ts: ShapeType,
    ta: ShapeType,
    m: &mut IndexedDataMapOfShapeListOfShape,
) {
    // OCCT L87-105: visit ancestors.
    for a_anc in topexp_explorer(brep, the_s, ta) {
        for a_exs in topexp_explorer(brep, &a_anc, ts) {
            let key = shape_key(&a_exs);
            let entry = match m.get_mut(&key) {
                Some(e) => e,
                None => {
                    m.insert(key, (a_exs.clone(), Vec::new()));
                    m.get_mut(&key).unwrap()
                }
            };
            entry.1.push(a_anc.clone());
        }
    }
    // OCCT L107-117: visit shapes not under ancestors.
    for a_ex in topexp_explorer_avoid(brep, the_s, ts, ta) {
        let key = shape_key(&a_ex);
        m.entry(key).or_insert((a_ex, Vec::new()));
    }
}

/// OCCT TopExp::MapShapesAndUniqueAncestors (TopExp.cxx L124-176): the
/// ancestor list skips duplicates (IsSame), keeping the first occurrence.
pub fn map_shapes_and_unique_ancestors(
    brep: &mut BRep,
    the_s: &Shape,
    ts: ShapeType,
    ta: ShapeType,
    m: &mut IndexedDataMapOfShapeListOfShape,
) {
    // OCCT L133-165: visit ancestors.
    for a_anc in topexp_explorer(brep, the_s, ta) {
        for a_exs in topexp_explorer(brep, &a_anc, ts) {
            let key = shape_key(&a_exs);
            let entry = match m.get_mut(&key) {
                Some(e) => e,
                None => {
                    m.insert(key, (a_exs.clone(), Vec::new()));
                    m.get_mut(&key).unwrap()
                }
            };
            // OCCT L154-162: check if anc already exists in the list.
            if !entry.1.iter().any(|a| occt_is_same_shape(a, &a_anc)) {
                entry.1.push(a_anc.clone());
            }
        }
    }
    // OCCT L167-175: visit shapes not under ancestors.
    for a_ex in topexp_explorer_avoid(brep, the_s, ts, ta) {
        let key = shape_key(&a_ex);
        m.entry(key).or_insert((a_ex, Vec::new()));
    }
}

/// OCCT TopoDS_Shape::IsSame — same TShape with the same Location value
/// (rcad: the location-table index, see the module map note).
pub fn occt_is_same_shape(a: &Shape, b: &Shape) -> bool {
    a.ptr_id() == b.ptr_id() && a.location == b.location
}

/// OCCT TopExp::FirstVertex(E, CumOri) (TopExp.cxx L182-195).
pub fn first_vertex(brep: &mut BRep, e: &Shape, cum_ori: bool) -> Shape {
    for it in iter_subshapes(brep, e, cum_ori, true) {
        if it.orientation == Orientation::Forward {
            return it;
        }
    }
    Shape::null()
}

/// OCCT TopExp::LastVertex(E, CumOri) (TopExp.cxx L198-211).
pub fn last_vertex(brep: &mut BRep, e: &Shape, cum_ori: bool) -> Shape {
    for it in iter_subshapes(brep, e, cum_ori, true) {
        if it.orientation == Orientation::Reversed {
            return it;
        }
    }
    Shape::null()
}

/// OCCT TopExp::Vertices(E, Vfirst, Vlast, CumOri) (TopExp.cxx L214-250).
pub fn vertices(brep: &mut BRep, e: &Shape, cum_ori: bool) -> (Shape, Shape) {
    let mut vfirst = Shape::null();
    let mut vlast = Shape::null();
    let mut is_first_defined = false;
    let mut is_last_defined = false;
    for it in iter_subshapes(brep, e, cum_ori, true) {
        if it.orientation == Orientation::Forward {
            vfirst = it;
            is_first_defined = true;
        } else if it.orientation == Orientation::Reversed {
            vlast = it;
            is_last_defined = true;
        }
    }
    if !is_first_defined {
        vfirst = Shape::null();
    }
    if !is_last_defined {
        vlast = Shape::null();
    }
    (vfirst, vlast)
}

/// OCCT TopExp::CommonVertex(E1, E2, V) (TopExp.cxx L316-334).
pub fn common_vertex(brep: &mut BRep, e1: &Shape, e2: &Shape) -> Option<Shape> {
    let (first_vertex1, last_vertex1) = vertices(brep, e1, false);
    let (first_vertex2, last_vertex2) = vertices(brep, e2, false);
    if occt_is_same_shape(&first_vertex1, &first_vertex2)
        || occt_is_same_shape(&first_vertex1, &last_vertex2)
    {
        return Some(first_vertex1);
    }
    if occt_is_same_shape(&last_vertex1, &first_vertex2)
        || occt_is_same_shape(&last_vertex1, &last_vertex2)
    {
        return Some(last_vertex1);
    }
    None
}

// ---------------------------------------------------------------------------
// BRep_Tool re-hosts (BRep_Tool.cxx) — the location-aware forms the trait
// surface does not carry.
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Surface(F, L) (BRep_Tool.cxx L228-243): the LOCAL surface
/// (the location is returned separately).
pub fn brep_tool_surface_loc(brep: &BRep, f: &Shape) -> (Option<Surface3>, u32) {
    match &*brep.tshapes[f.index] {
        TShape::Face(fd) => (fd.surface.clone(), fd.surface_location),
        _ => (None, 0),
    }
}

/// OCCT BRep_Tool::Curve(E, L, First, Last) (BRep_Tool.cxx L121-160): the
/// LOCAL curve with its location and range.
pub fn brep_tool_curve_loc(brep: &BRep, e: &Shape) -> (Option<Curve3>, u32, [f64; 2]) {
    match &*brep.tshapes[e.index] {
        TShape::Edge(ed) => (ed.curve.clone(), e.location, ed.range),
        _ => (None, 0, [0.0, 0.0]),
    }
}

/// OCCT BRep_Tool::Curve(E, First, Last) (BRep_Tool.cxx L162-195): the
/// location-applied curve; the parameters are rescaled by the transform
/// scale factor (Geom_Curve::TransformedParameter).
pub fn brep_tool_curve(brep: &BRep, e: &Shape) -> Option<(Curve3, f64, f64)> {
    let (c_opt, loc, range) = brep_tool_curve_loc(brep, e);
    let c = c_opt?;
    let (mut first, mut last) = (range[0], range[1]);
    let t = brep.get_location(loc);
    let c = if t != glam::DAffine3::IDENTITY {
        first = c.transformed_parameter(first);
        last = c.transformed_parameter(last);
        transform_curve(&c, &t)
    } else {
        c
    };
    Some((c, first, last))
}

/// OCCT BRep_Tool::Range(E, First, Last) (BRep_Tool.cxx L845-858) — the raw
/// 3D range of the edge.
pub fn brep_tool_range(brep: &BRep, e: &Shape) -> [f64; 2] {
    match &*brep.tshapes[e.index] {
        TShape::Edge(ed) => ed.range,
        _ => [0.0, 0.0],
    }
}

/// OCCT BRep_Tool::IsClosed(E, F) routed through the kernel trait
/// (BRep_Tool.cxx L795-841).
pub fn brep_tool_is_closed_edge_face(brep: &BRep, e: &Shape, f: &Shape) -> bool {
    brep.is_edge_closed_on_face(e, f)
}

// ---------------------------------------------------------------------------
// NCollection container helpers over the local aliases.
// ---------------------------------------------------------------------------

/// OCCT IndexedDataMap::FindFromKey with the "already inserted" guard —
/// returns the ancestor list clone (empty when the key is absent, the
/// OCCT not-found path where callers guard explicitly).
pub fn find_from_key<'a>(
    m: &'a IndexedDataMapOfShapeListOfShape,
    s: &Shape,
) -> Option<&'a (Shape, Vec<Shape>)> {
    m.get(&shape_key(s))
}

/// OCCT DataMap::Seek — an optional borrow of the value.
pub fn data_map_seek<'a>(
    m: &'a DataMapOfShapeMapOfShape,
    s: &Shape,
) -> Option<&'a (Shape, MapOfShape)> {
    m.get(&shape_key(s))
}

/// OCCT BRep_Tool::Degenerated(E) (BRep_Tool.cxx L1162-1171) — the edge
/// degenerated flag read.
pub fn is_edge_degenerated(brep: &BRep, e: &Shape) -> bool {
    matches!(&*brep.tshapes[e.index], TShape::Edge(ed) if ed.degenerated)
}
