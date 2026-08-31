//! Kernel-side helpers for the ShapeFix healing stack (TKBRep equivalents).
//!
//! These are the 1:1 translations of the small OCCT machinery the
//! ShapeBuild/ShapeFix packages call directly:
//! - `TopoDS_Shape::IsNull` -> [`shape_is_null`]
//! - `TopoDS_Shape::IsSame` / `IsPartner` -> [`occt_is_same`] / [`occt_is_partner`]
//! - `TopoDS_Iterator` -> [`iter_subshapes`]
//! - `TopExp_Explorer` -> [`topexp_explorer`]
//! - `TopoDS_Builder::Add` (the `Add` inherited by `BRep_Builder`) -> [`builder_add`]
//! - `BRep_Tool::IsClosed` -> [`brep_tool_is_closed`]
//!
//! Identity mapping note: OCCT compares `myLocation` by value while rcad
//! stores locations as indices into `BRep::locations`; the helpers therefore
//! compare the resolved `DAffine3` values, not the indices.

use rcad_kernel::topo::topods::{
    BRep, Orientation, Shape, ShapeType, TShape, tshape_flags,
};
use std::sync::Arc;

/// OCCT TopoDS_Shape::IsNull. rcad's null placeholder is built by
/// `Shape::null()` with `index == usize::MAX`.
pub fn shape_is_null(s: &Shape) -> bool {
    s.index == usize::MAX
}

/// OCCT TopoDS_Shape::IsPartner (TopoDS_Shape.hxx L263): same TShape only;
/// Locations and Orientations may differ. (rcad's `Shape::is_same` carries
/// this semantic despite the name.)
pub fn occt_is_partner(brep: &BRep, s1: &Shape, s2: &Shape) -> bool {
    let _ = brep;
    s1.ptr_id() == s2.ptr_id()
}

/// OCCT TopoDS_Shape::IsSame (TopoDS_Shape.hxx L268-271): same TShape with the
/// same Location value; Orientations may differ. (rcad's `Shape::is_partner`
 /// carries the ptr+location semantics despite the name; the location is
/// compared by resolved transform value, as OCCT compares the TopLoc value.)
pub fn occt_is_same(brep: &BRep, s1: &Shape, s2: &Shape) -> bool {
    s1.ptr_id() == s2.ptr_id() && location_values_equal(brep, s1.location, s2.location)
}

/// OCCT TopLoc_Location::operator== — the stored transforms are equal.
pub fn location_values_equal(brep: &BRep, l1: u32, l2: u32) -> bool {
    if l1 == l2 {
        return true;
    }
    brep.get_location(l1) == brep.get_location(l2)
}

/// OCCT TopoDS_Shape::Move(thePosition) (TopoDS_Shape.hxx L193-201):
/// `myLocation = thePosition * myLocation`. Registers the composed transform
/// in the BRep location table (identity stays index 0).
pub fn shape_move(brep: &mut BRep, s: &mut Shape, the_position: u32) {
    let new_loc = brep.get_location(the_position) * brep.get_location(s.location);
    s.location = brep.add_location(new_loc);
}

/// OCCT TopoDS_Iterator over the direct sub-shapes of `s`
/// (TopoDS_Iterator.cxx L28-70): children with the cumulative orientation
/// (`TopAbs::Compose`) and the cumulative location (`Move`) applied.
/// `cum_ori == false` keeps the stored child orientations ( FORWARD passed as
/// the parent orientation), `cum_loc == false` keeps the stored locations.
pub fn iter_subshapes(brep: &mut BRep, s: &Shape, cum_ori: bool, cum_loc: bool) -> Vec<Shape> {
    let my_orientation = if cum_ori { s.orientation } else { Orientation::Forward };
    let mut children = raw_subshapes(brep, s);
    for c in children.iter_mut() {
        c.orientation = my_orientation.compose(c.orientation);
        if cum_loc && s.location != 0 {
            // OCCT updateCurrentShape: myShape.Move(myLocation, false).
            let loc = brep.get_location(s.location);
            let composed = loc * brep.get_location(c.location);
            c.location = brep.add_location(composed);
        }
    }
    children
}

/// The raw children of a shape in the rcad TShape model (the equivalent of
/// OCCT `TShape->myShapes`). Architecture note: OCCT stores one generic list;
/// rcad stores typed fields. The typed-field order mirrors `add_tface` /
/// `add_twire` / `add_tshell` / `add_tsolid` construction order.
pub fn raw_subshapes(brep: &BRep, s: &Shape) -> Vec<Shape> {
    match s.data.as_ref() {
        TShape::Vertex(_) => Vec::new(),
        TShape::Edge(ed) => {
            let mut out = Vec::with_capacity(2);
            if !shape_is_null(&ed.first) {
                out.push(ed.first.clone());
            }
            if !shape_is_null(&ed.last) {
                out.push(ed.last.clone());
            }
            out
        }
        TShape::Wire(wd) => wd.edges.clone(),
        TShape::Face(fd) => {
            let mut out = Vec::with_capacity(2 + fd.inner_wires.len());
            if !shape_is_null(&fd.outer_wire) {
                out.push(fd.outer_wire.clone());
            }
            out.extend(fd.inner_wires.iter().cloned());
            out.extend(fd.internal_vertices.iter().cloned());
            out
        }
        TShape::Shell(sd) => sd.faces.clone(),
        TShape::Solid(sd) => {
            let mut out = Vec::with_capacity(sd.shells.len());
            out.extend(sd.shells.iter().cloned());
            out.extend(sd.internal_vertices.iter().cloned());
            out.extend(sd.internal_edges.iter().cloned());
            out
        }
        TShape::CompSolid(cd) => cd.clone(),
        TShape::Compound(cd) => cd.clone(),
    }
}

/// OCCT TopExp_Explorer(shape, to_find): a depth-first walk that collects every
/// occurrence of `to_find` and does NOT descend inside a found sub-shape. The
/// start shape itself is not tested. Orientations and locations are composed
/// cumulatively along the path (TopExp_Explorer.cxx L107-136, L172-279).
pub fn topexp_explorer(brep: &mut BRep, shape: &Shape, to_find: ShapeType) -> Vec<Shape> {
    let mut out = Vec::new();
    explore(brep, shape, to_find, shape.orientation, shape.location, &mut out);
    out
}

fn explore(
    brep: &mut BRep,
    current: &Shape,
    to_find: ShapeType,
    cum_ori: Orientation,
    cum_loc: u32,
    out: &mut Vec<Shape>,
) {
    for c in raw_subshapes(brep, current) {
        let ori = cum_ori.compose(c.orientation);
        let loc_val = brep.get_location(cum_loc) * brep.get_location(c.location);
        let loc = brep.add_location(loc_val);
        let child = Shape { orientation: ori, location: loc, ..c };
        if child.shape_type() == to_find {
            out.push(child);
        } else {
            explore(brep, &child, to_find, ori, loc, out);
        }
    }
}

/// OCCT TopoDS_Builder::Add (TopoDS_Builder.cxx L41-101) — the `Add` every
/// builder (including BRep_Builder) inherits: compatibility check, relative
/// orientation/location of the component, in-place append to the container's
/// sub-shape list.
///
/// rcad mapping of the container append: OCCT appends to one generic
/// `myShapes` list; rcad writes the typed field (`Compound` vec, `edges`,
/// `faces`, `shells`, face `outer_wire`/`inner_wires`/`internal_vertices`,
/// edge `first`/`last`). The first wire added to a face is the outer wire,
/// matching the `add_tface` construction order.
///
/// OCCT throws TopoDS_UnCompatibleShapes for an incompatible pair; rcad panics.
/// OCCT throws TopoDS_FrozenShape when the container is not Free; rcad does not
/// enforce the Free flag.
pub fn builder_add(brep: &mut BRep, a_shape: &Shape, a_component: &Shape) {
    // Compatibility table: aTb[componentType] & (1 << shapeType).
    let allowed: [u16; 9] = {
        let mut tb = [0u16; 9];
        let bit = |t: ShapeType| 1u16 << (t as u16);
        tb[ShapeType::Compound as usize] = bit(ShapeType::Compound);
        tb[ShapeType::CompSolid as usize] = bit(ShapeType::Compound);
        tb[ShapeType::Solid as usize] =
            bit(ShapeType::Compound) | bit(ShapeType::CompSolid);
        tb[ShapeType::Shell as usize] = bit(ShapeType::Compound) | bit(ShapeType::Solid);
        tb[ShapeType::Face as usize] = bit(ShapeType::Compound) | bit(ShapeType::Shell);
        tb[ShapeType::Wire as usize] = bit(ShapeType::Compound) | bit(ShapeType::Face);
        tb[ShapeType::Edge as usize] = bit(ShapeType::Compound)
            | bit(ShapeType::Solid)
            | bit(ShapeType::Wire);
        tb[ShapeType::Vertex as usize] = bit(ShapeType::Compound)
            | bit(ShapeType::Solid)
            | bit(ShapeType::Face)
            | bit(ShapeType::Edge);
        tb
    };
    let i_c = a_component.shape_type() as usize;
    let i_s = a_shape.shape_type() as usize;
    if (allowed[i_c] & (1 << i_s)) == 0 {
        panic!("TopoDS_UnCompatibleShapes: TopoDS_Builder::Add");
    }

    let mut a_child = a_component.clone();

    // Compute the relative Orientation.
    if a_shape.orientation == Orientation::Reversed {
        a_child.orientation = match a_child.orientation {
            Orientation::Forward => Orientation::Reversed,
            Orientation::Reversed => Orientation::Forward,
            other => other,
        };
    }

    // And the relative Location: aChild.Move(aLoc.Inverted()).
    if a_shape.location != 0 {
        let inv = brep.get_location(a_shape.location).inverse();
        let composed = inv * brep.get_location(a_child.location);
        a_child.location = brep.add_location(composed);
    }

    append_inplace(brep, a_shape, a_child);
}

/// The OCCT `TShape->myShapes.Append(aChild)` step, mapped onto rcad's typed
/// container fields. Mutates the TShape in place (shared-handle semantics, as
/// OCCT mutates the TShape through the handle) so every existing Shape handle
/// observes the added child.
fn append_inplace(brep: &mut BRep, container: &Shape, child: Shape) {
    // SAFETY: the healing chain builds containers sequentially right after
    // creation (BRep_Builder pattern); no aliased &TShape is alive. This is
    // the same in-place pattern as BRep::edge_mut_inplace.
    let ptr = Arc::as_ptr(&brep.tshapes[container.index]) as *mut TShape;
    let ts = unsafe { &mut *ptr };
    match ts {
        TShape::Compound(cd) => cd.push(child),
        TShape::CompSolid(cd) => cd.push(child),
        TShape::Solid(sd) => match child.shape_type() {
            ShapeType::Vertex => sd.internal_vertices.push(child),
            ShapeType::Edge => sd.internal_edges.push(child),
            _ => sd.shells.push(child),
        },
        TShape::Shell(sd) => sd.faces.push(child),
        TShape::Face(fd) => match child.shape_type() {
            ShapeType::Wire => {
                if shape_is_null(&fd.outer_wire) {
                    fd.outer_wire = child.clone();
                } else {
                    fd.inner_wires.push(child.clone());
                }
                fd.my_shapes.push(child);
            }
            ShapeType::Vertex => fd.internal_vertices.push(child),
            _ => {}
        },
        TShape::Wire(wd) => {
            wd.edges.push(child.clone());
            wd.my_shapes.push(child);
        }
        TShape::Edge(ed) => match child.orientation {
            Orientation::Forward => {
                if shape_is_null(&ed.first) {
                    ed.first = child;
                }
            }
            Orientation::Reversed => {
                if shape_is_null(&ed.last) {
                    ed.last = child;
                }
            }
            _ => {}
        },
        TShape::Vertex(_) => {}
    }
}

/// OCCT TopoDS_Shape::Closed(theFlag) — set/clear the CLOSED flag on the
/// TShape. Mutates in place (shared-handle semantics) so every existing Shape
/// handle observes the change, unlike BRep::set_flag (Arc::make_mut, which
/// splits identity when the Arc is shared).
pub fn set_flag_inplace(brep: &BRep, s: &Shape, flag: u16, on: bool) {
    // SAFETY: single-threaded healing chain; the flag write mirrors OCCT
    // mutating the TShape through the handle (same pattern as
    // BRep::edge_mut_inplace).
    let ptr = Arc::as_ptr(&brep.tshapes[s.index]) as *mut TShape;
    let ts = unsafe { &mut *ptr };
    let flags = match ts {
        TShape::Vertex(vd) => &mut vd.flags,
        TShape::Edge(ed) => &mut ed.flags,
        TShape::Wire(wd) => &mut wd.flags,
        TShape::Face(fd) => &mut fd.flags,
        TShape::Shell(sd) => &mut sd.flags,
        TShape::Solid(sd) => &mut sd.flags,
        _ => return,
    };
    if on {
        *flags |= flag;
    } else {
        *flags &= !flag;
    }
}

/// OCCT BRep_Tool::IsClosed(theShape) (BRep_Tool.cxx L1707-1761): SHELL counts
/// edge occurrences (degenerated and INTERNAL/EXTERNAL edges skipped), WIRE
/// counts vertex occurrences (INTERNAL/EXTERNAL vertices skipped), EDGE checks
/// its two vertices, everything else returns the CLOSED flag.
pub fn brep_tool_is_closed(brep: &BRep, the_shape: &Shape) -> bool {
    match the_shape.shape_type() {
        ShapeType::Shell => {
            let mut a_map: std::collections::HashSet<(u64, u32)> =
                std::collections::HashSet::new();
            let mut has_bound = false;
            // TopExp_Explorer over the FORWARD-oriented shell visits every
            // edge occurrence (shared edges once per face).
            let mut stack: Vec<Shape> = raw_subshapes(brep, the_shape);
            while let Some(s) = stack.pop() {
                match s.shape_type() {
                    ShapeType::Face => {
                        stack.extend(raw_subshapes(brep, &s));
                    }
                    ShapeType::Wire => {
                        for e in raw_subshapes(brep, &s) {
                            if e.shape_type() != ShapeType::Edge {
                                continue;
                            }
                            let degenerated =
                                matches!(e.data.as_ref(), TShape::Edge(ed) if ed.degenerated);
                            if degenerated
                                || e.orientation == Orientation::Internal
                                || e.orientation == Orientation::External
                            {
                                continue;
                            }
                            has_bound = true;
                            // TopTools_ShapeMapHasher: IsSame identity
                            // (TShape + Location). OCCT compares the location
                            // value; rcad keys by the location index, which is
                            // stable within one BRep's location table.
                            let key = (e.ptr_id(), e.location);
                            if !a_map.insert(key) {
                                a_map.remove(&key);
                            }
                        }
                    }
                    _ => {}
                }
            }
            has_bound && a_map.is_empty()
        }
        ShapeType::Wire => {
            let mut a_map: std::collections::HashSet<(u64, u32)> =
                std::collections::HashSet::new();
            let mut has_bound = false;
            for e in raw_subshapes(brep, the_shape) {
                if e.shape_type() != ShapeType::Edge {
                    continue;
                }
                if let TShape::Edge(ed) = e.data.as_ref() {
                    for sv in [&ed.first, &ed.last] {
                        let vori = match e.orientation {
                            Orientation::Reversed => match sv.orientation {
                                Orientation::Forward => Orientation::Reversed,
                                Orientation::Reversed => Orientation::Forward,
                                other => other,
                            },
                            _ => sv.orientation,
                        };
                        if vori == Orientation::Internal || vori == Orientation::External {
                            continue;
                        }
                        has_bound = true;
                        let key = (sv.ptr_id(), sv.location);
                        if !a_map.insert(key) {
                            a_map.remove(&key);
                        }
                    }
                }
            }
            has_bound && a_map.is_empty()
        }
        ShapeType::Edge => {
            if let TShape::Edge(ed) = the_shape.data.as_ref() {
                !shape_is_null(&ed.first)
                    && occt_is_same(brep, &ed.first, &ed.last)
            } else {
                false
            }
        }
        _ => brep.has_flag(the_shape.clone(), tshape_flags::CLOSED),
    }
}
