//! OCCT ShapeBuild_ReShape + its base BRepTools_ReShape
//! (TKShHealing ShapeBuild / TKBRep BRepTools packages).
//!
//! 1:1 translation of `ShapeBuild_ReShape.cxx` (L34-339) and
//! `BRepTools_ReShape.cxx` (L140-695). Rust has no inheritance, so the base
//! class state and methods are flattened into one struct; sections are marked
//! with the OCCT class they translate. The base class's own
//! `BRepTools_ReShape::Apply/applyImpl` (BRepTools_ReShape.cxx L439-602) is
//! fully overridden by the derived `ShapeBuild_ReShape::applyImpl` at OCCT
//! runtime (virtual dispatch), so only the derived version is translated.
//!
//! Every method takes the owning `BRep` pool: OCCT reaches TShape data
//! through the `TopoDS_Shape` handle against the global shape space, rcad
//! through the pool (architecture difference).
//!
//! Identity mapping notes:
//! - OCCT map keys use `TopTools_ShapeMapHasher` = `IsSame` (TShape +
//!   Location, orientation ignored). rcad keys are `(ptr_id, location)`;
//!   the location index stands for the location value within one BRep.
//! - OCCT `IsSame` (TShape+Location) and `IsPartner` (TShape only) map to
//!   [`brep_tool::occt_is_same`] / [`brep_tool::occt_is_partner`].
//! - OCCT `TopAbs_ShapeEnum` ranks COMPOUND=0 .. SHAPE=8; rcad's `ShapeType`
//!   is the reverse order, so OCCT rank comparisons go through
//!   [`occt_type_rank`].

use crate::shhealing::shape_build::brep_tool::{
    brep_tool_is_closed, builder_add, iter_subshapes, occt_is_partner, occt_is_same,
    set_flag_inplace, shape_is_null, topexp_explorer,
};
use crate::shhealing::shape_build::edge::ShapeBuildEdge;
use crate::shhealing::shape_extend::{ShapeExtendStatus, decode_status, encode_status};
use rcad_kernel::topo::topods::{BRep, Orientation, Shape, ShapeType};
use std::collections::{HashMap, HashSet};

/// OCCT TopAbs_ShapeEnum rank (COMPOUND=0 .. SHAPE=8); rcad's enum is the
/// reverse order so the rank is `8 - value`.
fn occt_type_rank(t: ShapeType) -> i32 {
    8 - (t as i32)
}

/// Map key with TopTools_ShapeMapHasher semantics (IsSame: TShape + Location;
/// orientation ignored).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct ShapeKey(u64, u32);

fn shape_key(s: &Shape) -> ShapeKey {
    ShapeKey(s.ptr_id(), s.location)
}

/// OCCT BRepTools_ReShape::TReplacementKind (BRepTools_ReShape.hxx L164-170).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TReplacementKind {
    Remove = 1,
    Modify = 2,
    MergeMain = 4,
    MergeOrdinary = 8,
}

/// OCCT BRepTools_ReShape::TReplacement (BRepTools_ReShape.hxx L201-237).
#[derive(Debug, Clone)]
struct TReplacement {
    my_result: Shape,
    my_kind: TReplacementKind,
}

impl TReplacement {
    /// OCCT TReplacement::Result: Merge_Ordinary replacements carry no result
    /// shape (the product is only reachable from the main part).
    fn result(&self) -> Shape {
        if self.my_kind != TReplacementKind::MergeOrdinary {
            self.my_result.clone()
        } else {
            Shape::null()
        }
    }
}

/// OCCT ShapeBuild_ReShape : BRepTools_ReShape.
#[derive(Debug, Clone)]
pub struct ShapeBuildReShape {
    /// OCCT `myShapeToReplacement`: maps each shape to its replacement. If a
    /// shape is not bound then the shape is replaced by itself.
    my_shape_to_replacement: HashMap<ShapeKey, TReplacement>,
    /// OCCT `myNewShapes`.
    my_new_shapes: HashSet<ShapeKey>,
    /// OCCT `myStatus` (ShapeExtend bit field; -1 = Apply not yet run).
    my_status: i32,
    /// OCCT `myConsiderLocation`.
    my_consider_location: bool,
}

impl Default for ShapeBuildReShape {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeBuildReShape {
    // -------------------------------------------------------------------
    // BRepTools_ReShape base-class section (BRepTools_ReShape.cxx)
    // -------------------------------------------------------------------

    /// OCCT BRepTools_ReShape() (L140-144): myStatus = -1, considerLocation
    /// off.
    pub fn new() -> Self {
        ShapeBuildReShape {
            my_shape_to_replacement: HashMap::new(),
            my_new_shapes: HashSet::new(),
            my_status: -1,
            my_consider_location: false,
        }
    }

    /// OCCT BRepTools_ReShape::Clear (L148-152): clears all substitution
    /// requests.
    pub fn clear(&mut self) {
        self.my_shape_to_replacement.clear();
        self.my_new_shapes.clear();
    }

    /// OCCT BRepTools_ReShape::Remove (L156-160): sets a request to remove a
    /// shape whatever the orientation.
    pub fn remove(&mut self, brep: &mut BRep, shape: &Shape) {
        let nulshape = Shape::null();
        self.replace_impl(brep, shape, &nulshape, TReplacementKind::Remove);
    }

    /// OCCT BRepTools_ReShape::Replace (BRepTools_ReShape.hxx L65-68): sets a
    /// request to replace a shape by a new one.
    pub fn replace(&mut self, brep: &mut BRep, shape: &Shape, newshape: &Shape) {
        self.replace_impl(brep, shape, newshape, TReplacementKind::Modify);
    }

    /// OCCT BRepTools_ReShape::Merge (BRepTools_ReShape.hxx L74-91): merges
    /// the parts to the single product; the first part is replaced by the
    /// product, the other parts are removed.
    pub fn merge(&mut self, brep: &mut BRep, parts: &[Shape], the_product: &Shape) {
        if let Some(first) = parts.first() {
            self.replace_impl(brep, first, the_product, TReplacementKind::MergeMain);
        }
        for part in parts.iter().skip(1) {
            self.replace_impl(brep, part, the_product, TReplacementKind::MergeOrdinary);
        }
    }

    /// OCCT BRepTools_ReShape::replace (L164-209): records the replacement
    /// after the reorientation rules:
    /// - REVERSED source: both shapes are reversed;
    /// - INTERNAL/EXTERNAL source: the new shape is oriented forward
    ///   (reversed) when its orientation is equal (not equal) to the source's,
    ///   and the source is oriented forward.
    fn replace_impl(
        &mut self,
        brep: &mut BRep,
        ashape: &Shape,
        anewshape: &Shape,
        the_kind: TReplacementKind,
    ) {
        let mut shape = ashape.clone();
        let mut newshape = anewshape.clone();
        if shape_is_null(&shape) || shape.is_equal(&newshape) {
            return;
        }

        if shape.orientation == Orientation::Reversed {
            reverse_orientation(&mut shape);
            reverse_orientation(&mut newshape);
        }
        // Protect against INTERNAL or EXTERNAL shape.
        else if shape.orientation == Orientation::Internal
            || shape.orientation == Orientation::External
        {
            newshape.orientation = if newshape.orientation == shape.orientation {
                Orientation::Forward
            } else {
                Orientation::Reversed
            };
            shape.orientation = Orientation::Forward;
        }

        if self.my_consider_location {
            // sln 29.11.01 Bug22: change location of 'newshape' in accordance
            // with location of 'shape':
            // newshape.Location(newshape.Location().Multiplied(
            //     shape.Location().Inverted())).
            let inv = brep.get_location(shape.location).inverse();
            let composed = inv * brep.get_location(newshape.location);
            newshape.location = brep.add_location(composed);
            shape.location = 0;
        }

        // Cycle handling: cycles in the replacement map are accepted at
        // insertion time; the DFS in-flight guard inside Apply breaks them at
        // traversal time (BRepTools_ReShape.cxx L203-206).
        self.my_shape_to_replacement.insert(
            shape_key(&shape),
            TReplacement {
                my_result: newshape.clone(),
                my_kind: the_kind,
            },
        );
        self.my_new_shapes.insert(shape_key(&newshape));
    }

    /// OCCT BRepTools_ReShape::IsRecorded (L213-226): tells if a shape is
    /// recorded for Replace/Remove.
    pub fn is_recorded(&self, ashape: &Shape) -> bool {
        if shape_is_null(ashape) {
            return false;
        }
        self.my_shape_to_replacement.contains_key(&shape_key(ashape))
    }

    /// OCCT BRepTools_ReShape::Value (L230-279): the new value for an
    /// individual shape. If not recorded, returns the original shape itself;
    /// if to be removed, returns a null shape; else the replacing item.
    pub fn value(&self, brep: &mut BRep, ashape: &Shape) -> Shape {
        if shape_is_null(ashape) {
            return Shape::null();
        }
        let mut shape = ashape.clone();
        if self.my_consider_location {
            shape.location = 0;
        }

        let mut res;
        let from_map;
        match self.my_shape_to_replacement.get(&shape_key(&shape)) {
            None => {
                res = shape.clone();
                from_map = false;
            }
            Some(replacement) => {
                res = replacement.result();
                if shape.orientation == Orientation::Reversed {
                    reverse_orientation(&mut res);
                }
                from_map = true;
            }
        }
        // For INTERNAL/EXTERNAL, since they are not fully supported, keep
        // orientation.
        if shape.orientation == Orientation::Internal || shape.orientation == Orientation::External
        {
            res.orientation = shape.orientation;
        }

        if self.my_consider_location {
            // sln 29.11.01 Bug22: recalculate location of the resulting shape
            // in accordance with whether the result is from the map or not.
            if from_map {
                let composed = brep.get_location(ashape.location) * brep.get_location(res.location);
                res.location = brep.add_location(composed);
            } else {
                res.location = ashape.location;
            }
        }

        res
    }

    /// OCCT BRepTools_ReShape::ValueLeaf (L283-318): follows the replacement
    /// chain to its leaf without descending into sub-shapes. Returns the final
    /// replacement, the original shape if not recorded, or a null shape when
    /// the chain terminates in a Remove.
    pub fn value_leaf(&self, brep: &mut BRep, the_shape: &Shape) -> Shape {
        if shape_is_null(the_shape) {
            return Shape::null();
        }
        // Visited TShapes (handle identity = ptr only), per OCCT's
        // NCollection_Map of TShape handles.
        let mut a_visited: HashSet<u64> = HashSet::new();
        let mut a_current = the_shape.clone();
        a_visited.insert(a_current.ptr_id());

        loop {
            let a_next = self.value(brep, &a_current);
            if shape_is_null(&a_next) {
                return a_next;
            }
            if occt_is_same(brep, &a_next, &a_current) {
                return a_next;
            }
            if !a_visited.insert(a_next.ptr_id()) {
                // Cycle in replacement data - return current best to avoid
                // looping.
                return a_next;
            }
            a_current = a_next;
        }
    }

    /// OCCT BRepTools_ReShape::Status (L322-387): complete substitution status
    /// for a shape. Returns `(res, newsh)`: `res` is 0 when not recorded
    /// (newsh = original shape), -1 when to be removed (newsh null), > 0 when
    /// to be replaced. With `last`, the final replacement is searched
    /// recursively via Apply.
    pub fn status(&mut self, brep: &mut BRep, ashape: &Shape, last: bool) -> (i32, Shape) {
        let mut res = 0;
        if shape_is_null(ashape) {
            return (0, Shape::null());
        }

        let mut shape = ashape.clone();
        let a_loc_sh = shape.location;
        if self.my_consider_location {
            shape.location = 0;
        }

        let mut newsh;
        match self.my_shape_to_replacement.get(&shape_key(&shape)) {
            None => {
                newsh = shape.clone();
                res = 0;
            }
            Some(replacement) => {
                newsh = replacement.my_result.clone();
                res = 1;
            }
        }
        if res > 0 {
            if shape_is_null(&newsh) {
                res = -1;
            } else if newsh.is_equal(&shape) {
                res = 0;
            } else if last
                && ((self.my_consider_location && !occt_is_partner(brep, &newsh, &shape))
                    || (!self.my_consider_location && !occt_is_same(brep, &newsh, &shape)))
            {
                // sln 29.11.01 Bug24: iterate to the final replacement.
                newsh = self.apply(brep, &shape, ShapeType::Shape);
                if shape_is_null(&newsh) {
                    res = -1;
                }
                if newsh.is_equal(&shape) {
                    res = 0;
                }
            }
        }
        if self.my_consider_location && !shape_is_null(&newsh) {
            let a_res_loc = if res > 0 && newsh.location != 0 {
                brep.get_location(a_loc_sh) * brep.get_location(newsh.location)
            } else {
                brep.get_location(a_loc_sh)
            };
            newsh.location = brep.add_location(a_res_loc);
        }
        (res, newsh)
    }

    /// OCCT BRepTools_ReShape::CopyVertex (L613-616): copy with a new
    /// tolerance, keeping the current position.
    pub fn copy_vertex(&mut self, brep: &mut BRep, the_v: &Shape, the_tol: f64) -> Shape {
        // OCCT BRep_Tool::Pnt(theV): the point in world coordinates.
        let the_new_pos = {
            let vd = brep.vertex(the_v.clone());
            brep.get_location(the_v.location).transform_point3(vd.point)
        };
        self.copy_vertex_at(brep, the_v, the_new_pos, the_tol)
    }

    /// OCCT BRepTools_ReShape::CopyVertex (L620-638): copy with a new position
    /// and tolerance; returns the modified copy when the original is not
    /// recorded, the modified original otherwise.
    pub fn copy_vertex_at(
        &mut self,
        brep: &mut BRep,
        the_v: &Shape,
        the_new_pos: glam::DVec3,
        the_tol: f64,
    ) -> Shape {
        let is_recorded = self.is_recorded(the_v);
        let a_vertex_copy = if is_recorded {
            let (_, newsh) = self.status(brep, the_v, false);
            newsh
        } else {
            brep.empty_copied(the_v)
        };

        // OCCT BRep_Builder::UpdateVertex(V, P, Tol) (BRep_Builder.cxx
        // L1203-1213): Pnt(P.Transformed(V.Location().Inverted())) +
        // UpdateTolerance(Tol) (max with the current tolerance).
        let a_new_tol = if the_tol > 0.0 {
            the_tol
        } else {
            brep.vertex(the_v.clone()).tolerance
        };
        let local = brep
            .get_location(a_vertex_copy.location)
            .inverse()
            .transform_point3(the_new_pos);
        let vd = brep.vertex_mut(a_vertex_copy.clone());
        vd.point = local;
        vd.tolerance = vd.tolerance.max(a_new_tol);

        if !is_recorded {
            self.replace(brep, the_v, &a_vertex_copy);
        }

        a_vertex_copy
    }

    /// OCCT BRepTools_ReShape::IsNewShape (L640-643): checks if the shape has
    /// been recorded by the reshaper as a value.
    pub fn is_new_shape(&self, the_shape: &Shape) -> bool {
        self.my_new_shapes.contains(&shape_key(the_shape))
    }

    // OCCT BRepTools_ReShape::History (L647-695) is not translated: it
    // depends on BRepTools_History, which has no rcad equivalent, and no
    // TKShHealing entry point (fixshape chain) calls it.

    // -------------------------------------------------------------------
    // ShapeBuild_ReShape derived-class section (ShapeBuild_ReShape.cxx)
    // -------------------------------------------------------------------

    /// OCCT ShapeBuild_ReShape::Apply(shape, until, buildmode) (L38-187).
    /// `buildmode` says how to rebuild a SOLID/SHELL when one of its
    /// sub-shapes has been changed:
    /// 0: at least one Replace or Remove -> COMPOUND, else as such;
    /// 1: at least one Remove (Replace ignored) -> COMPOUND;
    /// 2: Replace and Remove are both ignored.
    pub fn apply_buildmode(
        &mut self,
        brep: &mut BRep,
        shape: &Shape,
        until: ShapeType,
        buildmode: i32,
    ) -> Shape {
        if shape_is_null(shape) {
            return shape.clone();
        }
        let (stat, newsh) = self.status(brep, shape, false);
        if stat != 0 {
            return newsh;
        }

        let st = shape.shape_type();
        if st == until {
            return newsh; // Critere d arret.
        }

        if st == ShapeType::Compound || st == ShapeType::CompSolid {
            let mut modif = 0;
            let c = brep.add_tcompound(Vec::new());
            for sh in iter_subshapes(brep, shape, true, true) {
                let (stat, newsh) = self.status(brep, &sh, false);
                if stat != 0 {
                    modif = 1;
                }
                if stat >= 0 {
                    builder_add(brep, &c, &newsh);
                }
            }
            if modif == 0 {
                return shape.clone();
            }
            return c;
        }

        if st == ShapeType::Solid {
            let mut modif = 0;
            let c = brep.add_tcompound(Vec::new());
            let s = brep.add_tsolid(Vec::new());
            for sh in iter_subshapes(brep, shape, true, true) {
                let newsh = self.apply_buildmode(brep, &sh, until, buildmode);
                if shape_is_null(&newsh) {
                    modif = -1;
                } else if newsh.shape_type() != ShapeType::Shell {
                    let mut nbsub = 0;
                    for onesh in topexp_explorer(brep, &newsh, ShapeType::Shell) {
                        builder_add(brep, &s, &onesh);
                        nbsub += 1;
                    }
                    if nbsub == 0 {
                        modif = -1;
                    }
                    builder_add(brep, &c, &newsh); // c est tout
                } else {
                    if modif == 0 && !sh.is_equal(&newsh) {
                        modif = 1;
                    }
                    builder_add(brep, &c, &newsh);
                    builder_add(brep, &s, &newsh);
                }
            }

            if (modif < 0 && buildmode < 2) || (modif == 0 && buildmode < 1) {
                return c;
            }
            return s;
        }

        if st == ShapeType::Shell {
            let mut modif = 0;
            let c = brep.add_tcompound(Vec::new());
            let s = brep.add_tshell(Vec::new());
            for sh in iter_subshapes(brep, shape, true, true) {
                let newsh = self.apply_buildmode(brep, &sh, until, buildmode);
                if shape_is_null(&newsh) {
                    modif = -1;
                } else if newsh.shape_type() != ShapeType::Face {
                    let mut nbsub = 0;
                    for onef in topexp_explorer(brep, &newsh, ShapeType::Face) {
                        builder_add(brep, &s, &onef);
                        nbsub += 1;
                    }
                    if nbsub == 0 {
                        modif = -1;
                    }
                    builder_add(brep, &c, &newsh); // c est tout
                } else {
                    if modif == 0 && !sh.is_equal(&newsh) {
                        modif = 1;
                    }
                    builder_add(brep, &c, &newsh);
                    builder_add(brep, &s, &newsh);
                }
            }
            if (modif < 0 && buildmode < 2) || (modif == 0 && buildmode < 1) {
                return c;
            }
            // S.Closed(BRep_Tool::IsClosed(S)).
            let closed = brep_tool_is_closed(brep, &s);
            set_flag_inplace(brep, &s, rcad_kernel::topo::topods::tshape_flags::CLOSED, closed);
            return s;
        }
        println!("BRepTools_ReShape::Apply NOT YET IMPLEMENTED");
        shape.clone()
    }

    /// OCCT ShapeBuild_ReShape::Apply(shape, until) (L191-196) — applies the
    /// substitution requests to a shape, DFS with an in-flight cycle guard.
    pub fn apply(&mut self, brep: &mut BRep, the_shape: &Shape, the_until: ShapeType) -> Shape {
        let mut an_in_flight: HashSet<u64> = HashSet::new();
        self.apply_impl(brep, the_shape, the_until, &mut an_in_flight)
    }

    /// OCCT ShapeBuild_ReShape::applyImpl (L200-323).
    fn apply_impl(
        &mut self,
        brep: &mut BRep,
        the_shape: &Shape,
        the_until: ShapeType,
        the_in_flight: &mut HashSet<u64>,
    ) -> Shape {
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        if shape_is_null(the_shape) {
            return the_shape.clone();
        }

        // Apply direct replacement.
        let mut a_new_shape = self.value(brep, the_shape);

        // If shape removed, return NULL.
        if shape_is_null(&a_new_shape) {
            self.my_status = encode_status(ShapeExtendStatus::Done2);
            return a_new_shape;
        }

        // DFS cycle guard: if theShape is already being processed further up
        // the call stack, its replacement must be a compound that
        // transitively contains it. Return the direct replacement without
        // descending to break the cycle.
        if the_in_flight.contains(&the_shape.ptr_id()) {
            return a_new_shape;
        }

        // If shape was replaced, apply modifications to the result
        // recursively.
        let a_cons_loc = self.my_consider_location;
        if (a_cons_loc && !occt_is_partner(brep, &a_new_shape, the_shape))
            || (!a_cons_loc && !occt_is_same(brep, &a_new_shape, the_shape))
        {
            the_in_flight.insert(the_shape.ptr_id());
            let a_res = self.apply_impl(brep, &a_new_shape, the_until, the_in_flight);
            the_in_flight.remove(&the_shape.ptr_id());
            self.my_status |= encode_status(ShapeExtendStatus::Done1);
            return a_res;
        }

        let a_st = the_shape.shape_type();
        if occt_type_rank(a_st) >= occt_type_rank(the_until) {
            return a_new_shape; // Stop criterion.
        }
        if a_st == ShapeType::Vertex || a_st == ShapeType::Shape {
            return the_shape.clone();
        }

        let mut a_result = brep.empty_copied(the_shape);
        let an_orient = the_shape.orientation;
        a_result.orientation = Orientation::Forward; // Protect against INTERNAL or EXTERNAL shapes.
        let mut a_modif = false;
        let mut a_loc_status = self.my_status;

        // Apply recorded modifications to subshapes.
        the_in_flight.insert(the_shape.ptr_id());
        for a_sh in iter_subshapes(brep, the_shape, false, true) {
            a_new_shape = self.apply_impl(brep, &a_sh, the_until, the_in_flight);
            if !a_new_shape.is_equal(&a_sh) {
                if decode_status(self.my_status, ShapeExtendStatus::Done4) {
                    a_loc_status |= encode_status(ShapeExtendStatus::Done4);
                }
                a_modif = true;
            }
            if shape_is_null(&a_new_shape) {
                a_loc_status |= encode_status(ShapeExtendStatus::Done4);
                continue;
            }
            a_loc_status |= encode_status(ShapeExtendStatus::Done3);
            if a_st == ShapeType::Compound || a_new_shape.shape_type() == a_sh.shape_type() {
                // Fix for SAMTECH bug OCC322 about absence internal vertices
                // after sewing.
                builder_add(brep, &a_result, &a_new_shape);
                continue;
            }
            let mut a_nb_items = 0;
            // TopoDS_Iterator aSubIt(aNewShape): cumOri and cumLoc default on.
            for a_sub_sh in iter_subshapes(brep, &a_new_shape, true, true) {
                a_nb_items += 1;
                if a_sub_sh.shape_type() == a_sh.shape_type() {
                    builder_add(brep, &a_result, &a_sub_sh);
                } else {
                    a_loc_status |= encode_status(ShapeExtendStatus::Fail1);
                }
            }
            if a_nb_items == 0 {
                a_loc_status |= encode_status(ShapeExtendStatus::Fail1);
            }
        }
        the_in_flight.remove(&the_shape.ptr_id());
        if !a_modif {
            return the_shape.clone();
        }

        // Restore range on edge broken by EmptyCopied().
        if a_st == ShapeType::Edge {
            let an_sbe = ShapeBuildEdge;
            an_sbe.copy_ranges(brep, &a_result, the_shape, 0.0, 1.0);
        } else if a_st == ShapeType::Wire || a_st == ShapeType::Shell {
            // aResult.Closed(BRep_Tool::IsClosed(aResult)).
            let closed = brep_tool_is_closed(brep, &a_result);
            set_flag_inplace(
                brep,
                &a_result,
                rcad_kernel::topo::topods::tshape_flags::CLOSED,
                closed,
            );
        }
        a_result.orientation = an_orient;
        self.my_status = a_loc_status;

        let kind = if shape_is_null(&a_result) {
            TReplacementKind::Remove
        } else {
            TReplacementKind::Modify
        };
        self.replace_impl(brep, the_shape, &a_result, kind);

        a_result
    }

    /// OCCT ShapeBuild_ReShape::Status(shape, newsh, last) (L327-332): the
    /// base implementation.
    pub fn status_shape(&mut self, brep: &mut BRep, shape: &Shape, last: bool) -> (i32, Shape) {
        self.status(brep, shape, last)
    }

    /// OCCT ShapeBuild_ReShape::Status(ShapeExtend_Status) (L336-339): queries
    /// the status of the last call to Apply.
    /// OK: no (sub)shapes replaced or removed; DONE1: source shape replaced;
    /// DONE2: source shape removed; DONE3: some subshapes replaced;
    /// DONE4: some subshapes removed; FAIL1: some replacements not done
    /// because of bad type of subshape.
    pub fn status_flag(&self, the_status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status, the_status)
    }

    /// OCCT BRepTools_ReShape::ModeConsiderLocation accessor
    /// (BRepTools_ReShape.hxx L136).
    pub fn mode_consider_location(&mut self) -> &mut bool {
        &mut self.my_consider_location
    }
}

/// OCCT TopoDS_Shape::Reverse — flips Forward <-> Reversed.
fn reverse_orientation(s: &mut Shape) {
    s.orientation = match s.orientation {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        other => other,
    };
}

#[cfg(test)]
mod reshape_tests;
