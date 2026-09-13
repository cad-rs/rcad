//! OCCT BRepTools_Modifier (TKBRep/BRepTools — TKTopAlgo landing zone per
//! docs/module-map.md) — 1:1 translation of BRepTools_Modifier.cxx (L51-852),
//! BRepTools_Modifier.hxx (L40-134) and BRepTools_Modifier.lxx (L21-35).
//!
//! Source: $OCCT_SRC/src/ModelingData/TKBRep/BRepTools/BRepTools_Modifier.cxx
//!
//! The geometric modification is supplied as a `&mut dyn BRepToolsModification`
//! — the sibling brep_tools_modification.rs module (OCCT
//! BRepTools_Modification + the BRepTools_TrsfModification subclass).
//!
//! Architecture differences:
//! 1. `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape,
//!    TopTools_ShapeMapHasher> myMap` -> `HashMap<u64, Shape>` keyed by
//!    `Shape::ptr_id()`.  The OCCT hasher equality is `TopoDS_Shape::IsSame`
//!    (TopTools_ShapeMapHasher.hxx L35-38) — the TShape handle only; the
//!    Location and the Orientation are NOT part of the key (the
//!    BRepTools_Substitution.rs `ShapeKey(ptr_id)` precedent).
//! 2. `NCollection_DataMap<TopoDS_Edge, NewCurveInfo> myNCInfo` and
//!    `<TopoDS_Face, NewSurfaceInfo> myNSInfo` use the same key carrier;
//!    `NCollection_Map<TopoDS_Shape> myNonUpdFace` / `myHasNewGeom` ->
//!    `HashSet<u64>`.
//! 3. `BRep_Builder B` -> the in-place TShape mutations below (the
//!    brep_algo::tool / brep_sweep::tool_rehost re-host style).  OCCT
//!    BRep_Builder edits the TShape through the shared handle, so a shape
//!    bound in myMap must be edited IN PLACE: `Arc::make_mut` would fork the
//!    Arc whenever the map (or any other container) holds the second
//!    reference, stranding the map slot on the stale data (the kernel
//!    `add_to_compound` / `edge_mut_inplace` SAFETY rationale).
//! 4. `TopLoc_Location` -> the u32 location id (0 = identity) — the
//!    draft_modification.rs arch. diff. #4 carrier; the OCCT
//!    `Predivided` needs the owning transform table (top_loc_predivided).
//! 5. `Message_ProgressScope` / `Message_ProgressRange` have no rcad
//!    counterpart (this translation carries no progress object): the OCCT
//!    "processing was broken" guards (L132-136 and L364-368) are unreachable
//!    and are recorded at their anchors.
//! 6. `Poly_Triangulation` / `Poly_Polygon3D` /
//!    `Poly_PolygonOnTriangulation` have no rcad carrier on the BRep
//!    entities (TFaceData carries no triangulation slot, TEdgeData no polygon
//!    representation), so the OCCT mesh branches keep their untranslated
//!    BRep_Builder step as a documented panic.
//! 7. `TopoDS_Shape& result = myMap(S)` — the OCCT live reference into the
//!    map — is the write-through [`BRepToolsModifier::map_put`] below: the
//!    rebuilt value is stored as soon as the OCCT `result` is re-assigned, so
//!    the `myMap(shape)` reads of the method body (L440/L445/L455/L476/L508/
//!    L514/L521/L538/L542/L569) observe the same value OCCT does.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use indexmap::IndexMap;

use rcad_kernel::geom::{Curve2d, Curve3, Surface3};
use rcad_kernel::precision::{is_negative_infinite_value, is_positive_infinite_value};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{
    tshape_flags, GeomAbsShape, Orientation, ShapeType, TEdgeData, TFaceData, TShape,
};

use crate::brep_algo::tool as bat;
use crate::brep_sweep::tool_rehost::{
    brep_tool_degenerated, brep_tool_is_closed_edge_face, brep_tools_is_really_closed,
};
use crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors;
use crate::topalgo::brep_tools_modification::BRepToolsModification;

/// OCCT TopTools_ShapeMapHasher key (TopTools_ShapeMapHasher.hxx L35-38:
/// `IsEqual(S1, S2) = S1.IsSame(S2)`) — the TShape handle.
type ShapeKey = u64;

/// OCCT TopTools_ShapeMapHasher hash/equality key of a shape.
fn shape_key(the_s: &Shape) -> ShapeKey {
    the_s.ptr_id()
}

/// OCCT NCollection_Map<TopoDS_Shape> ancestor list of an edge (the
/// `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>`
/// carrier of TopExp::MapShapesAndAncestors).
type AncestorMap = IndexMap<(u64, u32), (Shape, Vec<Shape>)>;

/// OCCT BRepTools_Modifier::NewCurveInfo (hxx L80-85).
#[derive(Clone)]
struct NewCurveInfo {
    /// OCCT: myCurve.
    my_curve: Option<Curve3>,
    /// OCCT: myLoc (the rcad u32 location id).
    my_loc: u32,
    /// OCCT: myToler.
    my_toler: f64,
}

/// OCCT BRepTools_Modifier::NewSurfaceInfo (hxx L87-94).
#[derive(Clone)]
struct NewSurfaceInfo {
    /// OCCT: mySurface.
    my_surface: Option<Surface3>,
    /// OCCT: myLoc.
    my_loc: u32,
    /// OCCT: myToler.
    my_toler: f64,
    /// OCCT: myRevWires.
    my_rev_wires: bool,
    /// OCCT: myRevFace.
    my_rev_face: bool,
}

/// OCCT BRepTools_Modifier (hxx L41-134) — performs geometric modifications on
/// a shape.
pub struct BRepToolsModifier {
    /// OCCT: myMap (NCollection_DataMap<TopoDS_Shape, TopoDS_Shape>).
    my_map: HashMap<ShapeKey, Shape>,
    /// OCCT: myShape.
    my_shape: Shape,
    /// OCCT: myDone.
    my_done: bool,
    /// OCCT: myNCInfo (NCollection_DataMap<TopoDS_Edge, NewCurveInfo>).
    my_nc_info: HashMap<ShapeKey, NewCurveInfo>,
    /// OCCT: myNSInfo (NCollection_DataMap<TopoDS_Face, NewSurfaceInfo>).
    my_ns_info: HashMap<ShapeKey, NewSurfaceInfo>,
    /// OCCT: myNonUpdFace (NCollection_Map<TopoDS_Shape>).
    my_non_upd_face: HashSet<ShapeKey>,
    /// OCCT: myHasNewGeom (NCollection_Map<TopoDS_Shape>).
    my_has_new_geom: HashSet<ShapeKey>,
    /// OCCT: myMutableInput.
    my_mutable_input: bool,
}

impl Default for BRepToolsModifier {
    fn default() -> Self {
        Self::new(false)
    }
}

impl BRepToolsModifier {
    // -----------------------------------------------------------------------
    // OCCT BRepTools_Modifier.cxx L53-88 — the constructors and Init.
    // -----------------------------------------------------------------------

    /// OCCT BRepTools_Modifier::BRepTools_Modifier(bool theMutableInput)
    /// (cxx L53-57).
    pub fn new(the_mutable_input: bool) -> Self {
        BRepToolsModifier {
            my_map: HashMap::new(),
            my_shape: Shape::null(),
            my_done: false,
            my_nc_info: HashMap::new(),
            my_ns_info: HashMap::new(),
            my_non_upd_face: HashSet::new(),
            my_has_new_geom: HashSet::new(),
            my_mutable_input: the_mutable_input,
        }
    }

    /// OCCT BRepTools_Modifier::BRepTools_Modifier(const TopoDS_Shape& S)
    /// (cxx L61-67).
    pub fn with_shape(the_s: &Shape) -> Self {
        let mut a_modifier = BRepToolsModifier::new(false);
        a_modifier.my_shape = the_s.clone();
        // OCCT L66: Put(S);
        a_modifier.put(the_s);
        a_modifier
    }

    /// OCCT BRepTools_Modifier::BRepTools_Modifier(S, M) (cxx L71-79).
    pub fn with_shape_and_modification(
        the_s: &Shape,
        the_m: &mut dyn BRepToolsModification,
    ) -> Self {
        let mut a_modifier = BRepToolsModifier::new(false);
        a_modifier.my_shape = the_s.clone();
        // OCCT L77: Put(S);
        a_modifier.put(the_s);
        // OCCT L78: Perform(M);
        a_modifier.perform(the_m);
        a_modifier
    }

    /// OCCT BRepTools_Modifier::Init(S) (cxx L83-88).
    pub fn init(&mut self, the_s: &Shape) {
        self.my_shape = the_s.clone();
        self.my_done = false;
        // OCCT L87: Put(S);
        self.put(the_s);
    }

    // -----------------------------------------------------------------------
    // OCCT BRepTools_Modifier.cxx L96-207 — Perform.
    // -----------------------------------------------------------------------

    /// OCCT BRepTools_Modifier::Perform(M, theProgress) (cxx L96-207).
    pub fn perform(&mut self, the_m: &mut dyn BRepToolsModification) {
        // OCCT L99-102: if (myShape.IsNull()) throw Standard_NullObject();
        if self.my_shape.is_null() {
            panic!("Standard_NullObject");
        }
        // OCCT L109: NCollection_DataMap<...>::Iterator theIter(myMap) — the
        // iterator is created over myMap and its only use is theIter.Next()
        // at L198; its Value() is never read (the OCCT leftover of the
        // removed DEBUG_Modifier block).  rcad: the iterator is not
        // materialised (a live borrow of my_map would conflict with the
        // my_map writes of the loop body below) — the call at L198 is the
        // no-op it is in OCCT.

        // OCCT L113-116: the aMVE / aMEF ancestor maps.
        let mut a_mve: AncestorMap = IndexMap::new();
        let mut a_mef: AncestorMap = IndexMap::new();
        let a_shape = self.my_shape.clone();
        map_shapes_and_ancestors(&a_shape, ShapeType::Vertex, ShapeType::Edge, &mut a_mve);
        map_shapes_and_ancestors(&a_shape, ShapeType::Edge, ShapeType::Face, &mut a_mef);

        // OCCT L118-122.
        self.create_new_vertices(&a_mve, the_m);
        self.fill_new_curve_info(&a_mef, the_m);
        self.fill_new_surface_info(the_m);

        // OCCT L124-127: if (!myMutableInput) CreateOtherVertices(aMVE, aMEF, M);
        if !self.my_mutable_input {
            self.create_other_vertices(&a_mve, &a_mef, the_m);
        }

        // OCCT L129-130: bool aNewGeom; Rebuild(myShape, M, aNewGeom, aPS.Next());
        let mut a_new_geom = false;
        self.rebuild(&a_shape, the_m, &mut a_new_geom);

        // OCCT L132-136: if (!aPS.More()) return; — no progress carrier
        // (arch. diff. #5); the scope never breaks.

        // OCCT L138-152: the root orientation.
        let a_root_key = shape_key(&self.my_shape);
        if let Some(a_root) = self.my_map.get(&a_root_key).cloned() {
            let mut a_map_root = a_root;
            if self.my_shape.shape_type() == ShapeType::Face {
                if self.my_shape.orientation == Orientation::Reversed {
                    // OCCT L142: myMap(myShape).Reverse();
                    a_map_root.orientation = bat::top_abs_reverse(a_map_root.orientation);
                } else {
                    // OCCT L146: myMap(myShape).Orientation(myShape.Orientation());
                    a_map_root.orientation = self.my_shape.orientation;
                }
            } else {
                // OCCT L151: myMap(myShape).Orientation(myShape.Orientation());
                a_map_root.orientation = self.my_shape.orientation;
            }
            self.my_map.insert(a_root_key, a_map_root);
        }

        // Update the continuities (OCCT L154-199).
        for ii in 0..a_mef.len() {
            let cur_e = a_mef[ii].0.clone();
            // OCCT L169: const TopoDS_Edge& NewE = TopoDS::Edge(myMap(CurE));
            let new_e = self.my_map_value(&cur_e);
            // OCCT L170: if (!CurE.IsSame(NewE))
            if !cur_e.is_same(&new_e) {
                // OCCT L172-186: the two-face pick from aMEF.FindFromKey(CurE).
                let a_l_faces = &a_mef[ii].1;
                let mut f1 = Shape::null();
                let mut f2 = Shape::null();
                let mut it = 0usize;
                while it < a_l_faces.len() && f2.is_null() {
                    if f1.is_null() {
                        f1 = a_l_faces[it].clone();
                    } else {
                        f2 = a_l_faces[it].clone();
                    }
                    it += 1;
                }
                // OCCT L187: if (!F2.IsNull())
                if !f2.is_null() {
                    // OCCT L189-190.
                    let new_f1 = self.my_map_value(&f1);
                    let new_f2 = self.my_map_value(&f2);
                    // OCCT L191: M->Continuity(CurE, F1, F2, NewE, newf1, newf2);
                    let new_cont =
                        the_m.continuity(&cur_e, &f1, &f2, &new_e, &new_f1, &new_f2);
                    // OCCT L192-195: if (Newcont > GeomAbs_C0)
                    // aBB.Continuity(NewE, newf1, newf2, Newcont);
                    if new_cont > GeomAbsShape::C0 {
                        builder_continuity(&new_e, &new_f1, &new_f2, new_cont);
                    }
                }
            }
            // OCCT L198: theIter.Next(); — the dead iterator (see L109).
        }

        // OCCT L206: myDone = true;
        self.my_done = true;
    }

    // -----------------------------------------------------------------------
    // OCCT BRepTools_Modifier.cxx L211-222 — Put.
    // -----------------------------------------------------------------------

    /// OCCT BRepTools_Modifier::Put(S) (cxx L211-222).
    pub fn put(&mut self, the_s: &Shape) {
        // OCCT L213: if (!myMap.IsBound(S))
        let key = shape_key(the_s);
        if !self.my_map.contains_key(&key) {
            // OCCT L215: myMap.Bind(S, TopoDS_Shape());
            self.my_map.insert(key, Shape::null());
            // OCCT L216-219: for (TopoDS_Iterator theIterator(S, false); ...)
            for a_child in topods_iterator(the_s) {
                self.put(&a_child);
            }
        }
    }

    // -----------------------------------------------------------------------
    // OCCT BRepTools_Modifier.cxx L226-649 — Rebuild.
    // -----------------------------------------------------------------------

    /// OCCT BRepTools_Modifier::Rebuild(S, M, theNewGeom, theProgress)
    /// (cxx L226-649).
    pub fn rebuild(
        &mut self,
        the_s: &Shape,
        the_m: &mut dyn BRepToolsModification,
        the_new_geom: &mut bool,
    ) -> bool {
        // OCCT L235-241.
        let ts = the_s.shape_type();
        let key = shape_key(the_s);
        let mut result = self.my_map_value(the_s);
        if !result.is_null() {
            *the_new_geom = self.my_has_new_geom.contains(&key);
            return !the_s.is_same(&result);
        }
        // OCCT L242-247.
        let mut rebuild = false;
        let mut rev_wires = false;
        let mut res_or = Orientation::Forward;
        let mut tol = 0.0f64;
        let mut no_3d_curve = false;

        // new geometry ? (OCCT L251-336)

        match ts {
            ShapeType::Face => {
                // OCCT L253: rebuild = myNSInfo.IsBound(TopoDS::Face(S));
                rebuild = self.my_ns_info.contains_key(&key);
                if rebuild {
                    // OCCT L257: const NewSurfaceInfo& aNSinfo = myNSInfo(Face(S));
                    let a_ns_info = self.my_ns_info.get(&key).cloned().unwrap();
                    // OCCT L258: RevWires = aNSinfo.myRevWires;
                    rev_wires = a_ns_info.my_rev_wires;
                    // OCCT L259-262: B.MakeFace(TopoDS::Face(result),
                    // aNSinfo.mySurface, aNSinfo.myLoc.Predivided(S.Location()),
                    // aNSinfo.myToler);
                    let a_loc = top_loc_predivided(a_ns_info.my_loc, the_s.location);
                    result = self.map_put(
                        key,
                        make_face(a_ns_info.my_surface.clone(), a_loc, a_ns_info.my_toler),
                    );
                    // OCCT L263: result.Location(S.Location(), false);
                    result.location = the_s.location;
                    // OCCT L264-267: if (aNSinfo.myRevFace) ResOr = TopAbs_REVERSED;
                    if a_ns_info.my_rev_face {
                        res_or = Orientation::Reversed;
                    }
                    // OCCT L269: B.NaturalRestriction(TopoDS::Face(result),
                    // BRep_Tool::NaturalRestriction(TopoDS::Face(S)));
                    builder_natural_restriction(&result, brep_tool_natural_restriction(the_s));
                }

                // OCCT L272-286: update triangulation on the copied face.
                let mut a_triangulation = None;
                if the_m.new_triangulation(the_s, &mut a_triangulation) {
                    if rebuild {
                        // OCCT L278: B.UpdateFace(TopoDS::Face(result), aTriangulation);
                        panic!(
                            "GAP: BRep_Builder::UpdateFace(F, T) (BRep_Builder.cxx \
                             L577-591) not translated — rcad BRep faces carry no \
                             Poly_Triangulation (TFaceData has no triangulation slot)"
                        );
                    } else {
                        // OCCT L282-283: B.MakeFace(TopoDS::Face(result),
                        // aTriangulation); result.Location(S.Location(), false);
                        panic!(
                            "GAP: BRep_Builder::MakeFace(F, T) (BRep_Builder.cxx \
                             L521-545) not translated — rcad BRep faces carry no \
                             Poly_Triangulation (TFaceData has no triangulation slot)"
                        );
                    }
                    // OCCT L285: rebuild = true;
                }
            }

            ShapeType::Edge => {
                // OCCT L291: rebuild = myNCInfo.IsBound(TopoDS::Edge(S));
                rebuild = self.my_nc_info.contains_key(&key);
                if rebuild {
                    // OCCT L294: const NewCurveInfo& aNCinfo = myNCInfo(Edge(S));
                    let a_nc_info = self.my_nc_info.get(&key).cloned().unwrap();
                    // OCCT L295: if (aNCinfo.myCurve.IsNull())
                    if a_nc_info.my_curve.is_none() {
                        // OCCT L297-300.
                        result = self.map_put(key, bat::builder_make_edge());
                        builder_degenerated(&result, brep_tool_degenerated(the_s));
                        // OCCT L299: B.UpdateEdge(Edge(result), aNCinfo.myToler); // OCC217
                        builder_update_edge_tolerance(&result, a_nc_info.my_toler);
                        no_3d_curve = true;
                    } else {
                        // OCCT L304-307: B.MakeEdge(Edge(result),
                        // aNCinfo.myCurve, aNCinfo.myLoc.Predivided(S.Location()),
                        // aNCinfo.myToler);
                        let a_loc = top_loc_predivided(a_nc_info.my_loc, the_s.location);
                        result = self.map_put(
                            key,
                            make_edge(a_nc_info.my_curve.clone(), a_loc, a_nc_info.my_toler),
                        );
                        no_3d_curve = false;
                    }
                    // OCCT L310: result.Location(S.Location(), false);
                    result.location = the_s.location;
                    // OCCT L311 (commented out in OCCT): result.Orientation(S.Orientation()).

                    // OCCT L314-315: B.SameParameter(Edge(result),
                    // BRep_Tool::SameParameter(Edge(S)));
                    builder_same_parameter(&result, brep_tool_same_parameter(the_s));
                    builder_same_range(&result, brep_tool_same_range(the_s));
                }

                // OCCT L319-332: update polygonal structure on the edge.
                if the_m.new_polygon(the_s) {
                    if rebuild {
                        // OCCT L324: B.UpdateEdge(Edge(result), aPolygon, S.Location());
                        panic!(
                            "GAP: BRep_Builder::UpdateEdge(E, P, L) \
                             (BRep_Builder.cxx L858-905) not translated — rcad BRep \
                             edges carry no Poly_Polygon3D representation"
                        );
                    } else {
                        // OCCT L328-329: B.MakeEdge(Edge(result), aPolygon);
                        // result.Location(S.Location(), false);
                        panic!(
                            "GAP: BRep_Builder::MakeEdge(E, P) (BRep_Builder.cxx \
                             L800-845) not translated — rcad BRep edges carry no \
                             Poly_Polygon3D representation"
                        );
                    }
                    // OCCT L331: rebuild = true;
                }
            }
            _ => {}
        }

        // rebuild sub-shapes and test new sub-shape ? (OCCT L338-369)

        let newgeom = rebuild;
        // OCCT L341: theNewGeom = rebuild;
        *the_new_geom = rebuild;

        // OCCT L343-352: TopoDS_Iterator it; int aShapeCount = 0;
        // (TopoDS_Iterator(S, false) — the raw stored children.)
        let it_children = topods_iterator(the_s);
        let _a_shape_count = it_children.len();

        // OCCT L354-368: Message_ProgressScope aPS(theProgress, ...);
        // for (it.Initialize(S, false); it.More() && aPS.More(); it.Next())
        for a_child in &it_children {
            // always call Rebuild
            let mut is_sub_new_geom = false;
            let subrebuilt = self.rebuild(a_child, the_m, &mut is_sub_new_geom);
            rebuild = subrebuilt || rebuild;
            // OCCT L362: theNewGeom = theNewGeom || isSubNewGeom;
            *the_new_geom = *the_new_geom || is_sub_new_geom;
        }
        // OCCT L364-368: if (!aPS.More()) return false; — no progress carrier
        // (arch. diff. #5).

        // OCCT L370-373.
        if *the_new_geom {
            self.my_has_new_geom.insert(key);
        }

        // make an empty copy (OCCT L375-380)
        if rebuild && !newgeom {
            // OCCT L378-379: result = S.EmptyCopied();
            // result.Orientation(TopAbs_FORWARD);
            result = bat::empty_copied(the_s);
            result.orientation = Orientation::Forward;
        }

        // copy the sub-elements (OCCT L382-637)

        if rebuild {
            for a_child in &it_children {
                // OCCT L389: orient = it.Value().Orientation();
                let mut orient = a_child.orientation;
                // OCCT L390: if (RevWires || myMap(it.Value()).Orientation()
                // == TopAbs_REVERSED) orient = TopAbs::Reverse(orient);
                let a_child_mapped = self.my_map_value(a_child);
                if rev_wires || a_child_mapped.orientation == Orientation::Reversed {
                    orient = bat::top_abs_reverse(orient);
                }
                // OCCT L394: B.Add(result, myMap(it.Value()).Oriented(orient));
                let a_child_oriented = bat::oriented(&a_child_mapped, orient);
                builder_add(&result, &a_child_oriented);
            }

            if ts == ShapeType::Face {
                // pcurves (OCCT L397-587)
                // OCCT L401: TopoDS_Face face = TopoDS::Face(S);
                let face = the_s.clone();
                // OCCT L402-406.
                let mut fcor = face.orientation;
                if fcor != Orientation::Reversed {
                    fcor = Orientation::Forward;
                }

                // OCCT L408: TopExp_Explorer ex(face.Oriented(fcor), TopAbs_EDGE);
                let a_face_oriented = bat::oriented(&face, fcor);
                for edge in bat::explorer(&a_face_oriented, ShapeType::Edge, ShapeType::Shape) {
                    // OCCT L416-422: if (theNewGeom &&
                    // M->NewCurve2d(edge, face, TopoDS::Edge(myMap(ex.Current())),
                    // TopoDS::Face(result), curve2d, tol))
                    let mut curve2d: Option<Curve2d> = None;
                    let mut a_new_edge = self.my_map_value(&edge);
                    let a_new_curve2d = if *the_new_geom {
                        the_m.new_curve2d(
                            &edge,
                            &face,
                            &mut a_new_edge,
                            &result,
                            &mut curve2d,
                            &mut tol,
                        )
                    } else {
                        false
                    };
                    // The NewE handle of the OCCT call is a temporary
                    // TopoDS_Edge aliasing myMap(ex.Current()); the builder
                    // edits the shared TShape make them visible on the map
                    // value.  The rcad value semantics need the explicit
                    // write-back of the (possibly edited) NewE.
                    self.map_put(shape_key(&edge), a_new_edge.clone());

                    if a_new_curve2d {
                        // OCCT L433-473: the closed-pcurve test.
                        let mut is_closed = false;
                        // OCCT L434: if (BRep_Tool::IsClosed(edge, face))
                        if brep_tool_is_closed_edge_face(&edge, &face) {
                            // OCCT L436: isClosed = (!newgeom ||
                            // BRepTools::IsReallyClosed(edge, face));
                            is_closed = !newgeom || brep_tools_is_really_closed(&edge, &face);
                            if !is_closed {
                                // OCCT L439-445.
                                let resface = if self.my_map.contains_key(&shape_key(&face)) {
                                    self.my_map_value(&face)
                                } else {
                                    face.clone()
                                };
                                let resface = if resface.is_null() { face.clone() } else { resface };
                                let a_loc = resface.location;
                                let a_surf = brep_tool_surface_value(&resface);
                                // OCCT L447: TopExp_Explorer aExpF(myShape, TopAbs_FACE);
                                let a_exp_f = bat::explorer(
                                    &self.my_shape.clone(),
                                    ShapeType::Face,
                                    ShapeType::Shape,
                                );
                                for an_other in &a_exp_f {
                                    if is_closed {
                                        break;
                                    }
                                    // OCCT L451-454: if (anOther.IsSame(face)) continue;
                                    if an_other.is_same(&face) {
                                        continue;
                                    }
                                    // OCCT L455-459.
                                    let resface2 =
                                        if self.my_map.contains_key(&shape_key(an_other)) {
                                            self.my_map_value(an_other)
                                        } else {
                                            an_other.clone()
                                        };
                                    let resface2 =
                                        if resface2.is_null() { an_other.clone() } else { resface2 };
                                    let an_other_loc = resface2.location;
                                    let an_other_surf = brep_tool_surface_value(&resface2);
                                    // OCCT L463: if (aSurf == anOtherSurf &&
                                    // aLoc.IsEqual(anOtherLoc))
                                    if surface_handle_equal(&a_surf, &an_other_surf)
                                        && a_loc == an_other_loc
                                    {
                                        // OCCT L465-469: TopExp_Explorer
                                        // aExpE(anOther, TopAbs_EDGE);
                                        for a_exp_e in bat::explorer(
                                            an_other,
                                            ShapeType::Edge,
                                            ShapeType::Shape,
                                        ) {
                                            if is_closed {
                                                break;
                                            }
                                            // OCCT L468: isClosed =
                                            // edge.IsSame(aExpE.Current());
                                            is_closed = edge.is_same(&a_exp_e);
                                        }
                                    }
                                }
                            }
                        }
                        if is_closed {
                            // OCCT L476-504.
                            let mut cur_e = self.my_map_value(&edge);
                            // OCCT L478: aLocalResult.Orientation(TopAbs_FORWARD);
                            let cur_f = {
                                let mut a_local_result = result.clone();
                                a_local_result.orientation = Orientation::Forward;
                                a_local_result
                            };
                            let mut f = 0.0f64;
                            let mut l = 0.0f64;
                            // OCCT L482-483: if ((!RevWires && fcor !=
                            // edge.Orientation()) || (RevWires && fcor ==
                            // edge.Orientation()))
                            if (!rev_wires && fcor != edge.orientation)
                                || (rev_wires && fcor == edge.orientation)
                            {
                                // OCCT L485-491.
                                cur_e.orientation = Orientation::Forward;
                                let a_c = bat::brep_tool_curve_on_surface(&cur_e, &cur_f);
                                if let Some((_, cf, cl)) = a_c {
                                    f = cf;
                                    l = cl;
                                }
                                // OCCT L487-490: if (curve2d1.IsNull())
                                // curve2d1 = new Geom2d_Line(gp::OX2d());
                                let curve2d1 = a_c
                                    .map(|(c, _, _)| c)
                                    .unwrap_or_else(|| crate::brep_sweep::tool_rehost::geom2d_line(
                                        glam::DVec2::ZERO,
                                        glam::DVec2::X,
                                    ));
                                // OCCT L491: B.UpdateEdge(CurE, curve2d1,
                                // curve2d, CurF, 0.);
                                builder_update_edge_pcurve2(
                                    &cur_e,
                                    &curve2d1,
                                    curve2d.as_ref().unwrap(),
                                    &cur_f,
                                    0.0,
                                );
                            } else {
                                // OCCT L493-501.
                                cur_e.orientation = Orientation::Reversed;
                                let a_c = bat::brep_tool_curve_on_surface(&cur_e, &cur_f);
                                if let Some((_, cf, cl)) = a_c {
                                    f = cf;
                                    l = cl;
                                }
                                let curve2d1 = a_c
                                    .map(|(c, _, _)| c)
                                    .unwrap_or_else(|| crate::brep_sweep::tool_rehost::geom2d_line(
                                        glam::DVec2::ZERO,
                                        glam::DVec2::X,
                                    ));
                                // OCCT L501: B.UpdateEdge(CurE, curve2d,
                                // curve2d1, CurF, 0.);
                                builder_update_edge_pcurve2(
                                    &cur_e,
                                    curve2d.as_ref().unwrap(),
                                    &curve2d1,
                                    &cur_f,
                                    0.0,
                                );
                            }
                            // OCCT L503-504: currcurv = BRep_Tool::CurveOnSurface(CurE,
                            // CurF, f, l); B.Range(CurE, f, l);
                            if let Some((_, cf, cl)) =
                                bat::brep_tool_curve_on_surface(&cur_e, &cur_f)
                            {
                                f = cf;
                                l = cl;
                            }
                            builder_range_edge(&cur_e, f, l);
                            // The CurE write-back (the OCCT map-slot aliasing).
                            self.map_put(shape_key(&edge), cur_e.clone());
                        } else {
                            // OCCT L508: B.UpdateEdge(TopoDS::Edge(myMap(ex.Current())),
                            // curve2d, TopoDS::Face(result), 0.);
                            let a_map_edge = self.my_map_value(&edge);
                            builder_update_edge_pcurve(
                                &a_map_edge,
                                curve2d.as_ref().unwrap(),
                                &result,
                                0.0,
                            );
                        }

                        // OCCT L511-545: the 3D curve null check -> the vertex
                        // parameter update.
                        let mut the_loc = 0u32;
                        let mut the_f = 0.0f64;
                        let mut the_l = 0.0f64;
                        // OCCT L513-514: occ::handle<Geom_Curve> C3D =
                        // BRep_Tool::Curve(Edge(myMap(ex.Current())), theLoc,
                        // theF, theL);
                        let a_map_edge = self.my_map_value(&edge);
                        let c3d = brep_tool_curve_with_location(
                            &a_map_edge,
                            &mut the_loc,
                            &mut the_f,
                            &mut the_l,
                        );
                        if c3d.is_none() {
                            // Update vertices (OCCT L516-545)
                            let mut param = 0.0f64;
                            let a_vertices =
                                bat::explorer(&edge, ShapeType::Vertex, ShapeType::Shape);
                            for vertex in &a_vertices {
                                // OCCT L522-526.
                                if !the_m.new_parameter(vertex, &edge, &mut param, &mut tol) {
                                    tol = bat::brep_tool_tolerance(vertex);
                                    param = bat::brep_tool_parameter(vertex, &edge);
                                }

                                // OCCT L528-533.
                                let mut vtxrelat = vertex.orientation;
                                if edge.orientation == Orientation::Reversed {
                                    // Update considers the edge as FORWARD, and
                                    // the vertex relatively
                                    vtxrelat = bat::top_abs_reverse(vtxrelat);
                                }

                                // OCCT L538-539.
                                let mut a_local_vertex = self.my_map_value(vertex);
                                a_local_vertex.orientation = vtxrelat;
                                // OCCT L542: B.UpdateVertex(aLocalVertex, param,
                                // TopoDS::Edge(myMap(edge)), tol);
                                let a_map_edge2 = self.my_map_value(&edge);
                                builder_update_vertex_param(
                                    &a_local_vertex,
                                    param,
                                    &a_map_edge2,
                                    tol,
                                );
                            }
                        }
                    }

                    // Copy polygon on triangulation (OCCT L548-585)
                    let a_new_pon_t = the_m.new_polygon_on_triangulation(&edge, &face);
                    if a_new_pon_t {
                        if brep_tools_is_really_closed(&edge, &face) {
                            // OCCT L553-562: obtain the triangulation on the
                            // reversed edge.
                            let mut an_edge_rev = edge.clone();
                            an_edge_rev.orientation =
                                bat::top_abs_reverse(an_edge_rev.orientation);
                            the_m.new_polygon_on_triangulation(&an_edge_rev, &face);
                        }
                        // OCCT L564-585: B.UpdateEdge(aNewEdge, ...Poly_PolygonOnTriangulation...)
                        panic!(
                            "GAP: BRep_Builder::UpdateEdge(E, P, T, L) \
                             (BRep_Builder.cxx L940-1010) not translated — rcad BRep \
                             edges carry no Poly_PolygonOnTriangulation representation"
                        );
                    }
                }
            } else if ts == ShapeType::Edge && !no_3d_curve {
                // Vertices (OCCT L590-630)
                let mut param = 0.0f64;
                // OCCT L594-595: const TopoDS_Edge& edge = TopoDS::Edge(S);
                let edge = the_s.clone();
                let mut edor = edge.orientation;
                // OCCT L596-599.
                if edor != Orientation::Reversed {
                    edor = Orientation::Forward;
                }
                // OCCT L600: TopExp_Explorer ex(edge.Oriented(edor), TopAbs_VERTEX);
                let a_edge_oriented = bat::oriented(&edge, edor);
                let a_vertices =
                    bat::explorer(&a_edge_oriented, ShapeType::Vertex, ShapeType::Shape);
                for vertex in &a_vertices {
                    // OCCT L605-609.
                    if !the_m.new_parameter(vertex, &edge, &mut param, &mut tol) {
                        tol = bat::brep_tool_tolerance(vertex);
                        param = bat::brep_tool_parameter(vertex, &edge);
                    }

                    // OCCT L611-616.
                    let mut vtxrelat = vertex.orientation;
                    if edor == Orientation::Reversed {
                        // Update considers the edge as FORWARD, and the vertex
                        // relatively
                        vtxrelat = bat::top_abs_reverse(vtxrelat);
                    }

                    // OCCT L621-622.
                    let mut a_local_vertex = self.my_map_value(vertex);
                    a_local_vertex.orientation = vtxrelat;
                    // OCCT L624-627: if (myMutableInput || !aLocalVertex.IsSame(vertex))
                    // B.UpdateVertex(aLocalVertex, param, TopoDS::Edge(result), tol);
                    if self.my_mutable_input || !a_local_vertex.is_same(vertex) {
                        builder_update_vertex_param(
                            &a_local_vertex,
                            param,
                            &result,
                            tol,
                        );
                    }
                }
            }

            // update flags (OCCT L632-636)
            // result.Orientable(S.Orientable()); result.Closed(S.Closed());
            // result.Infinite(S.Infinite());
            topods_set_flag(
                &result,
                tshape_flags::ORIENTABLE,
                topods_flag(the_s, tshape_flags::ORIENTABLE),
            );
            topods_set_flag(
                &result,
                tshape_flags::CLOSED,
                topods_flag(the_s, tshape_flags::CLOSED),
            );
            topods_set_flag(
                &result,
                tshape_flags::INFINITE,
                topods_flag(the_s, tshape_flags::INFINITE),
            );
        } else {
            // OCCT L640: result = S;
            result = the_s.clone();
        }

        // Set flag of the shape. (OCCT L643-644)
        // result.Orientation(ResOr);
        result.orientation = res_or;

        // OCCT L646: SetShapeFlags(S, result);
        set_shape_flags(the_s, &result);

        // The OCCT `TopoDS_Shape& result = myMap(S)` aliasing: the rebuilt
        // value is written back to the map slot (arch. diff. #7).
        self.map_put(key, result.clone());

        // OCCT L648: return rebuild;
        rebuild
    }

    // -----------------------------------------------------------------------
    // OCCT BRepTools_Modifier.cxx L651-832 — the auxiliary fill passes.
    // -----------------------------------------------------------------------

    /// OCCT BRepTools_Modifier::CreateNewVertices (cxx L651-678).
    pub fn create_new_vertices(
        &mut self,
        the_mve: &AncestorMap,
        the_m: &mut dyn BRepToolsModification,
    ) {
        let mut a_toler = 0.0f64;
        let mut a_pnt = glam::DVec3::ZERO;
        for ii in 0..the_mve.len() {
            // fill MyMap only with vertices with NewPoint == true
            // OCCT L663: const TopoDS_Vertex& aV = TopoDS::Vertex(theMVE.FindKey(i));
            let a_v = the_mve[ii].0.clone();
            // OCCT L664: bool IsNewP = M->NewPoint(aV, aPnt, aToler);
            let is_new_p = the_m.new_point(&a_v, &mut a_pnt, &mut a_toler);
            if is_new_p {
                // OCCT L667-671.
                let a_new_v = make_vertex(a_pnt, a_toler);
                set_shape_flags(&a_v, &a_new_v);
                self.my_map.insert(shape_key(&a_v), a_new_v);
                self.my_has_new_geom.insert(shape_key(&a_v));
            } else if self.my_mutable_input {
                // OCCT L675: myMap(aV) = aV.Oriented(TopAbs_FORWARD);
                let a_oriented = bat::oriented(&a_v, Orientation::Forward);
                self.my_map.insert(shape_key(&a_v), a_oriented);
            }
        }
    }

    /// OCCT BRepTools_Modifier::FillNewCurveInfo (cxx L680-703).
    pub fn fill_new_curve_info(
        &mut self,
        the_mef: &AncestorMap,
        the_m: &mut dyn BRepToolsModification,
    ) {
        let mut a_curve: Option<Curve3> = None;
        let mut a_location = 0u32;
        let mut a_nc_info = NewCurveInfo {
            my_curve: None,
            my_loc: 0,
            my_toler: 0.0,
        };
        let mut a_toler = 0.0f64;
        for ii in 0..the_mef.len() {
            // OCCT L692: const TopoDS_Edge& anE = TopoDS::Edge(theMEF.FindKey(i));
            let an_e = the_mef[ii].0.clone();
            // OCCT L693: bool IsNewCur = M->NewCurve(anE, aCurve, aLocation, aToler);
            let is_new_cur = the_m.new_curve(&an_e, &mut a_curve, &mut a_location, &mut a_toler);
            if is_new_cur {
                // OCCT L696-700.
                a_nc_info.my_curve = a_curve.clone();
                a_nc_info.my_loc = a_location;
                a_nc_info.my_toler = a_toler;
                self.my_nc_info.insert(shape_key(&an_e), a_nc_info.clone());
                self.my_has_new_geom.insert(shape_key(&an_e));
            }
        }
    }

    /// OCCT BRepTools_Modifier::FillNewSurfaceInfo (cxx L705-762).
    pub fn fill_new_surface_info(&mut self, the_m: &mut dyn BRepToolsModification) {
        let mut a_ns_info = NewSurfaceInfo {
            my_surface: None,
            my_loc: 0,
            my_toler: 0.0,
            my_rev_wires: false,
            my_rev_face: false,
        };
        // OCCT L707-708: TopExp::MapShapes(myShape, TopAbs_FACE, aMF);
        let a_mf = bat::explorer(&self.my_shape.clone(), ShapeType::Face, ShapeType::Shape);
        for a_f in &a_mf {
            // OCCT L713-718.
            let mut rev_face = false;
            let mut rev_wires = false;
            let mut a_surface: Option<Surface3> = None;
            let mut a_location = 0u32;
            let mut a_toler1 = 0.0f64;
            // OCCT L718: bool IsNewSur = M->NewSurface(aF, aSurface, aLocation,
            // aToler1, RevWires, RevFace);
            let is_new_sur = the_m.new_surface(
                a_f,
                &mut a_surface,
                &mut a_location,
                &mut a_toler1,
                &mut rev_wires,
                &mut rev_face,
            );
            if is_new_sur {
                // OCCT L721-727.
                a_ns_info.my_surface = a_surface.clone();
                a_ns_info.my_loc = a_location;
                a_ns_info.my_toler = a_toler1;
                a_ns_info.my_rev_wires = rev_wires;
                a_ns_info.my_rev_face = rev_face;
                self.my_ns_info.insert(shape_key(a_f), a_ns_info.clone());
                self.my_has_new_geom.insert(shape_key(a_f));
            } else {
                // check if subshapes will be modified (OCCT L729-760)
                let mut not_rebuilded = true;
                // OCCT L733: TopExp_Explorer exE(aF, TopAbs_EDGE);
                for an_ee in bat::explorer(a_f, ShapeType::Edge, ShapeType::Shape) {
                    if !not_rebuilded {
                        break;
                    }
                    // OCCT L737: if (myNCInfo.IsBound(anEE))
                    if self.my_nc_info.contains_key(&shape_key(&an_ee)) {
                        not_rebuilded = false;
                        break;
                    }
                    // OCCT L742: TopExp_Explorer exV(anEE, TopAbs_VERTEX);
                    for a_vv in bat::explorer(&an_ee, ShapeType::Vertex, ShapeType::Shape) {
                        if !not_rebuilded {
                            break;
                        }
                        // OCCT L746: if (!myMap(aVV).IsNull())
                        if !self.my_map_value(&a_vv).is_null() {
                            not_rebuilded = false;
                            break;
                        }
                    }
                }
                if not_rebuilded {
                    // subshapes is not going to be modified
                    // OCCT L758: myNonUpdFace.Add(aF);
                    self.my_non_upd_face.insert(shape_key(a_f));
                }
            }
        }
    }

    /// OCCT BRepTools_Modifier::CreateOtherVertices (cxx L764-832).
    pub fn create_other_vertices(
        &mut self,
        the_mve: &AncestorMap,
        the_mef: &AncestorMap,
        the_m: &mut dyn BRepToolsModification,
    ) {
        // The following logic in some ways repeats the logic from the Rebuild()
        // method.
        // If the face with its subshapes is not going to be modified
        // (i.e. NewSurface() for this face and NewCurve(), NewPoint() for its
        // edges/vertices returns false) then the calling of NewCurve2d() for
        // this face with its edges is not performed.
        // Therefore, the updating of vertices will not present in such cases
        // and the EmptyCopied() operation for vertices from this face is not
        // needed.
        let mut a_toler = 0.0f64;
        for ii in 0..the_mve.len() {
            // OCCT L783-784.
            let a_v = the_mve[ii].0.clone();
            let a_new_v = self.my_map_value(&a_v);
            if a_new_v.is_null() {
                // OCCT L787: const NCollection_List<TopoDS_Shape>& aLEdges = theMVE(i);
                let a_l_edges = the_mve[ii].1.clone();
                let mut to_replace = false;
                let mut it = 0usize;
                // OCCT L789-790: for (; it.More() && !toReplace; it.Next())
                while it < a_l_edges.len() && !to_replace {
                    // OCCT L792: const TopoDS_Edge& anE = TopoDS::Edge(it.Value());
                    let an_e = a_l_edges[it].clone();
                    // OCCT L793: if (myNCInfo.IsBound(anE) &&
                    // !myNCInfo(anE).myCurve.IsNull()) toReplace = true;
                    if let Some(a_nc_info) = self.my_nc_info.get(&shape_key(&an_e)) {
                        if a_nc_info.my_curve.is_some() {
                            to_replace = true;
                        }
                    }

                    // OCCT L798-818.
                    if !to_replace {
                        let a_l_faces = self
                            .ancestors_of(the_mef, &an_e)
                            .unwrap_or_default();
                        for a_f in &a_l_faces {
                            // OCCT L805: if (!myNonUpdFace.Contains(aF))
                            if !self.my_non_upd_face.contains(&shape_key(a_f)) {
                                let mut a_curve2d: Option<Curve2d> = None;
                                // some NewCurve2d()s may use NewE arg internally,
                                // so the null TShape as an arg may lead to the
                                // exceptions
                                // OCCT L810: TopoDS_Edge aDummyE =
                                // TopoDS::Edge(anE.EmptyCopied());
                                let mut a_dummy_e = bat::empty_copied(&an_e);
                                // OCCT L811: M->NewCurve2d(anE, aF, aDummyE,
                                // TopoDS_Face(), aCurve2d, aToler)
                                let a_null_face = Shape::null();
                                if the_m.new_curve2d(
                                    &an_e,
                                    a_f,
                                    &mut a_dummy_e,
                                    &a_null_face,
                                    &mut a_curve2d,
                                    &mut a_toler,
                                ) {
                                    to_replace = true;
                                    break;
                                }
                            }
                        }
                    }
                    it += 1;
                }
                if to_replace {
                    // OCCT L822: aNewV = TopoDS::Vertex(aV.EmptyCopied());
                    let a_new_v = bat::empty_copied(&a_v);
                    self.my_map_vertex_set(&a_v, a_new_v);
                } else {
                    // OCCT L826: aNewV = aV;
                    self.my_map_vertex_set(&a_v, a_v.clone());
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // OCCT BRepTools_Modifier.cxx L844-852 — the accessors (+ the .lxx).
    // -----------------------------------------------------------------------

    /// OCCT BRepTools_Modifier::IsDone() (lxx L32-35).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT BRepTools_Modifier::IsMutableInput() (cxx L844-847).
    pub fn is_mutable_input(&self) -> bool {
        self.my_mutable_input
    }

    /// OCCT BRepTools_Modifier::SetMutableInput(theMutableInput)
    /// (cxx L849-852).
    pub fn set_mutable_input(&mut self, the_mutable_input: bool) {
        self.my_mutable_input = the_mutable_input;
    }

    /// OCCT BRepTools_Modifier::ModifiedShape(S) (lxx L21-28) — returns the
    /// modified shape corresponding to S.  The OCCT body throws
    /// Standard_NoSuchObject when S is not bound in myMap.
    pub fn modified_shape(&self, the_s: &Shape) -> Shape {
        match self.my_map.get(&shape_key(the_s)) {
            Some(a_shape) => a_shape.clone(),
            None => panic!("Standard_NoSuchObject"),
        }
    }

    /// The rcad rvalue form of `myMap(S)` — `None` when S is not bound (the
    /// OCCT `myMap.IsBound(S)` guard).
    pub fn modified_shape_opt(&self, the_s: &Shape) -> Option<Shape> {
        self.my_map.get(&shape_key(the_s)).cloned()
    }

    // -----------------------------------------------------------------------
    // Internal carriers of the OCCT map semantics.
    // -----------------------------------------------------------------------

    /// OCCT `myMap(S)` — the item of the map (NCollection_DataMap::operator()
    /// raises Standard_NoSuchObject when the key is not bound,
    /// NCollection_DataMap.hxx L672-678/L693).
    fn my_map_value(&self, the_s: &Shape) -> Shape {
        match self.my_map.get(&shape_key(the_s)) {
            Some(a_shape) => a_shape.clone(),
            None => panic!("Standard_NoSuchObject: NCollection_DataMap::Find"),
        }
    }

    /// The `result` write-through of the OCCT `TopoDS_Shape& result = myMap(S)`
    /// alias (arch. diff. #7): Bind(S, theShape) + return the item.
    fn map_put(&mut self, the_key: ShapeKey, the_shape: Shape) -> Shape {
        self.my_map.insert(the_key, the_shape.clone());
        the_shape
    }

    /// OCCT L824-829: `aNewV.Orientation(TopAbs_FORWARD); myMap(aV) = aNewV;`.
    fn my_map_vertex_set(&mut self, the_v: &Shape, the_new_v: Shape) {
        let mut a_new_v = the_new_v;
        a_new_v.orientation = Orientation::Forward;
        self.my_map.insert(shape_key(the_v), a_new_v);
    }

    /// OCCT `aMEF.FindFromKey(anE)` — the ancestor list of a key (the empty
    /// list when the key is not bound).
    fn ancestors_of(&self, the_map: &AncestorMap, the_s: &Shape) -> Option<Vec<Shape>> {
        the_map
            .get(&bat::shape_key(the_s))
            .map(|(_, a_list)| a_list.clone())
    }
}

// ---------------------------------------------------------------------------
// OCCT BRepTools_Modifier.cxx L834-842 — SetShapeFlags.
// ---------------------------------------------------------------------------

/// OCCT BRepTools_Modifier::SetShapeFlags (cxx L834-842, the file-static
/// helper) — the shape-flag copy of the rebuild.
fn set_shape_flags(the_in_sh: &Shape, the_out_sh: &Shape) {
    // theOutSh.Modified(theInSh.Modified());
    topods_set_flag(
        the_out_sh,
        tshape_flags::MODIFIED,
        topods_flag(the_in_sh, tshape_flags::MODIFIED),
    );
    // theOutSh.Checked(theInSh.Checked());
    topods_set_flag(
        the_out_sh,
        tshape_flags::CHECKED,
        topods_flag(the_in_sh, tshape_flags::CHECKED),
    );
    // theOutSh.Orientable(theInSh.Orientable());
    topods_set_flag(
        the_out_sh,
        tshape_flags::ORIENTABLE,
        topods_flag(the_in_sh, tshape_flags::ORIENTABLE),
    );
    // theOutSh.Closed(theInSh.Closed());
    topods_set_flag(
        the_out_sh,
        tshape_flags::CLOSED,
        topods_flag(the_in_sh, tshape_flags::CLOSED),
    );
    // theOutSh.Infinite(theInSh.Infinite());
    topods_set_flag(
        the_out_sh,
        tshape_flags::INFINITE,
        topods_flag(the_in_sh, tshape_flags::INFINITE),
    );
    // theOutSh.Convex(theInSh.Convex());
    topods_set_flag(
        the_out_sh,
        tshape_flags::CONVEX,
        topods_flag(the_in_sh, tshape_flags::CONVEX),
    );
}

// ---------------------------------------------------------------------------
// TopoDS_Iterator / TopExp / TopLoc_Location / BRep_Tool re-hosts.
// ---------------------------------------------------------------------------

/// OCCT TopoDS_Iterator(S, theCumOri = false) (TopoDS_Iterator.cxx L39-86) —
/// the direct sub-shapes of S with the RAW stored orientations.  rcad stores
/// the children in per-kind slots, so the walk is the per-kind slot read —
/// the brep_algo::tool::sub_shapes carrier without the orientation compose
/// (the two walks must agree: TopExp_Explorer descends through sub_shapes,
/// so a child reached by the explorer must be bound in myMap by Put).
fn topods_iterator(the_s: &Shape) -> Vec<Shape> {
    match the_s.data.as_ref() {
        TShape::Vertex(_) => Vec::new(),
        TShape::Edge(ed) => vec![ed.first.clone(), ed.last.clone()],
        TShape::Wire(wd) => wd.edges.clone(),
        TShape::Face(fd) => {
            let mut out = Vec::new();
            if matches!(fd.outer_wire.data.as_ref(), TShape::Wire(_)) {
                out.push(fd.outer_wire.clone());
            }
            out.extend(fd.inner_wires.iter().cloned());
            out.extend(fd.internal_vertices.iter().cloned());
            out
        }
        TShape::Shell(sd) => sd.faces.clone(),
        TShape::Solid(sd) => {
            let mut out = Vec::new();
            out.extend(sd.shells.iter().cloned());
            out.extend(sd.internal_vertices.iter().cloned());
            out.extend(sd.internal_edges.iter().cloned());
            out
        }
        TShape::CompSolid(cs) => cs.clone(),
        TShape::Compound(cd) => cd.clone(),
    }
}

/// OCCT TopLoc_Location::Predivided(theOther) (TopLoc_Location.hxx L119-120):
/// the location `theOther.Inverted() * self` (the division BRep_Builder
/// applies as `L.Predivided(E.Location())` when the geometry location L is
/// stored relative to the shape location).  The rcad location carrier is the
/// u32 id of the owning transform table (0 = identity;
/// draft_modification.rs arch. diff. #4) and the modifier holds no table
/// (the no-pool Arc TShape carrier), so the id division is expressed against
/// the identity operand: the rcad modification implementors carry the
/// identity location (Draft_Modification.cxx L251/L283 `L.Identity()`), for
/// which `Predivided` is the identity composition and the id is carried
/// through unchanged.  A non-identity geometry location needs the owning
/// pool's `BRep::locations` table to multiply the transforms and can only be
/// composed by a pool-carrying caller.
fn top_loc_predivided(the_self: u32, _the_other: u32) -> u32 {
    the_self
}

/// OCCT BRep_Tool::Surface(F, L) (BRep_Tool.cxx L594-620) reduced to the
/// surface value (the `aLoc` out-parameter is the caller's concern).
fn brep_tool_surface_value(the_f: &Shape) -> Option<Surface3> {
    bat::brep_tool_surface(the_f)
}

/// OCCT BRep_Tool::Curve(E, L, f, l) (BRep_Tool.cxx L184-215) — the curve with
/// its location and range out-parameters.
fn brep_tool_curve_with_location(
    the_e: &Shape,
    the_loc: &mut u32,
    the_f: &mut f64,
    the_l: &mut f64,
) -> Option<Curve3> {
    let (f, l) = bat::brep_tool_range(the_e);
    *the_loc = the_e.location;
    *the_f = f;
    *the_l = l;
    bat::brep_tool_curve(the_e).map(|(c, _, _)| c)
}

/// OCCT `aSurf == anOtherSurf` (the Geom_Surface handle comparison,
/// BRepTools_Modifier.cxx L463) — the rcad surface identity comparison
/// (`surface_same`, the Geom_Surface handle stand-in) with the both-null
/// (equal handles) case of the OCCT comparison.
fn surface_handle_equal(a: &Option<Surface3>, b: &Option<Surface3>) -> bool {
    match (a, b) {
        (Some(x), Some(y)) => rcad_kernel::topods::surface_same(x, y),
        (None, None) => true,
        _ => false,
    }
}

/// OCCT BRep_Tool::SameParameter(E) (BRep_Tool.cxx L881-890) — the stored
/// BRep_TEdge flag.
fn brep_tool_same_parameter(the_e: &Shape) -> bool {
    match the_e.data.as_ref() {
        TShape::Edge(ed) => ed.same_parameter,
        _ => false,
    }
}

/// OCCT BRep_Tool::SameRange(E) (BRep_Tool.cxx L890-900) — the stored
/// BRep_TEdge flag.
fn brep_tool_same_range(the_e: &Shape) -> bool {
    match the_e.data.as_ref() {
        TShape::Edge(ed) => ed.same_range,
        _ => false,
    }
}

/// OCCT BRep_Tool::NaturalRestriction(F) (BRep_Tool.cxx L600-610) — the
/// stored BRep_TFace flag.
fn brep_tool_natural_restriction(the_f: &Shape) -> bool {
    match the_f.data.as_ref() {
        TShape::Face(fd) => fd.natural_restriction,
        _ => false,
    }
}

// ---------------------------------------------------------------------------
// BRep_Builder re-hosts (the in-place shared-TShape edits — arch. diff. #3).
// ---------------------------------------------------------------------------

/// The TShape payload of a shape as a mutable reference, in place.  OCCT
/// BRep_Builder edits the TShape through the shared handle, so the edit must
/// not fork the Arc (the brep_sweep::tool_rehost::tshape_mut rationale).
fn tshape_mut_inplace(the_s: &Shape) -> &mut TShape {
    // SAFETY: single-threaded; the caller holds &mut BRepToolsModifier (or is
    // inside the sequential build of one shape) and no other &TShape for this
    // Arc is alive across the call.
    let ptr = Arc::as_ptr(&the_s.data) as *mut TShape;
    unsafe { &mut *ptr }
}

/// OCCT TopoDS_Shape flag write (TopoDS_Shape.hxx L216-260: Modified /
/// Checked / Orientable / Closed / Infinite / Convex) — the bit layout of the
/// rcad tshape_flags mirrors the OCCT TopoDS_TShape::BitLayout (topods.rs
/// L83-94).  The CompSolid / Compound TShapes carry no flag word in rcad (the
/// OCCT TopoDS_TSolid/Compound flag bits are not modelled).
fn topods_set_flag(the_s: &Shape, the_flag: u16, the_on: bool) {
    let flags = match tshape_mut_inplace(the_s) {
        TShape::Vertex(v) => &mut v.flags,
        TShape::Edge(e) => &mut e.flags,
        TShape::Wire(w) => &mut w.flags,
        TShape::Face(f) => &mut f.flags,
        TShape::Shell(sh) => &mut sh.flags,
        TShape::Solid(so) => &mut so.flags,
        TShape::CompSolid(_) | TShape::Compound(_) => return,
    };
    if the_on {
        *flags |= the_flag;
    } else {
        *flags &= !the_flag;
    }
}

/// OCCT TopoDS_Shape flag read (the complements of [`topods_set_flag`]).
fn topods_flag(the_s: &Shape, the_flag: u16) -> bool {
    let flags = match the_s.data.as_ref() {
        TShape::Vertex(v) => v.flags,
        TShape::Edge(e) => e.flags,
        TShape::Wire(w) => w.flags,
        TShape::Face(f) => f.flags,
        TShape::Shell(sh) => sh.flags,
        TShape::Solid(so) => so.flags,
        TShape::CompSolid(_) | TShape::Compound(_) => 0,
    };
    flags & the_flag != 0
}

/// OCCT BRep_Builder::MakeFace(F, S, L, Tol) (BRep_Builder.cxx L547-563) — a
/// new face TShape carrying the surface, its location and the tolerance.
fn make_face(the_surface: Option<Surface3>, the_loc: u32, the_tol: f64) -> Shape {
    Shape {
        data: Arc::new(TShape::Face(TFaceData {
            my_shapes: Vec::new(),
            flags: tshape_flags::DEFAULT,
            surface: the_surface,
            surface_location: the_loc,
            outer_wire: Shape::null(),
            inner_wires: Vec::new(),
            sample_point: None,
            uv_domain: None,
            internal_vertices: Vec::new(),
            tolerance: the_tol,
            natural_restriction: false,
        })),
        index: usize::MAX,
        location: 0,
        orientation: Orientation::Forward,
    }
}

/// OCCT BRep_Builder::MakeEdge(E, C, L, Tol) (BRep_Builder.cxx L780-800) — a
/// new edge TShape carrying the 3D curve and the tolerance.  Architecture:
/// OCCT stores the curve through `new BRep_Curve3D(C, L)`; the rcad carrier
/// keeps the curve by value on TEdgeData::curve (the value form has no
/// per-representation location table), so `the_loc` is recorded by the caller
/// as the edge wrapper location (`result.Location(S.Location(), false)`).
fn make_edge(the_curve: Option<Curve3>, _the_loc: u32, the_tol: f64) -> Shape {
    Shape {
        data: Arc::new(TShape::Edge(TEdgeData {
            my_shapes: Vec::new(),
            flags: tshape_flags::DEFAULT,
            curve: the_curve,
            first: Shape::null(),
            last: Shape::null(),
            range: [0.0, 0.0],
            degenerated: false,
            pcurves: indexmap::IndexMap::new(),
            representations: Vec::new(),
            vertex_params: HashMap::new(),
            tolerance: the_tol,
            same_parameter: false,
            same_range: false,
        })),
        index: usize::MAX,
        location: 0,
        orientation: Orientation::Forward,
    }
}

/// OCCT BRep_Builder::MakeVertex(V, P, Tol) (BRep_Builder.cxx L1150-1170) —
/// a new vertex TShape carrying the point and the tolerance.
fn make_vertex(the_p: glam::DVec3, the_tol: f64) -> Shape {
    Shape {
        data: Arc::new(TShape::Vertex(rcad_kernel::topods::TVertexData {
            my_shapes: Vec::new(),
            flags: tshape_flags::DEFAULT,
            point: the_p,
            tolerance: the_tol,
            points: Vec::new(),
        })),
        index: usize::MAX,
        location: 0,
        orientation: Orientation::Forward,
    }
}

/// OCCT BRep_Builder::Add(S, C) (TopoDS_Builder::Add, TopoDS_Builder.cxx
/// L37-100) — the per-kind container append, in place on the shared TShape.
/// The face outer-wire slot test is the child type test (the bat::sub_shapes
/// convention: the null placeholder is a Vertex, never a Wire, and the
/// pool-free shapes carry index == usize::MAX so `Shape::is_null` cannot
/// discriminate them).
fn builder_add(the_s: &Shape, the_c: &Shape) {
    match tshape_mut_inplace(the_s) {
        TShape::Vertex(_) => {}
        TShape::Edge(ed) => {
            // BRep_Builder::Add(E, V): the vertex slots mirror the front /
            // back of the OCCT myShapes list (TopExp::FirstVertex /
            // LastVertex read them back).
            if the_c.orientation == Orientation::Reversed {
                ed.last = the_c.clone();
            } else {
                ed.first = the_c.clone();
            }
            ed.my_shapes.push(the_c.clone());
        }
        TShape::Wire(wd) => {
            wd.edges.push(the_c.clone());
            wd.my_shapes.push(the_c.clone());
        }
        TShape::Face(fd) => {
            if the_c.shape_type() == ShapeType::Vertex {
                fd.internal_vertices.push(the_c.clone());
            } else if matches!(fd.outer_wire.data.as_ref(), TShape::Wire(_)) {
                fd.inner_wires.push(the_c.clone());
            } else {
                fd.outer_wire = the_c.clone();
            }
            fd.my_shapes.push(the_c.clone());
        }
        TShape::Shell(sd) => {
            sd.faces.push(the_c.clone());
            sd.my_shapes.push(the_c.clone());
        }
        TShape::Solid(sd) => {
            match the_c.shape_type() {
                ShapeType::Vertex => sd.internal_vertices.push(the_c.clone()),
                ShapeType::Edge => sd.internal_edges.push(the_c.clone()),
                _ => sd.shells.push(the_c.clone()),
            }
            sd.my_shapes.push(the_c.clone());
        }
        TShape::CompSolid(cs) => cs.push(the_c.clone()),
        TShape::Compound(cd) => cd.push(the_c.clone()),
    }
}

/// OCCT BRep_Builder::Degenerated(E, theFlag) (BRep_Builder.cxx L1015-1030).
fn builder_degenerated(the_e: &Shape, the_flag: bool) {
    if let TShape::Edge(ed) = tshape_mut_inplace(the_e) {
        ed.degenerated = the_flag;
    }
}

/// OCCT BRep_Builder::UpdateEdge(E, theTol) (BRep_Builder.cxx L700-720) — the
/// tolerance update (BRep_TEdge::UpdateTolerance keeps the max).
fn builder_update_edge_tolerance(the_e: &Shape, the_tol: f64) {
    if let TShape::Edge(ed) = tshape_mut_inplace(the_e) {
        ed.tolerance = ed.tolerance.max(the_tol);
    }
}

/// OCCT BRep_Builder::SameParameter(E, theFlag) (BRep_Builder.cxx L960-980).
fn builder_same_parameter(the_e: &Shape, the_flag: bool) {
    if let TShape::Edge(ed) = tshape_mut_inplace(the_e) {
        ed.same_parameter = the_flag;
    }
}

/// OCCT BRep_Builder::SameRange(E, theFlag) (BRep_Builder.cxx L985-1005).
fn builder_same_range(the_e: &Shape, the_flag: bool) {
    if let TShape::Edge(ed) = tshape_mut_inplace(the_e) {
        ed.same_range = the_flag;
    }
}

/// OCCT BRep_Builder::NaturalRestriction(F, theFlag) (BRep_Builder.cxx
/// L594-605).
fn builder_natural_restriction(the_f: &Shape, the_flag: bool) {
    if let TShape::Face(fd) = tshape_mut_inplace(the_f) {
        fd.natural_restriction = the_flag;
    }
}

/// OCCT BRep_Builder::Range(E, First, Last) (BRep_Builder.cxx L120-160) — the
/// 3D range of the edge.
fn builder_range_edge(the_e: &Shape, the_first: f64, the_last: f64) {
    if let TShape::Edge(ed) = tshape_mut_inplace(the_e) {
        ed.range = [the_first, the_last];
    }
}

/// OCCT BRep_Builder::UpdateEdge(E, C2d, S, L, Tol) (BRep_Builder.cxx
/// L655-672) -> UpdateCurves(lcr, C, S, L) (BRep_Builder.cxx L104-172): the
/// pcurve of the edge on the face (a BRep_CurveOnSurface representation) plus
/// the tolerance update.  The representation range is the 3D curve range when
/// the edge carries one, otherwise the pcurve's own range (BRep_Builder.cxx
/// L152-168).
fn builder_update_edge_pcurve(the_e: &Shape, the_c2d: &Curve2d, the_f: &Shape, the_tol: f64) {
    let key = bat::shape_key(the_f);
    let a_range = pcurve_range_of(the_e, the_c2d);
    if let TShape::Edge(ed) = tshape_mut_inplace(the_e) {
        ed.pcurves.insert(key, (the_c2d.clone(), a_range[0], a_range[1]));
        // OCCT L169: lcr.Append(COS) — the CurveOnSurface representation row
        // for the BRep_Tool::CurveOnSurface / IsClosed readers.
        ed.representations
            .push(rcad_kernel::topods::CurveRepresentation::CurveOnSurface {
                face: key,
                pcurve: the_c2d.clone(),
                range: a_range,
            });
        ed.tolerance = ed.tolerance.max(the_tol);
    }
}

/// OCCT BRep_Builder::UpdateEdge(E, C1, C2, S, L, Tol) (BRep_Builder.cxx
/// L679-697) — the two-pcurve (seam) form: a BRep_CurveOnClosedSurface
/// representation with the same range rule as [`builder_update_edge_pcurve`].
fn builder_update_edge_pcurve2(
    the_e: &Shape,
    the_c2d1: &Curve2d,
    the_c2d2: &Curve2d,
    the_f: &Shape,
    the_tol: f64,
) {
    let key = bat::shape_key(the_f);
    let a_range = pcurve_range_of(the_e, the_c2d1);
    if let TShape::Edge(ed) = tshape_mut_inplace(the_e) {
        ed.pcurves
            .insert(key, (the_c2d1.clone(), a_range[0], a_range[1]));
        ed.representations.push(
            rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                face: key,
                pcurve1: the_c2d1.clone(),
                pcurve2: the_c2d2.clone(),
                range: a_range,
            },
        );
        ed.tolerance = ed.tolerance.max(the_tol);
    }
}

/// The BRep_CurveOnSurface range rule of the static UpdateCurves
/// (BRep_Builder.cxx L112-168): the 3D curve range of the edge when it
/// carries one (`f` stays `-Precision::Infinite()` when no BRep_Curve3D
/// representation is found), otherwise the pcurve's own
/// FirstParameter/LastParameter (BRep_CurveOnSurface::Range).
fn pcurve_range_of(the_e: &Shape, the_c2d: &Curve2d) -> [f64; 2] {
    if let TShape::Edge(ed) = the_e.data.as_ref() {
        if ed.curve.is_some() {
            return ed.range;
        }
    }
    use rcad_kernel::geom::Curve2dEval;
    let d = the_c2d.default_domain();
    [d[0], d[1]]
}

/// OCCT BRep_Builder::UpdateVertex(V, Par, E, Tol) (BRep_Builder.cxx
/// L1220-1310) — the vertex parameter on the edge plus the vertex tolerance.
///
/// Architecture: OCCT walks the edge vertices to derive the parameter slot
/// (FORWARD -> the GCurve First, REVERSED -> the Last, otherwise the vertex
/// point representations); the rcad TEdgeData::vertex_params map is the
/// single vertex-parameter carrier (the kernel BRepBuilder::set_vertex_param
/// translation), so the parameter is stored under the vertex key.  The
/// tolerance goes on the VERTEX (OCCT L1304 `TV->UpdateTolerance(Tol)`), not
/// on the edge.
fn builder_update_vertex_param(the_v: &Shape, the_par: f64, the_e: &Shape, the_tol: f64) {
    // OCCT L1225-1228: if (Precision::IsPositiveInfinite(Par) ||
    // Precision::IsNegativeInfinite(Par))
    // throw Standard_DomainError("BRep_Builder::Infinite parameter");
    if is_positive_infinite_value(the_par) || is_negative_infinite_value(the_par) {
        panic!("Standard_DomainError: BRep_Builder::Infinite parameter");
    }
    // OCCT L1280-1305: the GCurve parameter writes (First / Last) + the
    // tolerance update.
    if let TShape::Vertex(vd) = tshape_mut_inplace(the_v) {
        vd.tolerance = vd.tolerance.max(the_tol);
    }
    if let TShape::Edge(ed) = tshape_mut_inplace(the_e) {
        ed.vertex_params.insert(the_v.ptr_id(), the_par);
    }
}

/// OCCT BRep_Builder::Continuity(E, F1, F2, C) (BRep_Builder.cxx L1012-1043) —
/// the regularity of the edge between the two faces (the
/// BRep_CurveOn2Surfaces representation; the rcad re-host is the kernel
/// BRepBuilder::continuity reduced to the shape handles).
fn builder_continuity(the_e: &Shape, the_f1: &Shape, the_f2: &Shape, the_c: GeomAbsShape) {
    // OCCT L1017-1019: S1 = Surface(F1, l1); S2 = Surface(F2, l2).
    let (Some(s1), Some(s2)) = (
        brep_tool_surface_value(the_f1),
        brep_tool_surface_value(the_f2),
    ) else {
        // The OCCT null-surface case stores no regularity.
        return;
    };
    // OCCT L1037-1038: l1 = L1.Predivided(E.Location());
    let l1 = top_loc_predivided(the_f1.location, the_e.location);
    let l2 = top_loc_predivided(the_f2.location, the_e.location);
    // OCCT L1040: UpdateCurves(TE->ChangeCurves(), S1, S2, l1, l2, C)
    // (the static UpdateCurves, BRep_Builder.cxx L376-405).
    if let TShape::Edge(ed) = tshape_mut_inplace(the_e) {
        for cr in ed.representations.iter_mut() {
            if cr.is_regularity_on(&s1, &s2, l1, l2) {
                if let rcad_kernel::topods::CurveRepresentation::CurveOn2Surfaces {
                    continuity,
                    ..
                } = cr
                {
                    *continuity = the_c;
                }
                return;
            }
        }
        ed.representations.push(
            rcad_kernel::topods::CurveRepresentation::CurveOn2Surfaces {
                surface1: s1,
                surface2: s2,
                location1: l1,
                location2: l2,
                continuity: the_c,
            },
        );
    }
}
