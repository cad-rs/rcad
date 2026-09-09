//! 1:1 translation of OCCT `ShapeUpgrade_UnifySameDomain`
//! (`TKShHealing/ShapeUpgrade/ShapeUpgrade_UnifySameDomain.hxx` + `.cxx`,
//! 4,687 LOC — W1-6 docket row `B-1:1`).
//!
//! The class unifies faces and edges of the shape which lie on the same
//! geometry.  Faces/edges are considered "same-domain" if a group of
//! neighbouring faces/edges lie on coincident surfaces/curves.
//!
//! Function-count equation (OCCT ShapeUpgrade_UnifySameDomain.cxx = rcad):
//! 32 file statics (`IsOnSingularity`, `IsUiso`, `IsLinear`, `SplitWire`,
//! `TrueValueOfOffset`, `GetFaceFromSeq`, `UpdateBoundaries`, `TryMakeLine`,
//! `RemoveEdgeFromMap`, `ComputeMinEdgeSize`, `FindCoordBounds`,
//! `getCurveParams`, `RelocatePCurvesToNewUorigin`, `InsertWiresIntoFaces`,
//! `FindCommonFace`, `FindClosestPoints`, `ReconstructMissedSeam`,
//! `SameSurf`, `TransformPCurves`, `AddPCurves`, `AddOrdinaryEdges`,
//! `getCylinder`, `ClearRts`, `GetNormalToSurface`, `IsSameDomain`,
//! `UpdateMapOfShapes`, `GlueEdgesWith3DCurves`, `IsMergingPossible`,
//! `GetLineEdgePoints`, `CheckSharedVertices`, `SetFixWireModes`,
//! `isSameSets`) = 32 rcad statics;
//! 17 class members (2 ctors, `Initialize`, `AllowInternalEdges`,
//! `SetSafeInputMode`, `KeepShape`, `KeepShapes`, `UnifyFaces`,
//! `IntUnifyFaces`, `UnifyEdges`, `Build`, `FillHistory`, `MergeEdges`,
//! `MergeSeq`, `MergeSubSeq`, `UnionPCurves`, static `generateSubSeq`,
//! struct `SubSequenceOfEdges`) = 17 rcad items.
//! Header inline accessors (`SetLinearTolerance`, `SetAngularTolerance`,
//! `Shape`, `History`) land here as direct hxx translations.
//!
//! Segment schedule (translation units follow the .cxx comment blocks):
//! - `statics_a` — cxx L93-533
//! - `statics_b` — cxx L536-1750
//! - `union_pcurves` — cxx L1754-2157
//! - `merge_sub_seq` — cxx L2164-2506 (+ `IsMergingPossible` L2508-2638,
//!   `GetLineEdgePoints` L2642-2675)
//! - `merge_edges` — cxx L2677-2948
//! - `unify_faces` — cxx L3047-3181
//! - `int_unify_faces` — cxx L3185-4320
//! - `unify_edges` — cxx L4324-4450
//! - `fill_history` — cxx L4475-4558
//! - `split_wire` — cxx L4560-4687 (the deferred `SplitWire` body)

mod fill_history;
mod gap_deps;
mod int_unify_faces;
mod merge_edges;
mod merge_sub_seq;
mod split_wire;
mod statics_a;
mod statics_b;
mod topexp;
mod unify_edges;
mod unify_faces;
mod union_pcurves;

use std::collections::HashMap;

use crate::bop::history::BRepToolsHistory;
use indexmap::IndexMap;
use rcad_kernel::geom::Plane;
use rcad_kernel::precision::{ANGULAR, CONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, Orientation, ShapeType};

use crate::shhealing::shape_build::reshape::ShapeBuildReShape;

/// Shape-key of the TopTools_ShapeMapHasher identity (TShape + Location).
/// Architecture note: OCCT compares the location by value; rcad keys by the
/// location-table index, stable within one BRep (the reshape.rs /
/// loc_ope_glued_shape.rs precedent).
#[inline]
pub(crate) fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT `NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>` /
/// `NCollection_IndexedMap` — the insertion-ordered rcad set; `Add`
/// returning "was not yet present" is expressed by
/// [`map_add`].
pub(crate) type MapOfShape = IndexMap<(u64, u32), Shape>;

/// OCCT `NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
/// TopTools_ShapeMapHasher>`.
pub(crate) type IndexedDataMapOfShapeListOfShape = IndexMap<(u64, u32), (Shape, Vec<Shape>)>;

/// OCCT `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, ...>` (myFaceNewFace).
pub(crate) type DataMapOfShapeShape = HashMap<(u64, u32), Shape>;

/// OCCT `DataMapOfFacePlane` (the hxx typedef) — `Geom_Plane` handle maps to
/// the rcad plane value.
pub(crate) type DataMapOfFacePlane = HashMap<(u64, u32), Plane>;

/// OCCT `DataMapOfShapeMapOfShape` (the hxx typedef).
pub(crate) type DataMapOfShapeMapOfShape = HashMap<(u64, u32), (Shape, MapOfShape)>;

/// OCCT `NCollection_Map::Add` — returns true when the shape was added
/// (not yet present), false when it was already contained.
pub(crate) fn map_add(m: &mut MapOfShape, s: &Shape) -> bool {
    m.insert(shape_key(s), s.clone()).is_none()
}

/// OCCT TopAbs::Reverse (TopAbs.hxx L80-91) on an orientation.
pub(crate) fn occt_reverse(o: Orientation) -> Orientation {
    o.compose(Orientation::Reversed)
}

/// OCCT `SubSequenceOfEdges` (cxx L2677-2681).
#[derive(Clone)]
pub(crate) struct SubSequenceOfEdges {
    /// OCCT SeqsEdges.
    pub seqs_edges: Vec<Shape>,
    /// OCCT UnionEdges.
    pub union_edges: Shape,
}

/// OCCT `ShapeUpgrade_UnifySameDomain` (the hxx class, L65-225).
pub struct ShapeUpgradeUnifySameDomain {
    /// OCCT myInitShape.
    pub(crate) my_init_shape: Shape,
    /// OCCT myLinTol.
    pub(crate) my_lin_tol: f64,
    /// OCCT myAngTol.
    pub(crate) my_ang_tol: f64,
    /// OCCT myUnifyFaces.
    pub(crate) my_unify_faces: bool,
    /// OCCT myUnifyEdges.
    pub(crate) my_unify_edges: bool,
    /// OCCT myConcatBSplines.
    pub(crate) my_concat_bsplines: bool,
    /// OCCT myAllowInternal.
    pub(crate) my_allow_internal: bool,
    /// OCCT mySafeInputMode.
    pub(crate) my_safe_input_mode: bool,
    /// OCCT myShape.
    pub(crate) my_shape: Shape,
    /// OCCT myContext (handle -> owned reshape context).
    pub(crate) my_context: ShapeBuildReShape,
    /// OCCT myKeepShapes.
    pub(crate) my_keep_shapes: MapOfShape,
    /// OCCT myFacePlaneMap.
    pub(crate) my_face_plane_map: DataMapOfFacePlane,
    /// OCCT myEFmap.
    pub(crate) my_ef_map: IndexedDataMapOfShapeListOfShape,
    /// OCCT myFaceNewFace.
    pub(crate) my_face_new_face: DataMapOfShapeShape,
    /// OCCT myHistory (the handle may be null in OCCT -> Option).
    pub(crate) my_history: Option<BRepToolsHistory>,
}

impl Default for ShapeUpgradeUnifySameDomain {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeUpgradeUnifySameDomain {
    // OCCT ShapeUpgrade_UnifySameDomain.cxx L2952-2963: the empty
    // constructor.
    pub fn new() -> Self {
        ShapeUpgradeUnifySameDomain {
            my_init_shape: Shape::null(),
            my_lin_tol: CONFUSION,
            my_ang_tol: ANGULAR,
            my_unify_faces: true,
            my_unify_edges: true,
            my_concat_bsplines: false,
            my_allow_internal: false,
            my_safe_input_mode: true,
            my_shape: Shape::null(),
            my_context: ShapeBuildReShape::new(),
            my_keep_shapes: MapOfShape::new(),
            my_face_plane_map: DataMapOfFacePlane::new(),
            my_ef_map: IndexedDataMapOfShapeListOfShape::new(),
            my_face_new_face: DataMapOfShapeShape::new(),
            my_history: Some(BRepToolsHistory::new()),
        }
    }

    // OCCT ShapeUpgrade_UnifySameDomain.cxx L2967-2983: the constructor
    // defining the input shape and the flags (no unification performed).
    pub fn with_shape(
        a_shape: &Shape,
        unify_edges: bool,
        unify_faces: bool,
        concat_bsplines: bool,
    ) -> Self {
        ShapeUpgradeUnifySameDomain {
            my_init_shape: a_shape.clone(),
            my_lin_tol: CONFUSION,
            my_ang_tol: ANGULAR,
            my_unify_faces: unify_faces,
            my_unify_edges: unify_edges,
            my_concat_bsplines: concat_bsplines,
            my_allow_internal: false,
            my_safe_input_mode: true,
            my_shape: a_shape.clone(),
            my_context: ShapeBuildReShape::new(),
            my_keep_shapes: MapOfShape::new(),
            my_face_plane_map: DataMapOfFacePlane::new(),
            my_ef_map: IndexedDataMapOfShapeListOfShape::new(),
            my_face_new_face: DataMapOfShapeShape::new(),
            my_history: Some(BRepToolsHistory::new()),
        }
    }

    // OCCT ShapeUpgrade_UnifySameDomain.cxx L2987-3004: Initialize.
    pub fn initialize(
        &mut self,
        brep: &mut BRep,
        a_shape: &Shape,
        unify_edges: bool,
        unify_faces: bool,
        concat_bsplines: bool,
    ) {
        self.my_init_shape = a_shape.clone();
        self.my_shape = a_shape.clone();
        self.my_unify_edges = unify_edges;
        self.my_unify_faces = unify_faces;
        self.my_concat_bsplines = concat_bsplines;

        self.my_context.clear();
        self.my_keep_shapes.clear();
        self.my_face_plane_map.clear();
        self.my_ef_map.clear();
        self.my_face_new_face.clear();
        // OCCT L3003: myHistory->Clear() — the rcad BRepToolsHistory has no
        // Clear yet (NEEDED EDIT IN bop/history.rs: add Clear, OCCT
        // BRepTools_History.cxx); re-creation is semantically identical for
        // the handle held from the constructor.
        self.my_history = Some(BRepToolsHistory::new());
    }

    // OCCT ShapeUpgrade_UnifySameDomain.cxx L3008-3011: AllowInternalEdges.
    pub fn allow_internal_edges(&mut self, the_value: bool) {
        self.my_allow_internal = the_value;
    }

    // OCCT hxx L118: SetSafeInputMode (cxx L3015-3018).
    pub fn set_safe_input_mode(&mut self, the_value: bool) {
        self.my_safe_input_mode = the_value;
    }

    // OCCT hxx L121-122: SetLinearTolerance (inline in the hxx).
    pub fn set_linear_tolerance(&mut self, the_value: f64) {
        self.my_lin_tol = the_value;
    }

    // OCCT hxx L126-129: SetAngularTolerance (inline in the hxx).
    pub fn set_angular_tolerance(&mut self, the_value: f64) {
        self.my_ang_tol = if the_value < ANGULAR {
            ANGULAR
        } else {
            the_value
        };
    }

    // OCCT hxx L135: Shape() — the resulting shape.
    pub fn shape(&self) -> &Shape {
        &self.my_shape
    }

    // OCCT hxx L138: History() — the collected history (None = the null
    // handle the user may have installed).
    pub fn history(&self) -> Option<&BRepToolsHistory> {
        self.my_history.as_ref()
    }

    // OCCT cxx L3022-3028: KeepShape.
    pub fn keep_shape(&mut self, the_shape: &Shape) {
        if the_shape.shape_type() == ShapeType::Edge || the_shape.shape_type() == ShapeType::Vertex
        {
            map_add(&mut self.my_keep_shapes, the_shape);
        }
    }

    // OCCT cxx L3032-3043: KeepShapes.
    pub fn keep_shapes(&mut self, the_shapes: &MapOfShape) {
        for (_, v) in the_shapes.iter() {
            if v.shape_type() == ShapeType::Edge || v.shape_type() == ShapeType::Vertex {
                map_add(&mut self.my_keep_shapes, v);
            }
        }
    }

    // OCCT ShapeUpgrade_UnifySameDomain.cxx L4454-4469: Build.
    pub fn build(&mut self, brep: &mut BRep) {
        // OCCT L4456: TopExp::MapShapesAndAncestors(myInitShape, EDGE, FACE,
        // myEFmap).
        self.my_ef_map = IndexedDataMapOfShapeListOfShape::new();
        topexp::map_shapes_and_ancestors(
            brep,
            &self.my_init_shape,
            ShapeType::Edge,
            ShapeType::Face,
            &mut self.my_ef_map,
        );

        if self.my_unify_faces {
            self.unify_faces(brep);
        }
        if self.my_unify_edges {
            self.unify_edges(brep);
        }

        // OCCT L4467-4468: Fill the history of modifications.
        self.fill_history(brep);
    }
}
