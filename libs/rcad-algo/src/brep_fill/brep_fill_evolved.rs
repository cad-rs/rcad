//! OCCT BRepFill_Evolved — 1:1 translation (part 1: class, constructors,
//! Perform, PrivatePerform, shared data carriers).
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepFill/
//!         BRepFill_Evolved.hxx (L40-199) +
//!         BRepFill_Evolved.cxx (L99-544, L548-618).
//!
//! The translation is split across sibling files (the <2000-line rule):
//! - this file: hxx class data + L99-103 (BRepFill_Confusion) + L172-182
//!   (EdgeVertices) + the constructors + Perform x2 + PrivatePerform +
//!   SetWork;
//! - brep_fill_evolved_b.rs: ElementaryPerform (L657-1263);
//! - brep_fill_evolved_c.rs: PlanarPerform .. ContinuityOnOffsetEdge
//!   (L1267-2544);
//! - brep_fill_evolved_d.rs: the remaining file statics.
//!
//! Architecture differences (referenced from the affected functions):
//! 1. `NCollection_DataMap<TopoDS_Shape, Item, TopTools_ShapeMapHasher>`
//!    maps to [`DataMapOfShapeItem`] (a Vec of (key, item) pairs keyed by
//!    the shape identity — the offset_ancestors.rs reduction; insertion
//!    order preserved, OCCT operator()/IsBound/Bind map to find/is_bound/
//!    bind).
//! 2. `TopLoc_Location` maps to a u32 index into the per-object location
//!    table [`LocationTable`] (index 0 = identity; the OCCT global TopLoc
//!    Datum table stand-in).  `Shape::Moved` / `Shape::Location(L)` map to
//!    [`location_shape_moved`] / [`location_shape_set`].
//! 3. The BRep_Builder operations map to the pool-free
//!    `brep_algo::tool` re-hosts (the shapes carry their TShapes as
//!    Arc handles — the OCCT global TShape registry stand-in).  The
//!    regularity / tolerance calls without a kernel table
//!    (`B.Continuity`, `BRepLib::UpdateTolerances`, `BRepLib::
//!    SameParameter`, `B.SameRange`) map to documented no-op stubs (the
//!    brep_lib.rs precedent).
//! 4. GAP carriers (dependencies outside the module, plan §0.6):
//!    - `BRepMAT2d_BisectingLocus` / `BRepMAT2d_Explorer` /
//!      `BRepMAT2d_LinkTopoBilo` (TKTopAlgo/BRepMAT2d — the wrapper
//!      classes are not translated; the underlying MAT / MAT2d / Bisector
//!      stacks are the real topalgo translations);
//!    - `BRepSweep_Prism` / `BRepSweep_Revol` (TKTopAlgo/BRepSweep);
//!    - `BRepTools_Quilt` (TKTopAlgo/BRepTools);
//!    - `Geom2dAPI_ExtremaCurveCurve` (TKGeomAlgo/Geom2dAPI — the
//!      Extrema_ExtCC2d engine is a GAP, the mat2d_mini_path.rs precedent);
//!    - `BndLib_Add2dCurve` (TKMath/BndLib);
//!    - `BRepLProp::Continuity` (TKTopAlgo/BRepLProp);
//!    - `BRepLib_FindSurface` / `BRepLib_MakeFace(Wire, OnlyPlane)` (the
//!      brep_fill_axe.rs carriers, reused);
//!    - `BRepFill_OffsetWire::PerformWithBiLo` panics inside the
//!      offset_wire.rs translation (its own GAP annotation);
//!    - the BRepFill_TrimSurfaceTool::Project MultiLine/ApproxSeewing chain
//!      (carried inside the trim_surface_tool.rs translation).
//!
//! first consumer: LocOpe_Pipe / LocOpe_RevolutionForm (the 3c batch) and
//! BRepOffsetAPI_MakeEvolved (2e).

use std::sync::Arc;

use glam::{DAffine3, DVec3};

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType, TShape};

use crate::brep_fill::brep_fill_pipe::top_exp_vertices;
use crate::brep_fill::brep_fill_trim_edge_tool::GeomAbsJoinType;
use crate::brep_fill::generator::{shape_key, shape_oriented, shape_reversed, ShapeKey};
use crate::brep_fill::offset_wire::MatSide;
use crate::topalgo::bisector::bisector_bisec::BisectorBisec;
use crate::topalgo::mat::{HandleMatArc, HandleMatBasicElt, HandleMatGraph, HandleMatNode};

// ---------------------------------------------------------------------------
// Shared data-map reductions (architecture difference #1)
// ---------------------------------------------------------------------------

/// OCCT NCollection_DataMap<TopoDS_Shape, Item, TopTools_ShapeHasher> —
/// the Vec-of-pairs reduction keyed by the shape identity (TShape pointer;
/// the TopTools_ShapeMapHasher location component is folded into the key
/// shape itself, the generator::shape_key precedent).
#[derive(Debug, Clone)]
pub(super) struct DataMapOfShapeItem<T> {
    entries: Vec<(ShapeKey, Shape, T)>,
}

impl<T> Default for DataMapOfShapeItem<T> {
    fn default() -> Self {
        DataMapOfShapeItem {
            entries: Vec::new(),
        }
    }
}

impl<T> DataMapOfShapeItem<T> {
    pub(super) fn new() -> Self {
        Self::default()
    }

    /// OCCT Clear().
    pub(super) fn clear(&mut self) {
        self.entries.clear();
    }

    /// OCCT IsBound(K).
    pub(super) fn is_bound(&self, k: &Shape) -> bool {
        let key = shape_key(k);
        self.entries.iter().any(|(kk, _, _)| *kk == key)
    }

    /// OCCT Bind(K, Item) — replace or insert.
    pub(super) fn bind(&mut self, k: &Shape, item: T) {
        let key = shape_key(k);
        for (kk, _, it) in self.entries.iter_mut() {
            if *kk == key {
                *it = item;
                return;
            }
        }
        self.entries.push((key, k.clone(), item));
    }

    /// OCCT operator()(K) — a missing key raises Standard_NoSuchObject.
    pub(super) fn find(&self, k: &Shape) -> &T {
        let key = shape_key(k);
        self.entries
            .iter()
            .find(|(kk, _, _)| *kk == key)
            .map(|(_, _, it)| it)
            .expect("Standard_NoSuchObject: NCollection_DataMap::operator()")
    }

    /// OCCT operator()(K) mutable.
    pub(super) fn find_mut(&mut self, k: &Shape) -> &mut T {
        let key = shape_key(k);
        self.entries
            .iter_mut()
            .find(|(kk, _, _)| *kk == key)
            .map(|(_, _, it)| it)
            .expect("Standard_NoSuchObject: NCollection_DataMap::operator()")
    }

    /// Iteration (OCCT Initialize/More/Next over the DataMap).
    pub(super) fn iter(&self) -> impl Iterator<Item = (&Shape, &T)> {
        self.entries.iter().map(|(_, ks, it)| (ks, it))
    }
}

/// OCCT NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>.
pub(super) type DataMapOfShapeListOfShape = DataMapOfShapeItem<Vec<Shape>>;

/// OCCT NCollection_DataMap<TopoDS_Shape, NCollection_DataMap<TopoDS_Shape,
/// NCollection_List<TopoDS_Shape>>> — the myMap type of BRepFill_Evolved.
pub(super) type EvolutionMap = DataMapOfShapeItem<DataMapOfShapeListOfShape>;

// ---------------------------------------------------------------------------
// The per-object TopLoc_Location table (architecture difference #2)
// ---------------------------------------------------------------------------

/// OCCT TopLoc_Identity — the identity location index.
pub(super) const LOCATION_IDENTITY: u32 = 0;

/// OCCT TopLoc_Location table — index 0 is the identity (the BRep
/// locations encoding; the per-object stand-in of the OCCT global table).
#[derive(Debug, Default, Clone)]
pub(super) struct LocationTable {
    transforms: Vec<DAffine3>,
}

impl LocationTable {
    /// OCCT TopLoc_Location(Trsf) — intern the transform (identity folds
    /// to index 0).
    pub(super) fn intern(&mut self, t: DAffine3) -> u32 {
        if t == DAffine3::IDENTITY {
            return LOCATION_IDENTITY;
        }
        if let Some(i) = self.transforms.iter().position(|m| *m == t) {
            return (i + 1) as u32;
        }
        self.transforms.push(t);
        self.transforms.len() as u32
    }

    /// OCCT TopLoc_Location::Transformation().
    pub(super) fn transform(&self, l: u32) -> DAffine3 {
        if l == LOCATION_IDENTITY {
            DAffine3::IDENTITY
        } else {
            self.transforms
                .get((l - 1) as usize)
                .copied()
                .unwrap_or(DAffine3::IDENTITY)
        }
    }

    /// OCCT TopLoc_Location::Inverted().
    pub(super) fn inverted(&mut self, l: u32) -> u32 {
        self.intern(self.transform(l).inverse())
    }

    /// OCCT TopLoc_Location::Multiplied(L).
    fn multiplied(&mut self, a: u32, b: u32) -> u32 {
        self.intern(self.transform(a) * self.transform(b))
    }
}

/// OCCT TopoDS_Shape::Moved(L) — a copy with L * Location composed on the
/// per-object table.
pub(super) fn location_shape_moved(table: &mut LocationTable, s: &Shape, l: u32) -> Shape {
    let composed = table.multiplied(l, s.location);
    s.clone().with_location(composed)
}

/// OCCT TopoDS_Shape::Location(L) — a copy with the location replaced.
pub(super) fn location_shape_set(_table: &mut LocationTable, s: &Shape, l: u32) -> Shape {
    s.clone().with_location(l)
}

// ---------------------------------------------------------------------------
// GAP carriers (architecture difference #4)
// ---------------------------------------------------------------------------

/// GAP: BRepMAT2d_Explorer (TKTopAlgo/BRepMAT2d — the wrapper class is not
/// translated; the underlying MAT2d stack is the real topalgo translation).
/// OCCT anchor: BRepMAT2d_Explorer.hxx L28-77.
#[derive(Debug, Default)]
pub(super) struct BRepMAT2dExplorerCarrier;

impl BRepMAT2dExplorerCarrier {
    /// OCCT BRepMAT2d_Explorer::BRepMAT2d_Explorer().
    pub(super) fn new() -> Self {
        BRepMAT2dExplorerCarrier
    }

    /// OCCT BRepMAT2d_Explorer::BRepMAT2d_Explorer(aFace).
    pub(super) fn new_with_face(_a_face: &Shape) -> Self {
        panic!(
            "GAP: BRepMAT2d_Explorer (TKTopAlgo/BRepMAT2d not translated) — see file header"
        )
    }

    /// OCCT BRepMAT2d_Explorer::Perform(aFace).
    pub(super) fn perform(&mut self, _a_face: &Shape) {
        panic!(
            "GAP: BRepMAT2d_Explorer (TKTopAlgo/BRepMAT2d not translated) — see file header"
        )
    }
}

/// GAP: BRepMAT2d_BisectingLocus (TKTopAlgo/BRepMAT2d — the wrapper class
/// is not translated; the underlying MAT / MAT2d / Bisector stacks are the
/// real topalgo translations).  OCCT anchor:
/// BRepMAT2d_BisectingLocus.hxx L60-133.
#[derive(Debug, Default)]
pub(super) struct BRepMAT2dBisectingLocusCarrier;

impl BRepMAT2dBisectingLocusCarrier {
    /// OCCT BRepMAT2d_BisectingLocus::BRepMAT2d_BisectingLocus().
    pub(super) fn new() -> Self {
        BRepMAT2dBisectingLocusCarrier
    }

    /// OCCT BRepMAT2d_BisectingLocus::Compute(anExplo, LineIndex, aSide).
    pub(super) fn compute(
        &mut self,
        _an_explo: &mut BRepMAT2dExplorerCarrier,
        _line_index: i32,
        _a_side: MatSide,
    ) {
        panic!(
            "GAP: BRepMAT2d_BisectingLocus (TKTopAlgo/BRepMAT2d not translated) — \
             see file header"
        )
    }

    /// OCCT BRepMAT2d_BisectingLocus::Graph() — the real MAT_Graph handle
    /// type (the topalgo translation).
    pub(super) fn graph(&self) -> HandleMatGraph {
        panic!(
            "GAP: BRepMAT2d_BisectingLocus (TKTopAlgo/BRepMAT2d not translated) — \
             see file header"
        )
    }

    /// OCCT BRepMAT2d_BisectingLocus::GeomBis(A, Reverse).
    pub(super) fn geom_bis(&self, _a: &HandleMatArc, _reverse: &mut bool) -> BisectorBisec {
        panic!(
            "GAP: BRepMAT2d_BisectingLocus (TKTopAlgo/BRepMAT2d not translated) — \
             see file header"
        )
    }
}

/// GAP: BRepMAT2d_LinkTopoBilo (TKTopAlgo/BRepMAT2d — the wrapper class is
/// not translated).  OCCT anchor: BRepMAT2d_LinkTopoBilo.hxx L34-71.
#[derive(Debug, Default)]
pub(super) struct BRepMAT2dLinkTopoBiloCarrier;

impl BRepMAT2dLinkTopoBiloCarrier {
    /// OCCT BRepMAT2d_LinkTopoBilo::BRepMAT2d_LinkTopoBilo().
    pub(super) fn new() -> Self {
        BRepMAT2dLinkTopoBiloCarrier
    }

    /// OCCT BRepMAT2d_LinkTopoBilo::BRepMAT2d_LinkTopoBilo(Explo, BiLo).
    pub(super) fn new_with_links(
        _explo: &BRepMAT2dExplorerCarrier,
        _bilo: &BRepMAT2dBisectingLocusCarrier,
    ) -> Self {
        panic!(
            "GAP: BRepMAT2d_LinkTopoBilo (TKTopAlgo/BRepMAT2d not translated) — \
             see file header"
        )
    }

    /// OCCT BRepMAT2d_LinkTopoBilo::Init(S).
    pub(super) fn init(&mut self, _s: &Shape) {
        panic!(
            "GAP: BRepMAT2d_LinkTopoBilo (TKTopAlgo/BRepMAT2d not translated) — \
             see file header"
        )
    }

    /// OCCT BRepMAT2d_LinkTopoBilo::More().
    pub(super) fn more(&self) -> bool {
        panic!(
            "GAP: BRepMAT2d_LinkTopoBilo (TKTopAlgo/BRepMAT2d not translated) — \
             see file header"
        )
    }

    /// OCCT BRepMAT2d_LinkTopoBilo::Next().
    pub(super) fn next(&mut self) {
        panic!(
            "GAP: BRepMAT2d_LinkTopoBilo (TKTopAlgo/BRepMAT2d not translated) — \
             see file header"
        )
    }

    /// OCCT BRepMAT2d_LinkTopoBilo::Value() — the real MAT_BasicElt handle
    /// type (the topalgo translation).
    pub(super) fn value(&self) -> HandleMatBasicElt {
        panic!(
            "GAP: BRepMAT2d_LinkTopoBilo (TKTopAlgo/BRepMAT2d not translated) — \
             see file header"
        )
    }

    /// OCCT BRepMAT2d_LinkTopoBilo::GeneratingShape(aBE).
    pub(super) fn generating_shape(&self, _a_be: &HandleMatBasicElt) -> Shape {
        panic!(
            "GAP: BRepMAT2d_LinkTopoBilo (TKTopAlgo/BRepMAT2d not translated) — \
             see file header"
        )
    }
}

/// GAP: BRepTools_Quilt (TKTopAlgo/BRepTools — not translated; the
/// bi_tgte / make_offset local carrier precedent).  OCCT anchor:
/// BRepTools_Quilt.hxx.
#[derive(Debug, Default)]
pub(super) struct BRepToolsQuiltCarrier;

impl BRepToolsQuiltCarrier {
    /// OCCT BRepTools_Quilt::BRepTools_Quilt().
    pub(super) fn new() -> Self {
        BRepToolsQuiltCarrier
    }

    /// OCCT BRepTools_Quilt::Bind(E1, E2).
    pub(super) fn bind(&mut self, _e1: &Shape, _e2: &Shape) {
        panic!("GAP: BRepTools_Quilt (TKTopAlgo/BRepTools not translated) — see file header")
    }

    /// OCCT BRepTools_Quilt::Add(S).
    pub(super) fn add(&mut self, _s: &Shape) {
        panic!("GAP: BRepTools_Quilt (TKTopAlgo/BRepTools not translated) — see file header")
    }

    /// OCCT BRepTools_Quilt::Shells().
    pub(super) fn shells(&self) -> Shape {
        panic!("GAP: BRepTools_Quilt (TKTopAlgo/BRepTools not translated) — see file header")
    }

    /// OCCT BRepTools_Quilt::IsCopied(S).
    pub(super) fn is_copied(&self, _s: &Shape) -> bool {
        panic!("GAP: BRepTools_Quilt (TKTopAlgo/BRepTools not translated) — see file header")
    }

    /// OCCT BRepTools_Quilt::Copy(S).
    pub(super) fn copy(&self, _s: &Shape) -> Shape {
        panic!("GAP: BRepTools_Quilt (TKTopAlgo/BRepTools not translated) — see file header")
    }
}

/// GAP: BRepSweep_Prism (TKTopAlgo/BRepSweep — not translated; the
/// feat::loc_ope_prism local carrier precedent).  OCCT anchor:
/// BRepSweep_Prism.hxx.
#[derive(Debug, Default)]
pub(super) struct BRepSweepPrismCarrier;

impl BRepSweepPrismCarrier {
    /// OCCT BRepSweep_Prism::BRepSweep_Prism(S, V, SulP).
    pub(super) fn new(_the_base: &Shape, _the_vec: DVec3, _sul_p: bool) -> Self {
        panic!("GAP: BRepSweep_Prism (TKTopAlgo/BRepSweep not translated) — see file header")
    }

    /// OCCT BRepSweep_Prism::Shape() — the whole sweep result.
    pub(super) fn shape(&self) -> Shape {
        panic!("GAP: BRepSweep_Prism (unreachable while the engine is a GAP)")
    }

    /// OCCT BRepSweep_Prism::LastShape() — the top face.
    pub(super) fn last_shape(&self) -> Shape {
        panic!("GAP: BRepSweep_Prism (unreachable while the engine is a GAP)")
    }

    /// OCCT BRepSweep_Prism::LastShape(S) — the generated last shape of an
    /// ancestor.
    pub(super) fn last_shape_of(&self, _the_s: &Shape) -> Shape {
        panic!("GAP: BRepSweep_Prism (unreachable while the engine is a GAP)")
    }

    /// OCCT BRepSweep_Prism::Shape(S) — the generated shape of an ancestor.
    pub(super) fn shape_of(&self, _the_s: &Shape) -> Shape {
        panic!("GAP: BRepSweep_Prism (unreachable while the engine is a GAP)")
    }
}

/// GAP: BRepSweep_Revol (TKTopAlgo/BRepSweep — not translated; the
/// feat::loc_ope_revol local carrier precedent).  OCCT anchor:
/// BRepSweep_Revol.hxx.
#[derive(Debug, Default)]
pub(super) struct BRepSweepRevolCarrier;

impl BRepSweepRevolCarrier {
    /// OCCT BRepSweep_Revol::BRepSweep_Revol(S, A, Pfin).
    pub(super) fn new(_the_base: &Shape, _the_axis: rcad_kernel::math::gp::Ax1, _pfin: bool) -> Self {
        panic!("GAP: BRepSweep_Revol (TKTopAlgo/BRepSweep not translated) — see file header")
    }

    /// OCCT BRepSweep_Revol::Shape(S) — the generated shape of an ancestor.
    pub(super) fn shape_of(&self, _the_s: &Shape) -> Shape {
        panic!("GAP: BRepSweep_Revol (unreachable while the engine is a GAP)")
    }
}

// ---------------------------------------------------------------------------
// Small re-hosts shared by the Evolved files
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Degenerated(E).
pub(super) fn brep_tool_degenerated(e: &Shape) -> bool {
    matches!(e.data.as_ref(), TShape::Edge(ed) if ed.degenerated)
}

/// OCCT BRep_Tool::Pnt(V).
pub(super) fn brep_tool_pnt(v: &Shape) -> DVec3 {
    match v.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => DVec3::ZERO,
    }
}

/// OCCT BRep_Tool::Surface(F, L) — the face surface.
pub(super) fn brep_tool_surface(f: &Shape) -> Option<rcad_kernel::geom::Surface3> {
    match f.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT BRep_Tool::Tolerance(F / V / E).
pub(super) fn brep_tool_tolerance(s: &Shape) -> f64 {
    match s.data.as_ref() {
        TShape::Face(fd) => fd.tolerance,
        TShape::Vertex(vd) => vd.tolerance,
        TShape::Edge(ed) => ed.tolerance,
        _ => 0.0,
    }
}

/// OCCT BRep_Builder::MakeCompound(C).
pub(super) fn builder_make_compound() -> Shape {
    crate::brep_algo::tool::builder_make_compound()
}

/// OCCT BRep_Builder::MakeWire(W).
pub(super) fn builder_make_wire() -> Shape {
    crate::brep_algo::tool::builder_make_wire()
}

/// OCCT BRep_Builder::MakeVertex(V).
pub(super) fn builder_make_vertex() -> Shape {
    crate::brep_algo::tool::builder_make_vertex()
}

/// OCCT BRep_Builder::MakeSolid(S).
pub(super) fn builder_make_solid() -> Shape {
    Shape {
        data: Arc::new(TShape::Solid(rcad_kernel::topods::TSolidData {
            my_shapes: Vec::new(),
            flags: rcad_kernel::topods::tshape_flags::DEFAULT,
            shells: Vec::new(),
            internal_vertices: Vec::new(),
            internal_edges: Vec::new(),
        })),
        index: usize::MAX,
        location: 0,
        orientation: Orientation::Forward,
    }
}

/// OCCT BRep_Builder::Add(Solid, Shell).
pub(super) fn builder_add_solid_shell(solid: &mut Shape, shell: &Shape) {
    if let TShape::Solid(sd) = Arc::make_mut(&mut solid.data) {
        sd.shells.push(shell.clone());
    }
}

/// OCCT BRep_Builder::Add(W, E).
pub(super) fn builder_add_wire_edge(w: &mut Shape, e: &Shape) {
    crate::brep_algo::tool::builder_add_wire_edge(w, e);
}

/// OCCT BRep_Builder::Add(C, S).
pub(super) fn builder_add_compound_shape(c: &mut Shape, s: &Shape) {
    crate::brep_algo::tool::builder_add_compound_shape(c, s);
}

/// OCCT BRep_Builder::Add(F, W) — the first wire is the outer wire.
pub(super) fn builder_add_face_wire(f: &mut Shape, w: &Shape) {
    crate::brep_algo::tool::builder_add_face_wire(f, w);
}

/// OCCT BRep_Builder::MakeFace(F, S, L, TolDegen) — the face over a
/// surface.
pub(super) fn builder_make_face_surface(
    s: &rcad_kernel::geom::Surface3,
    tol_degen: f64,
) -> Shape {
    Shape {
        data: Arc::new(TShape::Face(rcad_kernel::topods::TFaceData {
            my_shapes: Vec::new(),
            flags: rcad_kernel::topods::tshape_flags::DEFAULT,
            surface: Some(s.clone()),
            surface_location: 0,
            outer_wire: Shape::null(),
            inner_wires: Vec::new(),
            sample_point: None,
            uv_domain: None,
            internal_vertices: Vec::new(),
            tolerance: tol_degen,
            natural_restriction: false,
        })),
        index: usize::MAX,
        location: 0,
        orientation: Orientation::Forward,
    }
}

/// OCCT BRepLib_MakeWire(E) — the wire from a single edge.
pub(super) fn brep_lib_make_wire_from_edge(e: &Shape) -> Shape {
    let mut w = builder_make_wire();
    builder_add_wire_edge(&mut w, e);
    w
}

/// OCCT BRepLib::UpdateTolerances(S, WithShape) — no-op (the rcad kernel
/// keeps the tolerance invariants inline; the brep_lib.rs stub precedent).
pub(super) fn brep_lib_update_tolerances(_s: &Shape, _with_shape: bool) {}

/// OCCT BRepLib::SameParameter(E) — no-op (the brep_lib.rs stub precedent).
pub(super) fn brep_lib_same_parameter(_e: &Shape) {}

/// OCCT BRepLib::BuildCurves3d(S) — no-op (the rcad edge encodings carry
/// the 3d curve from construction; the offset_wire.rs BuildCurve3d GAP
/// note).
pub(super) fn brep_lib_build_curves3d(_s: &Shape) {}

/// OCCT BRep_Builder::Continuity(E, F1, F2, C) — the rcad edge TShape has
/// no regularity table yet (the kernel BRepBuilder::continuity is
/// pool-bound); documented no-op (the brep_lib.rs stub precedent).
pub(super) fn builder_continuity(_e: &Shape, _f1: &Shape, _f2: &Shape, _c: i32) {}

// ---------------------------------------------------------------------------
// File static (BRepFill_Evolved.cxx L99-103)
// ---------------------------------------------------------------------------

/// OCCT static BRepFill_Confusion (L99-103) — the package-local constant
/// function.
pub(super) fn brep_fill_confusion() -> f64 {
    let tol = 1.0e-6;
    tol
}

// ---------------------------------------------------------------------------
// BRepFill_Evolved (hxx L40-199)
// ---------------------------------------------------------------------------

/// OCCT BRepFill_Evolved (hxx L40-199) — constructs an evolved volume from
/// a spine (wire or face) and a profile (wire).
#[derive(Debug, Clone)]
pub struct BRepFillEvolved {
    pub(super) my_spine: Shape,               // OCCT: mySpine (TopoDS_Face)
    pub(super) my_profile: Shape,             // OCCT: myProfile (TopoDS_Wire)
    pub(super) my_shape: Shape,               // OCCT: myShape (TopoDS_Shape)
    pub(super) my_is_done: bool,              // OCCT: myIsDone
    pub(super) my_spine_type: bool,           // OCCT: mySpineType
    pub(super) my_join_type: GeomAbsJoinType, // OCCT: myJoinType
    pub(super) my_map: EvolutionMap,          // OCCT: myMap
    pub(super) my_top: Shape,                 // OCCT: myTop
    pub(super) my_bottom: Shape,              // OCCT: myBottom
    /// The per-object location table (architecture difference #2).
    pub(super) my_locations: LocationTable,
}

impl Default for BRepFillEvolved {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepFillEvolved {
    /// OCCT BRepFill_Evolved::BRepFill_Evolved() (cxx L186-190).
    pub fn new() -> Self {
        BRepFillEvolved {
            my_spine: Shape::null(),
            my_profile: Shape::null(),
            my_shape: Shape::null(),
            my_is_done: false,
            // OCCT L187-188: myIsDone(false), mySpineType(true).
            my_spine_type: true,
            // GeomAbs_JoinType default (GeomAbs_Arc) — hxx default arg.
            my_join_type: GeomAbsJoinType::Arc,
            my_map: EvolutionMap::new(),
            my_top: Shape::null(),
            my_bottom: Shape::null(),
            my_locations: LocationTable::default(),
        }
    }

    /// OCCT BRepFill_Evolved::BRepFill_Evolved(Spine, Profile, AxeProf,
    /// Join, Solid) — the wire spine (cxx L194-203).
    pub fn new_with_wire_spine(
        spine: &Shape,
        profile: &Shape,
        axe_prof: &rcad_kernel::math::gp::Ax3,
        join: GeomAbsJoinType,
        solid: bool,
    ) -> Self {
        let mut r = BRepFillEvolved::new();
        // OCCT L200: myIsDone(false) — the Default impl.
        r.perform_with_wire_spine(spine, profile, axe_prof, join, solid);
        r
    }

    /// OCCT BRepFill_Evolved::BRepFill_Evolved(Spine, Profile, AxeProf,
    /// Join, Solid) — the face spine (cxx L207-215).
    pub fn new_with_face_spine(
        spine: &Shape,
        profile: &Shape,
        axe_prof: &rcad_kernel::math::gp::Ax3,
        join: GeomAbsJoinType,
        solid: bool,
    ) -> Self {
        let mut r = BRepFillEvolved::new();
        // OCCT L212: myIsDone(false).
        r.perform_with_face_spine(spine, profile, axe_prof, join, solid);
        r
    }

    /// OCCT BRepFill_Evolved::IsDone() (cxx L2057-2060).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT BRepFill_Evolved::Shape() (cxx L2064-2067) — returns the
    /// generated shape.
    pub fn shape(&self) -> &Shape {
        &self.my_shape
    }

    /// OCCT BRepFill_Evolved::GeneratedShapes(SpineShape, ProfShape)
    /// (cxx L1797-1810) — returns the shapes created from a subshape of
    /// the spine and a subshape of the profile.
    pub fn generated_shapes(&self, spine_shape: &Shape, prof_shape: &Shape) -> Vec<Shape> {
        if self.my_map.is_bound(spine_shape) && self.my_map.find(spine_shape).is_bound(prof_shape)
        {
            self.my_map.find(spine_shape).find(prof_shape).clone()
        } else {
            // OCCT L1807-1809: the static Empty list.
            Vec::new()
        }
    }

    /// OCCT BRepFill_Evolved::JoinType() (cxx L2071-2074).
    pub fn join_type(&self) -> GeomAbsJoinType {
        self.my_join_type
    }

    /// OCCT BRepFill_Evolved::Top() (cxx L1783-1786) — returns the face
    /// Top if Solid is True in the constructor.
    pub fn top(&self) -> &Shape {
        &self.my_top
    }

    /// OCCT BRepFill_Evolved::Bottom() (cxx L1790-1793) — returns the face
    /// Bottom if Solid is True in the constructor.
    pub fn bottom(&self) -> &Shape {
        &self.my_bottom
    }

    /// OCCT BRepFill_Evolved::Generated() (cxx L1814-1821) — the myMap
    /// accessor.
    pub(super) fn generated(&mut self) -> &mut EvolutionMap {
        &mut self.my_map
    }

    /// OCCT BRepFill_Evolved::ChangeShape() (cxx L1959-1962).
    pub(super) fn change_shape(&mut self) -> &mut Shape {
        &mut self.my_shape
    }

    /// OCCT BRepFill_Evolved::SetWork(Sp, Pr) (cxx L614-618).
    pub(super) fn set_work(&mut self, sp: &Shape, pr: &Shape) {
        self.my_spine = sp.clone();
        self.my_profile = pr.clone();
    }

    // -------------------------------------------------------------------
    // Perform (cxx L309-330)
    // -------------------------------------------------------------------

    /// OCCT BRepFill_Evolved::Perform(Spine(Wire), Profile, AxeProf, Join,
    /// Solid) (cxx L309-318).
    pub fn perform_with_wire_spine(
        &mut self,
        spine: &Shape,
        profile: &Shape,
        axe_prof: &rcad_kernel::math::gp::Ax3,
        join: GeomAbsJoinType,
        solid: bool,
    ) {
        self.my_spine_type = false;
        // OCCT L316: TopoDS_Face aFace = BRepLib_MakeFace(Spine, true).
        let a_face =
            crate::brep_fill::brep_fill_axe::BRepLibMakeFaceWire::new(&mut rcad_kernel::topods::BRep::new(), spine, true);
        let a_face = a_face.face();
        self.private_perform(&a_face, profile, axe_prof, join, solid);
    }

    /// OCCT BRepFill_Evolved::Perform(Spine(Face), Profile, AxeProf, Join,
    /// Solid) (cxx L322-330).
    pub fn perform_with_face_spine(
        &mut self,
        spine: &Shape,
        profile: &Shape,
        axe_prof: &rcad_kernel::math::gp::Ax3,
        join: GeomAbsJoinType,
        solid: bool,
    ) {
        self.my_spine_type = true;
        self.private_perform(spine, profile, axe_prof, join, solid);
    }

    // -------------------------------------------------------------------
    // PrivatePerform (cxx L334-544)
    // -------------------------------------------------------------------

    /// OCCT BRepFill_Evolved::PrivatePerform(Spine, Profile, AxeProf, Join,
    /// Solid) (cxx L334-544).
    pub(super) fn private_perform(
        &mut self,
        spine: &Shape,
        profile: &Shape,
        axe_prof: &rcad_kernel::math::gp::Ax3,
        join: GeomAbsJoinType,
        solid: bool,
    ) {
        // OCCT L340-345.
        let a_local_shape = shape_oriented(spine, Orientation::Forward);
        self.my_spine = a_local_shape;
        let a_local_shape = shape_oriented(profile, Orientation::Forward);
        self.my_profile = a_local_shape;
        self.my_join_type = join;
        self.my_map.clear();

        if matches!(self.my_join_type, GeomAbsJoinType::Intersection) {
            // OCCT L349-352: throw Standard_NotImplemented.
            panic!("Standard_NotImplemented: BRepFill_Evolved::PrivatePerform");
        }

        // OCCT L354-356: WorkProf, WorkSpine, WPIte.
        let mut work_prof: Vec<Shape> = Vec::new();
        let mut work_spine = Shape::null();

        //-------------------------------------------------------------------
        // Positioning of mySpine and myProfil in the workspace.
        //-------------------------------------------------------------------
        // OCCT L361-367.
        let l_spine = self.find_location(&self.my_spine.clone());
        let t = crate::brep_fill::brep_fill_evolved_d::trsf_set_transformation_ax3(axe_prof);
        let l_profile = self.my_locations.intern(t);
        let init_ls = self.my_spine.location;
        let init_lp = self.my_profile.location;
        self.transform_init_work(l_spine, l_profile);

        //------------------------------------------------------------------
        // projection of the profile and cut of the spine.
        //------------------------------------------------------------------
        // OCCT L372-375.
        let mut map_prof: DataMapOfShapeItem<Shape> = DataMapOfShapeItem::new();
        let mut map_spine: DataMapOfShapeItem<Shape> = DataMapOfShapeItem::new();

        self.prepare_profile(&mut work_prof, &mut map_prof);
        self.prepare_spine(&mut work_spine, &mut map_spine);

        // OCCT L377-380.
        let tol = brep_fill_confusion();
        let mut ya_left = false;
        let mut ya_right = false;
        let mut sp = Shape::null();

        for w in work_prof.clone() {
            sp = w;
            if crate::brep_fill::brep_fill_evolved_d::side(&sp, tol) < 4 {
                ya_left = true;
            } else {
                ya_right = true;
            }
            if ya_left && ya_right {
                break;
            }
        }

        // OCCT L399-400: Face; Locus.
        let mut face = Shape::null();
        let mut locus = BRepMAT2dBisectingLocusCarrier::new();

        //----------------------------------------------------------
        // Initialisation of cut volevo.
        // For each part of the profile create a volevo added to CutVevo
        //----------------------------------------------------------
        // OCCT L406-411.
        let mut cut_vevo = BRepFillEvolved::new();
        let mut wp = builder_make_wire();

        for w in &work_prof {
            for e in wire_edges(w) {
                builder_add_wire_edge(&mut wp, &e);
            }
        }
        cut_vevo.set_work(&work_spine.clone(), &wp);

        // OCCT L422-423.
        let mut glue = BRepToolsQuiltCarrier::new();
        let mut c_side;

        //---------------------------------
        // Construction of vevos to the left.
        //---------------------------------
        if ya_left {
            //-----------------------------------------------------
            // Calculate the map of bisector locations at the left.
            // and links Topology -> base elements of the map.
            //-----------------------------------------------------
            // OCCT L434-436.
            let mut exp = BRepMAT2dExplorerCarrier::new_with_face(&work_spine.clone());
            locus.compute(&mut exp, 1, MatSide::Left);
            let mut link = BRepMAT2dLinkTopoBiloCarrier::new_with_links(&exp, &locus);

            for w in work_prof.clone() {
                sp = w;
                c_side = crate::brep_fill::brep_fill_evolved_d::side(&sp, tol);
                //-----------------------------------------------
                // Construction and adding of elementary volevo.
                //-----------------------------------------------
                let mut vevo = BRepFillEvolved::new();
                if c_side == 1 {
                    vevo.elementary_perform(&work_spine.clone(), &sp, &locus, &mut link, join);
                } else if c_side == 2 {
                    vevo.planar_perform(&work_spine.clone(), &sp, &locus, &mut link, join);
                } else if c_side == 3 {
                    vevo.vertical_perform(&work_spine.clone(), &sp, &locus, &mut link, join);
                }
                cut_vevo.add(&mut vevo, &sp, &mut glue);
            }
        }

        //---------------------------------
        // Construction of vevos to the right.
        //---------------------------------
        if ya_right {
            //-----------------------------------
            // Decomposition of the face into wires.
            //-----------------------------------
            // OCCT L470-507.
            for spine_wire in crate::brep_algo::tool::explorer(
                &work_spine.clone(),
                ShapeType::Wire,
                ShapeType::Shape,
            ) {
                //----------------------------------------------
                // Calculate the map to the right of the current wire.
                //----------------------------------------------
                // OCCT L476-480.
                let b = builder_make_face_surface(
                    &rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane::new(
                        DVec3::ZERO,
                        DVec3::Z,
                    )),
                    brep_tool_tolerance(&work_spine.clone()),
                );
                let a_local_shape_rev = shape_reversed(&spine_wire);
                face = b;
                builder_add_face_wire(&mut face, &a_local_shape_rev);
                let mut exp = BRepMAT2dExplorerCarrier::new_with_face(&face.clone());
                locus.compute(&mut exp, 1, MatSide::Left);
                let mut link = BRepMAT2dLinkTopoBiloCarrier::new_with_links(&exp, &locus);

                for w in work_prof.clone() {
                    sp = w;
                    c_side = crate::brep_fill::brep_fill_evolved_d::side(&sp, tol);
                    //-----------------------------------------------
                    // Construction and adding of an elementary volevo
                    //-----------------------------------------------
                    let mut vevo = BRepFillEvolved::new();
                    if c_side == 4 {
                        vevo.elementary_perform(&face.clone(), &sp, &locus, &mut link, join);
                    } else if c_side == 5 {
                        vevo.planar_perform(&face.clone(), &sp, &locus, &mut link, join);
                    } else if c_side == 6 {
                        vevo.vertical_perform(&face.clone(), &sp, &locus, &mut link, join);
                    }
                    cut_vevo.add(&mut vevo, &sp, &mut glue);
                }
            }
        }

        // OCCT L510-513.
        if solid {
            cut_vevo.add_top_and_bottom(&mut glue);
        }

        //-------------------------------------------------------------------------
        // Gluing of regularites on parallel edges generate4d by vertices of the
        // cut of the profile.
        //-------------------------------------------------------------------------
        // OCCT L519.
        cut_vevo.continuity_on_offset_edge(&work_prof);

        //-----------------------------------------------------------------
        // construction of the shape via the quilt, ie:
        // - sharing of topologies of elementary added volevos.
        // - Orientation of faces correspondingly to each other.
        //-----------------------------------------------------------------
        // OCCT L526-527.
        let scv = cut_vevo.change_shape();
        *scv = glue.shells();

        //------------------------------------------------------------------------
        // Transfer of the map of generated elements and of the shape of Cutvevo
        // in myMap and repositioning in the initial space.
        //------------------------------------------------------------------------
        // OCCT L532.
        let l_spine_inverted = self.my_locations.inverted(l_spine);
        self.transfert(
            &mut cut_vevo,
            &map_prof,
            &map_spine,
            l_spine_inverted,
            init_ls,
            init_lp,
        );

        // Orientation of the solid.
        // OCCT L535-538.
        if solid {
            self.make_solid();
        }

        //  modified by NIZHNY-EAP Mon Jan 24 11:26:48 2000 ___BEGIN___
        // OCCT L541.
        brep_lib_update_tolerances(&self.my_shape.clone(), false);
        //  modified by NIZHNY-EAP Mon Jan 24 11:26:50 2000 ___END___
        self.my_is_done = true;
    }
}

// ---------------------------------------------------------------------------
// File statics re-exported to the sibling Evolved files
// ---------------------------------------------------------------------------

/// OCCT static EdgeVertices (L172-182) — the oriented end vertices.
pub(super) fn edge_vertices(e: &Shape) -> (Shape, Shape) {
    top_exp_vertices(e)
}

/// OCCT BRepTools_WireExplorer(W) — the wire edges (the
/// brep_fill_pipe.rs reduction).
pub(super) fn wire_edges(w: &Shape) -> Vec<Shape> {
    match w.data.as_ref() {
        TShape::Wire(wd) => wd.edges.clone(),
        _ => Vec::new(),
    }
}

/// OCCT NCollection_DataMap<Handle(MAT_Node), DataMap<TopoDS_Shape,
/// TopoDS_Shape>> — the MapNodeVertex type (handle keys compare by pointer
/// identity).
#[derive(Clone, Default)]
pub(super) struct MapNodeVertex {
    entries: Vec<(HandleMatNode, DataMapOfShapeItem<Shape>)>,
}

impl MapNodeVertex {
    pub(super) fn new() -> Self {
        MapNodeVertex {
            entries: Vec::new(),
        }
    }

    fn same_node(a: &HandleMatNode, b: &HandleMatNode) -> bool {
        Arc::ptr_eq(a, b)
    }

    /// OCCT IsBound(aNode) && operator()(aNode).IsBound(ShapeOnNode).
    pub(super) fn node_shape_is_bound(&self, node: &HandleMatNode, shape: &Shape) -> bool {
        self.entries
            .iter()
            .find(|(n, _)| Self::same_node(n, node))
            .map(|(_, m)| m.is_bound(shape))
            .unwrap_or(false)
    }

    /// OCCT operator()(aNode)(ShapeOnNode).
    pub(super) fn node_shape_find(&self, node: &HandleMatNode, shape: &Shape) -> Shape {
        self.entries
            .iter()
            .find(|(n, _)| Self::same_node(n, node))
            .map(|(_, m)| m.find(shape).clone())
            .expect("MapNodeVertex(node)(shape)")
    }

    /// OCCT Bind / operator()(aNode).Bind(ShapeOnNode, VN) sequence.
    pub(super) fn node_shape_bind(&mut self, node: &HandleMatNode, shape: &Shape, vn: &Shape) {
        if let Some((_, m)) = self
            .entries
            .iter_mut()
            .find(|(n, _)| Self::same_node(n, node))
        {
            m.bind(shape, vn.clone());
        } else {
            let mut m = DataMapOfShapeItem::new();
            m.bind(shape, vn.clone());
            self.entries.push((node.clone(), m));
        }
    }
}

/// OCCT MAT_Node read helpers over the real topalgo handle types (the
/// Arc<RwLock<...>> projections).
pub(super) fn node_on_basic_elt(n: &HandleMatNode) -> bool {
    n.read().expect("MAT_Node").on_basic_elt()
}

pub(super) fn node_infinite(n: &HandleMatNode) -> bool {
    n.read().expect("MAT_Node").infinite()
}

pub(super) fn node_distance(n: &HandleMatNode) -> f64 {
    n.read().expect("MAT_Node").distance()
}

pub(super) fn arc_first_element(a: &HandleMatArc) -> HandleMatBasicElt {
    a.read()
        .expect("MAT_Arc")
        .first_element()
        .expect("MAT_Arc::FirstElement")
}

pub(super) fn arc_second_element(a: &HandleMatArc) -> HandleMatBasicElt {
    a.read()
        .expect("MAT_Arc")
        .second_element()
        .expect("MAT_Arc::SecondElement")
}

pub(super) fn arc_first_node(a: &HandleMatArc) -> HandleMatNode {
    a.read()
        .expect("MAT_Arc")
        .first_node()
        .expect("MAT_Arc::FirstNode")
}

pub(super) fn arc_second_node(a: &HandleMatArc) -> HandleMatNode {
    a.read()
        .expect("MAT_Arc")
        .second_node()
        .expect("MAT_Arc::SecondNode")
}

pub(super) fn basic_elt_start_arc(be: &HandleMatBasicElt) -> HandleMatArc {
    be.read()
        .expect("MAT_BasicElt")
        .start_arc()
        .expect("MAT_BasicElt::StartArc")
}
