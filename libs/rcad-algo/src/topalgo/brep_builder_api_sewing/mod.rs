//! OCCT BRepBuilderAPI_Sewing (TKTopAlgo/BRepBuilderAPI).
//!
//! 1:1 translation of `BRepBuilderAPI_Sewing.hxx` (L17-424, + the .lxx
//! inline accessors) and `BRepBuilderAPI_Sewing.cxx` (L1-5946).  The class
//! body lives in [`sewing`], split across child files by OCCT method groups
//! (each child is one `impl BRepBuilderAPISewing` block + its OCCT static
//! file-scope helpers).
//!
//! Architecture differences (the brep_tools_quilt.rs numbering style):
//! 1. `NCollection_IndexedDataMap<TopoDS_Shape, T, TopTools_ShapeMapHasher>`
//!    -> `IndexMap<ShapeKey, (Shape, T)>` (the key shape travels in the
//!    value; the insertion order stands for the OCCT index order).
//!    `NCollection_IndexedMap` / `NCollection_Map` -> `IndexMap<ShapeKey,
//!    Shape>` (the OCCT set iteration order is bucket-based and
//!    unspecified; the rcad order is the insertion order).
//!    `NCollection_DataMap<K, V>` -> `HashMap<ShapeKey, V>` (Bind replaces).
//!    `ShapeKey = (TShape ptr, Location)` is the TopTools_ShapeMapHasher
//!    identity (orientation ignored).
//! 2. `occ::handle<BRepTools_ReShape>` -> the flattened
//!    `shhealing::shape_build::reshape::ShapeBuildReShape` value member.
//!    OCCT `Apply/Replace/Remove/Clear/IsRecorded` map one-to-one; the rcad
//!    body takes the owning `BRep` pool, so the OCCT `const` methods that
//!    touch the reshape (`IsDegenerated`, `IsModifiedSubShape`,
//!    `ModifiedSubShape`, `Dump`, `CreateOutputInformations`) become `&mut
//!    self` (arch. diff. S2).
//! 3. `BRep_Builder` / the global shape space -> the class-local
//!    `my_brep: BRep` pool + its `BRepBuilder` in-place mutators (the
//!    bi_tgte_blended.rs arch. diff. #4/#19 carrier).  Every TShape the
//!    sewing creates is a slot of `my_brep`; the input shapes are
//!    materialized into the pool on `Add` (`BRep::import_shape_tree`), so
//!    the products read back as non-null through `Shape::is_null()`.
//! 4. `NCollection_Sequence` -> `Vec` (0-based with the explicit -1 offset
//!    at the access sites), `NCollection_Array1` (1-based) -> a 0-based
//!    `Vec` + the explicit offset (the bi_tgte_blended.rs arch. diff. #21
//!    convention).
//! 5. `Message_ProgressScope` / `OCCT_DEBUG` chronometers — omitted (rcad
//!    carries no progress scope; the debug output is debug-only).
//! 6. `NCollection_UBTree<int, Bnd_Box>` + BndBoxTreeSelector and
//!    `NCollection_CellFilter` + VertexInspector — local re-hosts in
//!    [`sewing::selectors`] keeping the OCCT Add/Fill/Select/Inspect call
//!    form; `Select` walks the filled entries in insertion order (the same
//!    documented reduction as the ShapeAnalysis box_bnd_tree.rs carrier).

#[path = "sewing.rs"]
pub mod sewing;

pub use sewing::BRepBuilderAPISewing;
