//! BRepGraphInc — topology graph incidence storage backend.
//!
//! OCCT TKBRep BRepGraphInc package.
//!
//! Provides the typed entity/ref/relation storage model that underpins
//! BRepGraph.  This layer is the lowest-level data model; the higher-level
//! Populate (TopoDS → Storage) and Reconstruct (Storage → TopoDS) pipelines
//! can be added later as separate modules that only depend on the store API.

pub mod id;
pub mod def;
pub mod ref_mod;
pub mod rel;
pub mod store;

pub use id::*;
pub use def::*;
pub use ref_mod::*;
pub use rel::*;
pub use store::*;
