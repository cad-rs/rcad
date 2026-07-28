//! OCCT TKBRep alignment crate.
//!
//! Contains the Rust translation of OCCT's TKBRep toolkit:
//!
//! | Module   | OCCT Package      | Description              |
//! |----------|-------------------|--------------------------|
//! | adaptor  | BRepAdaptor       | Edge/surface adaptors    |
//! | tools    | BRepTools         | BRep I/O, queries, transforms |
//! | graph    | BRepGraph         | Topology graph analysis  |

pub mod adaptor;
pub mod lprop;
pub mod tools;
pub mod graph;
