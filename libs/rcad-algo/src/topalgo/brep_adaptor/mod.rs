//! BRepAdaptor (TKBRep) — the shape adaptors.
//!
//! 1:1 translations:
//! - [`curve2d::BRepCurve2d`] — BRepAdaptor_Curve2d.hxx/.cxx (the pcurve
//!   adaptor of an edge on a face).

pub mod curve2d;

pub use curve2d::BRepCurve2d;
