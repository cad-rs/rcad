//! Adaptor2d (TKG2d) — 2D curve adaptors over the surface domain.
//!
//! 1:1 translations:
//! - [`line2d::Line2dAdaptor`] — Adaptor2d_Line2d.hxx/.cxx (the iso-line
//!   restriction adaptor of the TopolTool).  It implements
//!   [`crate::geomalgo::geom2d_int::Curve2dAdaptor`], the OCCT
//!   Adaptor2d_Curve2d abstract interface.

pub mod line2d;

pub use line2d::Line2dAdaptor;
