//! OCCT Geom2dHatch package (TKGeomAlgo/Geom2dHatch) plus the HatchGen
//! package (TKGeomAlgo/HatchGen) — the leaf layer under the Geom2dHatch_
//! Hatcher main engine (not translated yet; it is the Stage-3d hard
//! prerequisite of HLRTopoBRep_FaceIsoLiner).
//!
//! 1:1 translations:
//!   - [`hatch_gen`] — HatchGen domain types: IntersectionType, ErrorStatus,
//!     IntersectionPoint, PointOnElement, PointOnHatching, Domain.
//!   - [`element`] — Geom2dHatch_Element.hxx/.cxx.
//!   - [`elements`] — Geom2dHatch_Elements.hxx/.cxx (the classifier's face
//!     explorer: probing Segment/OtherSegment + wire/edge traversal).
//!   - [`intersector`] — Geom2dHatch_Intersector.hxx/.cxx/.lxx (the
//!     Geom2dInt_GInter-derived line/edge intersector + LocalGeometry).
//!   - [`fclass2d`] — Geom2dHatch_FClass2dOfClassifier.hxx/.cxx (the
//!     TopClass_Classifier2d instantiation).
//!   - [`classifier`] — Geom2dHatch_Classifier.hxx/.cxx (the
//!     TopClass_FaceClassifier instantiation).
//!   - [`hatching`] — Geom2dHatch_Hatching.hxx/.cxx.
//!
//! OCCT template notes: TopClass_Classifier2d.pxx / TopClass_FaceClassifier
//! .pxx are C++ template headers; following the landed BRepClass precedent
//! their bodies are inlined into the concrete rcad classes with the
//! Geom2dHatch template arguments.

pub mod classifier;
pub mod element;
pub mod elements;
pub mod fclass2d;
pub mod hatch_gen;
pub mod hatching;
pub mod intersector;
pub mod hatcher;    // Geom2dHatch_Hatcher (main engine: Trim / ComputeDomains / Domain)

#[cfg(test)]
mod tests;

pub use classifier::Classifier;
pub use element::HatchElement;
pub use elements::HatchElements;
pub use fclass2d::FClass2dOfClassifier;
pub use hatch_gen::{Domain as HatchGenDomain, ErrorStatus, IntersectionType, PointOnElement};
pub use hatching::Hatching;
pub use intersector::HatchIntersector;

// HatchGen_PointOnHatching is re-exported under its OCCT-full name too (the
// bare `PointOnHatching` name stays inside hatch_gen next to its siblings).
pub use hatch_gen::PointOnHatching as HatchGenPointOnHatching;
