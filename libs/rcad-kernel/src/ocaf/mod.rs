//! OCAF (Open CASCADE Application Framework) data model.
//!
//! Analogous to OCCT's OCAF packages: XCAFDoc_ShapeTool (assembly hierarchy),
//! XCAFDoc_ColorTool (appearance), XCAFDimTolObjects (GD&T),
//! XCAFNoteObjects (annotations), OCAF TopoNaming (persistent naming).

pub mod annotation;
pub mod appearance;
pub mod assembly;
pub mod dim_tol;
pub mod naming;
pub mod persistent_naming;
