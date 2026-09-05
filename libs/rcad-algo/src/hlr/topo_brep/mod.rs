// OCCT HLRTopoBRep (TKHLR) — the Data structure storing per-face
// interference lists (internal/outlines/isolines) and per-edge split lists.
// Members land as sibling modules.

pub mod data;
pub mod face_data;
pub mod v_data;

pub use data::Data;
pub use face_data::FaceData;
pub use v_data::VData;

pub mod ds_filler;      // HLRTopoBRep_DSFiller
pub mod face_iso_liner; // HLRTopoBRep_FaceIsoLiner
pub mod out_liner;      // HLRTopoBRep_OutLiner
