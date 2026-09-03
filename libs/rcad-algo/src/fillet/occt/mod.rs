//! OCCT-aligned TKFillet translation (see AGENTS.md alignment methodology).
//!
//! Layer structure mirrors OCCT:
//!   - `chfi_ds`          — TKFillet/ChFiDS data structures
//!   - `chfi3d`           — ChFi3d_Builder / ChFi3d_FilBuilder /
//!                          ChFi3d_ChBuilder
//!   - `brep_fillet_api`  — BRepFilletAPI_MakeFillet / MakeChamfer
//!
//! Pending numerical-core boundaries are marked in each file with the OCCT
//! source line references.

pub mod brep_fillet_api;
pub mod chfi3d;
pub mod chfi_ds;

pub use brep_fillet_api::{
    edges_of_wire, explore_edges, explore_faces, explore_solids, explore_wires,
    BRepFilletAPIMakeChamfer, BRepFilletAPIMakeFillet,
};
pub use chfi3d::{ChFi3dBuilder, ChFi3dChBuilder, ChFi3dFilBuilder};
pub use chfi_ds::{
    ChFi3dFilletShape, ChFiDSSpineHandle, ChFiDS_State, ChFiDS_ChamfMethod, ChFiDS_ChamfMode,
    ChFiDS_ErrorStatus, ChFiDSChamfSpine, ChFiDSFilSpine, ChFiDSSpine, ChFiDSStripe, ChFiDSMap,
    ChFiDSStripeMap, SharedStripe,
};
