//! TKBO — Boolean Operation algorithms (OCCT TKBO toolkit).
//!
//! | Submodule   | OCCT Package  | Description                        |
//! |-------------|---------------|------------------------------------|
//! | algo        | BOPAlgo       | PaveFiller, Builder                |
//! | ds          | BOPDS         | Data structures (DS, Pave, etc.)   |
//! | tools       | BOPTools      | BVH, box tree                      |
//! | algo_api    | BRepAlgoAPI   | High-level fuse/common/cut API     |
//! | int_tools   | IntTools      | Edge-edge, edge-face, face-face intersection |

pub mod algo;
pub mod ds;
pub mod tools;
pub mod brep_algo_api;
pub mod int_tools;
