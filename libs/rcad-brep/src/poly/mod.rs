//! OCCT Poly package — the TKBRep-facing triangulation module (plan A0).
//!
//! The Poly VALUE types (`Poly_Triangle`, `Poly_MeshPurpose`,
//! `Poly_TriangulationParameters`, `Poly_Triangulation`, `Poly_Polygon3D`,
//! `Poly_PolygonOnTriangulation`) are translated 1:1 in
//! [`rcad_kernel::poly`] — OCCT places TKMath/Poly BELOW TKBRep in its own
//! toolkit dependency graph (BRep_TFace carries `Poly_Triangulation`
//! handles), and the rcad crate DAG mirrors that direction
//! (rcad-brep depends on rcad-kernel). They are re-exported here so the
//! `rcad_brep::poly` path is the single home of the Poly layer.
//!
//! This module hosts the pieces that need BRep/pool context or are pure
//! algorithms on the triangulation:
//! - [`connect`]: `Poly_Connect` — the adjacency walker (1:1).
//!
//! Not translated (plan BREPMESH_RMSH_PLAN §1.4 non-goals):
//! `Poly_CoherentTriangulation` family, `Poly_MakeLoops`,
//! `Poly_MergeNodesTool`; `Poly_Polygon2D` (2D polygon-on-surface mount,
//! outside A0 scope).

pub mod connect;

pub use connect::PolyConnect;

// Re-export the kernel-side Poly value types so `rcad_brep::poly` mirrors the
// OCCT Poly package surface (the types live in rcad-kernel; see module docs).
pub use rcad_kernel::poly::{
    PolyMeshPurpose, PolyPolygon3D, PolyPolygonOnTriangulation, PolyTriangle,
    PolyTriangulation, PolyTriangulationParameters,
    POLY_MESH_PURPOSE_ACTIVE, POLY_MESH_PURPOSE_ANY_FALLBACK, POLY_MESH_PURPOSE_CALCULATION,
    POLY_MESH_PURPOSE_LOADED, POLY_MESH_PURPOSE_NONE, POLY_MESH_PURPOSE_PRESENTATION,
    POLY_MESH_PURPOSE_USER,
};
