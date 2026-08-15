//! Compatibility helpers migrated from the legacy `rcad-algorithms` crate.
//!
//! These modules back the OCCT boolean grid tests (`tests/occt`) that were
//! generated against the old crate. They depend only on `rcad-kernel`,
//! `rcad-brep` and `rcad-modeling` and use the current `topods::BRep` pool API.

pub mod bspline_edit;
pub mod tolerance;
pub mod brep_check;
pub mod shape_analysis;
pub mod brep_repair;
pub mod healing;
pub mod shape_custom;
pub mod features;
pub mod revolve;
pub mod extrude_profile;
pub mod geom_populate;
pub mod fillet;

mod topods_ext;
pub mod bool_ops_ext;
pub mod brep_algo;
pub mod brep_tools;

// Re-export the healing chain (used by rcad-step STEP export).
pub use healing::{
    HealingMode, HealingOptions, HealingReport, analyze_and_heal, analyze_wire_issues,
};
pub use shape_custom::restrict_to_bspline;
pub use features::{extrude_polygon_solid, revolve_polygon_solid};
pub use extrude_profile::{extrude_profile_solid, ProfileSegment};
pub use fillet::make_fillet_edge;

// Re-export BSpline edit helpers at the module top level.
pub use bspline_edit::{
    move_bspline2_point, move_bspline2_tangent, move_bspline3_point, move_bspline3_tangent,
};

// Re-export extraction helpers.
pub use bool_ops_ext::{extract_shells, extract_solids, n_ary_partition};
pub use brep_algo::total_edge_length;

// Re-export surface/volume/count helpers from rcad-kernel / rcad-brep with the
// names the generated OCCT tests use.
pub use rcad_kernel::base::gprop::surface::surface_area as total_surface_area;
pub use rcad_kernel::base::gprop::volume::volume as total_volume;
pub use rcad_brep::tools::{
    count_edges, count_faces, count_shells, count_vertices, count_wires,
};
