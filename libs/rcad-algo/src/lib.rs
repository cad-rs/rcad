//! rcad-algo: OCCT TKBO (boolean ops) + TKTopAlgo (topological algorithms)
//! + TKGeomAlgo (geometric algorithms).
//!
//! Depends only on rcad-kernel (TKMath + TKGeomBase) and rcad-brep (TKBRep).

pub mod bop;
pub mod geomalgo;
pub mod topalgo;
pub mod algo_ext;
pub mod helix;
pub mod fillet;
pub mod shhealing;

// Re-export boolean operation API at top level
pub use crate::bop::brep_algo_api::{boolean_op, boolean_op_with_retry, common, cut, fuse};
pub use crate::bop::algo::BooleanOpType;

// Re-export compatibility helpers (legacy rcad-algorithms surface used by the
// generated OCCT boolean grid tests).
pub use crate::algo_ext::{
    bool_ops_ext, brep_algo, brep_tools,
    count_edges, count_faces, count_shells, count_vertices, count_wires,
    extract_shells, extract_solids, extrude_polygon_solid, move_bspline2_point,
    move_bspline2_tangent, move_bspline3_point, move_bspline3_tangent, n_ary_partition,
    make_fillet_edge, restrict_to_bspline, revolve_polygon_solid, total_edge_length,
    total_surface_area,
    total_volume,
};
pub use crate::algo_ext::extrude_profile::{extrude_profile_solid, ProfileSegment};

// Re-export healing chain (STEP export needs it; legacy rcad-algorithms surface).
pub use crate::algo_ext::healing::{HealingMode, HealingOptions, HealingReport, analyze_and_heal};
pub use crate::algo_ext::brep_check::analyze_wire_issues;

// Re-export from rcad-kernel
pub use rcad_kernel::base::bnd_lib;
pub use rcad_kernel::base::extrema;
pub use rcad_kernel::math::math_poly::solve_quartic;
pub use rcad_kernel::math::lin::inverse_3x3;
pub use rcad_kernel::math::opt::golden_section_max;
