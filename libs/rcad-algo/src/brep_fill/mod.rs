//! OCCT TKBool/BRepFill package — sweep/loft classes consumed on demand
//! (plan D3/D5: per-class introduction, new top-level brep_fill/ dir).
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepFill/

pub mod compatible_wires;
pub mod compatible_wires_b;
pub mod generator;
pub mod generator_b;
pub mod offset_wire;
pub mod offset_wire_b;
pub mod brep_fill_draft;
pub mod brep_fill_trim_edge_tool;
pub mod brep_fill_pipe;
pub mod brep_fill_filling;
pub mod brep_fill_filling_b;
pub mod brep_fill_pipe_shell;
pub mod brep_fill_pipe_shell_b;
pub mod brep_fill_evolved;
pub mod brep_fill_evolved_b;
pub mod brep_fill_evolved_c;
pub mod brep_fill_evolved_d;
pub mod brep_fill_offset_ancestors;
pub mod brep_fill_trim_surface_tool;
pub mod brep_fill_axe;
pub mod brep_fill_location_law;
pub mod brep_fill_section_placement;
pub mod brep_fill_section_law;
pub mod brep_fill_shape_law;
pub mod brep_fill_edge3d_law;
pub mod brep_fill_nsections;
pub mod brep_fill_sweep;
