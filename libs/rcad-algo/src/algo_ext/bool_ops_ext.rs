// Boolean-ops extension: solid/shell extraction, n-ary partition (splitter).
// Migrated from the legacy rcad-algorithms crate (bool_ops_ext/mod.rs);
// boolean kernels are the current rcad-algo `common`/`cut` entry points.

use crate::bop::algo::BooleanError;
use glam::DVec3;

// Topods-native extraction helpers (same module scope via include):
// compact_brep_topods, extract_solids_topods, extract_shells_topods
#[path = "topods_ext.rs"]
mod topods_ext;
pub use topods_ext::{compact_brep_topods, extract_shells_topods, extract_solids_topods};

/// Extract each solid from a BRep as a separate self-contained BRep.
pub fn extract_solids(brep: &rcad_kernel::BRep) -> Vec<rcad_kernel::BRep> {
    extract_solids_topods(brep)
}

/// Extract each shell from a BRep as a separate self-contained BRep.
pub fn extract_shells(brep: &rcad_kernel::BRep) -> Vec<rcad_kernel::BRep> {
    extract_shells_topods(brep)
}

/// Partition objects by tools.
///
/// OCCT BRepAlgoAPI_Splitter / BOPAlgo_Splitter (BOPAlgo_Splitter.cxx L54-93):
/// objects and tools are combined into ONE PaveFiller argument list, and the
/// Builder runs the GF pipeline (BOPAlgo_Builder::PerformInternal1) with
/// myArguments = objects only.  BOPAlgo_Builder::BuildResult (BOPAlgo_Builder_1.cxx
/// L130-168) adds the split parts (images) of the OBJECTS into the result,
/// excluding the split parts of the Tools.  The result is a compound of the
/// object images; the cell boundaries (split faces/edges/vertices) are SHARED
/// TShape instances across adjacent cells (single PaveFiller run).
pub fn n_ary_partition(
    objects: &[rcad_kernel::BRep],
    tools: &[rcad_kernel::BRep],
) -> Result<Vec<rcad_kernel::BRep>, BooleanError> {
    // OCCT BRepAlgoAPI_Splitter::Build L42-46: check argument counts.
    if objects.is_empty() || (objects.len() + tools.len()) < 2 {
        return Err(BooleanError::TooFewArguments);
    }
    // OCCT BOPAlgo_Splitter::Perform: aLS = myArguments + myTools combined.
    let result = crate::bop::brep_algo_api::splitter(objects, tools)?;
    // OCCT bsplit: the result shape is a compound of the split parts of the
    // objects.  Return it as the single cell (the compound preserves the
    // shared boundary TShapes across cells).
    Ok(vec![result])
}

pub fn make_face_half_space(
    plane: &rcad_kernel::geom::Plane,
    bbox: &[DVec3; 2],
    normal_side: bool,
) -> rcad_kernel::topods::BRep {
    let [bmin, bmax] = *bbox;
    let diag = bmax - bmin;
    let margin = diag.length().max(1.0) * 2.0;
    let n = if normal_side { plane.normal } else { -plane.normal }.normalize();
    let abs = n.abs();
    let candidate = if abs.x <= abs.y && abs.x <= abs.z {
        DVec3::X
    } else if abs.y <= abs.z {
        DVec3::Y
    } else {
        DVec3::Z
    };
    let u = n.cross(candidate).normalize();
    let v = n.cross(u);
    let origin = if normal_side {
        plane.origin - u * (margin / 2.0) - v * (margin / 2.0)
    } else {
        plane.origin - u * (margin / 2.0) - v * (margin / 2.0) - n * margin
    };
    rcad_modeling::make_box_brep(origin, u, v, margin, margin, margin)
        .expect("make_face_half_space: box construction failed")
}
