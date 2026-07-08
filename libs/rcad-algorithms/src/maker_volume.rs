//! Build solids from reusable split cells.
//!
//! # OCCT Reference — BOPAlgo_MakerVolume
//!
//! This module provides functionality analogous to OCCT `BOPAlgo_MakerVolume`
//! (`BOPAlgo_MakerVolume.cxx` / `BOPAlgo_MakerVolume.hxx`).
//!
//! ## OCCT Pipeline
//!
//! 1. **`CollectFaces()`** — Collects all input faces from the argument shapes.
//! 2. **`MakeBox()`** — Creates a bounding box around all collected faces.
//! 3. **`BuildSolids()`** — Uses `BOPAlgo_BuilderSolid` to compute 3D regions
//!    from the space between the box and the collected faces.
//! 4. **`RemoveBox()`** — Removes the box from each solid to yield the result.
//! 5. **`FillInternalShapes()`** — Detects and fills internal voids in the result.
//!
//! ## RCAD Approach (this file)
//!
//! RCAD does **not** use the bounding-box/BuilderSolid pipeline. Instead it works
//! with **pre-split cell solids** and fuses them with `general_fuse`:
//!
//! | Step | OCCT BOPAlgo_MakerVolume | RCAD MakerVolume |
//! |------|--------------------------|------------------|
//! | Input | `CollectFaces`: flattens shapes into face set | Pre-split cell solids from caller |
//! | Partitioning | `MakeBox`: bounding box around everything | No bounding box |
//! | Solid building | `BOPAlgo_BuilderSolid` on box-faces | `general_fuse` pairwise |
//! | Box removal | `RemoveBox`: subtracts box from result | Not applicable |
//! | Internal voids | `FillInternalShapes` | Not implemented |
//!
//! Despite these differences, the **interface is equivalent**: callers register
//! pre-split cell solids and then assemble a final solid from a region mask, an
//! explicit cell index list, or a [`CellExpr`] boolean expression.

use std::collections::HashSet;

use rcad_kernel::topods;

use crate::{BooleanError, GeneralFuseHistory, general_fuse, general_fuse_with_history};

/// Error type for MakerVolume-style solid assembly.
#[derive(Debug)]
pub enum MakerVolumeError {
    /// No cells were provided.
    EmptyInput,
    /// The region selection resolved to no cells.
    EmptySelection,
    /// A region mask did not match the number of registered cells.
    RegionMaskLengthMismatch { expected: usize, got: usize },
    /// A requested cell index was out of bounds.
    InvalidCellIndex { index: usize, count: usize },
    /// Boolean assembly failed.
    Boolean(BooleanError),
    /// All cells are too small for the requested operation.
    CellsTooSmall,
}

impl std::fmt::Display for MakerVolumeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "empty input"),
            Self::EmptySelection => write!(f, "empty cell selection"),
            Self::RegionMaskLengthMismatch { expected, got } => {
                write!(f, "region mask length mismatch: expected {expected}, got {got}")
            }
            Self::InvalidCellIndex { index, count } => {
                write!(f, "invalid cell index {index}; available cells: 0..{count}")
            }
            Self::Boolean(source) => write!(f, "boolean assembly failed: {source}"),
            Self::CellsTooSmall => write!(f, "all cells are too small"),
        }
    }
}

impl std::error::Error for MakerVolumeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Boolean(source) => Some(source),
            _ => None,
        }
    }
}

impl From<BooleanError> for MakerVolumeError {
    fn from(value: BooleanError) -> Self {
        Self::Boolean(value)
    }
}

/// Report describing which cells were selected for a MakerVolume assembly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MakerVolumeSelection {
    /// Total number of input cells registered on the builder.
    pub input_cell_count: usize,
    /// Unique selected cell indices, in stable encounter order.
    pub selected_cell_indices: Vec<usize>,
}

/// Reusable solid assembler over precomputed split cells.
#[derive(Debug, Clone, Default)]
pub struct MakerVolume {
    cells: Vec<topods::BRep>,
}

impl MakerVolume {
    /// Create an empty MakerVolume builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a MakerVolume builder from precomputed cells.
    pub fn from_cells(cells: Vec<topods::BRep>) -> Self {
        Self { cells }
    }

    /// Add one cell and return its index.
    pub fn add_cell(&mut self, cell: topods::BRep) -> usize {
        self.cells.push(cell);
        self.cells.len() - 1
    }

    /// Number of registered cells.
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Build a solid from all registered cells.
    pub fn build_all(&self) -> Result<topods::BRep, MakerVolumeError> {
        let indices: Vec<usize> = (0..self.cells.len()).collect();
        self.build_from_indices(&indices)
    }

    /// Build a solid from a boolean region mask.
    pub fn build_from_region_mask(&self, region_mask: &[bool]) -> Result<topods::BRep, MakerVolumeError> {
        let selection = self.selection_from_region_mask(region_mask)?;
        self.build_from_indices(&selection.selected_cell_indices)
    }

    /// Build a solid and per-step history from a region mask.
    pub fn build_from_region_mask_with_history(
        &self,
        region_mask: &[bool],
    ) -> Result<(topods::BRep, GeneralFuseHistory), MakerVolumeError> {
        let selection = self.selection_from_region_mask(region_mask)?;
        self.build_from_indices_with_history(&selection.selected_cell_indices)
    }

    /// Build a solid from an explicit cell index list.
    pub fn build_from_indices(&self, indices: &[usize]) -> Result<topods::BRep, MakerVolumeError> {
        let parts = self.selected_cells(indices)?;
        // Bridge: general_fuse still uses old BRep
        let old_parts: Vec<_> = parts.iter().map(|p| rcad_kernel::BRep::from_topods(p)).collect();
        let old_result = general_fuse(&old_parts)?;
        Ok(old_result.to_topods())
    }

    /// Build a solid and per-step history from an explicit cell index list.
    pub fn build_from_indices_with_history(
        &self,
        indices: &[usize],
    ) -> Result<(topods::BRep, GeneralFuseHistory), MakerVolumeError> {
        let parts = self.selected_cells(indices)?;
        let old_parts: Vec<_> = parts.iter().map(|p| rcad_kernel::BRep::from_topods(p)).collect();
        let (old_result, hist) = general_fuse_with_history(&old_parts)?;
        Ok((old_result.to_topods(), hist))
    }

    /// Convert a region mask into a validated selection report.
    pub fn selection_from_region_mask(
        &self,
        region_mask: &[bool],
    ) -> Result<MakerVolumeSelection, MakerVolumeError> {
        if self.cells.is_empty() {
            return Err(MakerVolumeError::EmptyInput);
        }
        if region_mask.len() != self.cells.len() {
            return Err(MakerVolumeError::RegionMaskLengthMismatch {
                expected: self.cells.len(),
                got: region_mask.len(),
            });
        }

        let selected_cell_indices: Vec<usize> = region_mask
            .iter()
            .enumerate()
            .filter_map(|(index, enabled)| enabled.then_some(index))
            .collect();

        if selected_cell_indices.is_empty() {
            return Err(MakerVolumeError::EmptySelection);
        }

        Ok(MakerVolumeSelection {
            input_cell_count: self.cells.len(),
            selected_cell_indices,
        })
    }

    fn selected_cells(&self, indices: &[usize]) -> Result<Vec<topods::BRep>, MakerVolumeError> {
        if self.cells.is_empty() {
            return Err(MakerVolumeError::EmptyInput);
        }

        let unique_indices = unique_cell_indices(indices);
        if unique_indices.is_empty() {
            return Err(MakerVolumeError::EmptySelection);
        }

        let mut parts = Vec::with_capacity(unique_indices.len());
        for index in unique_indices {
            let cell = self
                .cells
                .get(index)
                .cloned()
                .ok_or(MakerVolumeError::InvalidCellIndex {
                    index,
                    count: self.cells.len(),
                })?;
            parts.push(cell);
        }
        Ok(parts)
    }
}

/// Convenience helper: assemble a solid from a region mask.
pub fn make_solid_from_region(
    cells: &[rcad_kernel::BRep],
    region_mask: &[bool],
) -> Result<rcad_kernel::BRep, MakerVolumeError> {
    let cells_t: Vec<_> = cells.iter().map(|c| c.to_topods()).collect();
    let result_t = MakerVolume::from_cells(cells_t).build_from_region_mask(region_mask)?;
    Ok(rcad_kernel::BRep::from_topods(&result_t))
}

/// Convenience helper: assemble a solid from a region mask and report history.
pub fn make_solid_from_region_with_history(
    cells: &[rcad_kernel::BRep],
    region_mask: &[bool],
) -> Result<(rcad_kernel::BRep, GeneralFuseHistory), MakerVolumeError> {
    let cells_t: Vec<_> = cells.iter().map(|c| c.to_topods()).collect();
    let (result_t, hist) = MakerVolume::from_cells(cells_t).build_from_region_mask_with_history(region_mask)?;
    Ok((rcad_kernel::BRep::from_topods(&result_t), hist))
}

/// Convenience helper: assemble a solid from explicit cell indices.
pub fn make_solid_from_cell_indices(
    cells: &[rcad_kernel::BRep],
    indices: &[usize],
) -> Result<rcad_kernel::BRep, MakerVolumeError> {
    let cells_t: Vec<_> = cells.iter().map(|c| c.to_topods()).collect();
    let result_t = MakerVolume::from_cells(cells_t).build_from_indices(indices)?;
    Ok(rcad_kernel::BRep::from_topods(&result_t))
}

fn unique_cell_indices(indices: &[usize]) -> Vec<usize> {
    let mut seen = HashSet::with_capacity(indices.len());
    let mut out = Vec::with_capacity(indices.len());
    for &index in indices {
        if seen.insert(index) {
            out.push(index);
        }
    }
    out
}

