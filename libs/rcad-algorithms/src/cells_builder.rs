//! Reusable split-cell expression builder.
//!
//! This provides a lightweight CellsBuilder API analogous to OCCT
//! `BOPAlgo_CellsBuilder`: callers register reusable cell solids and then
//! evaluate boolean expressions over those cells.

use crate::{BooleanError, BooleanOpType, boolean_op};
use rcad_kernel::BRep;
use rcad_kernel::topods;

/// Boolean expression over registered cells.
#[derive(Debug, Clone)]
pub enum CellExpr {
    /// Reference a registered cell by index.
    Cell(usize),
    /// Union of two expressions.
    Union(Box<CellExpr>, Box<CellExpr>),
    /// Intersection of two expressions.
    Intersection(Box<CellExpr>, Box<CellExpr>),
    /// Difference of two expressions: left - right.
    Difference(Box<CellExpr>, Box<CellExpr>),
    /// XOR: symmetric difference (A xor B = (A - B) ∪ (B - A))
    Xor(Box<CellExpr>, Box<CellExpr>),
}

/// Error type for [`CellsBuilder`].
#[derive(Debug)]
pub enum CellsBuilderError {
    /// Referenced cell index does not exist.
    InvalidCellIndex { index: usize, count: usize },
    /// Underlying boolean operation failed.
    Boolean(BooleanError),
}

impl std::fmt::Display for CellsBuilderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidCellIndex { index, count } => {
                write!(f, "invalid cell index {index}; available cells: 0..{count}")
            }
            Self::Boolean(e) => write!(f, "boolean operation failed: {e}"),
        }
    }
}

impl std::error::Error for CellsBuilderError {}

impl From<BooleanError> for CellsBuilderError {
    fn from(value: BooleanError) -> Self {
        Self::Boolean(value)
    }
}

/// Reusable cell container and expression evaluator.
#[derive(Debug, Clone, Default)]
pub struct CellsBuilder {
    cells: Vec<BRep>,
}

impl CellsBuilder {
    /// Create an empty builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a builder from precomputed cells.
    pub fn from_cells(cells: Vec<BRep>) -> Self {
        Self { cells }
    }

    /// Add one cell and return its index.
    pub fn add_cell(&mut self, cell: BRep) -> usize {
        self.cells.push(cell);
        self.cells.len() - 1
    }

    /// Number of registered cells.
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Evaluate a boolean expression over registered cells.
    pub fn evaluate(&self, expr: &CellExpr) -> Result<BRep, CellsBuilderError> {
        self.eval_rec(expr)
    }

    fn eval_rec(&self, expr: &CellExpr) -> Result<BRep, CellsBuilderError> {
        match expr {
            CellExpr::Cell(i) => self
                .cells
                .get(*i)
                .cloned()
                .ok_or(CellsBuilderError::InvalidCellIndex {
                    index: *i,
                    count: self.cells.len(),
                }),
            CellExpr::Union(a, b) => self.eval_bin(BooleanOpType::Union, a, b),
            CellExpr::Intersection(a, b) => self.eval_bin(BooleanOpType::Intersection, a, b),
            CellExpr::Difference(a, b) => self.eval_bin(BooleanOpType::Difference, a, b),
            CellExpr::Xor(a, b) => {
                // XOR: (A - B) ∪ (B - A)
                let a_min_b = self.eval_bin(BooleanOpType::Difference, a, b)?;
                let b_min_a = self.eval_bin(BooleanOpType::Difference, b, a)?;
                let t = boolean_op(BooleanOpType::Union, &a_min_b, &b_min_a)?; Ok(rcad_kernel::BRep::from_topods(&t))
            }
        }
    }

    fn eval_bin(
        &self,
        op: BooleanOpType,
        a: &CellExpr,
        b: &CellExpr,
    ) -> Result<BRep, CellsBuilderError> {
        let left = self.eval_rec(a)?;
        let right = self.eval_rec(b)?;
        let t = boolean_op(op, &left, &right)?; Ok(rcad_kernel::BRep::from_topods(&t))
    }
}


