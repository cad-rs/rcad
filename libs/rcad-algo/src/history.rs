//! History tracking for boolean operations.

use rcad_kernel::topods;

/// Source operand a sub-shape originated from (OCCT myImages concept).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FaceOrigin {
    FromA(usize),
    FromB(usize),
    Generated,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EdgeOrigin {
    FromA(usize),
    FromB(usize),
    SplitFromA(usize),
    SplitFromB(usize),
    Generated,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VertexOrigin {
    FromA(usize),
    FromB(usize),
    Intersection,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShellOrigin {
    ShapeA,
    ShapeB,
    Generated,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolidOrigin {
    ShapeA,
    ShapeB,
    Generated,
}

/// History status for a source shape's result entities.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HistoryStatus {
    Modified,
    Generated,
    Deleted,
    Preserved,
}

/// Tracks per-entity origin and modification status for a single boolean operation.
#[derive(Debug, Clone, Default)]
pub struct HistoryTracker;

/// Minimal history entry.
#[derive(Debug, Clone)]
pub struct SourceShapeEntry {
    pub result: Vec<topods::Shape>,
}

/// Minimal history.
#[derive(Debug, Clone, Default)]
pub struct BooleanHistory {
    pub source_history: Vec<SourceShapeEntry>,
    pub vertex_origins: Vec<VertexOrigin>,
    pub edge_origins: Vec<EdgeOrigin>,
    pub face_origins: Vec<FaceOrigin>,
}
