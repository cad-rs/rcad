//! OCCT HLRAlgo small poly data structures (TKHLR HLRAlgo package).
//!
//! 1:1 translations:
//! - `HLRAlgo_PolyMask.hxx` (L20-35) — the triangle/face flag bit masks.
//! - `HLRAlgo_TriangleData.hxx` (L26-31).
//! - `HLRAlgo_PolyInternalSegment.hxx` (L26-31).
//! - `HLRAlgo_PolyInternalNode.hxx` (L30-69).
//! - `HLRAlgo_PolyHidingData.hxx` (L26-71).

use glam::{DVec2, DVec3};

/// OCCT enum HLRAlgo_PolyMask (PolyMask.hxx L20-35).
pub mod poly_mask {
    pub const E_MSK_OUT_LIN1: i32 = 1;
    pub const E_MSK_OUT_LIN2: i32 = 2;
    pub const E_MSK_OUT_LIN3: i32 = 4;
    pub const E_MSK_GR_A_LIN1: i32 = 8;
    pub const E_MSK_GR_A_LIN2: i32 = 16;
    pub const E_MSK_GR_A_LIN3: i32 = 32;
    pub const F_MSK_BACK: i32 = 64;
    pub const F_MSK_SIDE: i32 = 128;
    pub const F_MSK_HIDING: i32 = 256;
    pub const F_MSK_FLAT: i32 = 512;
    pub const F_MSK_ON_OUT_L: i32 = 1024;
    pub const F_MSK_OR_BACK: i32 = 2048;
    pub const F_MSK_FR_BACK: i32 = 4096;
}

/// OCCT struct HLRAlgo_TriangleData (TriangleData.hxx L26-31).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct TriangleData {
    pub node1: i32,
    pub node2: i32,
    pub node3: i32,
    pub flags: i32,
}

/// OCCT struct HLRAlgo_PolyInternalSegment (PolyInternalSegment.hxx L26-31).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct PolyInternalSegment {
    pub lst_sg1: i32,
    pub lst_sg2: i32,
    pub nxt_sg1: i32,
    pub nxt_sg2: i32,
    pub conex1: i32,
    pub conex2: i32,
}

/// OCCT HLRAlgo_PolyInternalNode::NodeIndices (PolyInternalNode.hxx L33-36).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct NodeIndices {
    pub nd_sg: i32,
    pub flag: i32,
    pub edg1: i32,
    pub edg2: i32,
}

/// OCCT HLRAlgo_PolyInternalNode::NodeData (PolyInternalNode.hxx L38-50).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct NodeData {
    pub point: DVec3,
    pub normal: DVec3,
    pub uv: DVec2,
    pub pcu1: f64,
    pub pcu2: f64,
    pub scal: f64,
}

/// OCCT HLRAlgo_PolyInternalNode (PolyInternalNode.hxx L30-69) — to Update
/// OutLines.  The OCCT handle (Standard_Transient) maps to the plain struct.
#[derive(Debug, Clone, Copy, Default)]
pub struct PolyInternalNode {
    my_indices: NodeIndices,
    my_data: NodeData,
}

impl PolyInternalNode {
    /// OCCT HLRAlgo_PolyInternalNode() — hxx L52-58.
    pub fn new() -> Self {
        PolyInternalNode {
            my_indices: NodeIndices::default(),
            my_data: NodeData::default(),
        }
    }

    /// OCCT Indices() — hxx L60.
    pub fn indices(&mut self) -> &mut NodeIndices {
        &mut self.my_indices
    }

    /// OCCT Data() — hxx L62.
    pub fn data(&mut self) -> &mut NodeData {
        &mut self.my_data
    }

    /// Immutable accesses (const reads in the algorithms).
    pub fn indices_ref(&self) -> &NodeIndices {
        &self.my_indices
    }

    /// Immutable data access (see `indices_ref`).
    pub fn data_ref(&self) -> &NodeData {
        &self.my_data
    }
}

/// OCCT HLRAlgo_PolyHidingData::TriangleIndices (PolyHidingData.hxx L31-34).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct HidingTriangleIndices {
    pub index: i32,
    pub min: i32,
    pub max: i32,
}

/// OCCT HLRAlgo_PolyHidingData::PlaneT (PolyHidingData.hxx L36-45).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct HidingPlaneT {
    pub normal: DVec3,
    pub d: f64,
}

/// OCCT HLRAlgo_PolyHidingData (PolyHidingData.hxx L26-71) — data structure
/// of a set of hiding triangles.
#[derive(Debug, Clone, Copy, Default)]
pub struct PolyHidingData {
    my_indices: HidingTriangleIndices,
    my_plane: HidingPlaneT,
}

impl PolyHidingData {
    /// OCCT HLRAlgo_PolyHidingData() = default — hxx L47.
    pub fn new() -> Self {
        PolyHidingData::default()
    }

    /// OCCT Set(Index, Minim, Maxim, A, B, C, D) — hxx L49-62.
    #[allow(clippy::too_many_arguments)]
    pub fn set(&mut self, index: i32, minim: i32, maxim: i32, a: f64, b: f64, c: f64, d: f64) {
        self.my_indices.index = index;
        self.my_indices.min = minim;
        self.my_indices.max = maxim;
        self.my_plane.normal = DVec3::new(a, b, c);
        self.my_plane.d = d;
    }

    /// OCCT Indices() — hxx L64.
    pub fn indices(&mut self) -> &mut HidingTriangleIndices {
        &mut self.my_indices
    }

    /// OCCT Plane() — hxx L66.
    pub fn plane(&mut self) -> &mut HidingPlaneT {
        &mut self.my_plane
    }

    /// Immutable accesses (const reads in the algorithms).
    pub fn indices_ref(&self) -> &HidingTriangleIndices {
        &self.my_indices
    }

    /// Immutable plane access (see `indices_ref`).
    pub fn plane_ref(&self) -> &HidingPlaneT {
        &self.my_plane
    }
}
