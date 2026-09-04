//! OCCT HLRAlgo_EdgesBlock (TKHLR HLRAlgo package).
//!
//! 1:1 translation of `HLRAlgo_EdgesBlock.hxx` (L17-164) +
//! `HLRAlgo_EdgesBlock.cxx` (L192-197).

use rcad_kernel::topods::Orientation;

/// OCCT HLRAlgo_EdgesBlock::MinMaxIndices (hxx L45-80) — the packed 4x4
/// encoded bounding box used for shell-block culling.  The bit packing done
/// by `HLRAlgo::EncodeMinMax` / `DecodeMinMax` is reproduced bit-for-bit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MinMaxIndices {
    pub min: [i32; 8],
    pub max: [i32; 8],
}

impl MinMaxIndices {
    /// OCCT Minimize(theMinMaxIndices) — hxx L49-63.
    pub fn minimize(&mut self, the_min_max_indices: &MinMaxIndices) -> &mut Self {
        for a_i in 0..8 {
            if self.min[a_i] > the_min_max_indices.min[a_i] {
                self.min[a_i] = the_min_max_indices.min[a_i];
            }
            if self.max[a_i] > the_min_max_indices.max[a_i] {
                self.max[a_i] = the_min_max_indices.max[a_i];
            }
        }
        self
    }

    /// OCCT Maximize(theMinMaxIndices) — hxx L65-79.
    pub fn maximize(&mut self, the_min_max_indices: &MinMaxIndices) -> &mut Self {
        for a_i in 0..8 {
            if self.min[a_i] < the_min_max_indices.min[a_i] {
                self.min[a_i] = the_min_max_indices.min[a_i];
            }
            if self.max[a_i] < the_min_max_indices.max[a_i] {
                self.max[a_i] = the_min_max_indices.max[a_i];
            }
        }
        self
    }
}

impl Default for MinMaxIndices {
    fn default() -> Self {
        MinMaxIndices {
            min: [0; 8],
            max: [0; 8],
        }
    }
}

/// OCCT enum EMskFlags (hxx L149-156) — private flag bit masks.
const E_MASK_ORIENT: i32 = 15;
const E_MASK_OUTLINE: i32 = 16;
const E_MASK_INTERNAL: i32 = 32;
const E_MASK_DOUBLE: i32 = 64;
const E_MASK_ISOLINE: i32 = 128;

/// OCCT HLRAlgo_EdgesBlock — a set of edges structuring a wire: an array of
/// edge indices plus per-edge flag bits (orientation, outline, internal,
/// double, isoline) and the block min-max box.
#[derive(Debug, Clone)]
pub struct EdgesBlock {
    /// NCollection_Array1<int> myEdges (1, NbEdges) — stored 0-based.
    my_edges: Vec<i32>,
    my_flags: Vec<i32>,
    my_min_max: MinMaxIndices,
}

impl EdgesBlock {
    /// OCCT HLRAlgo_EdgesBlock(NbEdges) — cxx L192-196.
    pub fn new(nb_edges: usize) -> Self {
        EdgesBlock {
            my_edges: vec![0; nb_edges],
            my_flags: vec![0; nb_edges],
            my_min_max: MinMaxIndices::default(),
        }
    }

    /// OCCT NbEdges() — hxx L85 (Upper() of the 1-based array).
    pub fn nb_edges(&self) -> usize {
        self.my_edges.len()
    }

    /// OCCT Edge(I, EI) — hxx L87.
    pub fn set_edge(&mut self, i: usize, e_i: i32) {
        self.my_edges[i - 1] = e_i;
    }

    /// OCCT Edge(I) — hxx L89.
    pub fn edge(&self, i: usize) -> i32 {
        self.my_edges[i - 1]
    }

    /// OCCT Orientation(I, Or) — hxx L91-95.
    pub fn set_orientation(&mut self, i: usize, or_: Orientation) {
        self.my_flags[i - 1] &= !E_MASK_ORIENT;
        self.my_flags[i - 1] |= (or_ as i32) & E_MASK_ORIENT;
    }

    /// OCCT Orientation(I) — hxx L97-100.
    pub fn orientation(&self, i: usize) -> Orientation {
        let v = self.my_flags[i - 1] & E_MASK_ORIENT;
        match v {
            1 => Orientation::Reversed,
            2 => Orientation::Internal,
            3 => Orientation::External,
            _ => Orientation::Forward,
        }
    }

    /// OCCT OutLine(I) — hxx L102.
    pub fn out_line(&self, i: usize) -> bool {
        (self.my_flags[i - 1] & E_MASK_OUTLINE) != 0
    }

    /// OCCT OutLine(I, B) — hxx L104-110.
    pub fn set_out_line(&mut self, i: usize, b: bool) {
        if b {
            self.my_flags[i - 1] |= E_MASK_OUTLINE;
        } else {
            self.my_flags[i - 1] &= !E_MASK_OUTLINE;
        }
    }

    /// OCCT Internal(I) — hxx L112.
    pub fn internal(&self, i: usize) -> bool {
        (self.my_flags[i - 1] & E_MASK_INTERNAL) != 0
    }

    /// OCCT Internal(I, B) — hxx L114-120.
    pub fn set_internal(&mut self, i: usize, b: bool) {
        if b {
            self.my_flags[i - 1] |= E_MASK_INTERNAL;
        } else {
            self.my_flags[i - 1] &= !E_MASK_INTERNAL;
        }
    }

    /// OCCT Double(I) — hxx L122.
    pub fn double(&self, i: usize) -> bool {
        (self.my_flags[i - 1] & E_MASK_DOUBLE) != 0
    }

    /// OCCT Double(I, B) — hxx L124-130.
    pub fn set_double(&mut self, i: usize, b: bool) {
        if b {
            self.my_flags[i - 1] |= E_MASK_DOUBLE;
        } else {
            self.my_flags[i - 1] &= !E_MASK_DOUBLE;
        }
    }

    /// OCCT IsoLine(I) — hxx L132.
    pub fn iso_line(&self, i: usize) -> bool {
        (self.my_flags[i - 1] & E_MASK_ISOLINE) != 0
    }

    /// OCCT IsoLine(I, B) — hxx L134-140.
    pub fn set_iso_line(&mut self, i: usize, b: bool) {
        if b {
            self.my_flags[i - 1] |= E_MASK_ISOLINE;
        } else {
            self.my_flags[i - 1] &= !E_MASK_ISOLINE;
        }
    }

    /// OCCT UpdateMinMax(TotMinMax) — hxx L142.
    pub fn update_min_max(&mut self, tot_min_max: MinMaxIndices) {
        self.my_min_max = tot_min_max;
    }

    /// OCCT MinMax() — hxx L144.
    pub fn min_max(&mut self) -> &mut MinMaxIndices {
        &mut self.my_min_max
    }

    /// OCCT MinMax() const access (used where OCCT binds a const reference).
    pub fn min_max_ref(&self) -> &MinMaxIndices {
        &self.my_min_max
    }
}
