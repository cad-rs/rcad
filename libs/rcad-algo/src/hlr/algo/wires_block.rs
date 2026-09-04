//! OCCT HLRAlgo_WiresBlock (TKHLR HLRAlgo package).
//!
//! 1:1 translation of `HLRAlgo_WiresBlock.hxx` (L35-65, fully inline).

use super::edges_block::{EdgesBlock, MinMaxIndices};

/// OCCT HLRAlgo_WiresBlock — a set of Blocks used by the DataStructure to
/// structure the edges: an array of EdgesBlock handles plus the block
/// min-max box.  The OCCT handle array (`NCollection_Array1<handle>`,
/// 1-based) maps to a Vec of the block structs.
#[derive(Debug, Clone)]
pub struct WiresBlock {
    my_wires: Vec<EdgesBlock>,
    my_min_max: MinMaxIndices,
}

impl WiresBlock {
    /// OCCT HLRAlgo_WiresBlock(NbWires) — hxx L39-42.
    pub fn new(nb_wires: usize) -> Self {
        WiresBlock {
            my_wires: (0..nb_wires).map(|_| EdgesBlock::new(0)).collect(),
            my_min_max: MinMaxIndices::default(),
        }
    }

    /// OCCT NbWires() — hxx L44.
    pub fn nb_wires(&self) -> usize {
        self.my_wires.len()
    }

    /// OCCT Set(I, W) — hxx L46.
    pub fn set(&mut self, i: usize, w: EdgesBlock) {
        self.my_wires[i - 1] = w;
    }

    /// OCCT Wire(I) — hxx L48-51.
    pub fn wire(&mut self, i: usize) -> &mut EdgesBlock {
        &mut self.my_wires[i - 1]
    }

    /// OCCT UpdateMinMax(theMinMaxes) — hxx L53-56.
    pub fn update_min_max(&mut self, the_min_maxes: MinMaxIndices) {
        self.my_min_max = the_min_maxes;
    }

    /// OCCT MinMax() — hxx L58.
    pub fn min_max(&mut self) -> &mut MinMaxIndices {
        &mut self.my_min_max
    }
}
