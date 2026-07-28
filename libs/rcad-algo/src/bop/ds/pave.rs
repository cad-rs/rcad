// OCCT BOPDS_Pave + BOPDS_PaveBlock 1:1 translation.
//
// BOPDS_Pave.hxx      ?vertex-on-edge parametric point
// BOPDS_PaveBlock.hxx ?edge segment between two paves

use std::sync::{Arc, RwLock};

// ===
// BOPDS_Pave
// ===
#[derive(Debug, Clone, Copy)]
pub struct Pave {
    pub vertex_idx: usize,
    pub param: f64,
}

impl Pave {
    pub fn new(index: usize, parameter: f64) -> Self {
        Pave { vertex_idx: index, param: parameter }
    }
    pub fn index(&self) -> usize { self.vertex_idx }
    pub fn set_index(&mut self, idx: usize) { self.vertex_idx = idx; }
    pub fn parameter(&self) -> f64 { self.param }
    pub fn set_parameter(&mut self, p: f64) { self.param = p; }
    pub fn is_less(&self, other: &Pave) -> bool { self.param < other.param }
    pub fn is_equal(&self, other: &Pave) -> bool { self.vertex_idx == other.vertex_idx && self.param == other.param }
}

impl PartialEq for Pave {
    fn eq(&self, other: &Self) -> bool { self.is_equal(other) }
}
impl PartialOrd for Pave {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.param.partial_cmp(&other.param)
    }
}

// ===
// NO_EDGE sentinel ?section edge PaveBlocks have no original edge
// ===
pub const NO_EDGE: usize = usize::MAX;

// ===
// SharedPaveBlock ?Arc<RwLock<PaveBlock>> (Standard_Transient equivalent)
// ===
#[derive(Debug, Clone)]
pub struct SharedPB(pub Arc<RwLock<PaveBlock>>);

impl SharedPB {
    pub fn new(pb: PaveBlock) -> Self { SharedPB(Arc::new(RwLock::new(pb))) }
    pub fn read(&self) -> std::sync::RwLockReadGuard<'_, PaveBlock> { self.0.read().unwrap() }
    pub fn write(&self) -> std::sync::RwLockWriteGuard<'_, PaveBlock> { self.0.write().unwrap() }
}

// ===
// BOPDS_PaveBlock
// ===
#[derive(Debug, Clone)]
pub struct PaveBlock {
    // BOPDS_PaveBlock.hxx fields
    pub edge: usize,            // myEdge ?index of the edge
    pub original_edge: usize,   // myOriginalEdge
    pub pave1: Pave,            // myPave1
    pub pave2: Pave,            // myPave2
    pub ext_paves: Vec<Pave>,   // myExtPaves (NCollection_List)
    pub ts1: f64,               // myTS1 ?shrunk range start
    pub ts2: f64,               // myTS2 ?shrunk range end
    pub common_block_idx: Option<usize>, // rcad: link to common block (OCCT uses DataMap)
}

impl PaveBlock {
    pub fn new(edge_idx: usize, p1: Pave, p2: Pave) -> Self {
        PaveBlock {
            edge: edge_idx,
            original_edge: edge_idx,
            pave1: p1,
            pave2: p2,
            ext_paves: Vec::new(),
            ts1: 0.0,
            ts2: 0.0,
            common_block_idx: None,
        }
    }

    // -- Edge --
    pub fn set_edge(&mut self, e: usize) { self.edge = e; }
    pub fn edge(&self) -> usize { self.edge }
    pub fn has_edge(&self) -> bool { self.edge != usize::MAX }

    // -- OriginalEdge --
    pub fn set_original_edge(&mut self, e: usize) { self.original_edge = e; }
    pub fn original_edge(&self) -> usize { self.original_edge }
    pub fn is_split_edge(&self) -> bool { self.edge != self.original_edge }

    // -- Pave1/Pave2 --
    pub fn set_pave1(&mut self, p: Pave) { self.pave1 = p; }
    pub fn pave1(&self) -> &Pave { &self.pave1 }
    pub fn set_pave2(&mut self, p: Pave) { self.pave2 = p; }
    pub fn pave2(&self) -> &Pave { &self.pave2 }

    // -- Range --
    pub fn range(&self) -> (f64, f64) { (self.pave1.param, self.pave2.param) }

    // -- Indices --
    pub fn indices(&self) -> (usize, usize) { (self.pave1.vertex_idx, self.pave2.vertex_idx) }

    // -- HasSameBounds --
    pub fn has_same_bounds(&self, other: &PaveBlock) -> bool {
        self.pave1.vertex_idx == other.pave1.vertex_idx
            && self.pave2.vertex_idx == other.pave2.vertex_idx
    }

    // -- ExtPaves --
    pub fn is_to_update(&self) -> bool { !self.ext_paves.is_empty() }
    pub fn append_ext_pave(&mut self, p: Pave) { self.ext_paves.push(p); }
    pub fn append_ext_pave1(&mut self, p: Pave) { self.ext_paves.push(p); }
    pub fn remove_ext_pave(&mut self, vert_num: usize) {
        self.ext_paves.retain(|p| p.vertex_idx != vert_num);
    }
    pub fn ext_paves(&self) -> &[Pave] { &self.ext_paves }
    pub fn change_ext_paves(&mut self) -> &mut Vec<Pave> { &mut self.ext_paves }

    // -- ShrunkData --
    pub fn set_shrunk_data(&mut self, ts1: f64, ts2: f64, _is_splittable: bool) {
        self.ts1 = ts1;
        self.ts2 = ts2;
    }
    pub fn shrunk_data(&self) -> (f64, f64, bool) { (self.ts1, self.ts2, self.is_splittable()) }
    pub fn has_shrunk_data(&self) -> bool { self.ts1 != 0.0 || self.ts2 != 0.0 }
    pub fn is_splittable(&self) -> bool {
        let (t1, t2) = self.range();
        (t2 - t1).abs() > 1e-15
    }

    /// OCCT: PaveBlock on a curve edge with default paves ?not split (section edge).
    pub fn new_curve_block() -> Self {
        PaveBlock {
            edge: NO_EDGE, original_edge: NO_EDGE,
            pave1: Pave::new(0, 0.0), pave2: Pave::new(0, 0.0),
            ext_paves: Vec::new(), ts1: 0.0, ts2: 0.0, common_block_idx: None,
        }
    }

    pub fn contains_parameter(&self, prm: f64, tol: f64, ind: &mut usize) -> bool {
        for pv in &self.ext_paves {
            if (pv.param - prm).abs() <= tol {
                *ind = pv.vertex_idx;
                return true;
            }
        }
        false
    }
}

// ===
// BOPDS_PaveBlock::Update ?split into sub-blocks using ext paves
// ===
pub fn update_pave_block(pb: &PaveBlock, lp: &mut Vec<SharedPB>, flag: bool) {
    if !pb.is_to_update() { return; }
    let mut paves: Vec<Pave> = pb.ext_paves.clone();
    paves.sort_by(|a, b| a.param.partial_cmp(&b.param).unwrap_or(std::cmp::Ordering::Equal));
    if flag {
        paves.insert(0, pb.pave1);
        paves.push(pb.pave2);
    }
    for i in 0..paves.len().saturating_sub(1) {
        let p1 = paves[i];
        let p2 = paves[i + 1];
        if (p2.param - p1.param).abs() < 1e-15 { continue; }
        let new_pb = SharedPB::new(PaveBlock {
            edge: pb.edge,
            original_edge: pb.original_edge,
            pave1: p1,
            pave2: p2,
            ext_paves: Vec::new(),
            ts1: 0.0,
            ts2: 0.0,
            common_block_idx: pb.common_block_idx,
        });
        lp.push(new_pb);
    }
}

/// Pair of PaveBlocks connected to the same vertex.
/// OCCT: used in EE intersection for vertex connection map.
#[derive(Debug, Clone)]
pub struct CoupleOfPaveBlocks {
    pub pb1: SharedPB,
    pub pb2: SharedPB,
}

impl CoupleOfPaveBlocks {
    pub fn new(pb1: SharedPB, pb2: SharedPB) -> Self {
        CoupleOfPaveBlocks { pb1, pb2 }
    }
}
