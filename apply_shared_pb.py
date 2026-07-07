"""Apply Arc<RwLock<PaveBlock>> refactoring in a single deterministic pass.
Usage: python apply_shared_pb.py
"""
import re

# Step 1: pave.rs - add imports + SharedPB struct
with open('libs/rcad-algorithms/src/bopds/pave.rs', 'rb') as f:
    c = f.read().decode('utf-8')
c = c.replace('use std::collections::HashSet;', 'use std::collections::HashSet;\nuse std::sync::{Arc, RwLock};')
c = c.rstrip() + '''

/// OCCT-aligned: shared PaveBlock via `Arc<RwLock<PaveBlock>>`.
#[derive(Debug, Clone)]
pub struct SharedPB(pub Arc<RwLock<PaveBlock>>);

impl SharedPB {
    pub fn new(pb: PaveBlock) -> Self { SharedPB(Arc::new(RwLock::new(pb))) }
}
'''
with open('libs/rcad-algorithms/src/bopds/pave.rs', 'wb') as f:
    f.write(c.encode('utf-8'))
print('[1/4] pave.rs done')

# Step 2: types.rs - field types + IC methods
with open('libs/rcad-algorithms/src/bopds/ds/types.rs', 'rb') as f:
    c = f.read().decode('utf-8')

# Change all three pave_blocks field types (DSEdge, IntersectionCurve, DS)
# Each one is `pub pave_blocks: Vec<PaveBlock>,` -> SharedPB
c = c.replace('pub pave_blocks: Vec<PaveBlock>,', 'pub pave_blocks: Vec<crate::bopds::pave::SharedPB>,')

# Replace IntersectionCurve methods (instance -> static)
old = '''impl IntersectionCurve {
  /// =OCCT-aligned: BOPDS_Curve::InitPaveBlock1 (lxx:85-92).
  /// OCCT only pushes an empty PB to the list. PB vertices are set by
  /// PutPavesOnCurve (ext_paves) -> Update(false) (sub-PBs from ext_paves).
  pub fn init_pave_block1(&mut self) {
  if self.pave_blocks.is_empty() {
  self.pave_blocks.push(PaveBlock::new_curve_block());
  }
  }

  /// =OCCT-aligned: BOPDS_Curve::ChangePaveBlock1 (lxx:96-100).
  pub fn change_pave_block1(&mut self) -> Option<&mut PaveBlock> {
  self.pave_blocks.first_mut()
  }
}'''
new = '''impl IntersectionCurve {
  /// Creates a curve PaveBlock in the global pool.
  pub fn init_pave_block1(ds: &mut DS) -> usize {
  let idx = ds.pave_blocks.len();
  ds.pave_blocks.push(crate::bopds::pave::SharedPB::new(PaveBlock::new_curve_block()));
  idx
  }
}'''
c = c.replace(old, new)

# Also fix test code that constructs IntersectionCurve with pave_blocks: Vec::new()
c = c.replace('pave_blocks: Vec::new(),', 'pave_blocks: Vec::new(),')

with open('libs/rcad-algorithms/src/bopds/ds/types.rs', 'wb') as f:
    f.write(c.encode('utf-8'))
print('[2/4] types.rs done')

# Step 3: ds/mod.rs - SharedPB import + method updates
with open('libs/rcad-algorithms/src/bopds/ds/mod.rs', 'rb') as f:
    c = f.read().decode('utf-8')

c = c.replace('use super::pave::{Pave, PaveBlock, NO_EDGE};', 'use super::pave::{Pave, PaveBlock, SharedPB, NO_EDGE};')

# change_pave_blocks return type
c = c.replace('pub fn change_pave_blocks(&mut self, edge_idx: usize) -> &mut Vec<PaveBlock>',
              'pub fn change_pave_blocks(&mut self, edge_idx: usize) -> &mut Vec<SharedPB>')
# pave_blocks return type
c = c.replace('pub fn pave_blocks(&self, edge_idx: usize) -> &[PaveBlock]',
              'pub fn pave_blocks(&self, edge_idx: usize) -> &[SharedPB]')

# init_pave_blocks_for_edge: wrap in SharedPB
c = c.replace('let pb = PaveBlock::new(edge_idx, pv1, pv2);', 'let pb = SharedPB::new(PaveBlock::new(edge_idx, pv1, pv2));')

# Closed-edge append_ext_pave block (handle both old and new formatting)
old_close = '''  if sv == ev {
  // OCCT L479: aPB->AppendExtPave(aP1)  ?first endpoint
  self.edges[edge_idx].pave_blocks[0]
  .append_ext_pave(Pave { vertex_idx: sv, param: tr0 });
  // OCCT L481: aPB->AppendExtPave(aP2)  ?second endpoint
  // (fence dedups by vertex_idx, so second push is accepted
  //  because vertex_idx differs from the first push? No, same
  //  vertex  ?the OCCT fence check also uses vertex_idx.
  //  The second AppendExtPave is rejected by fence in both
  //  implementations; OCCT still writes it for form clarity.)
  self.edges[edge_idx].pave_blocks[0]
  .append_ext_pave(Pave { vertex_idx: sv, param: tr1 });
  }'''
new_close = '''  if sv == ev {
  // SharedPB: append via Arc<RwLock>
  let spb = &self.edges[edge_idx].pave_blocks[0];
  spb.0.write().unwrap().append_ext_pave(Pave { vertex_idx: sv, param: tr0 });
  spb.0.write().unwrap().append_ext_pave(Pave { vertex_idx: sv, param: tr1 });
  }'''
c = c.replace(old_close, new_close)

# allocate_pave_block: wrap in SharedPB
c = c.replace('self.pave_blocks.push(pb);', 'self.pave_blocks.push(SharedPB::new(pb));')

# CommonBlock accessors
old_cb = '''  ///  ?OCCT-aligned: BOPDS_DS::IsCommonBlock (hxx:188).
  /// Returns true if the PaveBlock belongs to a CommonBlock.
  pub fn is_common_block(&self, pb: &PaveBlock) -> bool {
  pb.common_block_idx.is_some()
  }

  ///  ?OCCT-aligned: BOPDS_DS::CommonBlock (hxx:192-193).
  /// Returns a reference to the CommonBlock for a PaveBlock.
  pub fn common_block(&self, pb: &PaveBlock) -> Option<&CommonBlock> {
  pb.common_block_idx.and_then(|idx| self.common_blocks.get(idx))
  }

  ///  ?OCCT-aligned: BOPDS_DS::CommonBlock (hxx:192-193) =mutable.
  pub fn common_block_mut(&mut self, pb: &PaveBlock) -> Option<&mut CommonBlock> {
  pb.common_block_idx.and_then(|idx| self.common_blocks.get_mut(idx))
  }

  ///  ?OCCT-aligned: BOPDS_DS::RealPaveBlock (BOPDS_DS.cxx L658-663).
  /// If the PaveBlock belongs to a CommonBlock, returns the edge index of
  /// the first PaveBlock in that block (the "real" edge). Otherwise returns
  /// the given PaveBlock's new_edge.
  pub fn real_pave_block_edge(&self, edge_idx: usize, pb: &PaveBlock) -> Option<usize> {
  let cb = self.common_block(pb)?;
  let first_pb_idx = cb.pave_blocks().first()?.0;
  self.edges.get(edge_idx)
  .and_then(|e| e.pave_blocks.get(first_pb_idx))
  .and_then(|pbr| pbr.new_edge)
  }'''
new_cb = '''  pub fn is_common_block(&self, pb: &SharedPB) -> bool {
  pb.0.read().unwrap().common_block_idx.is_some()
  }

  pub fn common_block(&self, pb: &SharedPB) -> Option<&CommonBlock> {
  pb.0.read().unwrap().common_block_idx.and_then(|idx| self.common_blocks.get(idx))
  }

  pub fn common_block_mut(&mut self, pb: &SharedPB) -> Option<&mut CommonBlock> {
  pb.0.read().unwrap().common_block_idx.and_then(|idx| self.common_blocks.get_mut(idx))
  }

  pub fn real_pave_block_edge(&self, edge_idx: usize, pb: &SharedPB) -> Option<usize> {
  let cb = self.common_block(pb)?;
  let first_pb_idx = cb.pave_blocks().first()?.0;
  let pbr = self.edges.get(edge_idx)?.pave_blocks.get(first_pb_idx)?;
  pbr.0.read().unwrap().new_edge
  }'''
c = c.replace(old_cb, new_cb)

# build_edge_images loop
old_bei = '''  for ei in 0..n_edges {
  let edge = &self.edges[ei];
  for pb in &edge.pave_blocks {
  let sub_ei = pb.new_edge.unwrap_or(ei);
  if sub_ei < self.edges.len() {
  self.my_images[ei].push(sub_ei);
  self.my_origins.push(ei);
  }
  }
  }'''
new_bei = '''  for ei in 0..n_edges {
  let edge = &self.edges[ei];
  for spb in &edge.pave_blocks {
  let pb = spb.0.read().unwrap();
  let sub_ei = pb.new_edge.unwrap_or(ei);
  if sub_ei < self.edges.len() {
  self.my_images[ei].push(sub_ei);
  self.my_origins.push(ei);
  }
  }
  }'''
c = c.replace(old_bei, new_bei)

# update_pave_blocks_with_sd_vertices loops
old_upd = '''  // Apply replacement to all PaveBlocks on edges.
  for edge in &mut self.edges {
  for pb in &mut edge.pave_blocks {
  if let Some(&rep) = replace.get(&pb.pave1.vertex_idx) {
  pb.pave1.vertex_idx = rep;
  }
  if let Some(&rep) = replace.get(&pb.pave2.vertex_idx) {
  pb.pave2.vertex_idx = rep;
  }
  }
  }
  // Apply replacement to global PaveBlocks pool.
  for pb in &mut self.pave_blocks {
  if let Some(&rep) = replace.get(&pb.pave1.vertex_idx) {
  pb.pave1.vertex_idx = rep;
  }
  if let Some(&rep) = replace.get(&pb.pave2.vertex_idx) {
  pb.pave2.vertex_idx = rep;
  }
  }'''
new_upd = '''  // Apply SD replacement through SharedPB.
  for edge in &mut self.edges {
  for spb in &mut edge.pave_blocks {
  let mut pb = spb.0.write().unwrap();
  if let Some(&rep) = replace.get(&pb.pave1.vertex_idx) {
  pb.pave1.vertex_idx = rep;
  }
  if let Some(&rep) = replace.get(&pb.pave2.vertex_idx) {
  pb.pave2.vertex_idx = rep;
  }
  }
  }
  // Apply to global pool (curve and orphan PBs).
  for spb in &mut self.pave_blocks {
  let mut pb = spb.0.write().unwrap();
  if let Some(&rep) = replace.get(&pb.pave1.vertex_idx) {
  pb.pave1.vertex_idx = rep;
  }
  if let Some(&rep) = replace.get(&pb.pave2.vertex_idx) {
  pb.pave2.vertex_idx = rep;
  }
  }'''
c = c.replace(old_upd, new_upd)

# refine_face_info_on closure (pb.pave1/pb.pave2)
c = c.replace(
    '  pave_blocks.get(pb_idx).map_or(false, |pb| {\n  pb.pave1.vertex_idx != pb.pave2.vertex_idx',
    '  pave_blocks.get(pb_idx).map_or(false, |pb| {\n  pb.0.read().unwrap().pave1.vertex_idx != pb.0.read().unwrap().pave2.vertex_idx'
)
# refine_face_info_in closure (on_pb/pb original_edge/pave1/pave2)
c = c.replace(
    '  on_pb.original_edge == pb.original_edge\n  && on_pb.pave1.vertex_idx == pb.pave1.vertex_idx\n  && on_pb.pave2.vertex_idx == pb.pave2.vertex_idx',
    '  on_pb.0.read().unwrap().original_edge == pb.0.read().unwrap().original_edge\n  && on_pb.0.read().unwrap().pave1.vertex_idx == pb.0.read().unwrap().pave1.vertex_idx\n  && on_pb.0.read().unwrap().pave2.vertex_idx == pb.0.read().unwrap().pave2.vertex_idx'
)

with open('libs/rcad-algorithms/src/bopds/ds/mod.rs', 'wb') as f:
    f.write(c.encode('utf-8'))
print('[3/4] ds/mod.rs done')

# Step 4: tools.rs - field access through SharedPB
with open('libs/rcad-algorithms/src/bopds/tools.rs', 'rb') as f:
    c = f.read().decode('utf-8')

# All pb.FIELD accesses in tools.rs are on SharedPB (from ds.pave_blocks[] or edge.pave_blocks[])
fields = ['pave1', 'pave2', 'original_edge', 'common_block_idx', 'new_edge']
for f in fields:
    c = re.sub(r'(\w+)\.' + f + r'\b(?!\s*=)', lambda m, ff=f: m.group(0) if '.0.read' in m.group(1) or '.0.write' in m.group(1) else m.group(1)+'.0.read().unwrap().'+ff, c)
    c = re.sub(r'(\w+)\.' + f + r'\s*=', lambda m, ff=f: m.group(0) if '.0.read' in m.group(1) or '.0.write' in m.group(1) else m.group(1)+'.0.write().unwrap().'+ff+' =', c)

with open('libs/rcad-algorithms/src/bopds/tools.rs', 'wb') as f:
    f.write(c.encode('utf-8'))
print('[4/4] tools.rs done')

print('ALL DONE - now run cargo check')
