"""Complete Arc<RwLock<PaveBlock>> refactoring.
Run: uv run python do_all.py
"""
import re, os

# ===== 1. pave.rs =====
with open('libs/rcad-algorithms/src/bopds/pave.rs', 'rb') as f: c = f.read().decode('utf-8')
c = c.replace('use std::collections::HashSet;', 'use std::collections::HashSet;\nuse std::sync::{Arc, RwLock};')
c = c.rstrip() + '\n\n/// OCCT-aligned: shared PaveBlock via `Arc<RwLock<PaveBlock>>`.\n#[derive(Debug, Clone)]\npub struct SharedPB(pub Arc<RwLock<PaveBlock>>);\n\nimpl SharedPB {\n    pub fn new(pb: PaveBlock) -> Self { SharedPB(Arc::new(RwLock::new(pb))) }\n}\n'
with open('libs/rcad-algorithms/src/bopds/pave.rs', 'wb') as f: f.write(c.encode('utf-8'))

# ===== 2. types.rs =====
with open('libs/rcad-algorithms/src/bopds/ds/types.rs', 'rb') as f: c = f.read().decode('utf-8')
c = c.replace('pub pave_blocks: Vec<PaveBlock>,', 'pub pave_blocks: Vec<crate::bopds::pave::SharedPB>,')
# IntersectionCurve methods
old = '''impl IntersectionCurve {\n  /// =OCCT-aligned: BOPDS_Curve::InitPaveBlock1 (lxx:85-92).\n  /// OCCT only pushes an empty PB to the list. PB vertices are set by\n  /// PutPavesOnCurve (ext_paves) -> Update(false) (sub-PBs from ext_paves).\n  pub fn init_pave_block1(&mut self) {\n  if self.pave_blocks.is_empty() {\n  self.pave_blocks.push(PaveBlock::new_curve_block());\n  }\n  }\n\n  /// =OCCT-aligned: BOPDS_Curve::ChangePaveBlock1 (lxx:96-100).\n  pub fn change_pave_block1(&mut self) -> Option<&mut PaveBlock> {\n  self.pave_blocks.first_mut()\n  }\n}'''
new = '''impl IntersectionCurve {\n  pub fn init_pave_block1(ds: &mut DS) -> usize {\n  let idx = ds.pave_blocks.len();\n  ds.pave_blocks.push(crate::bopds::pave::SharedPB::new(PaveBlock::new_curve_block()));\n  idx\n  }\n}'''
c = c.replace(old, new)
with open('libs/rcad-algorithms/src/bopds/ds/types.rs', 'wb') as f: f.write(c.encode('utf-8'))

# ===== 3. ds/mod.rs =====
with open('libs/rcad-algorithms/src/bopds/ds/mod.rs', 'rb') as f: c = f.read().decode('utf-8')
c = c.replace('use super::pave::{Pave, PaveBlock, NO_EDGE};', 'use super::pave::{Pave, PaveBlock, SharedPB, NO_EDGE};')
c = c.replace('pub fn change_pave_blocks(&mut self, edge_idx: usize) -> &mut Vec<PaveBlock>', 'pub fn change_pave_blocks(&mut self, edge_idx: usize) -> &mut Vec<SharedPB>')
c = c.replace('pub fn pave_blocks(&self, edge_idx: usize) -> &[PaveBlock]', 'pub fn pave_blocks(&self, edge_idx: usize) -> &[SharedPB]')
c = c.replace('let pb = PaveBlock::new(edge_idx, pv1, pv2);\n  // OCCT L469-471: ChangePaveBlocksPool', 'let pb = SharedPB::new(PaveBlock::new(edge_idx, pv1, pv2));\n  // OCCT L469-471: ChangePaveBlocksPool')
c = c.replace('self.pave_blocks.push(pb);', 'self.pave_blocks.push(SharedPB::new(pb));')
c = c.replace('pub fn is_common_block(&self, pb: &PaveBlock) -> bool {\n  pb.common_block_idx.is_some()\n  }\n\n  pub fn common_block(&self, pb: &PaveBlock) -> Option<&CommonBlock> {\n  pb.common_block_idx.and_then(|idx| self.common_blocks.get(idx))\n  }\n\n  pub fn common_block_mut(&mut self, pb: &PaveBlock) -> Option<&mut CommonBlock> {\n  pb.common_block_idx.and_then(|idx| self.common_blocks.get_mut(idx))\n  }\n\n  pub fn real_pave_block_edge(&self, edge_idx: usize, pb: &PaveBlock) -> Option<usize> {\n  let cb = self.common_block(pb)?;\n  let first_pb_idx = cb.pave_blocks().first()?.0;\n  self.edges.get(edge_idx)\n  .and_then(|e| e.pave_blocks.get(first_pb_idx))\n  .and_then(|pbr| pbr.new_edge)\n  }',
             'pub fn is_common_block(&self, pb: &SharedPB) -> bool {\n  pb.0.read().unwrap().common_block_idx.is_some()\n  }\n\n  pub fn common_block(&self, pb: &SharedPB) -> Option<&CommonBlock> {\n  pb.0.read().unwrap().common_block_idx.and_then(|idx| self.common_blocks.get(idx))\n  }\n\n  pub fn common_block_mut(&mut self, pb: &SharedPB) -> Option<&mut CommonBlock> {\n  pb.0.read().unwrap().common_block_idx.and_then(|idx| self.common_blocks.get_mut(idx))\n  }\n\n  pub fn real_pave_block_edge(&self, edge_idx: usize, pb: &SharedPB) -> Option<usize> {\n  let cb = self.common_block(pb)?;\n  let first_pb_idx = cb.pave_blocks().first()?.0;\n  let pbr = self.edges.get(edge_idx)?.pave_blocks.get(first_pb_idx)?;\n  pbr.0.read().unwrap().new_edge\n  }')

# Close edge append_ext_pave
c = c.replace('.append_ext_pave(Pave { vertex_idx: sv, param: tr0 });', '.0.write().unwrap().append_ext_pave(Pave { vertex_idx: sv, param: tr0 });')
c = c.replace('.append_ext_pave(Pave { vertex_idx: sv, param: tr1 });', '.0.write().unwrap().append_ext_pave(Pave { vertex_idx: sv, param: tr1 });')

# Remove old duplicate CommonBlock/RealPaveBlock methods (take &PaveBlock)
lines = c.split('\n')
new_lines = []
skip = False
for line in lines:
    if '///  ?OCCT-aligned: BOPDS_DS::CommonBlock (hxx:192-193)' in line:
        skip = True
    if not skip:
        new_lines.append(line)
    if skip and '.and_then(|pbr| pbr.new_edge)' in line:
        skip = False
c = '\n'.join(new_lines)

# build_edge_images loop
c = c.replace('for pb in &edge.pave_blocks {\n  let sub_ei = pb.new_edge.unwrap_or(ei);', 
              'for spb in &edge.pave_blocks {\n  let pbb = spb.0.read().unwrap();\n  let sub_ei = pbb.new_edge.unwrap_or(ei);')

# update_pave_blocks_with_sd_vertices
c = c.replace('for pb in &mut edge.pave_blocks {', 'for spb in &mut edge.pave_blocks {\n  let mut pb = spb.0.write().unwrap();')
c = c.replace('for pb in &mut self.pave_blocks {', 'for spb in &mut self.pave_blocks {\n  let mut pb = spb.0.write().unwrap();')
# Close the write guard scopes by adding an extra }
# Actually this won't work reliably. Let me use a different approach.
# Instead, I'll match and replace the exact blocks.

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

# refine_face_info closures
c = c.replace('  pave_blocks.get(pb_idx).map_or(false, |pb| {\n  pb.pave1.vertex_idx != pb.pave2.vertex_idx', 
              '  pave_blocks.get(pb_idx).map_or(false, |pb| {\n  pb.0.read().unwrap().pave1.vertex_idx != pb.0.read().unwrap().pave2.vertex_idx')
c = c.replace('  on_pb.original_edge == pb.original_edge\n  && on_pb.pave1.vertex_idx == pb.pave1.vertex_idx\n  && on_pb.pave2.vertex_idx == pb.pave2.vertex_idx',
              '  on_pb.0.read().unwrap().original_edge == pb.0.read().unwrap().original_edge\n  && on_pb.0.read().unwrap().pave1.vertex_idx == pb.0.read().unwrap().pave1.vertex_idx\n  && on_pb.0.read().unwrap().pave2.vertex_idx == pb.0.read().unwrap().pave2.vertex_idx')

with open('libs/rcad-algorithms/src/bopds/ds/mod.rs', 'wb') as f: f.write(c.encode('utf-8'))

print('Infra done. Now fixing consumer files...')

# ===== 4. ALL consumer files: add .0.read().unwrap() / .0.write().unwrap() =====
BASE = 'libs/rcad-algorithms/src'
FILES = [
    'bopds/tools.rs', 'brep_algo_api/section.rs', 'boptools/mod.rs',
    'builder/builder_utils.rs', 'builder/edge_builders.rs', 'builder/filler.rs',
    'bop_occt_union.rs',
    'pave_filler/intersection.rs', 'pave_filler/make_blocks.rs', 
    'pave_filler/interf.rs', 'pave_filler/paves.rs',
    'pave_filler/ff_intersect.rs', 'pave_filler/mod.rs',
]
ALL_PB_FIELDS = 'pave1|pave2|original_edge|new_edge|shrunk_range|is_splittable|common_block_idx|ext_paves|my_shrunk_box'
ALL_PB_METHODS = ['indices()', 'range()', 'is_to_update()', 'has_shrunk_data()', 'update(', 'append_ext_pave(','remove_ext_pave(','has_same_bounds(','contains_parameter(']

for fname in FILES:
    path = BASE+'/'+fname
    if not os.path.exists(path):
        continue
    with open(path, 'rb') as f: c = f.read().decode('utf-8')
    orig = c

    # 4a. Array-index patterns (100% safe)
    for ctx in ['ds.pave_blocks', 'self\\.ds\\.pave_blocks', 'edges\\[\\w+\\]\\.pave_blocks', 'self\\.ds\\.edges\\[\\w+\\]\\.pave_blocks']:
        # reads
        c = re.sub(r'(' + ctx + r'\[\w+\])\.(' + ALL_PB_FIELDS + r')\b(?!\s*=)',
                   lambda m: m.group(0) if '.0.read' in m.group(0) or '.0.write' in m.group(0) else m.group(1)+'.0.read().unwrap().'+m.group(2), c)
        # writes
        c = re.sub(r'(' + ctx + r'\[\w+\])\.(' + ALL_PB_FIELDS + r')\s*=',
                   lambda m: m.group(0) if '.0.read' in m.group(0) or '.0.write' in m.group(0) else m.group(1)+'.0.write().unwrap().'+m.group(2)+' =', c)
        # method calls
        for m in ALL_PB_METHODS:
            c = re.sub(r'(' + ctx + r'\[\w+\])\.' + re.escape(m),
                       lambda mx, mm=m: mx.group(0) if '.0.read' in mx.group(0) or '.0.write' in mx.group(0) else mx.group(1)+'.0.write().unwrap().'+mm, c)

    # 4b. ic.pave_blocks patterns (IntersectionCurve)
    for ctx in ['ic\\.pave_blocks', 'self\\.ds\\.intersection_curves\\[\\w+\\]\\.pave_blocks']:
        c = re.sub(r'(' + ctx + r'\[\w+\])\.(original_edge|pave1|pave2|new_edge)\b(?!\s*=)',
                   lambda m: m.group(0) if '.0.read' in m.group(0) else m.group(1)+'.0.read().unwrap().'+m.group(2), c)
        for m in ALL_PB_METHODS:
            c = re.sub(r'(' + ctx + r'\[\w+\])\.' + re.escape(m),
                       lambda mx, mm=m: mx.group(0) if '.0.read' in mx.group(0) else mx.group(1)+'.0.write().unwrap().'+mm, c)

    # 4c. Local SharedPB variables: pb.FIELD -> pb.0.read().unwrap().FIELD
    # This catches ALL remaining patterns. Over-replaced cases fixed in step 5.
    c = re.sub(r'(\w+)\.(' + ALL_PB_FIELDS + r')\b(?!\s*=)',
               lambda m: m.group(0) if '.0.read' in m.group(1) or '.0.write' in m.group(1) else m.group(1)+'.0.read().unwrap().'+m.group(2), c)
    c = re.sub(r'(\w+)\.(' + ALL_PB_FIELDS + r')\s*=',
               lambda m: m.group(0) if '.0.read' in m.group(1) or '.0.write' in m.group(1) else m.group(1)+'.0.write().unwrap().'+m.group(2)+' =', c)
    for m in ALL_PB_METHODS:
        c = re.sub(r'(\w+)\.' + re.escape(m),
                   lambda mx, mm=m: mx.group(0) if '.0.read' in mx.group(1) or '.0.write' in mx.group(1) else mx.group(1)+'.0.write().unwrap().'+mm, c)

    # 4d. init_pave_block1 / change_pave_block1 calls
    c = c.replace('.init_pave_block1()', '')
    c = c.replace('.change_pave_block1()', '.pave_blocks.first().map(|spb| spb.0.write().unwrap())')
    c = c.replace('IntersectionCurve::init_pave_block1(self.ds)', '')
    # Fix the empty leftover from removing init_pave_block1
    c = c.replace('let _ = ;', '')

    if c != orig:
        with open(path, 'wb') as f: f.write(c.encode('utf-8'))
        print(f'  {fname}')

# ===== 5. Fix over-replaced PaveBlock variables =====
# These are local variables that are genuine PaveBlock (not SharedPB)
OVER_FIXES = [
    # a_pb from update() result (Vec<PaveBlock>)
    ('pave_filler/make_blocks.rs', 'a_pb.0.read().unwrap().indices()', 'a_pb.indices()'),
    ('pave_filler/make_blocks.rs', 'a_pb.0.read().unwrap().range()', 'a_pb.range()'),
    ('pave_filler/make_blocks.rs', 'a_pb.0.write().unwrap().indices()', 'a_pb.indices()'),
    ('pave_filler/make_blocks.rs', 'a_pb.0.write().unwrap().range()', 'a_pb.range()'),
    # sub_pb from update() result
    ('pave_filler/mod.rs', 'sub_pb.0.write().unwrap()', 'sub_pb'),
    ('pave_filler/mod.rs', 'sub_pb.0.read().unwrap()', 'sub_pb'),
    ('pave_filler/mod.rs', 'pb_clone.0.write().unwrap()', 'pb_clone'),
    ('pave_filler/mod.rs', 'pb_clone.0.read().unwrap()', 'pb_clone'),
    # ShrunkRange variables
    ('pave_filler/paves.rs', 'sr_cell.0.write().unwrap()', 'sr_cell'),
    ('pave_filler/paves.rs', 'sr_cell.0.read().unwrap()', 'sr_cell'),
    ('pave_filler/paves.rs', 'p_cell.0.write().unwrap()', 'p_cell'),
    ('pave_filler/paves.rs', 'p_cell.0.read().unwrap()', 'p_cell'),
    # ic variable in ff_intersect (IntersectionCurve, not SharedPB)
    ('pave_filler/ff_intersect.rs', 'ic.0.read().unwrap()', 'ic'),
    ('pave_filler/ff_intersect.rs', 'ic.0.write().unwrap()', 'ic'),
    # pb in ff_intersect (genuine PaveBlock, not SharedPB)
    ('pave_filler/ff_intersect.rs', 'pb.0.write().unwrap().original_edge', 'pb.original_edge'),
    ('pave_filler/ff_intersect.rs', 'pb.0.read().unwrap().original_edge', 'pb.original_edge'),
    # sub_pb in ff_intersect
    ('pave_filler/ff_intersect.rs', 'sub_pb.0.write().unwrap()', 'sub_pb'),
]

for fname, old, new in OVER_FIXES:
    path = BASE+'/'+fname
    with open(path, 'rb') as f: c = f.read().decode('utf-8')
    if old in c:
        c = c.replace(old, new)
        with open(path, 'wb') as f: f.write(c.encode('utf-8'))
        print(f'  OVER-FIX {fname}: {old[:40]}')

print('DONE')
