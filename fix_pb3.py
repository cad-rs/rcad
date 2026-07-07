"""Simple field access fix for SharedPB.
Only does: IDENTIFIER.FIELD -> IDENTIFIER.0.read().unwrap().FIELD
and IDENTIFIER.FIELD = -> IDENTIFIER.0.write().unwrap().FIELD =
Skips path containing .0.read() or .0.write()
"""
import re, sys

FIELDS = 'original_edge|new_edge|is_splittable|ext_paves|common_block_idx|shrunk_range|my_shrunk_box|pave1|pave2'
METHODS = 'indices|range|has_shrunk_data|is_to_update|update|append_ext_pave|remove_ext_pave|set_shrunk_data|set_shrunk_data_with_box|has_same_bounds|contains_parameter'

def fix_file(path):
    with open(path, 'rb') as f:
        raw = f.read()
    try:
        text = raw.decode('utf-8')
    except:
        return False
    orig = text

    # Pattern 1: IDENTIFIER.FIELD = (write)
    text = re.sub(r'(\w+)\.(' + FIELDS + r')\s*=',
        lambda m: m.group(0) if '.0.read' in m.group(1) or '.0.write' in m.group(1) else m.group(1) + '.0.write().unwrap().' + m.group(2) + ' =',
        text)

    # Pattern 2: IDENTIFIER.FIELD (read)
    text = re.sub(r'(\w+)\.(' + FIELDS + r')\b(?!\s*=)',
        lambda m: m.group(0) if '.0.read' in m.group(1) or '.0.write' in m.group(1) else m.group(1) + '.0.read().unwrap().' + m.group(2),
        text)

    # Pattern 3: IDENTIFIER.method(
    for m in ('update(', 'append_ext_pave(', 'remove_ext_pave(', 'set_shrunk_data(', 'set_shrunk_data_with_box('):
        text = re.sub(r'(\w+)\.' + re.escape(m),
            lambda x, m=m: x.group(0) if '.0.read' in x.group(1) or '.0.write' in x.group(1) else x.group(1) + '.0.write().unwrap().' + m,
            text)
    for m in ('indices()', 'range()', 'has_shrunk_data()', 'is_to_update()', 'has_same_bounds(', 'contains_parameter('):
        text = re.sub(r'(\w+)\.' + re.escape(m),
            lambda x, m=m: x.group(0) if '.0.read' in x.group(1) or '.0.write' in x.group(1) else x.group(1) + '.0.read().unwrap().' + m,
            text)

    if text != orig:
        with open(path, 'wb') as f:
            f.write(text.encode('utf-8'))
        return True
    return False

FILES = [
    'bopds/tools.rs', 'brep_algo_api/section.rs', 'boptools/mod.rs',
    'builder/builder_utils.rs', 'builder/edge_builders.rs', 'builder/filler.rs',
    'bop_occt_union.rs', 'pave_filler/intersection.rs',
    'pave_filler/make_blocks.rs', 'pave_filler/interf.rs',
    'pave_filler/paves.rs', 'pave_filler/ff_intersect.rs',
    'pave_filler/mod.rs', 'pipeline_dump.rs',
]

BASE = r'rcad/libs/rcad-algorithms/src'
for f in FILES:
    path = f'{BASE}/{f}'
    if fix_file(path):
        print(f'Fixed: {f}')
