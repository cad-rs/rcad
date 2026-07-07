"""Fix SharedPB field access: pb.FIELD -> pb.0.read().unwrap().FIELD
and pb.FIELD = x -> pb.0.write().unwrap().FIELD = x

Only targets files with errors. Avoids double-wrapping.
"""
import re, sys

# All PaveBlock fields accessed through SharedPB
FIELDS = [
    'original_edge', 'new_edge', 'is_splittable', 'ext_paves',
    'common_block_idx', 'shrunk_range', 'my_shrunk_box',
    'pave1', 'pave2', 'curve', 'pcurve_on_a', 'pcurve_on_b',
]
# All PaveBlock methods called through SharedPB  
METHODS = [
    'indices', 'range', 'has_shrunk_data', 'is_to_update',
    'update', 'append_ext_pave', 'remove_ext_pave',
    'set_shrunk_data', 'set_shrunk_data_with_box',
    'has_same_bounds', 'contains_parameter',
]

FIELD_PAT = '|'.join(FIELDS)
METHOD_PAT = '|'.join(METHODS)

# Files to fix (from error output)
FIX_FILES = [
    'bopds/tools.rs', 'brep_algo_api/section.rs', 'boptools/mod.rs',
    'builder/builder_utils.rs', 'builder/edge_builders.rs', 'builder/filler.rs',
    'bop_occt_union.rs', 'pave_filler/intersection.rs',
    'pave_filler/make_blocks.rs', 'pave_filler/interf.rs',
    'pave_filler/paves.rs', 'pave_filler/ff_intersect.rs',
    'pave_filler/mod.rs', 'pipeline_dump.rs',
    'bopds/ds/mod.rs',  # remaining field errors in mod.rs
]

BASE = r'rcad/libs/rcad-algorithms/src'

def fix_file(path, is_ds_mod=False):
    with open(path, 'rb') as f:
        raw = f.read()
    try:
        text = raw.decode('utf-8')
    except:
        return False
    orig = text

    # 1. FIELD access: foo.FIELD -> foo.0.read().unwrap().FIELD
    # Skip if already wrapped with .0.read() or .0.write()
    text = re.sub(
        r'(\w+(?:\.\w+)*)\.(?:' + FIELD_PAT + r')(?!\s*=)(?!(?:\.0\.read|\.0\.write))',
        lambda m: m.group(0) if any(x in m.group(0) for x in ['.0.read()', '.0.write()', '.read().unwrap()', '.write().unwrap()', 'self.pave_blocks', 'ds.pave_blocks']) else m.group(1) + '.0.read().unwrap().' + m.group(2),
        text
    )
    # 2. FIELD = value: foo.FIELD = -> foo.0.write().unwrap().FIELD =
    text = re.sub(
        r'(\w+(?:\.\w+)*)\.(?:' + FIELD_PAT + r')\s*=',
        lambda m: m.group(0) if any(x in m.group(0) for x in ['.0.read()', '.0.write()', '.read().unwrap()', '.write().unwrap()']) else m.group(1) + '.0.write().unwrap().' + m.group(2) + ' =',
        text
    )
    # 3. Method calls
    text = re.sub(
        r'(\w+(?:\.\w+)*)\.(?:' + METHOD_PAT + r')\(',
        lambda m: m.group(0) if any(x in m.group(0) for x in ['.0.read()', '.0.write()']) else m.group(1) + '.0.write().unwrap().' + m.group(2) + '(',
        text
    )

    if text != orig:
        with open(path, 'wb') as f:
            f.write(text.encode('utf-8'))
        return True
    return False

if __name__ == '__main__':
    for f in FIX_FILES:
        path = f'{BASE}/{f}'
        if fix_file(path):
            print(f'Fixed: {f}')
        else:
            print(f'No changes: {f}')
