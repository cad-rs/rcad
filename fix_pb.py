"""Fix SharedPB field access: replace `pb.FIELD` with `pb.read().unwrap().FIELD`
for reads and `pb.FIELD =` with `pb.write().unwrap().FIELD =` for writes.

Handles all variable names referencing SharedPB instances across the codebase.
Runs on single files passed as arguments.
"""
import re, sys

# Fields on PaveBlock that need access via SharedPB
ACCESSED_FIELDS = [
    'original_edge', 'new_edge', 'is_splittable', 'ext_paves',
    'common_block_idx', 'shrunk_range', 'my_shrunk_box',
    'pave1', 'pave2', 'curve', 'pcurve_on_a', 'pcurve_on_b',
]
# Methods on PaveBlock
ACCESSED_METHODS = [
    'indices()', 'range()', 'has_shrunk_data()', 'is_to_update()',
    'update(', 'append_ext_pave(', 'remove_ext_pave(',
    'set_shrunk_data(', 'set_shrunk_data_with_box(',
    'has_same_bounds(', 'contains_parameter(',
]

def fix_file(filepath):
    with open(filepath, 'rb') as f:
        raw = f.read()
    try:
        text = raw.decode('utf-8')
    except:
        return False

    orig = text

    # Pattern: `something.FIELD` that's NOT already `.read().unwrap().FIELD`
    # and NOT `ds.pave_blocks` or `self.pave_blocks` (DS field)
    for f in ACCESSED_FIELDS:
        # Replace pb.FIELD with pb.read().unwrap().FIELD
        # But NOT when preceded by .read().unwrap() or .write().unwrap()
        # or when the left side is ds.pave_blocks[idx] (which is already SharedPB)
        text = re.sub(
            r'(?<!\.read\(\)\.unwrap\(\))(?<!\.write\(\)\.unwrap\(\))'
            r'(\.)' + re.escape(f) + r'\b',
            r'.read().unwrap()' + r'\1' + f,
            text
        )
        # Hmm that's wrong. Let me be more careful.
        pass

    # Actually, let me try a completely different approach.
    # Just handle the most common patterns:
    
    text = orig
    
    # Read pattern: pb.FIELD -> pb.read().unwrap().FIELD for reads
    for f in ACCESSED_FIELDS:
        # Don't touch if already has .read().unwrap() or .write().unwrap()
        # Match: VAR.FIELD (not after ds., self., face_info., etc)
        text = re.sub(
            r'(\w+)\.' + re.escape(f) + r'\b(?!\s*\()',
            lambda m: m.group(0) if '.read().unwrap()' in m.group(0) or '.write().unwrap()' in m.group(0) else m.group(1) + '.read().unwrap().' + f,
            text
        )
    
    # Write pattern: pb.FIELD = -> pb.write().unwrap().FIELD =
    for f in ACCESSED_FIELDS:
        text = re.sub(
            r'(\w+)\.' + re.escape(f) + r'\s*=',
            lambda m: m.group(0) if '.write().unwrap()' in m.group(0) or '.read().unwrap()' in m.group(0) else m.group(1) + '.write().unwrap().' + f + ' =',
            text
        )
    
    # Method calls: pb.method() -> pb.read().unwrap().method()
    for m in ACCESSED_METHODS:
        if m.endswith('('):
            # Method with args: pb.method(args)
            text = re.sub(
                r'(\w+)\.' + re.escape(m),
                lambda m2: m2.group(0) if '.read().unwrap()' in m2.group(0) or '.write().unwrap()' in m2.group(0) else m2.group(1) + '.write().unwrap().' + m,
                text
            )
        else:
            text = re.sub(
                r'(\w+)\.' + re.escape(m),
                lambda m2: m2.group(0) if '.read().unwrap()' in m2.group(0) or '.write().unwrap()' in m2.group(0) else m2.group(1) + '.read().unwrap().' + m,
                text
            )
    
    if text != orig:
        with open(filepath, 'wb') as f:
            f.write(text.encode('utf-8'))
        return True
    return False

if __name__ == '__main__':
    for f in sys.argv[1:]:
        if fix_file(f):
            print(f'Fixed: {f}')
