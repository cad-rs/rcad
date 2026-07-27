"""Fix variable-binding patterns in SAFE files only (no scoping issues).
Skips: filler_mod.rs, glue.rs, classify.rs (have var reuse across scopes).
"""
import re
import os

def join_multiline(text):
    return re.sub(r'(\w+)\s*\n\s*\.(\w+)', r'\1.\2', text)

EDGE_FIELDS = [
    ('start_vertex', 'edge_start_vertex_ds', False),
    ('end_vertex', 'edge_end_vertex_ds', False),
    ('t_range', 'edge_range', False),
    ('pave_blocks', 'edge_pave_blocks', False),
    ('paves', 'edge_paves', False),
    ('face_reps', 'edge_face_reps', False),
    ('geom_tol', 'edge_tolerance', False),
    ('is_internal', 'edge_is_internal', False),
    ('location', 'edge_location', False),
    ('origin', 'edge_origin', False),
    ('curve', None, True),  # special: needs unwrap
]

FACE_FIELDS = [
    ('normal', 'face_normal', False),
    ('origin', 'face_origin', False),
    ('boundary_edges', 'face_boundary_edges', False),
    ('boundary_verts', 'face_boundary_verts', False),
    ('boundary_edge_forwards', 'face_boundary_edge_forwards', False),
    ('inner_boundary_edges', 'face_inner_boundary', False),
    ('outer_wire_idx', 'face_outer_wire_idx', False),
    ('inner_wire_idxs', 'face_inner_wire_idxs', False),
    ('geom_tol', 'face_tolerance', False),
    ('location', 'face_location', False),
    ('natural_restriction', 'face_natural_restriction', False),
    ('source_face_idx', 'source_face_idx', False),
    ('source_shell_idx', 'source_shell_idx', False),
    ('source_solid_idx', 'source_solid_idx', False),
    ('source_compsolid_idx', 'source_compsolid_idx', False),
    ('face_info', 'face_info', False),
    ('uv_boundary', 'face_uv_boundary', False),
    ('surface', None, True),  # special: needs unwrap
]

def fix_file(fpath):
    with open(fpath, 'r', encoding='utf-8') as f:
        content = f.read()
    orig = content

    content = join_multiline(content)

    # Find edge bindings
    for m in list(re.finditer(r'^(\s*)let\s+(\w+)\s*=\s*&((?:self\.)?ds)\.edges\[([^\]]+)\]\s*;\s*$', content, re.MULTILINE)):
        var, prefix, idx = m.group(2), m.group(3), m.group(4)
        binding = m.group(0)
        pref = f'{prefix}.'
        content = content.replace(binding, '')

        for field, acc, special in EDGE_FIELDS:
            old = f'{var}.{field}'
            if special and field == 'curve':
                content = content.replace(f'&{old}', f'{pref}edge_curve({idx}).unwrap()')
                content = content.replace(f'{old}.clone()', f'{pref}edge_curve({idx}).cloned().unwrap()')
                content = content.replace(old, f'{pref}edge_curve({idx})')
            else:
                content = content.replace(old, f'{pref}{acc}({idx})')

    # Find face bindings
    for m in list(re.finditer(r'^(\s*)let\s+(\w+)\s*=\s*&((?:self\.)?ds)\.faces\[([^\]]+)\]\s*;\s*$', content, re.MULTILINE)):
        var, prefix, idx = m.group(2), m.group(3), m.group(4)
        binding = m.group(0)
        pref = f'{prefix}.'
        content = content.replace(binding, '')

        for field, acc, special in FACE_FIELDS:
            old = f'{var}.{field}'
            if special and field == 'surface':
                content = content.replace(f'&{old}', f'{pref}face_surface({idx}).unwrap()')
                content = content.replace(f'{old}.clone()', f'{pref}face_surface({idx}).cloned().unwrap()')
                content = content.replace(f'match &{old}', f'match {pref}face_surface({idx}).unwrap()')
                content = content.replace(f'match {old}', f'match {pref}face_surface({idx}).unwrap()')
                content = content.replace(f'{old}.default_domain()', f'{pref}face_surface({idx}).unwrap().default_domain()')
            else:
                content = content.replace(old, f'{pref}{acc}({idx})')

    if content != orig:
        with open(fpath, 'w', encoding='utf-8') as f:
            f.write(content)
        return True
    return False

# SAFE files only (no variable reuse across scopes)
FILES = [
    'bopalgo/pave_filler/ff_intersect.rs',
    'bopalgo/pave_filler/interf.rs',
    'bopalgo/pave_filler/intersection.rs',
    'bopalgo/pave_filler/make_blocks.rs',
    'bopalgo/pave_filler/mod.rs',
    'bopalgo/pave_filler/paves.rs',
    'bopalgo/pave_filler/posttreat.rs',
    'bopds/checker_si.rs',
    'bopds/shell_splitter.rs',
    'boptools/extra.rs',
    'boptools/mod.rs',
    'inttools/context.rs',
]

src = 'libs/rcad-algorithms/src'
fixed = 0
for rel in FILES:
    f = os.path.join(src, rel)
    if os.path.exists(f):
        if fix_file(f):
            print(f'  Fixed: {rel}')
            fixed += 1
print(f'\nFixed {fixed} files')
