import re

with open(r'C:\Users\lilu\works\rcad-pro\rcad\libs\rcad-algorithms\src\pave_filler\mod.rs', 'r', encoding='utf-8') as f:
    original = f.read()

targets = [
    'find_face_face_curve_indices',
    'sampled_face_boundary_points',
    'closest_point_on_boundary_samples',
    'snap_polyline_endpoints_to_face_boundaries',
    'check_seam_edge_shift',
    'reverse_seam_edge_shift',
    'intersect_face_face',
    'intersect_plane_plane_faces',
    'intersect_plane_sphere_faces',
    'intersect_sphere_sphere_faces',
    'intersect_sphere_cylinder_faces',
    'intersect_cylinder_cylinder_faces',
    'cylinder_face_v_range',
    'intersect_plane_cylinder_faces',
    'intersect_plane_cone_faces',
    'intersect_cylinder_cone_faces',
    'intersect_cone_cone_faces',
    'register_torus_intersection',
    'intersect_torus_plane_faces',
    'intersect_torus_sphere_faces',
    'intersect_torus_cylinder_faces',
    'intersect_torus_cone_faces',
    'intersect_torus_torus_faces',
    'intersect_sphere_cone_faces',
    'intersect_ff_by_numeric_intss',
    'intersect_ff_by_marching',
    'make_marching_pcurves_with_reapprox',
    'generate_surface_samples',
    'generate_surface_samples_grid',
    'estimate_step_size',
]

def extract_method(text, start):
    lines = text[start:].split('\n')
    depth = 0
    result_lines = []
    for line in lines:
        result_lines.append(line)
        depth += line.count('{') - line.count('}')
        if depth == 0 and len(result_lines) > 0:
            break
    return '\n'.join(result_lines)

# Phase 1: find all method bodies in original text
to_remove = []
for t in targets:
    pattern = r'^    fn ' + re.escape(t) + r'\b'
    for m in re.finditer(pattern, original, re.MULTILINE):
        start = m.start()
        body = extract_method(original, start)
        to_remove.append((start, start + len(body)))
        break

# Also find propagate_ic_vertices_to_shared_faces
func_pattern = r'^fn propagate_ic_vertices_to_shared_faces\b'
for m in re.finditer(func_pattern, original, re.MULTILINE):
    start = m.start()
    body = extract_method(original, start)
    to_remove.append((start, start + len(body)))
    break

# Sort in reverse order by start position
to_remove.sort(key=lambda x: x[0], reverse=True)

# Phase 2: remove from text, from end to start
text = original
for start, end in to_remove:
    text = text[:start] + text[end:]

# Phase 3: add pub mod ff;
helpers_pos = text.find('pub mod helpers;')
if helpers_pos >= 0:
    eol = text.find('\n', helpers_pos)
    text = text[:eol+1] + 'pub mod ff;\n' + text[eol+1:]

# Clean triple+ blank lines to double
text = re.sub(r'\n{4,}', '\n\n\n', text)

with open(r'C:\Users\lilu\works\rcad-pro\rcad\libs\rcad-algorithms\src\pave_filler\mod.rs', 'w', encoding='utf-8') as f:
    f.write(text)

print(f'Removed {len(to_remove)} method definitions')
print(f'New mod.rs size: {len(text)} bytes')
