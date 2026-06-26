import re

with open(r'C:\Users\lilu\works\rcad-pro\rcad\libs\rcad-algorithms\src\pave_filler\mod.rs', 'r', encoding='utf-8') as f:
    text = f.read()

# Target methods
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

# Find each method start position
method_starts = {}
for t in targets:
    pattern = r'^    fn ' + re.escape(t) + r'\b'
    for m in re.finditer(pattern, text, re.MULTILINE):
        method_starts[t] = m.start()
        break

# Extract each method body (balanced braces)
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

# Sort methods by position in file
sorted_methods = sorted(method_starts.items(), key=lambda x: x[1])

# Build the output
out_lines = []
out_lines.append('use glam::DVec3;')
out_lines.append('use rcad_kernel::geom::*;')
out_lines.append('use rcad_kernel::geom::CurveEval;')
out_lines.append('use crate::bopds::ds::{DS, DSEdge, DSRepOnFace, Interference, IntersectionCurve, ShapeOrigin};')
out_lines.append('use crate::bopds::pave::*;')
out_lines.append('use crate::bvh::Bvh;')
out_lines.append('use crate::tolerance::*;')
out_lines.append('use crate::inttools;')
out_lines.append('use crate::inttools::context::Context as IntToolsContext;')
out_lines.append('use crate::inttools::fclass2d::{FClass2d, State};')
out_lines.append('use crate::pave_filler::helpers::*;')
out_lines.append('use std::collections::HashSet;')
out_lines.append('')
out_lines.append("impl<'a> super::PaveFiller<'a> {")

for name, pos in sorted_methods:
    body = extract_method(text, pos)
    # Remove trailing newlines
    body = body.rstrip('\n')
    out_lines.append('\n' + body)

out_lines.append('}')
out_lines.append('')

# Also extract propagate_ic_vertices_to_shared_faces free function
func_pattern = r'^fn propagate_ic_vertices_to_shared_faces\b'
func_match = re.search(func_pattern, text, re.MULTILINE)
if func_match:
    body = extract_method(text, func_match.start())
    body = body.rstrip('\n')
    out_lines.append('\n' + body)

result = '\n'.join(out_lines)

with open(r'C:\Users\lilu\works\rcad-pro\rcad\libs\rcad-algorithms\src\pave_filler\ff.rs', 'w', encoding='utf-8') as f:
    f.write(result)

print(f'Extracted {len(sorted_methods)} methods to ff.rs')
print(f'File size: {len(result)} bytes')
