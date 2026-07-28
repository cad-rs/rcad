"""Fix BRepTool trait in topods.rs: ShapeRef -> &Shape for params, ShapeRef -> Shape for returns"""
with open(r'C:\Users\lilu\works\rcad-pro\rcad\libs\rcad-kernel\src\topods.rs', 'r', encoding='utf-8') as f:
    content = f.read()

# Map of old signatures to new signatures for the trait (lines ~1382-1486)
replacements = [
    # Return type changes (ShapeRef -> Shape)
    (' -> ShapeRef', ' -> Shape'),
    # Parameter changes (ShapeRef by value -> &Shape by ref)
    # Only in the trait definition and default implementations
    ('fn vertex_position(&self, v: ShapeRef)', 'fn vertex_position(&self, v: &Shape)'),
    ('fn vertex_tolerance(&self, v: ShapeRef)', 'fn vertex_tolerance(&self, v: &Shape)'),
    ('fn is_edge_degenerated(&self, e: ShapeRef)', 'fn is_edge_degenerated(&self, e: &Shape)'),
    ('fn edge_other_vertex(&self, edge: ShapeRef, v: ShapeRef)', 'fn edge_other_vertex(&self, edge: &Shape, v: &Shape)'),
    ('fn first_vertex(&self, edge: ShapeRef)', 'fn first_vertex(&self, edge: &Shape)'),
    ('fn last_vertex(&self, edge: ShapeRef)', 'fn last_vertex(&self, edge: &Shape)'),
    ('fn oriented_first_vertex(&self, edge: ShapeRef, orientation: Orientation)', 'fn oriented_first_vertex(&self, edge: &Shape, orientation: Orientation)'),
    ('fn parameter_on_edge(&self, vertex: ShapeRef, edge: ShapeRef, face: ShapeRef)', 'fn parameter_on_edge(&self, vertex: &Shape, edge: &Shape, face: &Shape)'),
    ('fn curve_on_surface(&self, edge: ShapeRef, face: ShapeRef)', 'fn curve_on_surface(&self, edge: &Shape, face: &Shape)'),
    ('fn face_surface(&self, face: ShapeRef)', 'fn face_surface(&self, face: &Shape)'),
    ('fn face_surface_world(&self, face: ShapeRef)', 'fn face_surface_world(&self, face: &Shape)'),
    ('fn edge_curve_world(&self, edge: ShapeRef)', 'fn edge_curve_world(&self, edge: &Shape)'),
    ('fn u_resolution(&self, face: ShapeRef,', 'fn u_resolution(&self, face: &Shape,'),
    ('fn v_resolution(&self, face: ShapeRef,', 'fn v_resolution(&self, face: &Shape,'),
    ('fn vertex_orientation(&self, _v: ShapeRef)', 'fn vertex_orientation(&self, _v: &Shape)'),
    ('fn tolerance(&self, s: ShapeRef)', 'fn tolerance(&self, s: &Shape)'),
    ('fn shape_type(&self, s: ShapeRef)', 'fn shape_type(&self, s: &Shape)'),
    ('fn has_flag(&self, s: ShapeRef,', 'fn has_flag(&self, s: &Shape,'),
    ('fn edge_data(&self, e: ShapeRef)', 'fn edge_data(&self, e: &Shape)'),
    ('fn face_data(&self, f: ShapeRef)', 'fn face_data(&self, f: &Shape)'),
    # Default methods
    ('fn is_closed(&self, s: ShapeRef)', 'fn is_closed(&self, s: &Shape)'),
    ('fn edge_same_parameter(&self, e: ShapeRef)', 'fn edge_same_parameter(&self, e: &Shape)'),
    ('fn edge_same_range(&self, e: ShapeRef)', 'fn edge_same_range(&self, e: &Shape)'),
    ('fn face_natural_restriction(&self, f: ShapeRef)', 'fn face_natural_restriction(&self, f: &Shape)'),
    ('fn edge_curve_data(&self, e: ShapeRef)', 'fn edge_curve_data(&self, e: &Shape)'),
    ('fn edge_range(&self, e: ShapeRef)', 'fn edge_range(&self, e: &Shape)'),
    # Default method bodies that use ShapeRef - they need & on self calls
    # is_edge_closed_on_face default
    ('fn is_edge_closed_on_face(&self, edge: &Shape, face: &Shape)', 'fn is_edge_closed_on_face(&self, edge: &Shape, face: &Shape)'),
    # curve_on_surface_second default
    ('fn curve_on_surface_second(\n        &self,\n        edge: &Shape,\n        face: &Shape,\n    )', 'fn curve_on_surface_second(\n        &self,\n        edge: &Shape,\n        face: &Shape,\n    )'),
]

for old, new in replacements:
    if old != new:  # skip no-ops
        content = content.replace(old, new)

# Now fix default method bodies that reference traits methods with old signatures
# The edge_same_parameter default calls self.edge_data(e) - e is now &Shape, still works
# The edge_data trait method takes &Shape now, so self.edge_data(e) where e: &Shape is fine

# Fix the BRep impl (starts at line ~1491)
# Use shape_to_ref bridge for internal BRep methods
# Replace method signatures in the impl block

with open(r'C:\Users\lilu\works\rcad-pro\rcad\libs\rcad-kernel\src\topods.rs', 'w', encoding='utf-8') as f:
    f.write(content)

print("Trait signatures updated")
