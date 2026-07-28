"""
Fix ALL ShapeRef -> &Shape/Shape changes in topods.rs and add bridge methods.
"""
path = r'C:\Users\lilu\works\rcad-pro\rcad\libs\rcad-kernel\src\topods.rs'
with open(path, 'r', encoding='utf-8') as f:
    content = f.read()

# Step 1: Fix import
content = content.replace('use crate::topo_shape::Shape;', 'pub use crate::topo_shape::Shape;')

# Step 2: Add shape_idx and shape_to_ref to BRep
# Find position to insert - before the first BRep accessor method (vertex/pnt)
insert_point = 'pub fn vertex(&self, r: ShapeRef) -> &TVertexData {'
if insert_point in content:
    idx = content.find(insert_point)
    bridge_methods = '''
    /// Resolve the flat tshape index for a Shape via pointer identity.
    pub fn shape_idx(&self, s: &Shape) -> usize {
        self.tshapes.iter().position(|ts| std::sync::Arc::as_ptr(ts) == std::sync::Arc::as_ptr(&s.data)).unwrap_or(0)
    }

    /// Convert a &Shape to a ShapeRef for use with the accessor methods below.
    pub fn shape_to_ref(&self, s: &Shape) -> ShapeRef {
        let idx = self.shape_idx(s);
        ShapeRef {
            ptr_id: s.ptr_id(),
            index: idx,
            orientation: s.orientation,
            location: s.location,
        }
    }

    /// Convert a ShapeRef back to a Shape in this BRep.
    pub fn ref_to_shape(&self, r: ShapeRef) -> Shape {
        Shape::new(self.tshapes[r.index].clone(), r.location, r.orientation)
    }

'''
    content = content[:idx] + bridge_methods + content[idx:]

# Step 3: Fix BRepTool TRAIT method signatures
trait_fixes = [
    ('fn vertex_position(&self, v: ShapeRef)', 'fn vertex_position(&self, v: &Shape)'),
    ('fn vertex_tolerance(&self, v: ShapeRef)', 'fn vertex_tolerance(&self, v: &Shape)'),
    ('fn is_edge_degenerated(&self, e: ShapeRef)', 'fn is_edge_degenerated(&self, e: &Shape)'),
    ('fn edge_other_vertex(&self, edge: ShapeRef, v: ShapeRef) -> ShapeRef',
     'fn edge_other_vertex(&self, edge: &Shape, v: &Shape) -> Shape'),
    ('fn first_vertex(&self, edge: ShapeRef) -> ShapeRef',
     'fn first_vertex(&self, edge: &Shape) -> Shape'),
    ('fn last_vertex(&self, edge: ShapeRef) -> ShapeRef',
     'fn last_vertex(&self, edge: &Shape) -> Shape'),
    ('fn oriented_first_vertex(&self, edge: ShapeRef, orientation: Orientation) -> ShapeRef',
     'fn oriented_first_vertex(&self, edge: &Shape, orientation: Orientation) -> Shape'),
    ('fn parameter_on_edge(&self, vertex: ShapeRef, edge: ShapeRef, face: ShapeRef)',
     'fn parameter_on_edge(&self, vertex: &Shape, edge: &Shape, face: &Shape)'),
    ('fn curve_on_surface(&self, edge: ShapeRef, face: ShapeRef)',
     'fn curve_on_surface(&self, edge: &Shape, face: &Shape)'),
    ('fn face_surface(&self, face: ShapeRef)', 'fn face_surface(&self, face: &Shape)'),
    ('fn face_surface_world(&self, face: ShapeRef)', 'fn face_surface_world(&self, face: &Shape)'),
    ('fn edge_curve_world(&self, edge: ShapeRef)', 'fn edge_curve_world(&self, edge: &Shape)'),
    ('fn u_resolution(&self, face: ShapeRef,', 'fn u_resolution(&self, face: &Shape,'),
    ('fn v_resolution(&self, face: ShapeRef,', 'fn v_resolution(&self, face: &Shape,'),
    ('fn vertex_orientation(&self, _v: ShapeRef)', 'fn vertex_orientation(&self, _v: &Shape)'),
    ('fn is_edge_closed_on_face(&self, edge: ShapeRef, face: ShapeRef)',
     'fn is_edge_closed_on_face(&self, edge: &Shape, face: &Shape)'),
    ('fn is_closed(&self, s: ShapeRef)', 'fn is_closed(&self, s: &Shape)'),
    ('fn edge_same_parameter(&self, e: ShapeRef)', 'fn edge_same_parameter(&self, e: &Shape)'),
    ('fn edge_same_range(&self, e: ShapeRef)', 'fn edge_same_range(&self, e: &Shape)'),
    ('fn face_natural_restriction(&self, f: ShapeRef)', 'fn face_natural_restriction(&self, f: &Shape)'),
    ('fn edge_curve_data(&self, e: ShapeRef)', 'fn edge_curve_data(&self, e: &Shape)'),
    ('fn edge_range(&self, e: ShapeRef)', 'fn edge_range(&self, e: &Shape)'),
    ('fn tolerance(&self, s: ShapeRef)', 'fn tolerance(&self, s: &Shape)'),
    ('fn shape_type(&self, s: ShapeRef)', 'fn shape_type(&self, s: &Shape)'),
    ('fn has_flag(&self, s: ShapeRef,', 'fn has_flag(&self, s: &Shape,'),
    ('fn edge_data(&self, e: ShapeRef)', 'fn edge_data(&self, e: &Shape)'),
    ('fn face_data(&self, f: ShapeRef)', 'fn face_data(&self, f: &Shape)'),
]

for old, new in trait_fixes:
    content = content.replace(old, new)

# Step 4: Fix BRep IMPL method signatures (use the same replacements)
for old, new in trait_fixes:
    content = content.replace(old, new)

# Step 5: Fix impl method bodies
# vertex_position
content = content.replace(
    'fn vertex_position(&self, v: &Shape) -> DVec3 {\n        let pt = self.vertex(v).point;\n        self.get_location(v.location).transform_point3(pt)\n    }',
    'fn vertex_position(&self, v: &Shape) -> DVec3 {\n        let sr = self.shape_to_ref(v);\n        let pt = self.vertex(sr).point;\n        self.get_location(sr.location).transform_point3(pt)\n    }'
)
# vertex_tolerance
content = content.replace(
    'fn vertex_tolerance(&self, v: &Shape) -> f64 {\n        self.vertex(v).tolerance\n    }',
    'fn vertex_tolerance(&self, v: &Shape) -> f64 {\n        self.vertex(self.shape_to_ref(v)).tolerance\n    }'
)
# is_edge_degenerated
content = content.replace(
    'fn is_edge_degenerated(&self, e: &Shape) -> bool {\n        self.edge(e).degenerated\n    }',
    'fn is_edge_degenerated(&self, e: &Shape) -> bool {\n        self.edge(self.shape_to_ref(e)).degenerated\n    }'
)
# edge_other_vertex
old_body = 'fn edge_other_vertex(&self, edge: &Shape, v: &Shape) -> Shape {\n        let ed = self.edge(edge);\n        if ed.first.index == v.index {\n            ed.last\n        } else {\n            ed.first\n        }\n    }'
new_body = 'fn edge_other_vertex(&self, edge: &Shape, v: &Shape) -> Shape {\n        let esr = self.shape_to_ref(edge);\n        let vsr = self.shape_to_ref(v);\n        let ed = self.edge(esr);\n        if ed.first.ptr_id() == vsr.ptr_id {\n            ed.last.clone()\n        } else {\n            ed.first.clone()\n        }\n    }'
content = content.replace(old_body, new_body)
# first_vertex
content = content.replace(
    'fn first_vertex(&self, edge: &Shape) -> Shape {\n        self.edge(edge).first\n    }',
    'fn first_vertex(&self, edge: &Shape) -> Shape {\n        self.edge(self.shape_to_ref(edge)).first.clone()\n    }'
)
# last_vertex
content = content.replace(
    'fn last_vertex(&self, edge: &Shape) -> Shape {\n        self.edge(edge).last\n    }',
    'fn last_vertex(&self, edge: &Shape) -> Shape {\n        self.edge(self.shape_to_ref(edge)).last.clone()\n    }'
)
# parameter_on_edge
old_body = 'fn parameter_on_edge(&self, vertex: &Shape, edge: &Shape, _face: &Shape) -> Option<f64> {\n        self.edge(edge).vertex_params.get(&vertex.index).copied()\n    }'
new_body = 'fn parameter_on_edge(&self, vertex: &Shape, edge: &Shape, _face: &Shape) -> Option<f64> {\n        let esr = self.shape_to_ref(edge);\n        let vsr = self.shape_to_ref(vertex);\n        self.edge(esr).vertex_params.get(&vsr.index).copied()\n    }'
content = content.replace(old_body, new_body)
# curve_on_surface
old_body = 'fn curve_on_surface(&self, edge: &Shape, face: &Shape) -> Option<&(Curve2d, f64, f64)> {\n        self.edge(edge).pcurves.get(&face.index)\n    }'
new_body = 'fn curve_on_surface(&self, edge: &Shape, face: &Shape) -> Option<&(Curve2d, f64, f64)> {\n        let esr = self.shape_to_ref(edge);\n        let fsr = self.shape_to_ref(face);\n        self.edge(esr).pcurves.get(&fsr.index)\n    }'
content = content.replace(old_body, new_body)
# face_surface
content = content.replace(
    'fn face_surface(&self, face: &Shape) -> Option<&Surface3> {\n        self.face(face).surface.as_ref()\n    }',
    'fn face_surface(&self, face: &Shape) -> Option<&Surface3> {\n        self.face(self.shape_to_ref(face)).surface.as_ref()\n    }'
)
# face_surface_world
old_body = 'fn face_surface_world(&self, face: &Shape) -> Option<Surface3> {\n        let fd = self.face(face);\n        let surface = fd.surface.as_ref()?.clone();\n        let loc = self.get_location(face.location);\n        if loc == glam::DAffine3::IDENTITY {\n            Some(surface)\n        } else {\n            Some(crate::geom::transform_surface(&surface, &loc))\n        }\n    }'
new_body = 'fn face_surface_world(&self, face: &Shape) -> Option<Surface3> {\n        let fsr = self.shape_to_ref(face);\n        let fd = self.face(fsr);\n        let surface = fd.surface.as_ref()?.clone();\n        let loc = self.get_location(fsr.location);\n        if loc == glam::DAffine3::IDENTITY {\n            Some(surface)\n        } else {\n            Some(crate::geom::transform_surface(&surface, &loc))\n        }\n    }'
content = content.replace(old_body, new_body)
# edge_curve_world
old_body = 'fn edge_curve_world(&self, edge: &Shape) -> Option<(Curve3, [f64; 2])> {\n        let ed = self.edge(edge);\n        let crv = ed.curve.as_ref()?.clone();\n        let loc = self.get_location(edge.location);\n        if loc == glam::DAffine3::IDENTITY {\n            Some((crv, ed.range))\n        } else {\n            Some((crate::geom::transform_curve(&crv, &loc), ed.range))\n        }\n    }'
new_body = 'fn edge_curve_world(&self, edge: &Shape) -> Option<(Curve3, [f64; 2])> {\n        let esr = self.shape_to_ref(edge);\n        let ed = self.edge(esr);\n        let crv = ed.curve.as_ref()?.clone();\n        let loc = self.get_location(esr.location);\n        if loc == glam::DAffine3::IDENTITY {\n            Some((crv, ed.range))\n        } else {\n            Some((crate::geom::transform_curve(&crv, &loc), ed.range))\n        }\n    }'
content = content.replace(old_body, new_body)
# u_resolution
old_body = 'fn u_resolution(&self, face: &Shape, tol3d: f64) -> f64 {\n        match self.face(face).surface.as_ref() {\n            Some(surf) => u_resolution_for_surface(surf, tol3d),\n            None => tol3d,\n        }\n    }'
new_body = 'fn u_resolution(&self, face: &Shape, tol3d: f64) -> f64 {\n        match self.face(self.shape_to_ref(face)).surface.as_ref() {\n            Some(surf) => u_resolution_for_surface(surf, tol3d),\n            None => tol3d,\n        }\n    }'
content = content.replace(old_body, new_body)
# v_resolution
old_body = 'fn v_resolution(&self, face: &Shape, tol3d: f64) -> f64 {\n        match self.face(face).surface.as_ref() {\n            Some(surf) => v_resolution_for_surface(surf, tol3d),\n            None => tol3d,\n        }\n    }'
new_body = 'fn v_resolution(&self, face: &Shape, tol3d: f64) -> f64 {\n        match self.face(self.shape_to_ref(face)).surface.as_ref() {\n            Some(surf) => v_resolution_for_surface(surf, tol3d),\n            None => tol3d,\n        }\n    }'
content = content.replace(old_body, new_body)
# tolerance
content = content.replace(
    'match &*self.tshapes[s.index] {',
    'match &*self.tshapes[self.shape_to_ref(s).index] {'
)
# has_flag
content = content.replace(
    'fn has_flag(&self, s: &Shape, flag: u16) -> bool {\n        let flags = match &*self.tshapes[s.index] {',
    'fn has_flag(&self, s: &Shape, flag: u16) -> bool {\n        let sr = self.shape_to_ref(s);\n        let flags = match &*self.tshapes[sr.index] {'
)
# edge_data
content = content.replace(
    'fn edge_data(&self, e: &Shape) -> Option<&TEdgeData> {\n        match &*self.tshapes[e.index] {',
    'fn edge_data(&self, e: &Shape) -> Option<&TEdgeData> {\n        let esr = self.shape_to_ref(e);\n        match &*self.tshapes[esr.index] {'
)
# face_data
content = content.replace(
    'fn face_data(&self, f: &Shape) -> Option<&TFaceData> {\n        match &*self.tshapes[f.index] {',
    'fn face_data(&self, f: &Shape) -> Option<&TFaceData> {\n        let fsr = self.shape_to_ref(f);\n        match &*self.tshapes[fsr.index] {'
)

with open(path, 'w', encoding='utf-8') as f:
    f.write(content)

print("All topods.rs changes applied")
