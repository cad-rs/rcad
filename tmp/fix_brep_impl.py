"""Fix BRep impl of BRepTool in topods.rs to match updated trait (ShapeRef -> &Shape)"""
with open(r'C:\Users\lilu\works\rcad-pro\rcad\libs\rcad-kernel\src\topods.rs', 'r', encoding='utf-8') as f:
    content = f.read()

# First fix method signatures in the impl block (starts with "impl BRepTool for BRep {")
# These are the same as the trait changes
impl_replacements = [
    ('fn vertex_position(&self, v: ShapeRef)', 'fn vertex_position(&self, v: &Shape)'),
    ('fn vertex_tolerance(&self, v: ShapeRef)', 'fn vertex_tolerance(&self, v: &Shape)'),
    ('fn is_edge_degenerated(&self, e: ShapeRef)', 'fn is_edge_degenerated(&self, e: &Shape)'),
    ('fn edge_other_vertex(&self, edge: ShapeRef, v: ShapeRef) -> ShapeRef', 'fn edge_other_vertex(&self, edge: &Shape, v: &Shape) -> Shape'),
    ('fn first_vertex(&self, edge: ShapeRef) -> ShapeRef', 'fn first_vertex(&self, edge: &Shape) -> Shape'),
    ('fn last_vertex(&self, edge: ShapeRef) -> ShapeRef', 'fn last_vertex(&self, edge: &Shape) -> Shape'),
    ('fn oriented_first_vertex(&self, edge: ShapeRef, orientation: Orientation) -> ShapeRef', 'fn oriented_first_vertex(&self, edge: &Shape, orientation: Orientation) -> Shape'),
    ('fn parameter_on_edge(&self, vertex: ShapeRef, edge: ShapeRef, _face: ShapeRef) -> Option<f64>', 'fn parameter_on_edge(&self, vertex: &Shape, edge: &Shape, _face: &Shape) -> Option<f64>'),
    ('fn curve_on_surface(&self, edge: ShapeRef, face: ShapeRef)', 'fn curve_on_surface(&self, edge: &Shape, face: &Shape)'),
    ('fn is_edge_closed_on_face(&self, edge: ShapeRef, face: ShapeRef) -> bool', 'fn is_edge_closed_on_face(&self, edge: &Shape, face: &Shape) -> bool'),
    ('fn curve_on_surface_second(\n        &self,\n        edge: ShapeRef,\n        face: ShapeRef,\n    )', 'fn curve_on_surface_second(\n        &self,\n        edge: &Shape,\n        face: &Shape,\n    )'),
    ('fn face_surface(&self, face: ShapeRef)', 'fn face_surface(&self, face: &Shape)'),
    ('fn face_surface_world(&self, face: ShapeRef)', 'fn face_surface_world(&self, face: &Shape)'),
    ('fn edge_curve_world(&self, edge: ShapeRef)', 'fn edge_curve_world(&self, edge: &Shape)'),
    ('fn u_resolution(&self, face: ShapeRef,', 'fn u_resolution(&self, face: &Shape,'),
    ('fn v_resolution(&self, face: ShapeRef,', 'fn v_resolution(&self, face: &Shape,'),
    ('fn tolerance(&self, s: ShapeRef)', 'fn tolerance(&self, s: &Shape)'),
    ('fn shape_type(&self, s: ShapeRef)', 'fn shape_type(&self, s: &Shape)'),
    ('fn has_flag(&self, s: ShapeRef,', 'fn has_flag(&self, s: &Shape,'),
    ('fn edge_data(&self, e: ShapeRef)', 'fn edge_data(&self, e: &Shape)'),
    ('fn face_data(&self, f: ShapeRef)', 'fn face_data(&self, f: &Shape)'),
]

for old, new in impl_replacements:
    if old != new:
        content = content.replace(old, new)

# Now fix method BODIES that use .index, .location, self.edge(v) etc.
# The methods need to use shape_to_ref to convert &Shape to ShapeRef for internal calls

# vertex_position body: use shape_to_ref
old_body = """fn vertex_position(&self, v: &Shape) -> DVec3 {
        let pt = self.vertex(v).point;
        self.get_location(v.location).transform_point3(pt)
    }"""
new_body = """fn vertex_position(&self, v: &Shape) -> DVec3 {
        let sr = self.shape_to_ref(v);
        let pt = self.vertex(sr).point;
        self.get_location(sr.location).transform_point3(pt)
    }"""
content = content.replace(old_body, new_body)

# vertex_tolerance
old_body = """fn vertex_tolerance(&self, v: &Shape) -> f64 {
        self.vertex(v).tolerance
    }"""
new_body = """fn vertex_tolerance(&self, v: &Shape) -> f64 {
        self.vertex(self.shape_to_ref(v)).tolerance
    }"""
content = content.replace(old_body, new_body)

# is_edge_degenerated
old_body = """fn is_edge_degenerated(&self, e: &Shape) -> bool {
        self.edge(e).degenerated
    }"""
new_body = """fn is_edge_degenerated(&self, e: &Shape) -> bool {
        self.edge(self.shape_to_ref(e)).degenerated
    }"""
content = content.replace(old_body, new_body)

# edge_other_vertex
old_body = """fn edge_other_vertex(&self, edge: &Shape, v: &Shape) -> Shape {
        let ed = self.edge(edge);
        if ed.first.index == v.index {
            ed.last
        } else {
            ed.first
        }
    }"""
new_body = """fn edge_other_vertex(&self, edge: &Shape, v: &Shape) -> Shape {
        let esr = self.shape_to_ref(edge);
        let vsr = self.shape_to_ref(v);
        let ed = self.edge(esr);
        if ed.first.ptr_id() == vsr.ptr_id {
            ed.last.clone()
        } else {
            ed.first.clone()
        }
    }"""
content = content.replace(old_body, new_body)

# first_vertex
old_body = """fn first_vertex(&self, edge: &Shape) -> Shape {
        self.edge(edge).first
    }"""
new_body = """fn first_vertex(&self, edge: &Shape) -> Shape {
        self.edge(self.shape_to_ref(edge)).first.clone()
    }"""
content = content.replace(old_body, new_body)

# last_vertex
old_body = """fn last_vertex(&self, edge: &Shape) -> Shape {
        self.edge(edge).last
    }"""
new_body = """fn last_vertex(&self, edge: &Shape) -> Shape {
        self.edge(self.shape_to_ref(edge)).last.clone()
    }"""
content = content.replace(old_body, new_body)

# oriented_first_vertex
old_body = """fn oriented_first_vertex(&self, edge: &Shape, orientation: Orientation) -> Shape {
        if orientation == Orientation::Reversed {
            self.last_vertex(edge)
        } else {
            self.first_vertex(edge)
        }
    }"""
new_body = """fn oriented_first_vertex(&self, edge: &Shape, orientation: Orientation) -> Shape {
        if orientation == Orientation::Reversed {
            self.last_vertex(edge)
        } else {
            self.first_vertex(edge)
        }
    }"""
content = content.replace(old_body, new_body)

# parameter_on_edge
old_body = """fn parameter_on_edge(&self, vertex: &Shape, edge: &Shape, _face: &Shape) -> Option<f64> {
        self.edge(edge).vertex_params.get(&vertex.index).copied()
    }"""
new_body = """fn parameter_on_edge(&self, vertex: &Shape, edge: &Shape, _face: &Shape) -> Option<f64> {
        let esr = self.shape_to_ref(edge);
        let vsr = self.shape_to_ref(vertex);
        self.edge(esr).vertex_params.get(&vsr.index).copied()
    }"""
content = content.replace(old_body, new_body)

# curve_on_surface
old_body = """fn curve_on_surface(&self, edge: &Shape, face: &Shape) -> Option<&(Curve2d, f64, f64)> {
        self.edge(edge).pcurves.get(&face.index)
    }"""
new_body = """fn curve_on_surface(&self, edge: &Shape, face: &Shape) -> Option<&(Curve2d, f64, f64)> {
        let esr = self.shape_to_ref(edge);
        let fsr = self.shape_to_ref(face);
        self.edge(esr).pcurves.get(&fsr.index)
    }"""
content = content.replace(old_body, new_body)

# is_edge_closed_on_face
old_body = """fn is_edge_closed_on_face(&self, edge: &Shape, face: &Shape) -> bool {
        let ed = self.edge(edge);
        ed.representations.iter().any(|r| matches!(r, CurveRepresentation::CurveOnClosedSurface { face: f, .. } if *f == face.index))
            || ed.pcurves.contains_key(&face.index)
                && ed.pcurves.contains_key(&(face.index + self.nb_faces()))
    }"""
new_body = """fn is_edge_closed_on_face(&self, edge: &Shape, face: &Shape) -> bool {
        let esr = self.shape_to_ref(edge);
        let fsr = self.shape_to_ref(face);
        let ed = self.edge(esr);
        ed.representations.iter().any(|r| matches!(r, CurveRepresentation::CurveOnClosedSurface { face: f, .. } if *f == fsr.index))
            || ed.pcurves.contains_key(&fsr.index)
                && ed.pcurves.contains_key(&(fsr.index + self.nb_faces()))
    }"""
content = content.replace(old_body, new_body)

# curve_on_surface_second
old_body = """fn curve_on_surface_second(
        &self,
        edge: &Shape,
        face: &Shape,
    ) -> Option<&(Curve2d, f64, f64)> {
        let shifted = face.index + self.nb_faces();
        self.edge(edge).pcurves.get(&shifted)
    }"""
new_body = """fn curve_on_surface_second(
        &self,
        edge: &Shape,
        face: &Shape,
    ) -> Option<&(Curve2d, f64, f64)> {
        let fsr = self.shape_to_ref(face);
        let shifted = fsr.index + self.nb_faces();
        self.edge(self.shape_to_ref(edge)).pcurves.get(&shifted)
    }"""
content = content.replace(old_body, new_body)

# face_surface
old_body = """fn face_surface(&self, face: &Shape) -> Option<&Surface3> {
        self.face(face).surface.as_ref()
    }"""
new_body = """fn face_surface(&self, face: &Shape) -> Option<&Surface3> {
        self.face(self.shape_to_ref(face)).surface.as_ref()
    }"""
content = content.replace(old_body, new_body)

# face_surface_world
old_body = """fn face_surface_world(&self, face: &Shape) -> Option<Surface3> {
        let fd = self.face(face);
        let surface = fd.surface.as_ref()?.clone();
        let loc = self.get_location(face.location);
        if loc == glam::DAffine3::IDENTITY {
            Some(surface)
        } else {
            Some(crate::geom::transform_surface(&surface, &loc))
        }
    }"""
new_body = """fn face_surface_world(&self, face: &Shape) -> Option<Surface3> {
        let fsr = self.shape_to_ref(face);
        let fd = self.face(fsr);
        let surface = fd.surface.as_ref()?.clone();
        let loc = self.get_location(fsr.location);
        if loc == glam::DAffine3::IDENTITY {
            Some(surface)
        } else {
            Some(crate::geom::transform_surface(&surface, &loc))
        }
    }"""
content = content.replace(old_body, new_body)

# edge_curve_world
old_body = """fn edge_curve_world(&self, edge: &Shape) -> Option<(Curve3, [f64; 2])> {
        let ed = self.edge(edge);
        let crv = ed.curve.as_ref()?.clone();
        let loc = self.get_location(edge.location);
        if loc == glam::DAffine3::IDENTITY {
            Some((crv, ed.range))
        } else {
            Some((crate::geom::transform_curve(&crv, &loc), ed.range))
        }
    }"""
new_body = """fn edge_curve_world(&self, edge: &Shape) -> Option<(Curve3, [f64; 2])> {
        let esr = self.shape_to_ref(edge);
        let ed = self.edge(esr);
        let crv = ed.curve.as_ref()?.clone();
        let loc = self.get_location(esr.location);
        if loc == glam::DAffine3::IDENTITY {
            Some((crv, ed.range))
        } else {
            Some((crate::geom::transform_curve(&crv, &loc), ed.range))
        }
    }"""
content = content.replace(old_body, new_body)

# u_resolution
old_body = """fn u_resolution(&self, face: &Shape, tol3d: f64) -> f64 {
        match self.face(face).surface.as_ref() {
            Some(surf) => u_resolution_for_surface(surf, tol3d),
            None => tol3d,
        }
    }"""
new_body = """fn u_resolution(&self, face: &Shape, tol3d: f64) -> f64 {
        match self.face(self.shape_to_ref(face)).surface.as_ref() {
            Some(surf) => u_resolution_for_surface(surf, tol3d),
            None => tol3d,
        }
    }"""
content = content.replace(old_body, new_body)

# v_resolution
old_body = """fn v_resolution(&self, face: &Shape, tol3d: f64) -> f64 {
        match self.face(face).surface.as_ref() {
            Some(surf) => v_resolution_for_surface(surf, tol3d),
            None => tol3d,
        }
    }"""
new_body = """fn v_resolution(&self, face: &Shape, tol3d: f64) -> f64 {
        match self.face(self.shape_to_ref(face)).surface.as_ref() {
            Some(surf) => v_resolution_for_surface(surf, tol3d),
            None => tol3d,
        }
    }"""
content = content.replace(old_body, new_body)

# tolerance
old_body = """fn tolerance(&self, s: &Shape) -> f64 {
        match &*self.tshapes[s.index] {"""
new_body = """fn tolerance(&self, s: &Shape) -> f64 {
        let sr = self.shape_to_ref(s);
        match &*self.tshapes[sr.index] {"""
content = content.replace(old_body, new_body)

# has_flag - need to fix the body
old_body = """fn has_flag(&self, s: &Shape, flag: u16) -> bool {
        let flags = match &*self.tshapes[s.index] {"""
new_body = """fn has_flag(&self, s: &Shape, flag: u16) -> bool {
        let sr = self.shape_to_ref(s);
        let flags = match &*self.tshapes[sr.index] {"""
content = content.replace(old_body, new_body)

# edge_data
old_body = """fn edge_data(&self, e: &Shape) -> Option<&TEdgeData> {
        match &*self.tshapes[e.index] {"""
new_body = """fn edge_data(&self, e: &Shape) -> Option<&TEdgeData> {
        let esr = self.shape_to_ref(e);
        match &*self.tshapes[esr.index] {"""
content = content.replace(old_body, new_body)

# face_data
old_body = """fn face_data(&self, f: &Shape) -> Option<&TFaceData> {
        match &*self.tshapes[f.index] {"""
new_body = """fn face_data(&self, f: &Shape) -> Option<&TFaceData> {
        let fsr = self.shape_to_ref(f);
        match &*self.tshapes[fsr.index] {"""
content = content.replace(old_body, new_body)

with open(r'C:\Users\lilu\works\rcad-pro\rcad\libs\rcad-kernel\src\topods.rs', 'w', encoding='utf-8') as f:
    f.write(content)

print("BRep impl updated")
