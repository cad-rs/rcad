use rcad_kernel::BRep;
use wgpu::util::DeviceExt;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SelectionMode {
    Face,
    Edge,
}

#[derive(Clone, Debug)]
pub struct SelectionState {
    pub mode: SelectionMode,
    pub additive_select: bool,
    pub selected_faces: Vec<usize>,
    pub selected_edges: Vec<usize>,
    pub hovered_face: Option<usize>,
    pub hovered_edge: Option<usize>,
    pub last_hover_pos: Option<(f32, f32)>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DisplayMode {
    SolidWithEdges,
    Solid,
    Wireframe,
    Transparent,
}

impl Default for DisplayMode {
    fn default() -> Self {
        Self::SolidWithEdges
    }
}

pub const DEFAULT_EDGE_PICK_RADIUS_PX: f32 = 8.0;

impl Default for SelectionState {
    fn default() -> Self {
        Self {
            mode: SelectionMode::Face,
            additive_select: false,
            selected_faces: Vec::new(),
            selected_edges: Vec::new(),
            hovered_face: None,
            hovered_edge: None,
            last_hover_pos: None,
        }
    }
}

impl SelectionState {
    pub fn set_mode(&mut self, mode: SelectionMode) {
        if self.mode != mode {
            self.mode = mode;
            self.clear_hover();
        }
    }

    pub fn click_at(
        &mut self,
        brep: &BRep,
        camera: &Camera,
        aspect: f32,
        viewport: [f32; 2],
        cursor: [f32; 2],
        edge_pick_radius_px: f32,
    ) {
        if cursor[0] < 0.0 || cursor[1] < 0.0 || cursor[0] > viewport[0] || cursor[1] > viewport[1]
        {
            return;
        }

        match self.mode {
            SelectionMode::Face => {
                let hit = pick_face(brep, camera, aspect, viewport, cursor);
                if self.additive_select {
                    if let Some(idx) = hit {
                        toggle_index(&mut self.selected_faces, idx);
                    }
                } else {
                    self.selected_faces.clear();
                    if let Some(idx) = hit {
                        self.selected_faces.push(idx);
                    }
                }
            }
            SelectionMode::Edge => {
                let hit = pick_edge(brep, camera, aspect, viewport, cursor, edge_pick_radius_px);
                if self.additive_select {
                    if let Some(idx) = hit {
                        toggle_index(&mut self.selected_edges, idx);
                    }
                } else {
                    self.selected_edges.clear();
                    if let Some(idx) = hit {
                        self.selected_edges.push(idx);
                    }
                }
            }
        }
    }

    pub fn hover_at(
        &mut self,
        brep: &BRep,
        camera: &Camera,
        aspect: f32,
        viewport: [f32; 2],
        cursor: [f32; 2],
        edge_pick_radius_px: f32,
    ) {
        if cursor[0] < 0.0 || cursor[1] < 0.0 || cursor[0] > viewport[0] || cursor[1] > viewport[1]
        {
            self.clear_hover();
            return;
        }

        self.last_hover_pos = Some((cursor[0], cursor[1]));
        match self.mode {
            SelectionMode::Face => {
                self.hovered_face = pick_face(brep, camera, aspect, viewport, cursor);
                self.hovered_edge = None;
            }
            SelectionMode::Edge => {
                self.hovered_edge =
                    pick_edge(brep, camera, aspect, viewport, cursor, edge_pick_radius_px);
                self.hovered_face = None;
            }
        }
    }

    pub fn clear_hover(&mut self) {
        self.hovered_face = None;
        self.hovered_edge = None;
        self.last_hover_pos = None;
    }

    pub fn highlighted_faces(&self) -> Vec<usize> {
        merged_indices(&self.selected_faces, self.hovered_face)
    }

    pub fn highlighted_edges(&self) -> Vec<usize> {
        merged_indices(&self.selected_edges, self.hovered_edge)
    }
}

fn toggle_index(list: &mut Vec<usize>, idx: usize) {
    if let Some(pos) = list.iter().position(|&v| v == idx) {
        list.swap_remove(pos);
    } else {
        list.push(idx);
    }
}

fn merged_indices(selected: &[usize], hovered: Option<usize>) -> Vec<usize> {
    let mut merged = selected.to_vec();
    if let Some(h) = hovered
        && !merged.contains(&h)
    {
        merged.push(h);
    }
    merged
}

pub fn pick_face(
    brep: &BRep,
    camera: &Camera,
    aspect: f32,
    viewport_size: [f32; 2],
    cursor_pos: [f32; 2],
) -> Option<usize> {
    let ray = screen_ray(camera, aspect, viewport_size, cursor_pos)?;
    let mut best: Option<(f32, usize)> = None;

    let mut face_idx = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                for tri in &face.triangles {
                    let a = to_vec3(brep.vertices.get(tri[0])?.point);
                    let b = to_vec3(brep.vertices.get(tri[1])?.point);
                    let c = to_vec3(brep.vertices.get(tri[2])?.point);
                    if let Some(t) = ray_triangle_intersection(ray.0, ray.1, a, b, c)
                        && t > 0.0
                    {
                        match best {
                            Some((best_t, _)) if t >= best_t => {}
                            _ => best = Some((t, face_idx)),
                        }
                    }
                }
                face_idx += 1;
            }
        }
    }

    best.map(|(_, idx)| idx)
}

pub fn pick_edge(
    brep: &BRep,
    camera: &Camera,
    aspect: f32,
    viewport_size: [f32; 2],
    cursor_pos: [f32; 2],
    max_distance_px: f32,
) -> Option<usize> {
    if viewport_size[0] <= 1.0 || viewport_size[1] <= 1.0 {
        return None;
    }

    let vp = glam::Mat4::from_cols_array_2d(&camera.build_view_projection_matrix(aspect.max(0.001)));
    let mut best: Option<(f32, f32, usize)> = None;

    for (idx, edge) in brep.edges.iter().enumerate() {
        let p0 = to_vec3(brep.vertices.get(edge.start)?.point);
        let p1 = to_vec3(brep.vertices.get(edge.end)?.point);

        let s0 = project_to_screen(vp, p0, viewport_size)?;
        let s1 = project_to_screen(vp, p1, viewport_size)?;
        let distance = point_segment_distance_2d(cursor_pos, [s0[0], s0[1]], [s1[0], s1[1]]);

        if distance > max_distance_px {
            continue;
        }

        let depth = (s0[2] + s1[2]) * 0.5;
        match best {
            Some((best_dist, best_depth, _))
                if distance > best_dist
                    || ((distance - best_dist).abs() < 1e-3 && depth >= best_depth) => {}
            _ => best = Some((distance, depth, idx)),
        }
    }

    best.map(|(_, _, idx)| idx)
}

#[derive(Clone, Debug)]
pub struct Mesh {
    pub vertices: Vec<[f32; 3]>,
    pub indices: Vec<u32>,
    pub line_indices: Vec<u32>,
}

#[repr(C)]
#[derive(Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
struct MaterialUniform {
    color: [f32; 4],
    flags: [f32; 4],
}

pub fn build_face_highlight_mesh(brep: &BRep, face_index: usize) -> Option<Mesh> {
    build_faces_highlight_mesh(brep, &[face_index])
}

pub fn build_faces_highlight_mesh(brep: &BRep, face_indices: &[usize]) -> Option<Mesh> {
    if face_indices.is_empty() {
        return None;
    }

    let selected: std::collections::HashSet<usize> = face_indices.iter().copied().collect();
    let mut current = 0usize;
    let mut indices: Vec<u32> = Vec::new();

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                if selected.contains(&current) {
                    for tri in &face.triangles {
                        indices.push(tri[0] as u32);
                        indices.push(tri[1] as u32);
                        indices.push(tri[2] as u32);
                    }
                }
                current += 1;
            }
        }
    }

    if indices.is_empty() {
        return None;
    }

    let vertices: Vec<[f32; 3]> = brep
        .vertices
        .iter()
        .map(|v| [v.point.x as f32, v.point.y as f32, v.point.z as f32])
        .collect();

    Some(Mesh {
        vertices,
        indices,
        line_indices: Vec::new(),
    })
}

pub fn build_edge_highlight_mesh(brep: &BRep, edge_index: usize) -> Option<Mesh> {
    build_edges_highlight_mesh(brep, &[edge_index])
}

pub fn build_edges_highlight_mesh(brep: &BRep, edge_indices: &[usize]) -> Option<Mesh> {
    if edge_indices.is_empty() {
        return None;
    }

    let mut indices = Vec::with_capacity(edge_indices.len() * 2);
    for &edge_index in edge_indices {
        let edge = brep.edges.get(edge_index)?;
        indices.push(edge.start as u32);
        indices.push(edge.end as u32);
    }

    let vertices: Vec<[f32; 3]> = brep
        .vertices
        .iter()
        .map(|v| [v.point.x as f32, v.point.y as f32, v.point.z as f32])
        .collect();

    Some(Mesh {
        vertices,
        indices,
        line_indices: Vec::new(),
    })
}

pub fn merge_meshes(meshes: &[&Mesh]) -> Option<Mesh> {
    if meshes.is_empty() {
        return None;
    }

    let total_vertices = meshes.iter().map(|mesh| mesh.vertices.len()).sum();
    let total_indices = meshes.iter().map(|mesh| mesh.indices.len()).sum();
    let total_line_indices = meshes.iter().map(|mesh| mesh.line_indices.len()).sum();

    if total_vertices == 0 || (total_indices == 0 && total_line_indices == 0) {
        return None;
    }

    let mut vertices = Vec::with_capacity(total_vertices);
    let mut indices = Vec::with_capacity(total_indices);
    let mut line_indices = Vec::with_capacity(total_line_indices);
    let mut vertex_offset = 0u32;

    for mesh in meshes {
        vertices.extend_from_slice(&mesh.vertices);
        indices.extend(mesh.indices.iter().map(|index| index + vertex_offset));
        line_indices.extend(mesh.line_indices.iter().map(|index| index + vertex_offset));
        vertex_offset += mesh.vertices.len() as u32;
    }

    Some(Mesh {
        vertices,
        indices,
        line_indices,
    })
}

pub struct Tessellator;

impl Tessellator {
    pub fn tessellate(brep: &BRep) -> Mesh {
        let flat_verts: Vec<[f32; 3]> = brep
            .vertices
            .iter()
            .map(|v| [v.point.x as f32, v.point.y as f32, v.point.z as f32])
            .collect();

        let mut indices: Vec<u32> = Vec::new();
        let mut line_indices: Vec<u32> = Vec::with_capacity(brep.edges.len() * 2);

        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    for tri in &face.triangles {
                        indices.push(tri[0] as u32);
                        indices.push(tri[1] as u32);
                        indices.push(tri[2] as u32);
                    }
                }
            }
        }

        for edge in &brep.edges {
            line_indices.push(edge.start as u32);
            line_indices.push(edge.end as u32);
        }

        Mesh {
            vertices: flat_verts,
            indices,
            line_indices,
        }
    }
}

#[repr(C)]
#[derive(Copy, Clone, bytemuck::Pod, bytemuck::Zeroable)]
pub struct CameraUniform {
    view_proj: [[f32; 4]; 4],
    eye_pos: [f32; 4],
    light_dir: [f32; 4],
}

#[derive(Clone, Copy, Debug)]
pub struct Camera {
    pub rot_x: f32,
    pub rot_y: f32,
    pub distance: f32,
    pub target: glam::Vec3,
}

impl Camera {
    pub fn new() -> Self {
        Self {
            rot_x: 0.4,
            rot_y: 0.5,
            distance: 3.0,
            target: glam::Vec3::ZERO,
        }
    }

    pub fn build_view_projection_matrix(&self, aspect: f32) -> [[f32; 4]; 4] {
        let eye = self.eye_position();
        let target = self.target;
        let up = glam::Vec3::Y;

        let view = glam::Mat4::look_at_rh(eye, target, up);
        let proj = glam::Mat4::perspective_rh(45.0_f32.to_radians(), aspect, 0.1, 100.0);

        (proj * view).to_cols_array_2d()
    }

    pub fn eye_position(&self) -> glam::Vec3 {
        self.target + glam::Vec3::new(
            self.distance * self.rot_y.cos() * self.rot_x.cos(),
            self.distance * self.rot_x.sin(),
            self.distance * self.rot_y.sin() * self.rot_x.cos(),
        )
    }

    pub fn pan_pixels(&mut self, dx: f32, dy: f32) {
        let eye = self.eye_position();
        let forward = (self.target - eye).normalize_or_zero();
        if forward.length_squared() <= 1e-8 {
            return;
        }

        let mut right = forward.cross(glam::Vec3::Y);
        if right.length_squared() <= 1e-8 {
            right = forward.cross(glam::Vec3::X);
        }
        right = right.normalize_or_zero();
        let up = right.cross(forward).normalize_or_zero();

        let scale = self.distance.max(0.1) * 0.0015;
        self.target += (-dx * right + dy * up) * scale;
    }
}

impl Default for Camera {
    fn default() -> Self {
        Self::new()
    }
}

fn to_vec3(v: glam::DVec3) -> glam::Vec3 {
    glam::Vec3::new(v.x as f32, v.y as f32, v.z as f32)
}

fn screen_ray(
    camera: &Camera,
    aspect: f32,
    viewport_size: [f32; 2],
    cursor_pos: [f32; 2],
) -> Option<(glam::Vec3, glam::Vec3)> {
    if viewport_size[0] <= 1.0 || viewport_size[1] <= 1.0 {
        return None;
    }

    let ndc_x = (2.0 * cursor_pos[0] / viewport_size[0]) - 1.0;
    let ndc_y = 1.0 - (2.0 * cursor_pos[1] / viewport_size[1]);

    let vp = glam::Mat4::from_cols_array_2d(&camera.build_view_projection_matrix(aspect.max(0.001)));
    let inv = vp.inverse();

    let near = inv * glam::Vec4::new(ndc_x, ndc_y, 0.0, 1.0);
    let far = inv * glam::Vec4::new(ndc_x, ndc_y, 1.0, 1.0);

    if near.w.abs() < 1e-6 || far.w.abs() < 1e-6 {
        return None;
    }

    let p0 = (near / near.w).truncate();
    let p1 = (far / far.w).truncate();
    let dir = (p1 - p0).normalize_or_zero();
    if dir.length_squared() <= 1e-8 {
        return None;
    }
    Some((p0, dir))
}

pub fn cursor_point_on_plane(
    camera: &Camera,
    aspect: f32,
    viewport_size: [f32; 2],
    cursor_pos: [f32; 2],
    plane_origin: glam::DVec3,
    plane_normal: glam::DVec3,
) -> Option<glam::DVec3> {
    let (ray_origin, ray_dir) = screen_ray(camera, aspect, viewport_size, cursor_pos)?;
    let plane_origin = glam::Vec3::new(
        plane_origin.x as f32,
        plane_origin.y as f32,
        plane_origin.z as f32,
    );
    let plane_normal = glam::Vec3::new(
        plane_normal.x as f32,
        plane_normal.y as f32,
        plane_normal.z as f32,
    )
    .normalize_or_zero();
    if plane_normal.length_squared() <= 1e-8 {
        return None;
    }

    let denom = plane_normal.dot(ray_dir);
    if denom.abs() <= 1e-6 {
        return None;
    }

    let t = plane_normal.dot(plane_origin - ray_origin) / denom;
    if !t.is_finite() || t < 0.0 {
        return None;
    }

    let point = ray_origin + ray_dir * t;
    Some(glam::DVec3::new(point.x as f64, point.y as f64, point.z as f64))
}

fn ray_triangle_intersection(
    ray_origin: glam::Vec3,
    ray_dir: glam::Vec3,
    v0: glam::Vec3,
    v1: glam::Vec3,
    v2: glam::Vec3,
) -> Option<f32> {
    let e1 = v1 - v0;
    let e2 = v2 - v0;
    let pvec = ray_dir.cross(e2);
    let det = e1.dot(pvec);
    if det.abs() < 1e-8 {
        return None;
    }
    let inv_det = 1.0 / det;
    let tvec = ray_origin - v0;
    let u = tvec.dot(pvec) * inv_det;
    if !(0.0..=1.0).contains(&u) {
        return None;
    }
    let qvec = tvec.cross(e1);
    let v = ray_dir.dot(qvec) * inv_det;
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    let t = e2.dot(qvec) * inv_det;
    if t.is_finite() {
        Some(t)
    } else {
        None
    }
}

fn project_to_screen(vp: glam::Mat4, p: glam::Vec3, viewport_size: [f32; 2]) -> Option<[f32; 3]> {
    let clip = vp * p.extend(1.0);
    if clip.w.abs() < 1e-6 {
        return None;
    }
    let ndc = (clip / clip.w).truncate();
    let x = (ndc.x + 1.0) * 0.5 * viewport_size[0];
    let y = (1.0 - ndc.y) * 0.5 * viewport_size[1];
    Some([x, y, ndc.z])
}

fn point_segment_distance_2d(p: [f32; 2], a: [f32; 2], b: [f32; 2]) -> f32 {
    let ab = [b[0] - a[0], b[1] - a[1]];
    let ap = [p[0] - a[0], p[1] - a[1]];
    let ab_len2 = ab[0] * ab[0] + ab[1] * ab[1];
    if ab_len2 <= 1e-8 {
        return ((p[0] - a[0]).powi(2) + (p[1] - a[1]).powi(2)).sqrt();
    }
    let t = ((ap[0] * ab[0] + ap[1] * ab[1]) / ab_len2).clamp(0.0, 1.0);
    let q = [a[0] + ab[0] * t, a[1] + ab[1] * t];
    ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt()
}

struct AxisBuffers {
    vertex_buffer: wgpu::Buffer,
    tri_index_buffer: wgpu::Buffer,
    tri_index_count: u32,
    line_index_buffer: wgpu::Buffer,
    line_index_count: u32,
}

fn build_axis_arrow_mesh(
    direction: glam::Vec3,
    shaft_length: f32,
    cone_radius: f32,
    cone_height: f32,
    segments: u32,
) -> (Vec<[f32; 3]>, Vec<u32>, Vec<u32>) {
    let dir = direction.normalize();

    // Build a local frame: dir is the axis, u and v are perpendicular
    let arbitrary = if dir.y.abs() < 0.9 {
        glam::Vec3::Y
    } else {
        glam::Vec3::X
    };
    let u = dir.cross(arbitrary).normalize();
    let v = dir.cross(u);

    // Shaft: two vertices (origin → shaft_length along dir)
    let shaft_end = dir * shaft_length;
    let mut vertices = vec![[0.0, 0.0, 0.0], shaft_end.to_array()];

    // Line indices for the shaft
    let line_indices = vec![0, 1];

    // Cone: base ring at shaft_end, tip at shaft_end + cone_height * dir
    let tip = shaft_end + dir * cone_height;
    let tip_idx = vertices.len() as u32;
    vertices.push(tip.to_array());

    let base_start = vertices.len() as u32;
    for i in 0..segments {
        let angle = 2.0 * std::f32::consts::PI * (i as f32) / (segments as f32);
        let point = shaft_end + u * (angle.cos() * cone_radius) + v * (angle.sin() * cone_radius);
        vertices.push(point.to_array());
    }

    // Triangle fan for cone
    let mut tri_indices = Vec::new();
    for i in 0..segments {
        let curr = base_start + i;
        let next = base_start + (i + 1) % segments;
        // Side face
        tri_indices.push(tip_idx);
        tri_indices.push(curr);
        tri_indices.push(next);
    }

    // Base cap center
    let base_center_idx = vertices.len() as u32;
    vertices.push(shaft_end.to_array());
    for i in 0..segments {
        let curr = base_start + i;
        let next = base_start + (i + 1) % segments;
        tri_indices.push(base_center_idx);
        tri_indices.push(next);
        tri_indices.push(curr);
    }

    (vertices, tri_indices, line_indices)
}

struct GridBuffers {
    vertex_buffer: wgpu::Buffer,
    major_index_buffer: wgpu::Buffer,
    major_index_count: u32,
    minor_index_buffer: wgpu::Buffer,
    minor_index_count: u32,
}

/// Build a grid on the XZ plane (Y=0). Returns (vertices, major_line_indices, minor_line_indices).
fn build_grid_mesh(
    extent: f32,
    major_spacing: f32,
    minor_spacing: f32,
) -> (Vec<[f32; 3]>, Vec<u32>, Vec<u32>) {
    let mut vertices = Vec::new();
    let mut major_indices = Vec::new();
    let mut minor_indices = Vec::new();

    // Generate lines along X (varying Z)
    let steps = ((extent * 2.0 / minor_spacing).round() as i32).max(1);
    for i in 0..=steps {
        let z = -extent + i as f32 * minor_spacing;
        // Snap to avoid floating point drift
        let z = (z / minor_spacing).round() * minor_spacing;
        if z.abs() > extent + 0.001 {
            continue;
        }

        let idx = vertices.len() as u32;
        vertices.push([-extent, 0.0, z]);
        vertices.push([extent, 0.0, z]);

        let is_major = (z / major_spacing).round() * major_spacing - z < minor_spacing * 0.1
            && z.abs() > 0.001; // skip origin (drawn by axes)
        let is_origin_line = z.abs() < 0.001;

        if is_origin_line {
            // Skip — the axis handles this
        } else if is_major {
            major_indices.push(idx);
            major_indices.push(idx + 1);
        } else {
            minor_indices.push(idx);
            minor_indices.push(idx + 1);
        }
    }

    // Generate lines along Z (varying X)
    for i in 0..=steps {
        let x = -extent + i as f32 * minor_spacing;
        let x = (x / minor_spacing).round() * minor_spacing;
        if x.abs() > extent + 0.001 {
            continue;
        }

        let idx = vertices.len() as u32;
        vertices.push([x, 0.0, -extent]);
        vertices.push([x, 0.0, extent]);

        let is_major = (x / major_spacing).round() * major_spacing - x < minor_spacing * 0.1
            && x.abs() > 0.001;
        let is_origin_line = x.abs() < 0.001;

        if is_origin_line {
            // Skip
        } else if is_major {
            major_indices.push(idx);
            major_indices.push(idx + 1);
        } else {
            minor_indices.push(idx);
            minor_indices.push(idx + 1);
        }
    }

    (vertices, major_indices, minor_indices)
}

pub struct WgpuRenderer {
    pipeline: wgpu::RenderPipeline,
    pipeline_depth: wgpu::RenderPipeline,
    pipeline_line: wgpu::RenderPipeline,
    pipeline_line_depth: wgpu::RenderPipeline,
    camera_buffer: wgpu::Buffer,
    camera_bind_group: wgpu::BindGroup,
    material_bind_group: wgpu::BindGroup,
    material_face_highlight_bind_group: wgpu::BindGroup,
    material_edge_highlight_bind_group: wgpu::BindGroup,
    vertex_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    index_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    index_count: std::sync::Mutex<u32>,
    line_index_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    line_index_count: std::sync::Mutex<u32>,
    highlight_face_vertex_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    highlight_face_index_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    highlight_face_index_count: std::sync::Mutex<u32>,
    highlight_edge_vertex_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    highlight_edge_index_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    highlight_edge_index_count: std::sync::Mutex<u32>,
    depth_texture: std::sync::Mutex<Option<wgpu::Texture>>,
    depth_view: std::sync::Mutex<Option<wgpu::TextureView>>,
    depth_size: std::sync::Mutex<(u32, u32)>,
    axes_buffers: [AxisBuffers; 3],
    axes_material_bind_groups: [wgpu::BindGroup; 3],
    show_axes: std::sync::Mutex<bool>,
    display_mode: std::sync::Mutex<DisplayMode>,
    material_transparent_bind_group: wgpu::BindGroup,
    material_buffer: wgpu::Buffer,
    grid: GridBuffers,
    grid_major_material_bind_group: wgpu::BindGroup,
    grid_minor_material_bind_group: wgpu::BindGroup,
    show_grid: std::sync::Mutex<bool>,
    light_dir: std::sync::Mutex<glam::Vec3>,
}

unsafe impl Send for WgpuRenderer {}
unsafe impl Sync for WgpuRenderer {}

impl WgpuRenderer {
    pub fn default_clear_color() -> wgpu::Color {
        wgpu::Color {
            r: 0.07,
            g: 0.07,
            b: 0.11,
            a: 1.0,
        }
    }

    fn create_pipeline(
        device: &wgpu::Device,
        shader: &wgpu::ShaderModule,
        surface_format: wgpu::TextureFormat,
        layout: &wgpu::PipelineLayout,
        topology: wgpu::PrimitiveTopology,
        with_depth: bool,
    ) -> wgpu::RenderPipeline {
        device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some(if with_depth {
                "Render Pipeline (Depth)"
            } else {
                "Render Pipeline"
            }),
            layout: Some(layout),
            cache: None,
            vertex: wgpu::VertexState {
                module: shader,
                entry_point: Some("vs_main"),
                compilation_options: Default::default(),
                buffers: &[wgpu::VertexBufferLayout {
                    array_stride: std::mem::size_of::<[f32; 3]>() as wgpu::BufferAddress,
                    step_mode: wgpu::VertexStepMode::Vertex,
                    attributes: &wgpu::vertex_attr_array![0 => Float32x3],
                }],
            },
            fragment: Some(wgpu::FragmentState {
                module: shader,
                entry_point: Some("fs_main"),
                compilation_options: Default::default(),
                targets: &[Some(wgpu::ColorTargetState {
                    format: surface_format,
                    blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive: wgpu::PrimitiveState {
                topology,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: None,
                polygon_mode: wgpu::PolygonMode::Fill,
                unclipped_depth: false,
                conservative: false,
            },
            depth_stencil: if with_depth {
                Some(wgpu::DepthStencilState {
                    format: wgpu::TextureFormat::Depth32Float,
                    depth_write_enabled: true,
                    depth_compare: wgpu::CompareFunction::LessEqual,
                    stencil: wgpu::StencilState::default(),
                    bias: wgpu::DepthBiasState::default(),
                })
            } else {
                None
            },
            multisample: wgpu::MultisampleState {
                count: 1,
                mask: !0,
                alpha_to_coverage_enabled: false,
            },
            multiview: None,
        })
    }

    pub fn new(device: &wgpu::Device, surface_format: wgpu::TextureFormat) -> Self {
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("RCAD Shader"),
            source: wgpu::ShaderSource::Wgsl(include_str!("shader.wgsl").into()),
        });

        let camera_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Camera Buffer"),
            contents: bytemuck::cast_slice(&[CameraUniform {
                view_proj: glam::Mat4::IDENTITY.to_cols_array_2d(),
                eye_pos: [0.0, 0.0, 3.0, 1.0],
                light_dir: [0.45, 0.85, 0.35, 0.0],
            }]),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        let camera_bind_group_layout =
            device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("Camera Bind Group Layout"),
                entries: &[wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                }],
            });

        let camera_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Camera Bind Group"),
            layout: &camera_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: camera_buffer.as_entire_binding(),
            }],
        });

        let material_bind_group_layout =
            device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("Material Bind Group Layout"),
                entries: &[wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                }],
            });

        let material_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Material Buffer"),
            contents: bytemuck::cast_slice(&[MaterialUniform {
                color: [0.18, 0.64, 0.96, 1.0],
                flags: [0.0, 0.0, 0.0, 0.0],
            }]),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let material_transparent_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Transparent Material Buffer"),
            contents: bytemuck::cast_slice(&[MaterialUniform {
                color: [0.18, 0.64, 0.96, 0.3],
                flags: [0.0, 0.0, 0.0, 0.0],
            }]),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let material_face_highlight_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Face Highlight Material Buffer"),
            contents: bytemuck::cast_slice(&[MaterialUniform {
                color: [1.0, 0.45, 0.05, 0.45],
                flags: [1.0, 0.0, 0.0, 0.0],
            }]),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let material_edge_highlight_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Edge Highlight Material Buffer"),
            contents: bytemuck::cast_slice(&[MaterialUniform {
                color: [1.0, 0.95, 0.1, 1.0],
                flags: [1.0, 0.0, 0.0, 0.0],
            }]),
            usage: wgpu::BufferUsages::UNIFORM,
        });

        let material_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Material Bind Group"),
            layout: &material_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: material_buffer.as_entire_binding(),
            }],
        });
        let material_face_highlight_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Face Highlight Material Bind Group"),
            layout: &material_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: material_face_highlight_buffer.as_entire_binding(),
            }],
        });
        let material_edge_highlight_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Edge Highlight Material Bind Group"),
            layout: &material_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: material_edge_highlight_buffer.as_entire_binding(),
            }],
        });
        let material_transparent_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Transparent Material Bind Group"),
            layout: &material_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: material_transparent_buffer.as_entire_binding(),
            }],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Render Pipeline Layout"),
            bind_group_layouts: &[&camera_bind_group_layout, &material_bind_group_layout],
            push_constant_ranges: &[],
        });

        let pipeline = Self::create_pipeline(
            device,
            &shader,
            surface_format,
            &pipeline_layout,
            wgpu::PrimitiveTopology::TriangleList,
            false,
        );
        let pipeline_depth = Self::create_pipeline(
            device,
            &shader,
            surface_format,
            &pipeline_layout,
            wgpu::PrimitiveTopology::TriangleList,
            true,
        );
        let pipeline_line = Self::create_pipeline(
            device,
            &shader,
            surface_format,
            &pipeline_layout,
            wgpu::PrimitiveTopology::LineList,
            false,
        );
        let pipeline_line_depth = Self::create_pipeline(
            device,
            &shader,
            surface_format,
            &pipeline_layout,
            wgpu::PrimitiveTopology::LineList,
            true,
        );

        // Build background grid
        let grid_major_material_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Grid Major Material Buffer"),
            contents: bytemuck::cast_slice(&[MaterialUniform {
                color: [0.35, 0.35, 0.35, 0.5],
                flags: [1.0, 0.0, 0.0, 0.0], // unlit
            }]),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let grid_minor_material_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Grid Minor Material Buffer"),
            contents: bytemuck::cast_slice(&[MaterialUniform {
                color: [0.25, 0.25, 0.25, 0.3],
                flags: [1.0, 0.0, 0.0, 0.0], // unlit
            }]),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let grid_major_material_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Grid Major Material Bind Group"),
            layout: &material_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: grid_major_material_buffer.as_entire_binding(),
            }],
        });
        let grid_minor_material_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("Grid Minor Material Bind Group"),
            layout: &material_bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: grid_minor_material_buffer.as_entire_binding(),
            }],
        });

        let (grid_verts, grid_major_idx, grid_minor_idx) = build_grid_mesh(5.0, 1.0, 0.2);
        let grid = GridBuffers {
            vertex_buffer: device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Grid Vertex Buffer"),
                contents: bytemuck::cast_slice(&grid_verts),
                usage: wgpu::BufferUsages::VERTEX,
            }),
            major_index_buffer: device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Grid Major Index Buffer"),
                contents: bytemuck::cast_slice(&grid_major_idx),
                usage: wgpu::BufferUsages::INDEX,
            }),
            major_index_count: grid_major_idx.len() as u32,
            minor_index_buffer: device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("Grid Minor Index Buffer"),
                contents: bytemuck::cast_slice(&grid_minor_idx),
                usage: wgpu::BufferUsages::INDEX,
            }),
            minor_index_count: grid_minor_idx.len() as u32,
        };

        // Build axis arrows (X=red, Y=green, Z=blue)
        let axis_colors: [[f32; 4]; 3] = [
            [1.0, 0.2, 0.2, 1.0],  // X — red
            [0.2, 1.0, 0.2, 1.0],  // Y — green
            [0.3, 0.5, 1.0, 1.0],  // Z — blue
        ];
        let axis_dirs = [glam::Vec3::X, glam::Vec3::Y, glam::Vec3::Z];
        let axis_names = ["X", "Y", "Z"];

        let mut axes_material_bind_groups_vec = Vec::with_capacity(3);
        let mut axes_buffers_vec = Vec::with_capacity(3);

        for i in 0..3 {
            let mat_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(&format!("Axis {} Material Buffer", axis_names[i])),
                contents: bytemuck::cast_slice(&[MaterialUniform {
                    color: axis_colors[i],
                    flags: [1.0, 0.0, 0.0, 0.0], // unlit
                }]),
                usage: wgpu::BufferUsages::UNIFORM,
            });
            let bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some(&format!("Axis {} Material Bind Group", axis_names[i])),
                layout: &material_bind_group_layout,
                entries: &[wgpu::BindGroupEntry {
                    binding: 0,
                    resource: mat_buf.as_entire_binding(),
                }],
            });
            axes_material_bind_groups_vec.push(bg);

            let (verts, tri_idx, line_idx) =
                build_axis_arrow_mesh(axis_dirs[i], 1.0, 0.03, 0.1, 8);

            let vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(&format!("Axis {} Vertex Buffer", axis_names[i])),
                contents: bytemuck::cast_slice(&verts),
                usage: wgpu::BufferUsages::VERTEX,
            });
            let tri_index_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(&format!("Axis {} Tri Index Buffer", axis_names[i])),
                contents: bytemuck::cast_slice(&tri_idx),
                usage: wgpu::BufferUsages::INDEX,
            });
            let line_index_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(&format!("Axis {} Line Index Buffer", axis_names[i])),
                contents: bytemuck::cast_slice(&line_idx),
                usage: wgpu::BufferUsages::INDEX,
            });

            axes_buffers_vec.push(AxisBuffers {
                vertex_buffer,
                tri_index_buffer,
                tri_index_count: tri_idx.len() as u32,
                line_index_buffer,
                line_index_count: line_idx.len() as u32,
            });
        }

        // Convert Vecs to fixed-size arrays
        let axes_material_bind_groups = {
            let mut iter = axes_material_bind_groups_vec.into_iter();
            [
                iter.next().unwrap(),
                iter.next().unwrap(),
                iter.next().unwrap(),
            ]
        };
        let axes_buffers = {
            let mut iter = axes_buffers_vec.into_iter();
            [
                iter.next().unwrap(),
                iter.next().unwrap(),
                iter.next().unwrap(),
            ]
        };

        Self {
            pipeline,
            pipeline_depth,
            pipeline_line,
            pipeline_line_depth,
            camera_buffer,
            camera_bind_group,
            material_bind_group,
            material_face_highlight_bind_group,
            material_edge_highlight_bind_group,
            vertex_buffer: std::sync::Mutex::new(None),
            index_buffer: std::sync::Mutex::new(None),
            index_count: std::sync::Mutex::new(0),
            line_index_buffer: std::sync::Mutex::new(None),
            line_index_count: std::sync::Mutex::new(0),
            highlight_face_vertex_buffer: std::sync::Mutex::new(None),
            highlight_face_index_buffer: std::sync::Mutex::new(None),
            highlight_face_index_count: std::sync::Mutex::new(0),
            highlight_edge_vertex_buffer: std::sync::Mutex::new(None),
            highlight_edge_index_buffer: std::sync::Mutex::new(None),
            highlight_edge_index_count: std::sync::Mutex::new(0),
            depth_texture: std::sync::Mutex::new(None),
            depth_view: std::sync::Mutex::new(None),
            depth_size: std::sync::Mutex::new((0, 0)),
            axes_buffers,
            axes_material_bind_groups,
            show_axes: std::sync::Mutex::new(true),
            display_mode: std::sync::Mutex::new(DisplayMode::default()),
            material_transparent_bind_group,
            material_buffer,
            grid,
            grid_major_material_bind_group,
            grid_minor_material_bind_group,
            show_grid: std::sync::Mutex::new(true),
            light_dir: std::sync::Mutex::new(glam::Vec3::new(0.45, 0.85, 0.35)),
        }
    }

    pub fn ensure_depth_texture(&self, device: &wgpu::Device, width: u32, height: u32) {
        let width = width.max(1);
        let height = height.max(1);

        {
            let size = self.depth_size.lock().unwrap();
            let has_view = self.depth_view.lock().unwrap().is_some();
            if has_view && *size == (width, height) {
                return;
            }
        }

        let texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("RCAD Depth Texture"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Depth32Float,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        });
        let view = texture.create_view(&wgpu::TextureViewDescriptor::default());

        *self.depth_texture.lock().unwrap() = Some(texture);
        *self.depth_view.lock().unwrap() = Some(view);
        *self.depth_size.lock().unwrap() = (width, height);
    }

    pub fn prepare_scene(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        mesh: &Mesh,
        camera: &Camera,
        aspect: f32,
    ) {
        self.upload_mesh(device, mesh);
        self.update_camera(queue, camera, aspect.max(0.001));
    }

    pub fn prepare_scene_with_depth(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        mesh: &Mesh,
        camera: &Camera,
        aspect: f32,
        depth_size: (u32, u32),
    ) {
        self.ensure_depth_texture(device, depth_size.0, depth_size.1);
        self.prepare_scene(device, queue, mesh, camera, aspect);
    }

    pub fn upload_mesh(&self, device: &wgpu::Device, mesh: &Mesh) {
        *self.vertex_buffer.lock().unwrap() = Some(device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Vertex Buffer"),
            contents: bytemuck::cast_slice(&mesh.vertices),
            usage: wgpu::BufferUsages::VERTEX,
        }));

        *self.index_buffer.lock().unwrap() = Some(device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Index Buffer"),
            contents: bytemuck::cast_slice(&mesh.indices),
            usage: wgpu::BufferUsages::INDEX,
        }));

        *self.index_count.lock().unwrap() = mesh.indices.len() as u32;

        if mesh.line_indices.is_empty() {
            *self.line_index_buffer.lock().unwrap() = None;
            *self.line_index_count.lock().unwrap() = 0;
        } else {
            *self.line_index_buffer.lock().unwrap() = Some(device.create_buffer_init(
                &wgpu::util::BufferInitDescriptor {
                    label: Some("Line Index Buffer"),
                    contents: bytemuck::cast_slice(&mesh.line_indices),
                    usage: wgpu::BufferUsages::INDEX,
                },
            ));
            *self.line_index_count.lock().unwrap() = mesh.line_indices.len() as u32;
        }
    }

    pub fn upload_highlights(
        &self,
        device: &wgpu::Device,
        face_mesh: Option<&Mesh>,
        edge_mesh: Option<&Mesh>,
    ) {
        if let Some(mesh) = face_mesh {
            *self.highlight_face_vertex_buffer.lock().unwrap() = Some(device.create_buffer_init(
                &wgpu::util::BufferInitDescriptor {
                    label: Some("Highlight Face Vertex Buffer"),
                    contents: bytemuck::cast_slice(&mesh.vertices),
                    usage: wgpu::BufferUsages::VERTEX,
                },
            ));
            *self.highlight_face_index_buffer.lock().unwrap() = Some(device.create_buffer_init(
                &wgpu::util::BufferInitDescriptor {
                    label: Some("Highlight Face Index Buffer"),
                    contents: bytemuck::cast_slice(&mesh.indices),
                    usage: wgpu::BufferUsages::INDEX,
                },
            ));
            *self.highlight_face_index_count.lock().unwrap() = mesh.indices.len() as u32;
        } else {
            *self.highlight_face_vertex_buffer.lock().unwrap() = None;
            *self.highlight_face_index_buffer.lock().unwrap() = None;
            *self.highlight_face_index_count.lock().unwrap() = 0;
        }

        if let Some(mesh) = edge_mesh {
            *self.highlight_edge_vertex_buffer.lock().unwrap() = Some(device.create_buffer_init(
                &wgpu::util::BufferInitDescriptor {
                    label: Some("Highlight Edge Vertex Buffer"),
                    contents: bytemuck::cast_slice(&mesh.vertices),
                    usage: wgpu::BufferUsages::VERTEX,
                },
            ));
            *self.highlight_edge_index_buffer.lock().unwrap() = Some(device.create_buffer_init(
                &wgpu::util::BufferInitDescriptor {
                    label: Some("Highlight Edge Index Buffer"),
                    contents: bytemuck::cast_slice(&mesh.indices),
                    usage: wgpu::BufferUsages::INDEX,
                },
            ));
            *self.highlight_edge_index_count.lock().unwrap() = mesh.indices.len() as u32;
        } else {
            *self.highlight_edge_vertex_buffer.lock().unwrap() = None;
            *self.highlight_edge_index_buffer.lock().unwrap() = None;
            *self.highlight_edge_index_count.lock().unwrap() = 0;
        }
    }

    pub fn update_camera(&self, queue: &wgpu::Queue, camera: &Camera, aspect: f32) {
        let eye = camera.eye_position();
        let ld = *self.light_dir.lock().unwrap();
        let uniform = CameraUniform {
            view_proj: camera.build_view_projection_matrix(aspect),
            eye_pos: [eye.x, eye.y, eye.z, 1.0],
            light_dir: [ld.x, ld.y, ld.z, 0.0],
        };
        queue.write_buffer(&self.camera_buffer, 0, bytemuck::cast_slice(&[uniform]));
    }

    pub fn draw_in_render_pass(
        &self,
        render_pass: &mut wgpu::RenderPass<'_>,
        use_depth_pipeline: bool,
    ) {
        let mode = *self.display_mode.lock().unwrap();

        // Draw grid first (behind everything)
        if *self.show_grid.lock().unwrap() {
            self.draw_grid_in_render_pass(render_pass, use_depth_pipeline);
        }

        let vb_guard = self.vertex_buffer.lock().unwrap();
        let ib_guard = self.index_buffer.lock().unwrap();
        let count = *self.index_count.lock().unwrap();
        let lib_guard = self.line_index_buffer.lock().unwrap();
        let lcount = *self.line_index_count.lock().unwrap();

        // Draw model based on display mode
        let draw_triangles = matches!(mode, DisplayMode::Solid | DisplayMode::SolidWithEdges | DisplayMode::Transparent);
        let draw_wireframe = matches!(mode, DisplayMode::Wireframe | DisplayMode::SolidWithEdges | DisplayMode::Transparent);

        // In transparent mode, draw wireframe first so it's behind the translucent surface
        if mode == DisplayMode::Transparent && draw_wireframe && lcount > 0 {
            if let (Some(vb), Some(lib)) = (vb_guard.as_ref(), lib_guard.as_ref()) {
                if use_depth_pipeline {
                    render_pass.set_pipeline(&self.pipeline_line_depth);
                } else {
                    render_pass.set_pipeline(&self.pipeline_line);
                }
                render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
                render_pass.set_bind_group(1, &self.material_bind_group, &[]);
                render_pass.set_vertex_buffer(0, vb.slice(..));
                render_pass.set_index_buffer(lib.slice(..), wgpu::IndexFormat::Uint32);
                render_pass.draw_indexed(0..lcount, 0, 0..1);
            }
        }

        // Draw triangles
        if draw_triangles {
            if let (Some(vb), Some(ib)) = (vb_guard.as_ref(), ib_guard.as_ref()) {
                if use_depth_pipeline {
                    render_pass.set_pipeline(&self.pipeline_depth);
                } else {
                    render_pass.set_pipeline(&self.pipeline);
                }
                render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
                let mat = if mode == DisplayMode::Transparent {
                    &self.material_transparent_bind_group
                } else {
                    &self.material_bind_group
                };
                render_pass.set_bind_group(1, mat, &[]);
                render_pass.set_vertex_buffer(0, vb.slice(..));
                render_pass.set_index_buffer(ib.slice(..), wgpu::IndexFormat::Uint32);
                render_pass.draw_indexed(0..count, 0, 0..1);
            }
        }

        // Draw wireframe (non-transparent modes)
        if draw_wireframe && mode != DisplayMode::Transparent && lcount > 0 {
            if let (Some(vb), Some(lib)) = (vb_guard.as_ref(), lib_guard.as_ref()) {
                if use_depth_pipeline {
                    render_pass.set_pipeline(&self.pipeline_line_depth);
                } else {
                    render_pass.set_pipeline(&self.pipeline_line);
                }
                render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
                render_pass.set_bind_group(1, &self.material_bind_group, &[]);
                render_pass.set_vertex_buffer(0, vb.slice(..));
                render_pass.set_index_buffer(lib.slice(..), wgpu::IndexFormat::Uint32);
                render_pass.draw_indexed(0..lcount, 0, 0..1);
            }
        }

        // Draw face highlights
        let hvb_guard = self.highlight_face_vertex_buffer.lock().unwrap();
        let hib_guard = self.highlight_face_index_buffer.lock().unwrap();
        let hcount = *self.highlight_face_index_count.lock().unwrap();
        if hcount > 0
            && let (Some(vb), Some(ib)) = (hvb_guard.as_ref(), hib_guard.as_ref())
        {
            if use_depth_pipeline {
                render_pass.set_pipeline(&self.pipeline_depth);
            } else {
                render_pass.set_pipeline(&self.pipeline);
            }
            render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
            render_pass.set_bind_group(1, &self.material_face_highlight_bind_group, &[]);
            render_pass.set_vertex_buffer(0, vb.slice(..));
            render_pass.set_index_buffer(ib.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..hcount, 0, 0..1);
        }

        // Draw edge highlights
        let evb_guard = self.highlight_edge_vertex_buffer.lock().unwrap();
        let eib_guard = self.highlight_edge_index_buffer.lock().unwrap();
        let ecount = *self.highlight_edge_index_count.lock().unwrap();
        if ecount > 0
            && let (Some(vb), Some(ib)) = (evb_guard.as_ref(), eib_guard.as_ref())
        {
            if use_depth_pipeline {
                render_pass.set_pipeline(&self.pipeline_line_depth);
            } else {
                render_pass.set_pipeline(&self.pipeline_line);
            }
            render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
            render_pass.set_bind_group(1, &self.material_edge_highlight_bind_group, &[]);
            render_pass.set_vertex_buffer(0, vb.slice(..));
            render_pass.set_index_buffer(ib.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..ecount, 0, 0..1);
        }

        // Draw coordinate axes
        if *self.show_axes.lock().unwrap() {
            self.draw_axes_in_render_pass(render_pass, use_depth_pipeline);
        }
    }

    fn draw_axes_in_render_pass(
        &self,
        render_pass: &mut wgpu::RenderPass<'_>,
        use_depth_pipeline: bool,
    ) {
        for i in 0..3 {
            let axis = &self.axes_buffers[i];

            // Draw cone (triangles)
            if use_depth_pipeline {
                render_pass.set_pipeline(&self.pipeline_depth);
            } else {
                render_pass.set_pipeline(&self.pipeline);
            }
            render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
            render_pass.set_bind_group(1, &self.axes_material_bind_groups[i], &[]);
            render_pass.set_vertex_buffer(0, axis.vertex_buffer.slice(..));
            render_pass.set_index_buffer(axis.tri_index_buffer.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..axis.tri_index_count, 0, 0..1);

            // Draw shaft (line)
            if use_depth_pipeline {
                render_pass.set_pipeline(&self.pipeline_line_depth);
            } else {
                render_pass.set_pipeline(&self.pipeline_line);
            }
            render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
            render_pass.set_bind_group(1, &self.axes_material_bind_groups[i], &[]);
            render_pass.set_vertex_buffer(0, axis.vertex_buffer.slice(..));
            render_pass.set_index_buffer(axis.line_index_buffer.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..axis.line_index_count, 0, 0..1);
        }
    }

    fn draw_grid_in_render_pass(
        &self,
        render_pass: &mut wgpu::RenderPass<'_>,
        use_depth_pipeline: bool,
    ) {
        if use_depth_pipeline {
            render_pass.set_pipeline(&self.pipeline_line_depth);
        } else {
            render_pass.set_pipeline(&self.pipeline_line);
        }
        render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
        render_pass.set_vertex_buffer(0, self.grid.vertex_buffer.slice(..));

        // Draw minor lines first (thinner/dimmer)
        if self.grid.minor_index_count > 0 {
            render_pass.set_bind_group(1, &self.grid_minor_material_bind_group, &[]);
            render_pass.set_index_buffer(self.grid.minor_index_buffer.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..self.grid.minor_index_count, 0, 0..1);
        }

        // Draw major lines on top
        if self.grid.major_index_count > 0 {
            render_pass.set_bind_group(1, &self.grid_major_material_bind_group, &[]);
            render_pass.set_index_buffer(self.grid.major_index_buffer.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..self.grid.major_index_count, 0, 0..1);
        }
    }

    pub fn set_show_axes(&self, show: bool) {
        *self.show_axes.lock().unwrap() = show;
    }

    pub fn show_axes(&self) -> bool {
        *self.show_axes.lock().unwrap()
    }

    pub fn set_display_mode(&self, mode: DisplayMode) {
        *self.display_mode.lock().unwrap() = mode;
    }

    pub fn display_mode(&self) -> DisplayMode {
        *self.display_mode.lock().unwrap()
    }

    pub fn set_show_grid(&self, show: bool) {
        *self.show_grid.lock().unwrap() = show;
    }

    pub fn show_grid(&self) -> bool {
        *self.show_grid.lock().unwrap()
    }

    pub fn set_model_color(&self, queue: &wgpu::Queue, color: [f32; 4]) {
        let uniform = MaterialUniform {
            color,
            flags: [0.0, 0.0, 0.0, 0.0],
        };
        queue.write_buffer(&self.material_buffer, 0, bytemuck::cast_slice(&[uniform]));
    }

    pub fn set_model_unlit(&self, queue: &wgpu::Queue, color: [f32; 4], unlit: bool) {
        let uniform = MaterialUniform {
            color,
            flags: [if unlit { 1.0 } else { 0.0 }, 0.0, 0.0, 0.0],
        };
        queue.write_buffer(&self.material_buffer, 0, bytemuck::cast_slice(&[uniform]));
    }

    pub fn set_light_direction(&self, dir: glam::Vec3) {
        *self.light_dir.lock().unwrap() = dir;
    }

    pub fn light_direction(&self) -> glam::Vec3 {
        *self.light_dir.lock().unwrap()
    }

    /// Set light direction to match the camera eye direction (headlight mode).
    pub fn set_headlight(&self, camera: &Camera) {
        let eye = camera.eye_position();
        let dir = (eye - camera.target).normalize_or_zero();
        if dir.length_squared() > 1e-6 {
            *self.light_dir.lock().unwrap() = dir;
        }
    }

    pub fn render(
        &self,
        view: &wgpu::TextureView,
        encoder: &mut wgpu::CommandEncoder,
        clear_color: wgpu::Color,
        clip_bounds: Option<(u32, u32, u32, u32)>,
    ) {
        let use_depth = clip_bounds.is_some();
        let depth_view_guard = self.depth_view.lock().unwrap();
        let depth_attachment = if use_depth {
            depth_view_guard
                .as_ref()
                .map(|depth_view| wgpu::RenderPassDepthStencilAttachment {
                    view: depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Discard,
                    }),
                    stencil_ops: None,
                })
        } else {
            None
        };

        let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("Render Pass"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: if clip_bounds.is_some() {
                        wgpu::LoadOp::Load
                    } else {
                        wgpu::LoadOp::Clear(clear_color)
                    },
                    store: wgpu::StoreOp::Store,
                },
                depth_slice: None,
            })],
            depth_stencil_attachment: depth_attachment,
            timestamp_writes: None,
            occlusion_query_set: None,
        });

        let use_depth_pipeline = use_depth && depth_view_guard.is_some();

        if let Some((x, y, width, height)) = clip_bounds {
            if width > 0 && height > 0 {
                render_pass.set_scissor_rect(x, y, width.max(1), height.max(1));
            }
        }

        self.draw_in_render_pass(&mut render_pass, use_depth_pipeline);
    }

    pub fn render_with_defaults(
        &self,
        view: &wgpu::TextureView,
        encoder: &mut wgpu::CommandEncoder,
        clip_bounds: Option<(u32, u32, u32, u32)>,
    ) {
        self.render(view, encoder, Self::default_clear_color(), clip_bounds);
    }

    /// Render the current scene to an offscreen texture and return it as an RGBA image.
    pub fn screenshot(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        camera: &Camera,
        mesh: &Mesh,
        width: u32,
        height: u32,
    ) -> image::RgbaImage {
        let width = width.max(1);
        let height = height.max(1);
        let aspect = width as f32 / height as f32;

        // Prepare scene data
        self.upload_mesh(device, mesh);
        self.update_camera(queue, camera, aspect);

        // Create offscreen color texture
        let color_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Screenshot Color Texture"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Rgba8UnormSrgb,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });
        let color_view = color_texture.create_view(&wgpu::TextureViewDescriptor::default());

        // Create offscreen depth texture
        let depth_texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Screenshot Depth Texture"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Depth32Float,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        });
        let depth_view = depth_texture.create_view(&wgpu::TextureViewDescriptor::default());

        // Render
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("Screenshot Encoder"),
        });

        {
            let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("Screenshot Render Pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &color_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(Self::default_clear_color()),
                        store: wgpu::StoreOp::Store,
                    },
                    depth_slice: None,
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Discard,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            self.draw_in_render_pass(&mut render_pass, true);
        }

        // Copy texture to staging buffer
        let bytes_per_pixel = 4u32;
        // wgpu requires rows to be aligned to 256 bytes
        let unpadded_bytes_per_row = width * bytes_per_pixel;
        let padded_bytes_per_row =
            (unpadded_bytes_per_row + 255) & !255;
        let buffer_size = (padded_bytes_per_row * height) as u64;

        let staging_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Screenshot Staging Buffer"),
            size: buffer_size,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        encoder.copy_texture_to_buffer(
            wgpu::TexelCopyTextureInfo {
                texture: &color_texture,
                mip_level: 0,
                origin: wgpu::Origin3d::ZERO,
                aspect: wgpu::TextureAspect::All,
            },
            wgpu::TexelCopyBufferInfo {
                buffer: &staging_buffer,
                layout: wgpu::TexelCopyBufferLayout {
                    offset: 0,
                    bytes_per_row: Some(padded_bytes_per_row),
                    rows_per_image: Some(height),
                },
            },
            wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
        );

        queue.submit(std::iter::once(encoder.finish()));

        // Map and read back
        let buffer_slice = staging_buffer.slice(..);
        let (sender, receiver) = std::sync::mpsc::channel();
        buffer_slice.map_async(wgpu::MapMode::Read, move |result| {
            sender.send(result).unwrap();
        });
        device.poll(wgpu::PollType::wait_indefinitely()).unwrap();
        receiver.recv().unwrap().unwrap();

        let data = buffer_slice.get_mapped_range();
        let mut pixels = Vec::with_capacity((width * height * bytes_per_pixel) as usize);
        for row in 0..height {
            let start = (row * padded_bytes_per_row) as usize;
            let end = start + (width * bytes_per_pixel) as usize;
            pixels.extend_from_slice(&data[start..end]);
        }
        drop(data);
        staging_buffer.unmap();

        image::RgbaImage::from_raw(width, height, pixels)
            .expect("pixel buffer size matches image dimensions")
    }

    /// Render the scene and save to a PNG file.
    pub fn screenshot_to_file(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        camera: &Camera,
        mesh: &Mesh,
        width: u32,
        height: u32,
        path: &std::path::Path,
    ) -> Result<(), String> {
        let img = self.screenshot(device, queue, camera, mesh, width, height);
        img.save(path).map_err(|e| e.to_string())
    }
}
