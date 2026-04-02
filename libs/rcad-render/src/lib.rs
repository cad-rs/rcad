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
        let uniform = CameraUniform {
            view_proj: camera.build_view_projection_matrix(aspect),
            eye_pos: [eye.x, eye.y, eye.z, 1.0],
        };
        queue.write_buffer(&self.camera_buffer, 0, bytemuck::cast_slice(&[uniform]));
    }

    pub fn draw_in_render_pass(
        &self,
        render_pass: &mut wgpu::RenderPass<'_>,
        use_depth_pipeline: bool,
    ) {
        if use_depth_pipeline {
            render_pass.set_pipeline(&self.pipeline_depth);
        } else {
            render_pass.set_pipeline(&self.pipeline);
        }
        render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
        render_pass.set_bind_group(1, &self.material_bind_group, &[]);

        let vb_guard = self.vertex_buffer.lock().unwrap();
        let ib_guard = self.index_buffer.lock().unwrap();
        let count = *self.index_count.lock().unwrap();

        if let (Some(vb), Some(ib)) = (vb_guard.as_ref(), ib_guard.as_ref()) {
            render_pass.set_vertex_buffer(0, vb.slice(..));
            render_pass.set_index_buffer(ib.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..count, 0, 0..1);
        }

        let lib_guard = self.line_index_buffer.lock().unwrap();
        let lcount = *self.line_index_count.lock().unwrap();
        if lcount > 0
            && let (Some(vb), Some(lib)) = (vb_guard.as_ref(), lib_guard.as_ref())
        {
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

            if use_depth_pipeline {
                render_pass.set_pipeline(&self.pipeline_depth);
            } else {
                render_pass.set_pipeline(&self.pipeline);
            }
            render_pass.set_bind_group(0, &self.camera_bind_group, &[]);
            render_pass.set_bind_group(1, &self.material_bind_group, &[]);
        }

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
}
