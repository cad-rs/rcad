use rcad_kernel::BRep;
use wgpu::util::DeviceExt;

#[derive(Clone, Debug)]
pub struct Mesh {
    pub vertices: Vec<[f32; 3]>,
    pub indices: Vec<u32>,
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

        Mesh {
            vertices: flat_verts,
            indices,
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
}

impl Camera {
    pub fn new() -> Self {
        Self {
            rot_x: 0.4,
            rot_y: 0.5,
            distance: 3.0,
        }
    }

    pub fn build_view_projection_matrix(&self, aspect: f32) -> [[f32; 4]; 4] {
        let eye = self.eye_position();
        let target = glam::Vec3::ZERO;
        let up = glam::Vec3::Y;

        let view = glam::Mat4::look_at_rh(eye, target, up);
        let proj = glam::Mat4::perspective_rh(45.0_f32.to_radians(), aspect, 0.1, 100.0);

        (proj * view).to_cols_array_2d()
    }

    pub fn eye_position(&self) -> glam::Vec3 {
        glam::Vec3::new(
            self.distance * self.rot_y.cos() * self.rot_x.cos(),
            self.distance * self.rot_x.sin(),
            self.distance * self.rot_y.sin() * self.rot_x.cos(),
        )
    }
}

impl Default for Camera {
    fn default() -> Self {
        Self::new()
    }
}

pub struct WgpuRenderer {
    pipeline: wgpu::RenderPipeline,
    pipeline_depth: wgpu::RenderPipeline,
    camera_buffer: wgpu::Buffer,
    camera_bind_group: wgpu::BindGroup,
    vertex_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    index_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    index_count: std::sync::Mutex<u32>,
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
                    blend: Some(wgpu::BlendState::REPLACE),
                    write_mask: wgpu::ColorWrites::ALL,
                })],
            }),
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
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

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Render Pipeline Layout"),
            bind_group_layouts: &[&camera_bind_group_layout],
            push_constant_ranges: &[],
        });

        let pipeline =
            Self::create_pipeline(device, &shader, surface_format, &pipeline_layout, false);
        let pipeline_depth =
            Self::create_pipeline(device, &shader, surface_format, &pipeline_layout, true);

        Self {
            pipeline,
            pipeline_depth,
            camera_buffer,
            camera_bind_group,
            vertex_buffer: std::sync::Mutex::new(None),
            index_buffer: std::sync::Mutex::new(None),
            index_count: std::sync::Mutex::new(0),
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

        let vb_guard = self.vertex_buffer.lock().unwrap();
        let ib_guard = self.index_buffer.lock().unwrap();
        let count = *self.index_count.lock().unwrap();

        if let (Some(vb), Some(ib)) = (vb_guard.as_ref(), ib_guard.as_ref()) {
            render_pass.set_vertex_buffer(0, vb.slice(..));
            render_pass.set_index_buffer(ib.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..count, 0, 0..1);
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
