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
        let eye = glam::Vec3::new(
            self.distance * self.rot_y.cos() * self.rot_x.cos(),
            self.distance * self.rot_x.sin(),
            self.distance * self.rot_y.sin() * self.rot_x.cos(),
        );
        let target = glam::Vec3::ZERO;
        let up = glam::Vec3::Y;

        let view = glam::Mat4::look_at_rh(eye, target, up);
        let proj = glam::Mat4::perspective_rh(45.0_f32.to_radians(), aspect, 0.1, 100.0);

        (proj * view).to_cols_array_2d()
    }
}

impl Default for Camera {
    fn default() -> Self {
        Self::new()
    }
}

pub struct WgpuRenderer {
    pub pipeline: wgpu::RenderPipeline,
    pub camera_buffer: wgpu::Buffer,
    pub camera_bind_group: wgpu::BindGroup,
    pub vertex_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    pub index_buffer: std::sync::Mutex<Option<wgpu::Buffer>>,
    pub index_count: std::sync::Mutex<u32>,
}

unsafe impl Send for WgpuRenderer {}
unsafe impl Sync for WgpuRenderer {}

impl WgpuRenderer {
    pub fn new(device: &wgpu::Device, surface_format: wgpu::TextureFormat) -> Self {
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("RCAD Shader"),
            source: wgpu::ShaderSource::Wgsl(include_str!("shader.wgsl").into()),
        });

        let camera_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Camera Buffer"),
            contents: bytemuck::cast_slice(&[CameraUniform {
                view_proj: glam::Mat4::IDENTITY.to_cols_array_2d(),
            }]),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        let camera_bind_group_layout =
            device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("Camera Bind Group Layout"),
                entries: &[wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::VERTEX,
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

        let pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("Render Pipeline"),
            layout: Some(&pipeline_layout),
            cache: None,
            vertex: wgpu::VertexState {
                module: &shader,
                entry_point: Some("vs_main"),
                compilation_options: Default::default(),
                buffers: &[wgpu::VertexBufferLayout {
                    array_stride: std::mem::size_of::<[f32; 3]>() as wgpu::BufferAddress,
                    step_mode: wgpu::VertexStepMode::Vertex,
                    attributes: &wgpu::vertex_attr_array![0 => Float32x3],
                }],
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader,
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
            depth_stencil: None,
            multisample: wgpu::MultisampleState {
                count: 1,
                mask: !0,
                alpha_to_coverage_enabled: false,
            },
            multiview: None,
        });

        Self {
            pipeline,
            camera_buffer,
            camera_bind_group,
            vertex_buffer: std::sync::Mutex::new(None),
            index_buffer: std::sync::Mutex::new(None),
            index_count: std::sync::Mutex::new(0),
        }
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
        let uniform = CameraUniform {
            view_proj: camera.build_view_projection_matrix(aspect),
        };
        queue.write_buffer(&self.camera_buffer, 0, bytemuck::cast_slice(&[uniform]));
    }

    pub fn render(
        &self,
        view: &wgpu::TextureView,
        encoder: &mut wgpu::CommandEncoder,
        clear_color: wgpu::Color,
    ) {
        let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("Render Pass"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: wgpu::LoadOp::Clear(clear_color),
                    store: wgpu::StoreOp::Store,
                },
                depth_slice: None,
            })],
            depth_stencil_attachment: None,
            timestamp_writes: None,
            occlusion_query_set: None,
        });

        render_pass.set_pipeline(&self.pipeline);
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
}
