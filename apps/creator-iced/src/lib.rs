use iced::widget::{column, row, text};
use iced::{Element, Length, Task};
use rcad_kernel::BRep;
use rcad_render::{Camera, Mesh, Tessellator, WgpuRenderer};

const SAMPLE_STEP: &str = include_str!("../../../assets/box.step");

// ─── App ─────────────────────────────────────────────────────────────────────

pub struct RCadApp {
    brep: BRep,
    mesh: Mesh,
    camera: Camera,
}

#[derive(Debug, Clone)]
pub enum Message {
    RotateCamera(f32, f32),
}

impl RCadApp {
    pub fn new() -> (Self, Task<Message>) {
        let brep = rcad_step::StepReader::parse_string(SAMPLE_STEP)
            .unwrap_or_else(|_| BRep::create_box(1.0, 1.0, 1.0));
        let mesh = Tessellator::tessellate(&brep);

        (
            Self {
                brep,
                mesh,
                camera: Camera::new(),
            },
            Task::none(),
        )
    }

    pub fn update(&mut self, msg: Message) -> Task<Message> {
        match msg {
            Message::RotateCamera(dx, dy) => {
                self.camera.rot_y += dx * 0.01;
                self.camera.rot_x += dy * 0.01;
            }
        }
        Task::none()
    }

    pub fn view(&self) -> Element<'_, Message> {
        let face_count: usize = self
            .brep
            .solids
            .iter()
            .flat_map(|s| &s.shells)
            .map(|sh| sh.faces.len())
            .sum();

        let info = column![
            text("RCAD · iced").size(20),
            text("─────────────────"),
            text(format!("Vertices : {}", self.brep.vertices.len())),
            text(format!("Edges    : {}", self.brep.edges.len())),
            text(format!("Faces    : {}", face_count)),
            text(format!("Triangles: {}", self.mesh.indices.len() / 3)),
            text("─────────────────"),
            text("Drag to rotate"),
        ]
        .spacing(4)
        .padding(12)
        .width(Length::Fixed(180.0));

        let viewport: Element<'_, Message> = iced::widget::shader(Scene {
            mesh: &self.mesh,
            camera: &self.camera,
        })
        .into();

        row![info, viewport].into()
    }
}

// ─── Shader Integration ──────────────────────────────────────────────────────

struct Scene<'a> {
    mesh: &'a Mesh,
    camera: &'a Camera,
}

#[derive(Default)]
struct SceneState {}

pub struct RCadPipeline {
    renderer: WgpuRenderer,
}

impl iced::widget::shader::Pipeline for RCadPipeline {
    fn new(
        device: &iced::wgpu::Device,
        _queue: &iced::wgpu::Queue,
        format: iced::wgpu::TextureFormat,
    ) -> Self {
        Self {
            renderer: WgpuRenderer::new(
                unsafe { std::mem::transmute(device) },
                unsafe { std::mem::transmute(format) },
            ),
        }
    }
}

impl<'a> iced::widget::shader::Program<Message> for Scene<'a> {
    type State = SceneState;
    type Primitive = Primitive;

    fn update(
        &self,
        _state: &mut Self::State,
        _event: &iced::Event,
        _bounds: iced::Rectangle,
        _cursor: iced::mouse::Cursor,
    ) -> Option<iced::widget::shader::Action<Message>> {
        None
    }

    fn draw(
        &self,
        _state: &Self::State,
        _cursor: iced::mouse::Cursor,
        _bounds: iced::Rectangle,
    ) -> Self::Primitive {
        Primitive {
            mesh: self.mesh.clone(),
            camera: *self.camera,
            aspect: _bounds.width / _bounds.height,
        }
    }

    fn mouse_interaction(
        &self,
        _state: &Self::State,
        _bounds: iced::Rectangle,
        _cursor: iced::mouse::Cursor,
    ) -> iced::mouse::Interaction {
        iced::mouse::Interaction::default()
    }
}

#[derive(Debug, Clone)]
struct Primitive {
    mesh: Mesh,
    camera: Camera,
    aspect: f32,
}

impl iced::widget::shader::Primitive for Primitive {
    type Pipeline = RCadPipeline;

    fn prepare(
        &self,
        pipeline: &mut Self::Pipeline,
        device: &iced::wgpu::Device,
        queue: &iced::wgpu::Queue,
        _bounds: &iced::Rectangle,
        _viewport: &iced::advanced::graphics::Viewport,
    ) {
        pipeline.renderer.update_camera(
            unsafe { std::mem::transmute(queue) },
            &self.camera,
            self.aspect,
        );
        pipeline.renderer.upload_mesh(unsafe { std::mem::transmute(device) }, &self.mesh);
    }

    fn render(
        &self,
        pipeline: &Self::Pipeline,
        encoder: &mut iced::wgpu::CommandEncoder,
        target: &iced::wgpu::TextureView,
        _clip_bounds: &iced::Rectangle<u32>,
    ) {
        let clear_color = iced::Color::from_rgb(0.07, 0.07, 0.11);
        pipeline.renderer.render(
            unsafe { std::mem::transmute(target) },
            unsafe { std::mem::transmute(encoder) },
            wgpu::Color {
                r: clear_color.r as f64,
                g: clear_color.g as f64,
                b: clear_color.b as f64,
                a: clear_color.a as f64,
            },
        );
    }
}

impl Default for RCadApp {
    fn default() -> Self {
        Self::new().0
    }
}

// ─── Native entry ────────────────────────────────────────────────────────────

pub fn run_native() -> iced::Result {
    iced::application(RCadApp::new, RCadApp::update, RCadApp::view)
        .title("RCAD Creator · iced")
        .run()
}

// ─── WASM entry ──────────────────────────────────────────────────────────────

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen(start)]
pub fn start() {
    console_error_panic_hook::set_once();
    run_native().expect("iced failed to start");
}
