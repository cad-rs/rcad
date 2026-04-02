use iced::widget::{button, checkbox, column, container, row, text};
use iced::{Element, Length, Task};
use rcad_kernel::BRep;
use rcad_render::{Camera, Mesh, Tessellator, WgpuRenderer};

const SAMPLE_STEP: &str = include_str!("../../../assets/box.step");

// ─── App ─────────────────────────────────────────────────────────────────────

pub struct RCadApp {
    brep: BRep,
    mesh: Mesh,
    camera: Camera,
    auto_rotate: bool,
}

#[derive(Debug, Clone)]
pub enum Message {
    RotateCamera(f32, f32),
    ZoomCamera(f32),
    ToggleAutoRotate(bool),
    ResetCamera,
    Tick(f32),
}

impl RCadApp {
    pub fn new(step_content: Option<String>) -> (Self, Task<Message>) {
        let brep = if let Some(content) = step_content {
            rcad_step::StepReader::parse_string(&content)
        } else {
            rcad_step::StepReader::parse_string(SAMPLE_STEP)
        }
        .unwrap_or_else(|_| BRep::create_box(1.0, 1.0, 1.0));
        let mesh = Tessellator::tessellate(&brep);

        (
            Self {
                brep,
                mesh,
                camera: Camera::new(),
                auto_rotate: true,
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
            Message::ZoomCamera(delta) => {
                self.camera.distance -= delta * 0.01 * self.camera.distance;
                self.camera.distance = self.camera.distance.clamp(1.0, 50.0);
            }
            Message::ToggleAutoRotate(on) => {
                self.auto_rotate = on;
            }
            Message::ResetCamera => {
                self.camera = Camera::new();
            }
            Message::Tick(dt) => {
                if self.auto_rotate {
                    self.camera.rot_y += dt * 0.6;
                }
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

        let info = container(column![
            text("RCAD  ·  iced").size(20),
            text("─────────────────"),
            text(format!("Vertices : {}", self.brep.vertices.len())),
            text(format!("Edges    : {}", self.brep.edges.len())),
            text(format!("Faces    : {}", face_count)),
            text(format!("Triangles: {}", self.mesh.indices.len() / 3)),
            text("─────────────────"),
            checkbox(self.auto_rotate).label("Auto-rotate").on_toggle(Message::ToggleAutoRotate),
            text("Drag to rotate"),
            button("Reset Camera").on_press(Message::ResetCamera),
        ]
        .spacing(8))
        .padding(12)
        .width(Length::Fixed(180.0))
        .height(Length::Fill);
        // .style(|_| container::Style {
        //     background: Some(Color::from_rgb(0.1, 0.1, 0.15).into()),
        //     border: Border {
        //         width: 0.0,
        //         color: Color::TRANSPARENT,
        //         radius: 0.0.into(),
        //     },
        //     ..Default::default()
        // });

        let viewport = container(iced::widget::shader(Scene {
            mesh: &self.mesh,
            camera: &self.camera,
        })
        .width(Length::Fill)
        .height(Length::Fill));


        row![info, viewport]
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }
}

// ─── Shader Integration ──────────────────────────────────────────────────────

struct Scene<'a> {
    mesh: &'a Mesh,
    camera: &'a Camera,
}

#[derive(Default)]
struct SceneState {
    is_dragging: bool,
    last_cursor_position: Option<iced::Point>,
}

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
        state: &mut Self::State,
        event: &iced::Event,
        _bounds: iced::Rectangle,
        cursor: iced::mouse::Cursor,
    ) -> Option<iced::widget::shader::Action<Message>> {
        match event {
            iced::Event::Mouse(iced::mouse::Event::ButtonPressed(iced::mouse::Button::Left)) => {
                if cursor.is_over(_bounds) {
                    state.is_dragging = true;
                    state.last_cursor_position = cursor.position_in(_bounds);
                }
            }
            iced::Event::Mouse(iced::mouse::Event::ButtonReleased(iced::mouse::Button::Left)) => {
                state.is_dragging = false;
            }
            iced::Event::Mouse(iced::mouse::Event::WheelScrolled { delta }) => {
                if cursor.is_over(_bounds) {
                    let y = match delta {
                        iced::mouse::ScrollDelta::Lines { y, .. } => *y * 20.0,
                        iced::mouse::ScrollDelta::Pixels { y, .. } => *y,
                    };
                    return Some(iced::widget::shader::Action::publish(Message::ZoomCamera(y)));
                }
            }
            iced::Event::Mouse(iced::mouse::Event::CursorMoved { .. }) => {
                if state.is_dragging {
                    if let Some(current_pos) = cursor.position_in(_bounds) {
                        if let Some(last_pos) = state.last_cursor_position {
                            let dx = current_pos.x - last_pos.x;
                            let dy = current_pos.y - last_pos.y;

                            state.last_cursor_position = Some(current_pos);

                            return Some(iced::widget::shader::Action::publish(
                                Message::RotateCamera(dx * 0.8, dy * 0.8),
                            ));
                        }
                        state.last_cursor_position = Some(current_pos);
                    }
                }
            }
            _ => {}
        }
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
        state: &Self::State,
        _bounds: iced::Rectangle,
        _cursor: iced::mouse::Cursor,
    ) -> iced::mouse::Interaction {
        if state.is_dragging {
            iced::mouse::Interaction::Grabbing
        } else if _cursor.is_over(_bounds) {
            iced::mouse::Interaction::Grab
        } else {
            iced::mouse::Interaction::default()
        }
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
        clip_bounds: &iced::Rectangle<u32>,
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
            Some((
                clip_bounds.x,
                clip_bounds.y,
                clip_bounds.width,
                clip_bounds.height,
            )),
        );
    }
}

impl Default for RCadApp {
    fn default() -> Self {
        RCadApp::new(None).0
    }
}

// ─── Native entry ────────────────────────────────────────────────────────────

pub fn run_native(step_content: Option<String>) -> iced::Result {
    iced::application(move || RCadApp::new(step_content.clone()), RCadApp::update, RCadApp::view)
        .title("RCAD Creator · iced")
        .window(iced::window::Settings {
            size: iced::Size::new(900.0, 600.0),
            ..Default::default()
        })
        .subscription(|_| {
            iced::time::every(std::time::Duration::from_millis(16))
                .map(|_| Message::Tick(0.016))
        })
        .run()
}

// ─── WASM entry ──────────────────────────────────────────────────────────────

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen(start)]
pub fn start() {
    console_error_panic_hook::set_once();
    run_native(None).expect("iced failed to start");
}
