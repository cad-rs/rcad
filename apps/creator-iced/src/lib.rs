use iced::widget::{button, checkbox, column, container, row, text};
use iced::{Element, Length, Task};
use rcad_kernel::BRep;
use rcad_render::{
    build_edges_highlight_mesh, build_faces_highlight_mesh, Camera, Mesh, SelectionMode,
    SelectionState, Tessellator, WgpuRenderer, DEFAULT_EDGE_PICK_RADIUS_PX,
};
use rcad_step::writer::{ExportSelection, StepWriter};

const SAMPLE_STEP: &str = include_str!("../../../assets/box.step");

// ─── App ─────────────────────────────────────────────────────────────────────

pub struct RCadApp {
    brep: BRep,
    mesh: Mesh,
    camera: Camera,
    selection: SelectionState,
    export_status: Option<String>,
}

#[derive(Debug, Clone)]
pub enum Message {
    RotateCamera(f32, f32),
    PanCamera(f32, f32),
    ZoomCamera(f32),
    SelectAt(f32, f32, f32, f32),
    HoverAt(f32, f32, f32, f32),
    ClearHover,
    SetSelectionMode(SelectionMode),
    SetAdditiveSelect(bool),
    ExportStep,
    ResetCamera,
}

impl RCadApp {
    pub fn new(step_content: Option<String>) -> (Self, Task<Message>) {
        let parse_result = if let Some(content) = step_content {
            rcad_step::StepReader::parse_string(&content)
        } else {
            rcad_step::StepReader::parse_string(SAMPLE_STEP)
        };

        let brep = match parse_result {
            Ok(brep) => {
                eprintln!(
                    "[rcad-step][iced] parsed STEP: vertices={}, edges={}, solids={}",
                    brep.vertices.len(),
                    brep.edges.len(),
                    brep.solids.len()
                );
                brep
            }
            Err(err) => {
                eprintln!("[rcad-step][iced] parse failed, fallback to box: {err}");
                BRep::create_box(1.0, 1.0, 1.0)
            }
        };
        let mesh = Tessellator::tessellate(&brep);

        (
            Self {
                brep,
                mesh,
                camera: Camera::new(),
                selection: SelectionState::default(),
                export_status: None,
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
            Message::PanCamera(dx, dy) => {
                self.camera.pan_pixels(dx, dy);
            }
            Message::ZoomCamera(delta) => {
                self.camera.distance -= delta * 0.01 * self.camera.distance;
                self.camera.distance = self.camera.distance.clamp(1.0, 50.0);
            }
            Message::SelectAt(x, y, w, h) => {
                let viewport = [w.max(1.0), h.max(1.0)];
                let cursor = [x, y];
                let aspect = viewport[0] / viewport[1];
                self.selection
                    .click_at(
                        &self.brep,
                        &self.camera,
                        aspect,
                        viewport,
                        cursor,
                        DEFAULT_EDGE_PICK_RADIUS_PX,
                    );
            }
            Message::HoverAt(x, y, w, h) => {
                let viewport = [w.max(1.0), h.max(1.0)];
                let cursor = [x, y];
                let aspect = viewport[0] / viewport[1];
                self.selection
                    .hover_at(
                        &self.brep,
                        &self.camera,
                        aspect,
                        viewport,
                        cursor,
                        DEFAULT_EDGE_PICK_RADIUS_PX,
                    );
            }
            Message::ClearHover => {
                self.selection.clear_hover();
            }
            Message::SetSelectionMode(mode) => {
                self.selection.set_mode(mode);
            }
            Message::SetAdditiveSelect(v) => {
                self.selection.additive_select = v;
            }
            Message::ExportStep => {
                self.export_status = Some(match export_step_file(&self.brep, &self.selection) {
                    Ok(path) => format!("Exported: {}", path),
                    Err(err) => format!("Export failed: {}", err),
                });
            }
            Message::ResetCamera => {
                self.camera = Camera::new();
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
            text(format!("Curves   : {}", self.brep.geom.curves.len())),
            text(format!("Surfaces : {}", self.brep.geom.surfaces.len())),
            text("─────────────────"),
            text("Selection Mode"),
            row![
                button("Select Face").on_press(Message::SetSelectionMode(SelectionMode::Face)),
                button("Select Edge").on_press(Message::SetSelectionMode(SelectionMode::Edge))
            ]
            .spacing(6),
            checkbox(self.selection.additive_select)
                .label("Additive Select")
                .on_toggle(Message::SetAdditiveSelect),
            text(format!(
                "Selected Faces: {}",
                self.selection.selected_faces.len()
            )),
            text(format!(
                "Selected Edges: {}",
                self.selection.selected_edges.len()
            )),
            text(format!(
                "Hover Face: {}",
                self.selection.hovered_face
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string())
            )),
            text(format!(
                "Hover Edge: {}",
                self.selection.hovered_edge
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string())
            )),
            text(format!(
                "Hover Pos: {}",
                self.selection.last_hover_pos
                    .map(|(x, y)| format!("{x:.1}, {y:.1}"))
                    .unwrap_or_else(|| "-".to_string())
            )),
            text("Click to select, toggle Additive Select for multi-select"),
            text("Left drag: rotate, Middle drag: pan"),
            button("Export STEP").on_press(Message::ExportStep),
            text(self.export_status.clone().unwrap_or_default()),
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
            brep: &self.brep,
            mesh: &self.mesh,
            camera: &self.camera,
            selected_faces: self.selection.highlighted_faces(),
            selected_edges: self.selection.highlighted_edges(),
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
    brep: &'a BRep,
    mesh: &'a Mesh,
    camera: &'a Camera,
    selected_faces: Vec<usize>,
    selected_edges: Vec<usize>,
}

#[cfg(not(target_arch = "wasm32"))]
fn export_step_file(brep: &BRep, selection: &SelectionState) -> Result<String, String> {
    let step = StepWriter::write_string(
        brep,
        ExportSelection {
            selected_faces: &selection.selected_faces,
            selected_edges: &selection.selected_edges,
        },
    );
    let path = std::env::current_dir()
        .map_err(|e| e.to_string())?
        .join("rcad_export.step");
    std::fs::write(&path, step).map_err(|e| e.to_string())?;
    Ok(path.display().to_string())
}

#[cfg(target_arch = "wasm32")]
fn export_step_file(_brep: &BRep, _selection: &SelectionState) -> Result<String, String> {
    Err("STEP export is only available in the native app".to_string())
}

#[derive(Default)]
struct SceneState {
    is_rotating: bool,
    is_panning: bool,
    rotate_drag_distance: f32,
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
                    state.is_rotating = true;
                    state.rotate_drag_distance = 0.0;
                    state.last_cursor_position = local_cursor_position(cursor, _bounds);
                }
            }
            iced::Event::Mouse(iced::mouse::Event::ButtonPressed(iced::mouse::Button::Middle)) => {
                if cursor.is_over(_bounds) {
                    state.is_panning = true;
                    state.last_cursor_position = local_cursor_position(cursor, _bounds);
                }
            }
            iced::Event::Mouse(iced::mouse::Event::ButtonReleased(iced::mouse::Button::Left)) => {
                if state.is_rotating
                    && state.rotate_drag_distance < 3.0
                    && let Some(pos) = local_cursor_position(cursor, _bounds)
                {
                    state.is_rotating = false;
                    state.last_cursor_position = Some(pos);
                    return Some(iced::widget::shader::Action::publish(Message::SelectAt(
                        pos.x,
                        pos.y,
                        _bounds.width,
                        _bounds.height,
                    )));
                }
                state.is_rotating = false;
            }
            iced::Event::Mouse(iced::mouse::Event::ButtonReleased(iced::mouse::Button::Middle)) => {
                state.is_panning = false;
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
                if state.is_rotating {
                    if let Some(current_pos) = local_cursor_position(cursor, _bounds) {
                        if let Some(last_pos) = state.last_cursor_position {
                            let dx = current_pos.x - last_pos.x;
                            let dy = current_pos.y - last_pos.y;
                            state.rotate_drag_distance += (dx * dx + dy * dy).sqrt();

                            state.last_cursor_position = Some(current_pos);

                            return Some(iced::widget::shader::Action::publish(
                                Message::RotateCamera(dx * 0.8, dy * 0.8),
                            ));
                        }
                        state.last_cursor_position = Some(current_pos);
                    }
                } else if state.is_panning {
                    if let Some(current_pos) = local_cursor_position(cursor, _bounds) {
                        if let Some(last_pos) = state.last_cursor_position {
                            let dx = current_pos.x - last_pos.x;
                            let dy = current_pos.y - last_pos.y;
                            state.last_cursor_position = Some(current_pos);
                            return Some(iced::widget::shader::Action::publish(
                                Message::PanCamera(dx, dy),
                            ));
                        }
                        state.last_cursor_position = Some(current_pos);
                    }
                } else if let Some(current_pos) = local_cursor_position(cursor, _bounds) {
                    return Some(iced::widget::shader::Action::publish(Message::HoverAt(
                        current_pos.x,
                        current_pos.y,
                        _bounds.width,
                        _bounds.height,
                    )));
                } else {
                    return Some(iced::widget::shader::Action::publish(Message::ClearHover));
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
            brep: self.brep.clone(),
            mesh: self.mesh.clone(),
            camera: *self.camera,
            aspect: _bounds.width / _bounds.height,
            selected_faces: self.selected_faces.clone(),
            selected_edges: self.selected_edges.clone(),
        }
    }

    fn mouse_interaction(
        &self,
        state: &Self::State,
        _bounds: iced::Rectangle,
        _cursor: iced::mouse::Cursor,
    ) -> iced::mouse::Interaction {
        if state.is_rotating || state.is_panning {
            iced::mouse::Interaction::Grabbing
        } else if _cursor.is_over(_bounds) {
            iced::mouse::Interaction::Grab
        } else {
            iced::mouse::Interaction::default()
        }
    }
}

fn local_cursor_position(cursor: iced::mouse::Cursor, bounds: iced::Rectangle) -> Option<iced::Point> {
    let absolute = cursor.position()?;
    if absolute.x < bounds.x
        || absolute.y < bounds.y
        || absolute.x > bounds.x + bounds.width
        || absolute.y > bounds.y + bounds.height
    {
        return None;
    }
    Some(iced::Point::new(absolute.x - bounds.x, absolute.y - bounds.y))
}

#[derive(Debug, Clone)]
struct Primitive {
    brep: BRep,
    mesh: Mesh,
    camera: Camera,
    aspect: f32,
    selected_faces: Vec<usize>,
    selected_edges: Vec<usize>,
}

impl iced::widget::shader::Primitive for Primitive {
    type Pipeline = RCadPipeline;

    fn prepare(
        &self,
        pipeline: &mut Self::Pipeline,
        device: &iced::wgpu::Device,
        queue: &iced::wgpu::Queue,
        _bounds: &iced::Rectangle,
        viewport: &iced::advanced::graphics::Viewport,
    ) {
        let face_mesh = build_faces_highlight_mesh(&self.brep, &self.selected_faces);
        let edge_mesh = build_edges_highlight_mesh(&self.brep, &self.selected_edges);
        pipeline.renderer.upload_highlights(
            unsafe { std::mem::transmute(device) },
            face_mesh.as_ref(),
            edge_mesh.as_ref(),
        );

        let physical_size = viewport.physical_size();
        pipeline.renderer.prepare_scene_with_depth(
            unsafe { std::mem::transmute(device) },
            unsafe { std::mem::transmute(queue) },
            &self.mesh,
            &self.camera,
            self.aspect,
            (physical_size.width, physical_size.height),
        );
    }

    fn render(
        &self,
        pipeline: &Self::Pipeline,
        encoder: &mut iced::wgpu::CommandEncoder,
        target: &iced::wgpu::TextureView,
        clip_bounds: &iced::Rectangle<u32>,
    ) {
        pipeline.renderer.render_with_defaults(
            unsafe { std::mem::transmute(target) },
            unsafe { std::mem::transmute(encoder) },
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
