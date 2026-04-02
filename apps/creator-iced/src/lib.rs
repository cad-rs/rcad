use iced::widget::{button, checkbox, column, container, row, text};
use iced::{Element, Length, Task};
use rcad_kernel::BRep;
use rcad_render::{
    build_edges_highlight_mesh, build_faces_highlight_mesh, merge_meshes, Camera, Mesh,
    SelectionMode, SelectionState, Tessellator, WgpuRenderer,
};
use rcad_scene::{append_brep, CreationController, Tool};
use rcad_step::writer::{ExportSelection, StepWriter};

const SAMPLE_STEP: &str = include_str!("../../../assets/box.step");

pub struct RCadApp {
    brep: BRep,
    mesh: Mesh,
    camera: Camera,
    selection: SelectionState,
    creation: CreationController,
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
    SetTool(Tool),
    CancelCommand,
    ConfirmCommand,
    UndoLastStep,
    GrowSelectedFaces,
    GrowSelectedEdges,
    SelectFaceBoundaryEdges,
    SelectEdgeIncidentFaces,
    ClearSelection,
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
                rcad_modeling::make_box_brep(
                    glam::DVec3::ZERO,
                    glam::DVec3::X,
                    glam::DVec3::Y,
                    1.0,
                    1.0,
                    1.0,
                )
                .expect("default fallback box should be valid")
            }
        };
        let mesh = Tessellator::tessellate(&brep);

        (
            Self {
                brep,
                mesh,
                camera: Camera::new(),
                selection: SelectionState::default(),
                creation: CreationController::default(),
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
                self.handle_primary_click([x, y], [w.max(1.0), h.max(1.0)]);
            }
            Message::HoverAt(x, y, w, h) => {
                self.handle_pointer_move([x, y], [w.max(1.0), h.max(1.0)]);
            }
            Message::ClearHover => {
                self.creation
                    .clear_hover_if_selection_tool(&mut self.selection);
            }
            Message::SetSelectionMode(mode) => {
                self.creation.set_selection_mode(mode, &mut self.selection);
            }
            Message::SetAdditiveSelect(v) => {
                self.creation.set_additive_select(&mut self.selection, v);
            }
            Message::SetTool(tool) => {
                self.creation.set_tool(tool, &mut self.selection);
            }
            Message::CancelCommand => {
                self.creation.cancel_active_command();
            }
            Message::ConfirmCommand => {
                self.commit_active_command();
            }
            Message::UndoLastStep => {
                self.creation.undo_last_step();
            }
            Message::GrowSelectedFaces => {
                self.creation
                    .grow_selected_faces(&self.brep, &mut self.selection);
            }
            Message::GrowSelectedEdges => {
                self.creation
                    .grow_selected_edges(&self.brep, &mut self.selection);
            }
            Message::SelectFaceBoundaryEdges => {
                self.creation
                    .select_face_boundary_edges(&self.brep, &mut self.selection);
            }
            Message::SelectEdgeIncidentFaces => {
                self.creation
                    .select_edge_incident_faces(&self.brep, &mut self.selection);
            }
            Message::ClearSelection => {
                self.creation.clear_selection(&mut self.selection);
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

        let toolbar = container(
            row![
                text("Tools"),
                button("Select Face").on_press(Message::SetTool(Tool::SelectFace)),
                button("Select Edge").on_press(Message::SetTool(Tool::SelectEdge)),
                button("Box").on_press(Message::SetTool(Tool::Box)),
                button("Sphere").on_press(Message::SetTool(Tool::Sphere)),
                text(self.creation.command_hint()),
            ]
            .spacing(8),
        )
        .padding(8)
        .width(Length::Fill);

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
            text("Topology Nav"),
            row![
                button("Grow Faces").on_press(Message::GrowSelectedFaces),
                button("Grow Edges").on_press(Message::GrowSelectedEdges)
            ]
            .spacing(6),
            row![
                button("Face -> Edges").on_press(Message::SelectFaceBoundaryEdges),
                button("Edge -> Faces").on_press(Message::SelectEdgeIncidentFaces)
            ]
            .spacing(6),
            button("Clear Selection").on_press(Message::ClearSelection),
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
                self.selection
                    .hovered_face
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string())
            )),
            text(format!(
                "Hover Edge: {}",
                self.selection
                    .hovered_edge
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string())
            )),
            text(format!(
                "Hover Pos: {}",
                self.selection
                    .last_hover_pos
                    .map(|(x, y)| format!("{x:.1}, {y:.1}"))
                    .unwrap_or_else(|| "-".to_string())
            )),
            text(format!("Active Tool: {}", self.creation.tool_name())),
            text("Left click: select/create, Alt+Left drag: rotate"),
            text("Middle drag: pan, Wheel: zoom, Esc: cancel, Enter: finish"),
            button("Export STEP").on_press(Message::ExportStep),
            text(self.export_status.clone().unwrap_or_default()),
            button("Reset Camera").on_press(Message::ResetCamera),
        ]
        .spacing(8))
        .padding(12)
        .width(Length::Fixed(220.0))
        .height(Length::Fill);

        let viewport = container(
            iced::widget::shader(Scene {
                brep: &self.brep,
                mesh: &self.mesh,
                camera: &self.camera,
                selected_faces: self.selection.highlighted_faces(),
                selected_edges: self.selection.highlighted_edges(),
                preview_mesh: self
                    .creation
                    .preview_brep(self.camera.distance)
                    .map(|brep| Tessellator::tessellate(&brep)),
            })
            .width(Length::Fill)
            .height(Length::Fill),
        );

        column![toolbar, row![info, viewport].width(Length::Fill).height(Length::Fill)]
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }

    fn handle_primary_click(&mut self, cursor: [f32; 2], viewport: [f32; 2]) {
        if let Some(new_brep) = self.creation.handle_primary_click(
            &self.brep,
            &self.camera,
            &mut self.selection,
            cursor,
            viewport,
        ) {
            append_brep(&mut self.brep, new_brep);
            self.mesh = Tessellator::tessellate(&self.brep);
        }
    }

    fn handle_pointer_move(&mut self, cursor: [f32; 2], viewport: [f32; 2]) {
        self.creation.handle_pointer_move(
            &self.brep,
            &self.camera,
            &mut self.selection,
            cursor,
            viewport,
        );
    }

    fn commit_active_command(&mut self) {
        if let Some(new_brep) = self.creation.confirm_active_command(&self.camera) {
            append_brep(&mut self.brep, new_brep);
            self.mesh = Tessellator::tessellate(&self.brep);
        }
    }
}

struct Scene<'a> {
    brep: &'a BRep,
    mesh: &'a Mesh,
    camera: &'a Camera,
    selected_faces: Vec<usize>,
    selected_edges: Vec<usize>,
    preview_mesh: Option<Mesh>,
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
    primary_pressed: bool,
    alt_pressed: bool,
    is_rotating: bool,
    is_panning: bool,
    primary_drag_distance: f32,
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
        bounds: iced::Rectangle,
        cursor: iced::mouse::Cursor,
    ) -> Option<iced::widget::shader::Action<Message>> {
        match event {
            iced::Event::Keyboard(iced::keyboard::Event::ModifiersChanged(modifiers)) => {
                state.alt_pressed = modifiers.alt();
            }
            iced::Event::Keyboard(iced::keyboard::Event::KeyPressed { key, .. }) => match key {
                iced::keyboard::Key::Named(iced::keyboard::key::Named::Escape) => {
                    return Some(iced::widget::shader::Action::publish(Message::CancelCommand));
                }
                iced::keyboard::Key::Named(iced::keyboard::key::Named::Enter) => {
                    return Some(iced::widget::shader::Action::publish(Message::ConfirmCommand));
                }
                iced::keyboard::Key::Named(iced::keyboard::key::Named::Backspace) => {
                    return Some(iced::widget::shader::Action::publish(Message::UndoLastStep));
                }
                _ => {}
            },
            iced::Event::Mouse(iced::mouse::Event::ButtonPressed(iced::mouse::Button::Left)) => {
                if cursor.is_over(bounds) {
                    state.primary_pressed = true;
                    state.is_rotating = state.alt_pressed;
                    state.primary_drag_distance = 0.0;
                    state.last_cursor_position = local_cursor_position(cursor, bounds);
                }
            }
            iced::Event::Mouse(iced::mouse::Event::ButtonPressed(iced::mouse::Button::Middle)) => {
                if cursor.is_over(bounds) {
                    state.is_panning = true;
                    state.last_cursor_position = local_cursor_position(cursor, bounds);
                }
            }
            iced::Event::Mouse(iced::mouse::Event::ButtonReleased(iced::mouse::Button::Left)) => {
                if state.primary_pressed
                    && state.primary_drag_distance < 3.0
                    && let Some(pos) = local_cursor_position(cursor, bounds)
                {
                    state.primary_pressed = false;
                    state.is_rotating = false;
                    state.last_cursor_position = Some(pos);
                    return Some(iced::widget::shader::Action::publish(Message::SelectAt(
                        pos.x,
                        pos.y,
                        bounds.width,
                        bounds.height,
                    )));
                }
                state.primary_pressed = false;
                state.is_rotating = false;
            }
            iced::Event::Mouse(iced::mouse::Event::ButtonReleased(iced::mouse::Button::Middle)) => {
                state.is_panning = false;
            }
            iced::Event::Mouse(iced::mouse::Event::WheelScrolled { delta }) => {
                if cursor.is_over(bounds) {
                    let y = match delta {
                        iced::mouse::ScrollDelta::Lines { y, .. } => *y * 20.0,
                        iced::mouse::ScrollDelta::Pixels { y, .. } => *y,
                    };
                    return Some(iced::widget::shader::Action::publish(Message::ZoomCamera(y)));
                }
            }
            iced::Event::Mouse(iced::mouse::Event::CursorMoved { .. }) => {
                if state.is_rotating {
                    if let Some(current_pos) = local_cursor_position(cursor, bounds) {
                        if let Some(last_pos) = state.last_cursor_position {
                            let dx = current_pos.x - last_pos.x;
                            let dy = current_pos.y - last_pos.y;
                            state.primary_drag_distance += (dx * dx + dy * dy).sqrt();
                            state.last_cursor_position = Some(current_pos);
                            return Some(iced::widget::shader::Action::publish(
                                Message::RotateCamera(dx * 0.8, dy * 0.8),
                            ));
                        }
                        state.last_cursor_position = Some(current_pos);
                    }
                } else if state.primary_pressed {
                    if let Some(current_pos) = local_cursor_position(cursor, bounds) {
                        if let Some(last_pos) = state.last_cursor_position {
                            let dx = current_pos.x - last_pos.x;
                            let dy = current_pos.y - last_pos.y;
                            state.primary_drag_distance += (dx * dx + dy * dy).sqrt();
                        }
                        state.last_cursor_position = Some(current_pos);
                    }
                } else if state.is_panning {
                    if let Some(current_pos) = local_cursor_position(cursor, bounds) {
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
                } else if let Some(current_pos) = local_cursor_position(cursor, bounds) {
                    return Some(iced::widget::shader::Action::publish(Message::HoverAt(
                        current_pos.x,
                        current_pos.y,
                        bounds.width,
                        bounds.height,
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
        bounds: iced::Rectangle,
    ) -> Self::Primitive {
        Primitive {
            brep: self.brep.clone(),
            mesh: self.mesh.clone(),
            camera: *self.camera,
            aspect: bounds.width / bounds.height,
            selected_faces: self.selected_faces.clone(),
            selected_edges: self.selected_edges.clone(),
            preview_mesh: self.preview_mesh.clone(),
        }
    }

    fn mouse_interaction(
        &self,
        state: &Self::State,
        bounds: iced::Rectangle,
        cursor: iced::mouse::Cursor,
    ) -> iced::mouse::Interaction {
        if state.is_rotating || state.is_panning {
            iced::mouse::Interaction::Grabbing
        } else if cursor.is_over(bounds) {
            iced::mouse::Interaction::Crosshair
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
    preview_mesh: Option<Mesh>,
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
        let selected_face_mesh = build_faces_highlight_mesh(&self.brep, &self.selected_faces);
        let mut face_parts: Vec<&Mesh> = Vec::new();
        if let Some(mesh) = selected_face_mesh.as_ref() {
            face_parts.push(mesh);
        }
        if let Some(mesh) = self.preview_mesh.as_ref() {
            face_parts.push(mesh);
        }
        let face_mesh = merge_meshes(&face_parts);
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

pub fn run_native(step_content: Option<String>) -> iced::Result {
    iced::application(move || RCadApp::new(step_content.clone()), RCadApp::update, RCadApp::view)
        .title("RCAD Creator · iced")
        .window(iced::window::Settings {
            size: iced::Size::new(900.0, 600.0),
            ..Default::default()
        })
        .run()
}

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen(start)]
pub fn start() {
    console_error_panic_hook::set_once();
    run_native(None).expect("iced failed to start");
}
