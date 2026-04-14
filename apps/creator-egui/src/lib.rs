use eframe::{egui, egui_wgpu};
use glam::DVec3;
use rcad_kernel::BRep;
use rcad_modeling::make_box_brep;
use rcad_render::{
    AxisGizmoHit, Camera, Mesh, SelectionMode, SelectionState, TessellationOptions,
    Tessellator, WgpuRenderer, axis_gizmo_hit_test, build_edges_highlight_mesh,
    build_faces_highlight_mesh, merge_meshes,
};
use rcad_scene::{CreationController, Tool, WorkPlane, append_brep};
use rcad_step::writer::{ExportSelection, StepWriter};

const SAMPLE_STEP: &str = include_str!("../../../assets/box.step");

// ─── App ─────────────────────────────────────────────────────────────────────

pub struct RCadApp {
    brep: BRep,
    mesh: Mesh,
    camera: Camera,
    has_renderer: bool,
    selection: SelectionState,
    creation: CreationController,
    export_status: Option<String>,
}

impl RCadApp {
    pub fn new(cc: &eframe::CreationContext<'_>, step_content: Option<String>) -> Self {
        let parse_result = if let Some(content) = step_content {
            rcad_step::StepReader::parse_string(&content)
        } else {
            rcad_step::StepReader::parse_string(SAMPLE_STEP)
        };

        let mut brep = match parse_result {
            Ok(brep) => {
                eprintln!(
                    "[rcad-step][egui] parsed STEP: vertices={}, edges={}, solids={}",
                    brep.vertices.len(),
                    brep.edges.len(),
                    brep.solids.len()
                );
                brep
            }
            Err(err) => {
                eprintln!("[rcad-step][egui] parse failed, fallback to box: {err}");
                make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
                    .expect("default fallback box should be valid")
            }
        };
        let tess_opts = TessellationOptions::default();
        let mesh = Tessellator::tessellate_with_options(&mut brep, &tess_opts);

        // Initialize wgpu renderer
        let mut has_renderer = false;
        if let Some(rs) = &cc.wgpu_render_state {
            // In egui 0.33 / wgpu 27, we use the device and queue from the render state
            let renderer = WgpuRenderer::new(&rs.device, rs.target_format);
            renderer.upload_mesh(&rs.device, &mesh);
            rs.renderer.write().callback_resources.insert(renderer);
            has_renderer = true;
        }

        Self {
            brep,
            mesh,
            camera: Camera::new(),
            has_renderer,
            selection: SelectionState::default(),
            creation: CreationController::default(),
            export_status: None,
        }
    }
}

// ─── Wgpu Callback ───────────────────────────────────────────────────────────

struct RenderCallback {
    brep: BRep,
    camera: Camera,
    aspect: f32,
    mesh: Mesh,
    viewport_origin_px: [u32; 2],
    viewport_size_px: [u32; 2],
    selected_faces: Vec<usize>,
    selected_edges: Vec<usize>,
    preview_mesh: Option<Mesh>,
}

impl egui_wgpu::CallbackTrait for RenderCallback {
    fn prepare(
        &self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        _screen_descriptor: &egui_wgpu::ScreenDescriptor,
        _egui_encoder: &mut wgpu::CommandEncoder,
        callback_resources: &mut egui_wgpu::CallbackResources,
    ) -> Vec<wgpu::CommandBuffer> {
        let Some(renderer) = callback_resources.get::<WgpuRenderer>() else {
            return Vec::new();
        };
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
        renderer.upload_highlights(device, face_mesh.as_ref(), edge_mesh.as_ref());
        renderer.prepare_scene(device, queue, &self.mesh, &self.camera, self.aspect);
        Vec::new()
    }

    fn paint(
        &self,
        _info: egui::PaintCallbackInfo,
        render_pass: &mut wgpu::RenderPass<'static>,
        callback_resources: &egui_wgpu::CallbackResources,
    ) {
        let Some(renderer) = callback_resources.get::<WgpuRenderer>() else {
            return;
        };
        renderer.draw_in_render_pass(render_pass, false);
        renderer.draw_axis_gizmo_in_render_pass(
            render_pass,
            self.viewport_origin_px,
            self.viewport_size_px,
        );
    }
}

// ─── eframe::App ─────────────────────────────────────────────────────────────

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::JsCast as _;

impl eframe::App for RCadApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Handle scroll zoom
        let scroll = ctx.input(|i| i.smooth_scroll_delta.y);
        if scroll != 0.0 {
            self.camera.distance -= scroll * 0.005 * self.camera.distance;
            self.camera.distance = self.camera.distance.clamp(1.0, 50.0);
            ctx.request_repaint();
        }

        if ctx.input(|i| i.key_pressed(egui::Key::Escape)) {
            self.creation.cancel_active_command();
        }
        if ctx.input(|i| i.key_pressed(egui::Key::Enter)) {
            self.commit_active_command();
        }
        if ctx.input(|i| i.key_pressed(egui::Key::Backspace)) {
            self.creation.undo_last_step();
        }

        egui::TopBottomPanel::top("toolbar").show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label("Tools");
                if ui
                    .selectable_label(
                        self.creation.active_tool() == Tool::SelectFace,
                        "Select Face",
                    )
                    .clicked()
                {
                    self.set_tool(Tool::SelectFace);
                }
                if ui
                    .selectable_label(
                        self.creation.active_tool() == Tool::SelectEdge,
                        "Select Edge",
                    )
                    .clicked()
                {
                    self.set_tool(Tool::SelectEdge);
                }
                if ui
                    .selectable_label(self.creation.active_tool() == Tool::Box, "Box")
                    .clicked()
                {
                    self.set_tool(Tool::Box);
                }
                if ui
                    .selectable_label(self.creation.active_tool() == Tool::Sphere, "Sphere")
                    .clicked()
                {
                    self.set_tool(Tool::Sphere);
                }
                if ui
                    .selectable_label(self.creation.active_tool() == Tool::Cylinder, "Cylinder")
                    .clicked()
                {
                    self.set_tool(Tool::Cylinder);
                }
                if ui
                    .selectable_label(self.creation.active_tool() == Tool::Cone, "Cone")
                    .clicked()
                {
                    self.set_tool(Tool::Cone);
                }
                if ui
                    .selectable_label(self.creation.active_tool() == Tool::Torus, "Torus")
                    .clicked()
                {
                    self.set_tool(Tool::Torus);
                }
                ui.separator();
                ui.label("Work Plane:");
                for plane in [WorkPlane::XY, WorkPlane::XZ, WorkPlane::YZ] {
                    if ui
                        .selectable_label(self.creation.work_plane() == plane, plane.label())
                        .clicked()
                    {
                        self.creation.set_work_plane(plane);
                    }
                }
                ui.separator();
                ui.label(self.creation.command_hint());
            });
        });

        egui::SidePanel::left("info")
            .min_width(180.0)
            .show(ctx, |ui| {
                ui.heading("RCAD  ·  egui");
                ui.separator();
                ui.label(format!("Vertices : {}", self.brep.vertices.len()));
                ui.label(format!("Edges    : {}", self.brep.edges.len()));
                let face_count: usize = self
                    .brep
                    .solids
                    .iter()
                    .flat_map(|s| &s.shells)
                    .map(|sh| sh.faces.len())
                    .sum();
                ui.label(format!("Faces    : {}", face_count));
                ui.label(format!("Triangles: {}", self.mesh.indices.len() / 3));
                ui.label(format!("Curves   : {}", self.brep.geom.curves.len()));
                ui.label(format!("Surfaces : {}", self.brep.geom.surfaces.len()));
                ui.separator();
                ui.label("Selection Mode");
                ui.horizontal(|ui| {
                    if ui.button("Select Face").clicked() {
                        self.creation
                            .set_selection_mode(SelectionMode::Face, &mut self.selection);
                    }
                    if ui.button("Select Edge").clicked() {
                        self.creation
                            .set_selection_mode(SelectionMode::Edge, &mut self.selection);
                    }
                });
                let mut additive_select = self.selection.additive_select;
                if ui
                    .checkbox(&mut additive_select, "Additive Select")
                    .changed()
                {
                    self.creation
                        .set_additive_select(&mut self.selection, additive_select);
                }
                ui.separator();
                ui.label("Topology Nav");
                ui.horizontal_wrapped(|ui| {
                    if ui.button("Grow Faces").clicked() {
                        self.creation
                            .grow_selected_faces(&self.brep, &mut self.selection);
                    }
                    if ui.button("Grow Edges").clicked() {
                        self.creation
                            .grow_selected_edges(&self.brep, &mut self.selection);
                    }
                });
                ui.horizontal_wrapped(|ui| {
                    if ui.button("Face -> Edges").clicked() {
                        self.creation
                            .select_face_boundary_edges(&self.brep, &mut self.selection);
                    }
                    if ui.button("Edge -> Faces").clicked() {
                        self.creation
                            .select_edge_incident_faces(&self.brep, &mut self.selection);
                    }
                });
                if ui.button("Clear Selection").clicked() {
                    self.creation.clear_selection(&mut self.selection);
                }
                ui.label(format!(
                    "Selected Faces: {}",
                    self.selection.selected_faces.len()
                ));
                ui.label(format!(
                    "Selected Edges: {}",
                    self.selection.selected_edges.len()
                ));
                ui.label(format!(
                    "Hover Face: {}",
                    self.selection
                        .hovered_face
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "-".to_string())
                ));
                ui.label(format!(
                    "Hover Edge: {}",
                    self.selection
                        .hovered_edge
                        .map(|v| v.to_string())
                        .unwrap_or_else(|| "-".to_string())
                ));
                ui.label(format!(
                    "Hover Pos: {}",
                    self.selection
                        .last_hover_pos
                        .map(|(x, y)| format!("{x:.1}, {y:.1}"))
                        .unwrap_or_else(|| "-".to_string())
                ));
                ui.label(format!("Active Tool: {}", self.creation.tool_name()));
                ui.label("Left click: select/create, Alt+Left drag: rotate");
                ui.label("Middle drag: pan, Wheel: zoom, Esc: cancel, Enter: finish");

                if ui.button("Export STEP").clicked() {
                    self.export_status =
                        Some(match export_step_file(&self.brep, &self.selection) {
                            Ok(path) => format!("Exported: {}", path),
                            Err(err) => format!("Export failed: {}", err),
                        });
                }
                if let Some(status) = &self.export_status {
                    ui.label(status);
                }

                if ui.button("Reset Camera").clicked() {
                    self.camera = Camera::new();
                }
            });

        // ── Topology Tree Panel ───────────────────────────────────────────────
        {
            // Pre-compute global face indices so collapsed headers don't desync counts.
            let mut fi_counter = 0usize;
            let solid_face_starts: Vec<Vec<usize>> = self
                .brep
                .solids
                .iter()
                .map(|solid| {
                    solid
                        .shells
                        .iter()
                        .map(|shell| {
                            let start = fi_counter;
                            fi_counter += shell.faces.len();
                            start
                        })
                        .collect()
                })
                .collect();

            let brep = &self.brep;
            let sel_faces: &[usize] = &self.selection.selected_faces;
            let sel_edges: &[usize] = &self.selection.selected_edges;
            let mut face_toggle: Option<usize> = None;
            let mut edge_toggle: Option<usize> = None;

            egui::SidePanel::right("topo_tree")
                .min_width(220.0)
                .resizable(true)
                .show(ctx, |ui| {
                    ui.heading("Topology Tree");
                    ui.separator();
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        for (si, (solid, shell_fi_starts)) in
                            brep.solids.iter().zip(solid_face_starts.iter()).enumerate()
                        {
                            egui::CollapsingHeader::new(format!("⬡ Body {si}"))
                                .id_salt(("body", si))
                                .default_open(si == 0)
                                .show(ui, |ui| {
                                    for (shi, (shell, &fi_start)) in
                                        solid.shells.iter().zip(shell_fi_starts.iter()).enumerate()
                                    {
                                        egui::CollapsingHeader::new(format!("Shell {shi}"))
                                            .id_salt(("shell", si, shi))
                                            .default_open(shi == 0)
                                            .show(ui, |ui| {
                                                for (lfi, face) in shell.faces.iter().enumerate() {
                                                    let fi = fi_start + lfi;
                                                    let is_fsel = sel_faces.contains(&fi);
                                                    let face_color = if is_fsel {
                                                        egui::Color32::GOLD
                                                    } else {
                                                        ui.visuals().text_color()
                                                    };
                                                    let face_text =
                                                        egui::RichText::new(format!("▲ Face {fi}"))
                                                            .color(face_color);
                                                    let fr = egui::CollapsingHeader::new(face_text)
                                                        .id_salt(("face", fi))
                                                        .show(ui, |ui| {
                                                            ui.label(
                                                                egui::RichText::new("Outer Wire")
                                                                    .small()
                                                                    .italics(),
                                                            );
                                                            for ei in &face.outer_wire.edges {
                                                                show_edge_item(
                                                                    ui,
                                                                    brep,
                                                                    ei.idx,
                                                                    sel_edges,
                                                                    &mut edge_toggle,
                                                                );
                                                            }
                                                            for (iwi, iw) in
                                                                face.inner_wires.iter().enumerate()
                                                            {
                                                                ui.label(
                                                                    egui::RichText::new(format!(
                                                                        "Inner Wire {iwi}"
                                                                    ))
                                                                    .small()
                                                                    .italics(),
                                                                );
                                                                for ei in &iw.edges {
                                                                    show_edge_item(
                                                                        ui,
                                                                        brep,
                                                                        ei.idx,
                                                                        sel_edges,
                                                                        &mut edge_toggle,
                                                                    );
                                                                }
                                                            }
                                                        });
                                                    if fr.header_response.clicked() {
                                                        face_toggle = Some(fi);
                                                    }
                                                }
                                            });
                                    }
                                });
                        }
                        // ── Flat vertex list ─────────────────────────────────
                        ui.separator();
                        egui::CollapsingHeader::new(format!("Vertices ({})", brep.vertices.len()))
                            .id_salt("all_vertices")
                            .default_open(false)
                            .show(ui, |ui| {
                                for (vi, v) in brep.vertices.iter().enumerate() {
                                    ui.label(
                                        egui::RichText::new(format!(
                                            "v{vi}: ({:.3}, {:.3}, {:.3})",
                                            v.point.x, v.point.y, v.point.z
                                        ))
                                        .monospace(),
                                    );
                                }
                            });
                    });
                });

            // Apply topology tree selection toggles
            if let Some(fi) = face_toggle {
                if self.selection.selected_faces.contains(&fi) {
                    self.selection.selected_faces.retain(|&f| f != fi);
                } else {
                    self.selection.selected_faces.push(fi);
                }
            }
            if let Some(ei) = edge_toggle {
                if self.selection.selected_edges.contains(&ei) {
                    self.selection.selected_edges.retain(|&e| e != ei);
                } else {
                    self.selection.selected_edges.push(ei);
                }
            }
        }

        egui::CentralPanel::default().show(ctx, |ui| {
            let size = ui.available_size();
            let (rect, response) = ui.allocate_exact_size(size, egui::Sense::click_and_drag());
            let pixels_per_point = ctx.pixels_per_point();
            let viewport_origin_px = [
                (rect.min.x * pixels_per_point).round().max(0.0) as u32,
                (rect.min.y * pixels_per_point).round().max(0.0) as u32,
            ];
            let viewport_size_px = [
                (rect.width() * pixels_per_point).round().max(1.0) as u32,
                (rect.height() * pixels_per_point).round().max(1.0) as u32,
            ];

            if response.dragged() {
                let delta = ui.input(|i| i.pointer.delta());
                let pan_with_middle =
                    ui.input(|i| i.pointer.button_down(egui::PointerButton::Middle));
                let rotate_with_alt = ui.input(|i| {
                    i.modifiers.alt && i.pointer.button_down(egui::PointerButton::Primary)
                });
                if pan_with_middle {
                    self.camera.pan_pixels(delta.x, delta.y);
                    ctx.request_repaint();
                } else if rotate_with_alt {
                    self.camera.rot_y += delta.x * 0.008;
                    self.camera.rot_x += delta.y * 0.008;
                    ctx.request_repaint();
                }
            }

            if response.clicked()
                && let Some(pointer) = response.interact_pointer_pos()
            {
                let local = pointer - rect.min;
                let alt_down = ui.input(|i| i.modifiers.alt);
                if !alt_down {
                    let pointer_px = [pointer.x * pixels_per_point, pointer.y * pixels_per_point];
                    if self.has_renderer
                        && let Some(hit) = axis_gizmo_hit_test(
                            &self.camera,
                            viewport_origin_px,
                            viewport_size_px,
                            pointer_px,
                        )
                    {
                        self.apply_axis_gizmo_hit(hit);
                        ctx.request_repaint();
                        return;
                    }
                    self.handle_primary_click([local.x, local.y], [rect.width(), rect.height()]);
                }
            }

            if let Some(pointer) = ui.input(|i| i.pointer.hover_pos()) {
                if rect.contains(pointer) {
                    let local = pointer - rect.min;
                    self.handle_pointer_move([local.x, local.y], [rect.width(), rect.height()]);
                    if self.creation.is_command_active() {
                        ctx.request_repaint();
                    }
                } else if matches!(
                    self.creation.active_tool(),
                    Tool::SelectFace | Tool::SelectEdge
                ) {
                    self.creation
                        .clear_hover_if_selection_tool(&mut self.selection);
                }
            } else if !response.hovered() {
                self.creation
                    .clear_hover_if_selection_tool(&mut self.selection);
            }

            if self.has_renderer {
                let aspect = rect.width() / rect.height();
                let preview_mesh = self
                    .creation
                    .preview_brep(self.camera.distance)
                    .map(|brep| Tessellator::tessellate(&brep));
                let cb = RenderCallback {
                    brep: self.brep.clone(),
                    camera: self.camera,
                    aspect,
                    mesh: self.mesh.clone(),
                    viewport_origin_px,
                    viewport_size_px,
                    selected_faces: self.selection.highlighted_faces(),
                    selected_edges: self.selection.highlighted_edges(),
                    preview_mesh,
                };

                ui.painter()
                    .add(egui_wgpu::Callback::new_paint_callback(rect, cb));
            } else {
                ui.centered_and_justified(|ui| {
                    ui.label("wgpu not available");
                });
            }
        });
    }
}

impl RCadApp {
    fn set_tool(&mut self, tool: Tool) {
        self.creation.set_tool(tool, &mut self.selection);
    }

    fn apply_axis_gizmo_hit(&mut self, hit: AxisGizmoHit) {
        match hit {
            AxisGizmoHit::X => self.camera.set_view_direction(glam::Vec3::X),
            AxisGizmoHit::Y => self.camera.set_view_direction(glam::Vec3::Y),
            AxisGizmoHit::Z => self.camera.set_view_direction(glam::Vec3::Z),
            AxisGizmoHit::Center => self.camera.set_isometric_view(),
        }
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
            let tess_opts = TessellationOptions::default();
            self.mesh = Tessellator::tessellate_with_options(&mut self.brep, &tess_opts);
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
            let tess_opts = TessellationOptions::default();
            self.mesh = Tessellator::tessellate_with_options(&mut self.brep, &tess_opts);
        }
    }
}

// ─── Topology tree helpers ────────────────────────────────────────────────────

/// Render a selectable edge row with a hover tooltip that shows vertex coordinates.
fn show_edge_item(
    ui: &mut egui::Ui,
    brep: &BRep,
    ei: usize,
    sel_edges: &[usize],
    edge_toggle: &mut Option<usize>,
) {
    let Some(edge) = brep.edges.get(ei) else {
        return;
    };
    let is_esel = sel_edges.contains(&ei);
    let color = if is_esel {
        egui::Color32::from_rgb(100, 180, 255)
    } else {
        ui.visuals().weak_text_color()
    };
    let label = egui::RichText::new(format!("─ Edge {ei}: v{} → v{}", edge.start, edge.end))
        .color(color)
        .monospace();
    let resp = ui.selectable_label(is_esel, label);
    if resp.clicked() {
        *edge_toggle = Some(ei);
    }
    resp.on_hover_ui(|ui| {
        if let Some(vs) = brep.vertices.get(edge.start) {
            ui.label(format!(
                "v{}: ({:.3}, {:.3}, {:.3})",
                edge.start, vs.point.x, vs.point.y, vs.point.z
            ));
        }
        if let Some(ve) = brep.vertices.get(edge.end) {
            ui.label(format!(
                "v{}: ({:.3}, {:.3}, {:.3})",
                edge.end, ve.point.x, ve.point.y, ve.point.z
            ));
        }
    });
}

// ─── Native entry ────────────────────────────────────────────────────────────

#[cfg(not(target_arch = "wasm32"))]
pub fn run_native(step_content: Option<String>) {
    let opts = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("RCAD Creator · egui")
            .with_inner_size([900.0, 600.0]),
        renderer: eframe::Renderer::Wgpu,
        ..Default::default()
    };
    eframe::run_native(
        "RCAD Creator (egui)",
        opts,
        Box::new(move |cc| Ok(Box::new(RCadApp::new(cc, step_content)))),
    )
    .expect("eframe failed");
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

// ─── WASM entry ──────────────────────────────────────────────────────────────

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen(start)]
pub fn start() {
    console_error_panic_hook::set_once();

    let web_options = eframe::WebOptions::default();

    wasm_bindgen_futures::spawn_local(async move {
        let document = web_sys::window().unwrap().document().unwrap();
        let canvas = document
            .get_element_by_id("main_canvas")
            .unwrap()
            .dyn_into::<web_sys::HtmlCanvasElement>()
            .unwrap();

        eframe::WebRunner::new()
            .start(
                canvas,
                web_options,
                Box::new(|cc| Ok(Box::new(RCadApp::new(cc, None)))),
            )
            .await
            .expect("eframe WebRunner failed");
    });
}
