use eframe::{egui, egui_wgpu};
use rcad_kernel::BRep;
use rcad_render::{
    build_edges_highlight_mesh, build_faces_highlight_mesh, Camera, Mesh,
    SelectionMode, SelectionState, Tessellator, WgpuRenderer, DEFAULT_EDGE_PICK_RADIUS_PX,
};

const SAMPLE_STEP: &str = include_str!("../../../assets/box.step");

// ─── App ─────────────────────────────────────────────────────────────────────

pub struct RCadApp {
    brep: BRep,
    mesh: Mesh,
    camera: Camera,
    has_renderer: bool,
    selection: SelectionState,
}

impl RCadApp {
    pub fn new(cc: &eframe::CreationContext<'_>, step_content: Option<String>) -> Self {
        let parse_result = if let Some(content) = step_content {
            rcad_step::StepReader::parse_string(&content)
        } else {
            rcad_step::StepReader::parse_string(SAMPLE_STEP)
        };

        let brep = match parse_result {
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
                BRep::create_box(1.0, 1.0, 1.0)
            }
        };
        let mesh = Tessellator::tessellate(&brep);

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
        }
    }
}

// ─── Wgpu Callback ───────────────────────────────────────────────────────────

struct RenderCallback {
    brep: BRep,
    camera: Camera,
    aspect: f32,
    mesh: Mesh,
    selected_faces: Vec<usize>,
    selected_edges: Vec<usize>,
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
        let face_mesh = build_faces_highlight_mesh(&self.brep, &self.selected_faces);
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

        egui::SidePanel::left("info").min_width(180.0).show(ctx, |ui| {
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
                if ui
                    .button("Select Face")
                    .clicked()
                {
                    self.selection.set_mode(SelectionMode::Face);
                }
                if ui
                    .button("Select Edge")
                    .clicked()
                {
                    self.selection.set_mode(SelectionMode::Edge);
                }
            });
            ui.checkbox(&mut self.selection.additive_select, "Additive Select");
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
            ui.label("Click to select, toggle Additive Select for multi-select");
            ui.label("Left drag: rotate, Middle drag: pan");

            if ui.button("Reset Camera").clicked() {
                self.camera = Camera::new();
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            let size = ui.available_size();
            let (rect, response) = ui.allocate_exact_size(size, egui::Sense::drag());

            if response.dragged() {
                let delta = ui.input(|i| i.pointer.delta());
                let pan_with_middle = ui.input(|i| i.pointer.button_down(egui::PointerButton::Middle));
                if pan_with_middle {
                    self.camera.pan_pixels(delta.x, delta.y);
                } else {
                    self.camera.rot_y += delta.x * 0.008;
                    self.camera.rot_x += delta.y * 0.008;
                }
            }

            if response.clicked()
                && let Some(pointer) = response.interact_pointer_pos()
            {
                let local = pointer - rect.min;
                let viewport = [rect.width(), rect.height()];
                let cursor = [local.x, local.y];
                let aspect = rect.width() / rect.height().max(1.0);
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

            if response.hovered()
                && !response.dragged()
                && let Some(pointer) = response.interact_pointer_pos()
            {
                let local = pointer - rect.min;
                let viewport = [rect.width(), rect.height()];
                let cursor = [local.x, local.y];
                let aspect = rect.width() / rect.height().max(1.0);
                self.selection
                    .hover_at(
                        &self.brep,
                        &self.camera,
                        aspect,
                        viewport,
                        cursor,
                        DEFAULT_EDGE_PICK_RADIUS_PX,
                    );
            } else if !response.hovered() {
                self.selection.clear_hover();
            }

            if self.has_renderer {
                let aspect = rect.width() / rect.height();
                let cb = RenderCallback {
                    brep: self.brep.clone(),
                    camera: self.camera,
                    aspect,
                    mesh: self.mesh.clone(),
                    selected_faces: self.selection.highlighted_faces(),
                    selected_edges: self.selection.highlighted_edges(),
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
