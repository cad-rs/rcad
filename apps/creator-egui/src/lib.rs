use eframe::{egui, egui_wgpu};
use rcad_kernel::BRep;
use rcad_render::{Camera, Mesh, Tessellator, WgpuRenderer};

const SAMPLE_STEP: &str = include_str!("../../../assets/box.step");

// ─── App ─────────────────────────────────────────────────────────────────────

pub struct RCadApp {
    brep: BRep,
    mesh: Mesh,
    camera: Camera,
    has_renderer: bool,
    auto_rotate: bool,
}

impl RCadApp {
    pub fn new(cc: &eframe::CreationContext<'_>, step_content: Option<String>) -> Self {
        let brep = if let Some(content) = step_content {
            rcad_step::StepReader::parse_string(&content)
        } else {
            rcad_step::StepReader::parse_string(SAMPLE_STEP)
        }
        .unwrap_or_else(|_| BRep::create_box(1.0, 1.0, 1.0));
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
            auto_rotate: true,
        }
    }
}

// ─── Wgpu Callback ───────────────────────────────────────────────────────────

struct RenderCallback {
    camera: Camera,
    aspect: f32,
    mesh: Mesh,
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

        if self.auto_rotate {
            self.camera.rot_y += ctx.input(|i| i.unstable_dt) * 0.6;
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
            ui.separator();
            ui.checkbox(&mut self.auto_rotate, "Auto-rotate");
            ui.label("Drag to rotate");

            if ui.button("Reset Camera").clicked() {
                self.camera = Camera::new();
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            let size = ui.available_size();
            let (rect, response) = ui.allocate_exact_size(size, egui::Sense::drag());

            if response.dragged() {
                let delta = response.drag_delta();
                self.camera.rot_y += delta.x * 0.008;
                self.camera.rot_x += delta.y * 0.008;
            }

            if self.has_renderer {
                let aspect = rect.width() / rect.height();
                let cb = RenderCallback {
                    camera: self.camera,
                    aspect,
                    mesh: self.mesh.clone(),
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
