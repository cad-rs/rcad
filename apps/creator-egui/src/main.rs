fn main() {
    #[cfg(not(target_arch = "wasm32"))]
    creator_egui::run_native();
}
