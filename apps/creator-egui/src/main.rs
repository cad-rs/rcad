fn main() {
    #[cfg(not(target_arch = "wasm32"))]
    {
        let step_content = std::env::args().nth(1).and_then(|path| {
            std::fs::read_to_string(path).ok()
        });
        creator_egui::run_native(step_content);
    }
}
