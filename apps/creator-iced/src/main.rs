fn main() -> iced::Result {
    #[cfg(not(target_arch = "wasm32"))]
    {
        let step_content = std::env::args().nth(1).and_then(|path| {
            std::fs::read_to_string(path).ok()
        });
        return creator_iced::run_native(step_content);
    }

    #[cfg(target_arch = "wasm32")]
    Ok(())
}
