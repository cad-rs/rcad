fn main() -> iced::Result {
    #[cfg(not(target_arch = "wasm32"))]
    return creator_iced::run_native();

    #[cfg(target_arch = "wasm32")]
    Ok(())
}
