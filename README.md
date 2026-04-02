# RCAD

A generic CAD engine written in pure Rust, targeting feature parity with Open CASCADE Technology (OCCT). Compiles to both native and WebAssembly.

## Current Status (April 2026)

- STEP import can load real sample models (for example `assets/hfss.step`) instead of a fixed placeholder box.
- `rcad-kernel` separates analytic geometry (`geom`) and connectivity topology (`topology`).
- `rcad-render` centralizes picking, selection state, and highlight rendering.
- `creator-egui` and `creator-iced` are interaction-aligned:
    - Left drag: rotate
    - Middle drag: pan
    - Mouse wheel: zoom
    - Click: select face/edge (with additive multi-select)

## Workspace Layout

```
rcad/
├── libs/
│   ├── rcad-kernel/      # B-Rep geometry and topology primitives
│   ├── rcad-algorithms/  # Boolean ops, sweeps, fillets
│   ├── rcad-step/        # STEP (ISO 10303) import / export
│   └── rcad-render/      # wgpu tessellation and mesh pipeline
└── apps/
    ├── creator-egui/     # egui-based modelling app
    └── creator-iced/     # iced-based modelling app
```

## Prerequisites

| Tool | Install |
|------|---------|
| Rust (stable, 1.85+) | `rustup update stable` |
| wasm32 target | `rustup target add wasm32-unknown-unknown` |
| Trunk | `cargo install trunk` |
| wasm-bindgen CLI (must match crate pin `=0.2.114`) | `cargo install wasm-bindgen-cli --version 0.2.114` |

## Run in the browser (Trunk / WASM)

Each app is a self-contained Trunk project. `cd` into the app directory and run `trunk serve`.

### creator-egui

```bash
cd apps/creator-egui
trunk serve --port 8080 --open
```

Then open <http://localhost:8080>.

- Left panel shows model info and selection state.
- Right panel renders the loaded model with face/edge highlight overlays.
- Controls: Left drag rotate, Middle drag pan, Wheel zoom, Click select.

### creator-iced

```bash
cd apps/creator-iced
trunk serve --port 8081 --open
```

Then open <http://localhost:8081>.

- Left panel shows model info and selection state.
- Right panel renders the loaded model with face/edge highlight overlays.
- Controls: Left drag rotate, Middle drag pan, Wheel zoom, Click select.

> **Note — release builds**: `trunk build --release` runs `wasm-opt` which is downloaded from GitHub. If the network blocks GitHub, build in dev mode (omit `--release`) or install `wasm-opt` manually and add it to `PATH`.

## Run natively

```bash
# egui app
cargo run -p creator-egui

# iced app
cargo run -p creator-iced

# load a specific STEP file
cargo run -p creator-egui -- assets/hfss.step
cargo run -p creator-iced -- assets/hfss.step
```

## Check / test all crates

```bash
cargo check --workspace
cargo test  --workspace
```

## WASM-specific notes

- `wasm-bindgen` is **pinned** to `=0.2.114` in the root `Cargo.toml` to match the installed CLI version. If you upgrade the CLI, update the pin to match.
- `creator-iced` uses **only the `tiny-skia` renderer** on WASM (`default-features = false, features = ["canvas", "tiny-skia"]`). Enabling `wgpu` alongside `tiny-skia` causes a canvas context conflict in browsers (WebGL context vs `CanvasRenderingContext2d` on the same canvas).
- `creator-egui` uses `eframe`'s `glow` backend (WebGL 2) on WASM.
