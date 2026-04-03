# RCAD

A generic CAD engine written in pure Rust, targeting feature parity with Open CASCADE Technology (OCCT). Compiles to both native and WebAssembly.

## Current Status (April 2026)

Four development phases (A–D) of the OCCT parity roadmap are complete.

**Phase A — Geometry/topology foundations**
- `CurveEval` / `SurfaceEval` traits with `point_at`, `tangent_at`, `normal_at` for all analytic types.
- Edge parameter ranges `[t1, t2]` in `GeomStore.edge_curve_range`.
- `WireEdge { idx, forward }` orientation flags on every Wire edge.
- `BSplineCurve3` / `BSplineSurface` data types with de Boor evaluation.
- `BRep::bounding_box()` for axis-aligned bounds.

**Phase B — Modeling API**
- `make_edge` / `make_wire` / `make_face` / `make_solid` (analogous to `BRepBuilderAPI_Make*`).
- `extrude(profile, direction, distance)` — linear prism.
- `revolve(profile, axis, angle)` — solid of revolution.

**Phase C — Algorithms and analysis**
- `surface_area`, `volume`, `centroid` global properties.
- `BRepCheck` shape validity checker.
- `section(brep, plane)` → cross-section polylines.

**Phase D — Data exchange and visualization**
- STEP colored export: per-face and solid-level (`COLOUR_RGB` → `STYLED_ITEM` chain).
- STEP assembly: multi-BRep with `PRODUCT` / `NEXT_ASSEMBLY_USAGE_OCCURRENCE` hierarchy.
- `B_SPLINE_CURVE_WITH_KNOTS` read/write in STEP.
- HLR (Hidden-Line Removal): ray-triangle occlusion testing, SVG output via `hlr_to_svg`.

**Core infrastructure (ongoing)**
- STEP import: LINE / CIRCLE / ELLIPSE / B-SPLINE curves; PLANE / CYL / SPHERE / CONE / TORUS surfaces; PCurve / SURFACE_CURVE chains; GEOMETRIC_CURVE_SET.
- `rcad-kernel` separates analytic geometry (`geom`) and connectivity topology (`topology`).
- `rcad-render` centralizes picking, selection state, and highlight rendering.
- `rcad-scene` centralizes creation command state machines (Box/Sphere flow, preview, confirm/cancel/undo).
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
│   ├── rcad-modeling/    # Builder-style analytic geometry and primitive creation
│   ├── rcad-algorithms/  # Boolean ops, sweeps, fillets
│   ├── rcad-step/        # STEP (ISO 10303) import / export
│   ├── rcad-render/      # wgpu tessellation and mesh pipeline
│   └── rcad-scene/       # scene/command interaction logic shared by apps
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
