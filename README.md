# RCAD

A generic CAD engine written in pure Rust, targeting feature parity with Open CASCADE Technology (OCCT). Compiles to both native and WebAssembly.

## Current Status (April 2026)

**Geometry kernel (`rcad-kernel`)**
- 3D curves: `Line3`, `Circle3`, `Ellipse3`, `Hyperbola3`, `Parabola3`, `BSplineCurve3`, `BezierCurve3`, `OffsetCurve3`
- 3D surfaces: `Plane`, `CylindricalSurface`, `SphericalSurface`, `ConicalSurface`, `ToroidalSurface`, `BSplineSurface`, `BezierSurface`, `LinearExtrusionSurface`, `RevolutionSurface`, `OffsetSurface`, `TrimmedSurface`
- 2D PCurves: `Line2d`, `Circle2d`, `Ellipse2d`, `BSplineCurve2`, `BezierCurve2`
- `CurveEval` / `SurfaceEval` traits; `BRep::bounding_box`; `BRep::apply_transform` / `transformed`
- Per-entity tolerance (`vertex/edge/face_tolerance`); `edge_same_parameter` / `edge_same_range` flags; `edge_curve_range`; `face_surface_range`; `curve2d_range`
- Curvature: `principal_curvatures`, `gaussian_curvature`, `mean_curvature`; arc length (Gauss-Legendre); inertia tensor
- Point projection: `closest_point_on_curve` / `closest_point_on_surface` (analytic + Newton-Raphson)
- Curve-curve extrema: `extrema_curve_curve`
- B-spline fitting: `interpolate_points`, `approximate_points`
- NURBS interop: `curve_to_bspline` / `surface_to_bspline` (exact for Line/Circle/Ellipse/Plane/Cylinder/Sphere/Bezier; sampling fallback for others); rational tensor-product `BSplineSurface::point_at`
- Curve/surface trim & extend: `trim_curve` (Boehm insertion), `extend_curve_to_point`, `extend_curve_by_length`, `trim_surface`, `extend_bspline_surface`

**Modeling API (`rcad-modeling`)**
- Primitives: `Box`, `Sphere`, `Cylinder`, `Cone`, `Torus`
- Sweeps: `extrude`, `revolve`, `sweep_pipe`, `sweep_pipe_variable`, `loft`
- Blending: `fillet_edge`, `chamfer_edge`, `fillet_edges`, `corner_blend`
- Shell sewing: `sew_shells`

**Algorithms (`rcad-algorithms`)**
- Boolean ops: `boolean_op(Union|Intersection|Difference)` + `*_with_history`; full support for planar solids; curved solids (Cylinder/Sphere) partially supported
- B-Rep repair: `merge_close_vertices`, `remove_degenerate_faces`, `recompute_face_normals`, `fix_wire_orientation`, `repair`
- Face imprinting: `imprint_brep`; gap/overlap detection: `detect_gaps_overlaps`
- Section: `section_polylines`, `section_curves` (analytic Circle/Ellipse/Line for Plane/Sphere/Cylinder/Cone faces)
- Shape distance: `min_distance`; topology query: `edge_adjacent_faces`, `face_edges`, `vertex_adjacent_edges`
- Surface-surface intersection: `intersect_surfaces` (analytic pairs + marching fallback)
- HLR: `hlr`, `hlr_to_svg`

**Data exchange (`rcad-step`)**
- STEP AP203/AP214 read/write for all geometry and topology types above
- Color import (`parse_string_with_color`) and export (`write_string_colored`)

**Rendering (`rcad-render`)**
- wgpu-based renderer with Blinn-Phong lighting (configurable light direction, headlight mode)
- Display modes: SolidWithEdges, Solid, Wireframe, Transparent
- Coordinate axes visualization (RGB arrows with cone heads)
- Background grid (XZ plane, major/minor lines)
- Face/edge picking (ray-cast) and selection highlighting
- Per-object color (`set_model_color`), screenshot export (`screenshot_to_file`)

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
