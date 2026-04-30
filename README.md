# RCAD

A parametric solid modeling kernel implemented in Rust, designed for CAD/CAE applications, electromagnetic simulation preprocessing, and engineering software.

## Features

- **Analytic-first geometry** - Exact analytic surfaces (planes, cylinders, spheres, cones, tori, B-splines) for modeling, triangulation only for rendering
- **B-Rep topology** - Full boundary representation with vertices, edges, wires, faces, shells, and solids
- **Boolean operations** - Union, intersection, difference with robust handling
- **STEP I/O** - ISO 10303 AP203/AP214 import and export
- **Layout I/O** - GDS-II (`rcad-gds`) and OASIS (`rcad-oas`) for mask/IC layout exchange
- **IGES** - B-Rep read/write (`rcad-iges`)
- **Rendering** - GPU-accelerated visualization via wgpu (WebGPU API)
- **Cross-platform** - Native and WebAssembly targets

## Architecture

```
[Applications]      creator-egui | creator-iced
--------------------------------------------
[Scene Commands]    rcad-scene (tool states, workflows)
--------------------------------------------
[Viewers]           rcad-render (wgpu)
--------------------------------------------
[Algorithms]        rcad-algorithms (Boolean, Fillet, Sweeps)
--------------------------------------------
[Kernel Core]       rcad-kernel (B-Rep, Topology, Geometry)
--------------------------------------------
[Data Sources]      rcad-step (ISO 10303), rcad-gds / rcad-oas (layout), rcad-iges
```

## Crates

| Crate | Description |
|-------|-------------|
| `rcad-kernel` | Geometry primitives, topology, B-Rep, transformations |
| `rcad-modeling` | Modeling API: sweeps, fillets, offsets, history tree |
| `rcad-algorithms` | Boolean operations, face imprint, cross-sections, HLR |
| `rcad-step` | STEP AP203/AP214 read/write, assembly I/O |
| `rcad-gds` | GDS-II layout read/write |
| `rcad-oas` | OASIS layout read/write |
| `rcad-iges` | IGES B-Rep read/write |
| `rcad-render` | wgpu rendering, picking, HLR stroke lines |
| `rcad-scene` | Scene management, tool states |
| `rcad-py` | Python bindings (PyO3): B-rep primitives, booleans, STEP/IGES |

## Building

```bash
# Native build
cargo build --release

# WebAssembly build
cargo build --target wasm32-unknown-unknown --release
```

### Python (uv + maturin)

The `rcad-py` crate publishes the `rcad` package. From `libs/rcad-py`:

```bash
uv sync --group dev
uv run maturin develop   # editable install with native extension
uv run python -c "import rcad; print(rcad.BRep.sphere((0,0,0), 1.0).volume())"
```

Requires a Rust toolchain (`cargo`) on `PATH`. The extension uses the stable **abi3** binary interface for Python 3.9+.

`BRep` exposes primitives and I/O (`read_step`, `read_step_with_metadata` for geometry plus a JSON-derived metadata `dict`, `write_step`, IGES), booleans (`union`, `intersection`, `difference`), sweeps (`extrude`, `revolve`, `loft`, `sweep_pipe`, `pipe_sweep_wire`, `sweep_wire_linear`), offset (`offset_solid`), projection (`project_wire_to_face`), local ops (`fillet_edge`, `fillet_edges`, `chamfer_edge`, `repair`), rigid transforms (`translate`, `rotate_axis_angle`, `scale_uniform`), and mass-properties / counts (`volume`, `surface_area`, `signed_volume`, `centroid`, `bounding_box`, `inertia_tensor`, `face_count`, …).

## Applications

- `creator-egui` - Desktop application using egui
- `creator-iced` - Desktop application using Iced

## License

MIT License
