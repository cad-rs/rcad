# RCAD

A parametric solid modeling kernel implemented in Rust, designed for CAD/CAE applications, electromagnetic simulation preprocessing, and engineering software.

## Features

- **Analytic-first geometry** - Exact analytic surfaces (planes, cylinders, spheres, cones, tori, B-splines) for modeling, triangulation only for rendering
- **B-Rep topology** - Full boundary representation with vertices, edges, wires, faces, shells, and solids
- **Boolean operations** - Union, intersection, difference with robust handling
- **STEP I/O** - ISO 10303 AP203/AP214 import and export
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
[Data Sources]      rcad-step (ISO 10303)
```

## Crates

| Crate | Description |
|-------|-------------|
| `rcad-kernel` | Geometry primitives, topology, B-Rep, transformations |
| `rcad-modeling` | Modeling API: sweeps, fillets, offsets, history tree |
| `rcad-algorithms` | Boolean operations, face imprint, cross-sections, HLR |
| `rcad-step` | STEP AP203/AP214 read/write, assembly I/O |
| `rcad-iges` | IGES format support |
| `rcad-render` | wgpu rendering, picking, HLR stroke lines |
| `rcad-scene` | Scene management, tool states |

## Building

```bash
# Native build
cargo build --release

# WebAssembly build
cargo build --target wasm32-unknown-unknown --release
```

## Applications

- `creator-egui` - Desktop application using egui
- `creator-iced` - Desktop application using Iced

## License

MIT License
