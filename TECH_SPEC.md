# RCAD Technical Specification (TECH_SPEC.md)

## 1. Project Overview
RCAD is a generic, high-performance CAD engine written entirely in Rust. It aims to provide a modern alternative to Open CASCADE Technology (OCCT) with a focus on memory safety, concurrency, and WebAssembly compatibility.

## 2. Core Technical Requirements
- **Language**: Rust (Edition 2024 preferred).
- **Architecture**: Cargo Workspace for modularity.
- **Target Platforms**: Native (Windows/Linux/macOS) and WebAssembly (via `wasm-bindgen` and `trunk`).
- **Rendering**: `wgpu` (WebGPU/Vulkan/Metal/DX12 abstraction).
- **Data Interchange**: Support for STEP (ISO 10303) import/export.
- **Kernel Architecture**: B-Rep (Boundary Representation).

## 3. Workspace Structure
- `libs/rcad-kernel`: Core geometry and topology engine.
- `libs/rcad-algorithms`: Boolean operations (CSG), Filleting, and Sweeping.
- `libs/rcad-step`: STEP file parser and generator.
- `libs/rcad-render`: wgpu-based visualization pipeline.
- `apps/creator-egui`: Validation app using `egui` and `eframe`.
- `apps/creator-iced`: Validation app using `iced`.

## 4. Feature Alignment (Targeting OCCT Parity)
### 4.1 Geometry (Geom)
- Primitive curves (Line, Circle, Ellipse, B-Spline, Bezier).
- Primitive surfaces (Plane, Cylinder, Sphere, Torus, B-Spline surfaces).
### 4.2 Topology (TopoDS)
- Vertex, Edge, Wire, Face, Shell, Solid, Compound.
- Connectivity graph management.
### 4.3 Modeling Algorithms
- Boolean operations (Union, Intersection, Difference).
- Filleting and Chamfering.
- Sweeping (Extrude, Revolve, Pipe).

## 5. Rendering Pipeline
- **Tessellation**: Fast B-Rep to Mesh conversion for visualization.
- **Shaders**: WGSL for PBR (Physically Based Rendering) and edge highlighting.
- **Interaction**: Ray-casting based picking and manipulation.

## 6. WASM Integration
- All libraries must be `no_std` compatible where possible or at least `wasm32-unknown-unknown` compatible.
- Use `wasm-bindgen` for JS interop.
- Use `trunk` for building and serving the web apps.

## 7. Performance and Safety
- **Memory Management**: Rust's ownership model replaces OCCT's `Handle` (smart pointers) with `Arc` and `Weak` for reference counting in cyclic graphs.
- **Concurrency**: Parallelize long-running algorithms (e.g., Boolean operations) using `rayon` for native targets and a compatible strategy for WASM (Web Workers).
- **Accuracy**: Double precision (`f64`) for all geometric calculations; configurable tolerances for topological consistency.

## 8. Development Roadmap
- **Phase 1 (MVP)**: Kernel primitives, Basic STEP reader, egui integration.
- **Phase 2 (Intermediate)**: Topology graph management, Tessellation for wgpu, basic Boolean operations.
- **Phase 3 (Advanced)**: Full STEP export, iced integration, complex modeling algorithms (fillets/sweeps).