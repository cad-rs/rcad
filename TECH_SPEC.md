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
- `libs/rcad-modeling`: OCCT-aligned modeling entry points for analytic geometry and primitive creation.
- `libs/rcad-algorithms`: Boolean operations (CSG), Filleting, and Sweeping.
- `libs/rcad-step`: STEP file parser and generator.
- `libs/rcad-render`: wgpu-based visualization pipeline.
- `libs/rcad-scene`: shared scene-level command logic (tool state, creation workflow, command preview).
- `apps/creator-egui`: Validation app using `egui` and `eframe`.
- `apps/creator-iced`: Validation app using `iced`.

## 4. Feature Alignment (Targeting OCCT Parity)
### 4.1 Geometry (Geom)
- Implemented: Line, Circle, Ellipse.
- Implemented: Plane, Cylinder, Sphere, Cone, Torus.
- Planned: B-Spline / Bezier curve and surface families.
### 4.2 Topology (TopoDS)
- Implemented: Vertex, Edge, Wire, Face, Shell, Solid.
- Planned: Compound and richer topology graph services.
### 4.3 Modeling Algorithms
- Boolean operations (Union, Intersection, Difference).
- Filleting and Chamfering.
- Sweeping (Extrude, Revolve, Pipe).

## 5. Rendering Pipeline
- **Tessellation**: Fast B-Rep to Mesh conversion for visualization.
- **Shaders**: WGSL with base material and explicit face/edge highlight overlays.
- **Interaction**:
	- Ray-casting based face picking.
	- Screen-space edge picking.
	- Shared `SelectionState` for mode, additive select, and hover state.
	- Camera orbit + zoom + middle-mouse pan.

## 6. WASM Integration
- All libraries must be `no_std` compatible where possible or at least `wasm32-unknown-unknown` compatible.
- Use `wasm-bindgen` for JS interop.
- Use `trunk` for building and serving the web apps.

## 7. Performance and Safety
- **Memory Management**: Rust's ownership model replaces OCCT's `Handle` (smart pointers) with `Arc` and `Weak` for reference counting in cyclic graphs.
- **Concurrency**: Parallelize long-running algorithms (e.g., Boolean operations) using `rayon` for native targets and a compatible strategy for WASM (Web Workers).
- **Accuracy**: Double precision (`f64`) for all geometric calculations; configurable tolerances for topological consistency.

## 8. Development Roadmap
- **Phase 1 (Completed)**: Kernel primitives, basic STEP reader, egui integration.
- **Phase 2 (Completed)**: iced integration, shared renderer pipeline, picking/highlight, dual-app interaction alignment.
- **Phase 3 (In Progress)**: richer topology/geometry entities, robust STEP coverage, modeling algorithms (booleans/fillets/sweeps).
- **Phase 4 (Planned)**: full STEP export and advanced CAD command system.

## 9. Modeling API Principles (Mandatory)

The following rules are mandatory for current and future modeling APIs. They are intended to keep RCAD close to OCCT-style construction patterns while remaining idiomatic in Rust.

### 9.1 Single Modeling Entry Layer
- Public geometry and primitive creation helpers MUST live in `libs/rcad-modeling`.
- `rcad-kernel` remains the storage/model layer for `Curve3`, `Surface3`, `PrimitiveSolid`, and `BRep`, but SHOULD NOT become the main user-facing construction API.

### 9.2 OCCT-Aligned Constructor Style
- User-facing modeling helpers SHOULD follow OCCT-style direct construction functions rather than fluent builder object chains.
- For Rust code in RCAD, that means public creation APIs SHOULD be exposed as plain functions in `libs/rcad-modeling/src/builder.rs`.
- Example direction: `line(origin, direction)`, `plane(origin, normal)`, `cylinder_brep(center, axis, ref_dir, radius, height, segments)`.

### 9.3 Functional API Over Builder Structs
- Do NOT introduce public fluent `*Builder` structs for standard analytic geometry and primitive creation unless there is a concrete need that cannot be expressed cleanly with functions.
- Validation SHOULD still happen inside `rcad-modeling`, with errors returned explicitly through `Result<..., BuildError>` or an equivalent typed error.

### 9.4 Separation of Responsibilities
- `rcad-modeling` defines how users create geometry.
- `rcad-kernel` defines what geometry and topology are.
- `rcad-kernel` primitive tessellation helpers such as `BRep::create_*` are implementation details and SHOULD NOT be used directly outside the kernel.
- `rcad-step` maps STEP entities to and from RCAD types.
- `rcad-scene` owns shared command state machines and creation workflow logic used by app crates.
- `rcad-render` and app crates MUST consume modeling results through public APIs rather than reimplementing geometry construction logic.

## 10. Rendering Coding Principles (Mandatory)

The following rules are mandatory for all current and future rendering work. They are designed to keep behavior consistent between [apps/creator-egui](apps/creator-egui) and [apps/creator-iced](apps/creator-iced).

### 10.1 Single Ownership of Rendering Logic
- All wgpu rendering logic MUST be implemented in [libs/rcad-render](libs/rcad-render).
- App crates MUST NOT implement or duplicate low-level rendering steps such as:
	- Pipeline creation
	- Bind group and uniform layout setup
	- Vertex/index draw calls
	- Depth attachment management
	- Default clear color selection

### 10.2 Allowed Responsibilities in App Layer
- App crates MAY handle only framework integration concerns:
	- UI layout and widgets
	- User input mapping (drag, zoom, toggle states)
	- Runtime wiring for egui/iced callback or shader program hooks
- App crates MUST call renderer APIs exposed by [libs/rcad-render/src/lib.rs](libs/rcad-render/src/lib.rs) for scene preparation and drawing.

### 10.3 API-First Evolution Rule
- When a new rendering feature is required (e.g., shadows, clipping planes, anti-aliasing policy, tone mapping), it MUST be added to [libs/rcad-render](libs/rcad-render) first.
- App crates MUST consume that feature only through public renderer APIs.
- Direct access to internal renderer fields from app crates is forbidden; renderer internals should remain encapsulated.

### 10.4 Cross-App Visual Consistency
- Default render behavior (camera update path, lighting model, depth-test policy, clear color, picking behavior) MUST be defined centrally in [libs/rcad-render](libs/rcad-render).
- [apps/creator-egui](apps/creator-egui) and [apps/creator-iced](apps/creator-iced) MUST use the same rcad-render defaults unless there is an explicitly documented exception.
- Any intentional visual difference between apps MUST be documented in this specification before merge.

### 10.5 Definition of Done for Rendering Changes
- A rendering-related change is complete only if all conditions are satisfied:
	- Code changes are centralized in [libs/rcad-render](libs/rcad-render) unless strictly UI integration.
	- Both app targets compile successfully:
		- `cargo check -p creator-egui`
		- `cargo check -p creator-iced`
	- If web path is affected, web build must pass for relevant app with trunk.
	- The change does not reintroduce duplicated rendering code in app crates.