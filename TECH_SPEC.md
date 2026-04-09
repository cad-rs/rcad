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
- Implemented: Line, Circle, Ellipse (3D).
- Implemented: Plane, Cylinder, Sphere, Cone, Torus.
- Implemented: Line2d, Circle2d (2D parameter-space curves for PCurves).
- Planned: B-Spline / Bezier curve and surface families (Curve3 and Curve2d variants).
### 4.2 Topology (TopoDS)
- Implemented: Vertex, Edge, Wire, Face, Shell, Solid.
- Implemented: PCurve (`GeomStore.edge_pcurves`) — parameter-space curve binding per edge per adjacent face surface. Analogous to OCCT `BRep_CurveOnSurface`.
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
- **Error Handling**: All public APIs in `rcad-algorithms`, `rcad-modeling`, and `rcad-step` MUST return `Result<T, E>` with typed error enums. `unwrap`/`panic` in production code paths is forbidden; use `let-else`, indexed access with prior length checks, or `unwrap_or` with documented fallback semantics.
  - `BooleanError`: `EmptyInput | MissingGeometry | DegenerateResult | NumericalFailure | EmptyCollection`
  - `BuildError`: `NonFiniteValue | NonPositiveValue | ZeroVector | ParallelVectors | DegenerateGeometry | InvalidIndex`
  - `StepError`: `Io | InvalidFormat | MissingEntity | EmptyResult`

## 8. Development Roadmap
- **Phase 1 (Completed)**: Kernel primitives, basic STEP reader, egui integration.
- **Phase 2 (Completed)**: iced integration, shared renderer pipeline, picking/highlight, dual-app interaction alignment.
- **Phase 3 (In Progress)**: richer topology/geometry entities, robust STEP coverage, modeling algorithms (booleans/fillets/sweeps).
  - Error handling hardening: `unwrap`/`panic` removed from boolean op critical paths; structured `BooleanError` variants.
  - Integration tests: `rcad-algorithms/tests/`, `rcad-kernel/tests/`, `rcad-modeling/tests/`, `rcad-scene/tests/`, `rcad-step/tests/`, `rcad-render/tests/`.
  - Benchmark coverage: fillet, loft, sweep, box-cylinder booleans added to `rcad-algorithms/benches/`.
- **Phase 4 (Planned)**: full STEP export and advanced CAD command system.

## 9. Modeling API Principles (Mandatory)

The following rules are mandatory for current and future modeling APIs. They are intended to keep RCAD close to OCCT-style construction patterns while remaining idiomatic in Rust.

### 9.0 Analytic-First Invariant (HIGHEST PRIORITY)

RCAD is a CAD/CAE engine. The internal geometric model must be **exact and analytic** at every layer. The following invariants are non-negotiable:

**Invariant A — Triangles are rendering metadata, not geometry.**
`Face.triangles` is an optional rendering cache. It carries no geometric meaning. Modeling, Boolean, STEP, and algorithm code MUST treat it as absent/irrelevant. A `Face` with empty `triangles` is a fully valid face.

**Invariant B — Every analytic face must be backed by a `GeomStore` surface.**
After any modeling operation that produces or modifies a `BRep`, every face that has an analytic surface (Plane, Cylinder, Sphere, Cone, Torus, or future B-Spline) MUST have a corresponding entry in `GeomStore.face_surface` pointing to the correct `Surface3` value. A face with `geom.face_surface[face_idx] == None` is incomplete.

**Invariant C — Every analytic edge must be backed by a `GeomStore` curve.**
Likewise every edge with a known analytic curve (Line, Circle, Ellipse, or future B-Spline) MUST have a corresponding entry in `GeomStore.edge_curve`. An edge with no curve entry may be tolerated only for degenerate or seam edges that genuinely have no analytic form.

**Invariant D — Primitives are never triangle soups.**
`BRep::create_sphere`, `create_cylinder`, `create_cone`, and `create_torus` MUST build proper analytic BRep topology: named vertices at feature points, edges with analytic curves, faces with edge loops and `GeomStore` surface bindings. Using `from_triangle_soup` for these shapes is a bug and must be fixed before any shape reaches the export or algorithm layers.

**Invariant E — Triangulation lives exclusively in `rcad-render`.**
The render pipeline is the only place allowed to produce triangle meshes from analytic surfaces. The `rcad-algorithms`, `rcad-modeling`, `rcad-step`, and `rcad-scene` crates MUST NOT call tessellation routines.

**Correctness test:** If a primitive solid exports from `rcad-step::StepWriter` as `ADVANCED_FACE` with the correct analytic surface type (SPHERICAL_SURFACE, CYLINDRICAL_SURFACE, CONICAL_SURFACE, TOROIDAL_SURFACE, PLANE), the modeling layer is correct for that shape. Triangle-face fallback in STEP output is a red flag that Invariant B or D is violated.


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