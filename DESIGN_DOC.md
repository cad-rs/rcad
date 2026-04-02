# RCAD Design and Development Document (DESIGN_DOC.md)

## 1. System Architecture Diagram
RCAD uses a layered architecture to separate low-level geometry calculations from high-level application logic.
```
[Applications]      creator-egui | creator-iced
--------------------------------------------
[Scene Commands]    rcad-scene (tool states, creation workflows)
--------------------------------------------
[Viewers]           rcad-render (wgpu)
--------------------------------------------
[Algorithms]        rcad-algorithms (BooleanOps, FilletOps, Sweeps)
--------------------------------------------
[Kernel Core]       rcad-kernel (B-Rep, Topology, Geometry)
--------------------------------------------
[Data Sources]      rcad-step (ISO 10303)
```

## 2. Core Kernel Design (rcad-kernel)

### 2.1 Geometry Primitives
- Implemented in `libs/rcad-kernel/src/geom.rs`.
- Uses `glam::DVec3` for double-precision geometry coordinates.
- Current analytic geometry coverage:
  - Curves: `Line`, `Circle`, `Ellipse`
  - Surfaces: `Plane`, `Cylinder`, `Sphere`, `Cone`, `Torus`
  - Primitive solids: `Box`, `Sphere`, `Cylinder`, `Cone`, `Torus`

### 2.2 Topological Structures
- Implemented in `libs/rcad-kernel/src/topology.rs`.
- `Vertex`: Represents a point in 3D space.
- `Edge`: Bounded portion of a curve (two vertices).
- `Wire`: Sequential collection of edges.
- `Face`: Bounded portion of a surface (one or more wires).
- `Shell`: Connected collection of faces.
- `Solid`: Bounded volume (one or more shells).

### 2.3 Topological Data Storage
- `BRep` stores topology arrays (`vertices`, `edges`, `solids`) plus geometric bindings (`geom: GeomStore`).
- `GeomStore` keeps curve/surface pools and mapping arrays from edges/faces to analytic geometry.

## 2.4 Modeling Entry Layer (rcad-modeling)
- Implemented in `libs/rcad-modeling/src/builder.rs`.
- Provides user-facing construction helpers for analytic geometry and primitive solids.
- API direction is intentionally aligned with OCCT constructor style:
  - Prefer direct public functions over fluent builder structs.
  - Keep validation at the modeling layer and return typed errors.
  - Return `Curve3`, `Surface3`, `PrimitiveSolid`, or `BRep` depending on the construction helper.

## 3. Visualization Pipeline (rcad-render)
- **Tessellation**: Converts B-Rep triangles to renderable mesh buffers.
- **Picking**:
  - Face picking by screen ray vs triangle intersection.
  - Edge picking by projected screen-space segment distance.
- **Selection State**:
  - `SelectionState` centralizes mode (`Face`/`Edge`), additive select, hover, and highlighted sets.
- **Wgpu Rendering**:
  - Main mesh pass + face highlight overlay + edge highlight overlay.
  - Shared renderer API used by both app frontends.
- **Camera Interaction**:
  - Orbit rotation, wheel zoom, and middle-mouse pan (`Camera::pan_pixels`).

## 3.1 Scene Command Layer (rcad-scene)
- Shared command state machine for creation tools (`SelectFace`, `SelectEdge`, `Box`, `Sphere`).
- Shared command lifecycle actions:
  - pointer click/move handling
  - preview generation
  - confirm/cancel/undo
- Shared BRep append utility used by creator apps after command confirmation.

## 4. STEP Importer/Exporter (rcad-step)
- **Parser**: Hand-written STEP Part 21 parser for core entities.
- **Mapping**: Converts common entities (point/direction/line/circle/surface + topology entities) into internal `BRep`.
- **Fallback behavior**: When shell/face topology is missing but points exist, importer falls back to a bbox solid for viewability.

## 5. Development Workflow
1. **Kernel updates** in `rcad-kernel` for type definitions and storage layout.
2. **Modeling API updates** in `rcad-modeling` for user-facing geometry construction.
3. **STEP mapping updates** in `rcad-step` with tests against sample assets.
4. **Rendering and interaction updates** in `rcad-render` first (API-first rule).
5. **Frontend wiring only** in `creator-egui` / `creator-iced`.
6. **Validation** with `cargo check` for both apps and target libs.

## 6. Project Directory Structure
```
rcad/
├── Cargo.toml          # Workspace root
├── libs/
│   ├── rcad-kernel/    # Primitives, Topology
│   ├── rcad-modeling/  # User-facing geometry construction helpers
│   ├── rcad-algorithms/# Boolean operations, Sweeps
│   ├── rcad-step/      # STEP Parser/Writer
│   ├── rcad-render/    # wgpu Rendering Engine
│   └── rcad-scene/     # Shared scene command workflows
├── apps/
│   ├── creator-egui/   # egui Modeling App
│   └── creator-iced/   # iced Modeling App
├── assets/             # Example STEP files, Shaders
└── scripts/            # Build/deploy scripts
```