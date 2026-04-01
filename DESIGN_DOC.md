# RCAD Design and Development Document (DESIGN_DOC.md)

## 1. System Architecture Diagram
RCAD uses a layered architecture to separate low-level geometry calculations from high-level application logic.
```
[Applications]      creator-egui | creator-iced
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
- `Point3D`, `Vector3D`, `Direction3D` based on `glam` or `nalgebra` for performance.
- `Curve` and `Surface` traits for unified evaluation, derivatives, and point inversion.

### 2.2 Topological Structures
- `Vertex`: Represents a point in 3D space.
- `Edge`: Bounded portion of a curve (two vertices).
- `Wire`: Sequential collection of edges.
- `Face`: Bounded portion of a surface (one or more wires).
- `Shell`: Connected collection of faces.
- `Solid`: Bounded volume (one or more shells).

### 2.3 Topological Data Storage
- Use `SlotMap` or similar for stable indices in a topological graph.
- `BRep` structure holds all topological and geometric data.

## 3. Visualization Pipeline (rcad-render)
- **Tessellation**: Converts `Face` surfaces to `Mesh` (Triangle soup) using adaptive sampling.
- **Wgpu Integration**: 
  - `RenderContext` initializes adapter, device, and queue.
  - `ShapePipeline` handles drawing solid bodies, wireframes, and isolated vertices.
  - `Uniforms` for camera/transformation data.

## 4. STEP Importer/Exporter (rcad-step)
- **Parser**: Hand-written or PEG-based for STEP physical file (Part 21).
- **Mapping**: Converts STEP entities (`CARTESIAN_POINT`, `B_SPLINE_CURVE_WITH_KNOTS`, `ADVANCED_FACE`) into internal `rcad-kernel` structures.
- **Accuracy**: Maintains exact numeric precision from the source file.

## 5. Development Workflow
1. **Core First**: Implement `rcad-kernel` primitives (`Point3D`, `Vector3D`, `Curve`, `Surface`).
2. **STEP Integration**: Basic STEP reader to load static geometries.
3. **Visualization**: Build `rcad-render` to show primitives using `wgpu`.
4. **Topology**: Implement `Vertex`, `Edge`, `Wire`, `Face` and their connectivity.
5. **Validation Apps**: Integrate `egui` and `iced` via `trunk` for web/native targets.

## 6. Project Directory Structure
```
rcad/
├── Cargo.toml          # Workspace root
├── libs/
│   ├── rcad-kernel/    # Primitives, Topology
│   ├── rcad-algorithms/# Boolean operations, Sweeps
│   ├── rcad-step/      # STEP Parser/Writer
│   └── rcad-render/    # wgpu Rendering Engine
├── apps/
│   ├── creator-egui/   # egui Modeling App
│   └── creator-iced/   # iced Modeling App
├── assets/             # Example STEP files, Shaders
└── scripts/            # Build/deploy scripts
```