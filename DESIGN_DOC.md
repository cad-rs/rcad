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

### 2.0 Analytic-First Modeling Principle (MANDATORY)

RCAD is a CAD/CAE engine. Its internal model of geometry must be **exact and analytic at all times**. Triangulation is a rendering artifact, not a modeling artifact.

**Rules that apply to every layer of the codebase:**

1. **The authoritative shape is analytic.** Every `Face` in a `BRep` that has an analytic backing surface MUST have that surface stored in `GeomStore.surfaces` and indexed via `GeomStore.face_surface`. Every `Edge` with an analytic curve MUST have it in `GeomStore.curves` / `GeomStore.edge_curve`.

2. **Triangles are rendering-only.** `Face.triangles` exists exclusively to give the render pipeline pre-computed mesh data. It has NO modeling significance. No modeling, Boolean, STEP, or algorithm code may depend on `Face.triangles` being populated; it is always considered optional.

3. **Primitives are not triangle soups.** `BRep::create_sphere`, `create_cylinder`, `create_cone`, and `create_torus` MUST produce analytically correct BReps with proper edge/wire topology and populated `GeomStore` entries. Using `from_triangle_soup` for these shapes is a bug.

4. **Triangulation happens in `rcad-render` only.** The render pipeline tessellates analytic surfaces on demand. App or algorithm code MUST NOT call tessellation routines as part of modeling.

5. **STEP export is the correctness test.** If a shape exports from `rcad-step::StepWriter` as `ADVANCED_FACE` with its proper analytic surface type (SPHERICAL_SURFACE, CYLINDRICAL_SURFACE, etc.) rather than as triangle faces, the modeling layer is correct.

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
- `Face`: Bounded portion of a surface (one or more wires). The `triangles` field is **rendering metadata only** and must not influence modeling logic.
- `Shell`: Connected collection of faces.
- `Solid`: Bounded volume (one or more shells).

### 2.3 Topological Data Storage
- `BRep` stores topology arrays (`vertices`, `edges`, `solids`) plus geometric bindings (`geom: GeomStore`).
- `GeomStore` keeps curve/surface pools and mapping arrays from edges/faces to analytic geometry.
- **`GeomStore` is the source of truth for shape.** A `BRep` without populated `GeomStore` entries is incomplete and must not leave `rcad-modeling` in that state.

### 2.4 PCurve (Parameter-Space Curve)

A **PCurve** is the image of a 3D edge in the 2D parameter domain (u, v) of an adjacent surface. This concept mirrors OCCT's `BRep_CurveOnSurface` and STEP's `PCURVE` / `SURFACE_CURVE` entities.

```
Edge
 ├── 3D curve  (Curve3)   — position in world space
 └── PCurve(s) per adjacent face:
      └── Curve2d on Surface parameter domain (u, v)
```

**Storage:**
- `GeomStore.curve2ds: Vec<Curve2d>` — pool of 2D analytic curves (`Line2d`, `Circle2d`)
- `GeomStore.edge_pcurves: Vec<Vec<PCurve>>` — per-edge list of `PCurve { surface_idx, curve2d_idx }`
- Seam edges on closed surfaces (sphere, cylinder, torus) have **two** PCurves — one for each boundary side

**STEP mapping:**
```
EDGE_CURVE → SURFACE_CURVE(#3d_curve, (#pcurve1, #pcurve2)) → EDGE_CURVE
PCURVE('', #surface, DEFINITIONAL_REPRESENTATION(...#2d_curve...))
```

PCurves are required for full OCCT/CAE interoperability. Edges without PCurves fall back to 3D-curve-only STEP export, which is valid but loses parametric surface information.

## 2.4 Modeling Entry Layer (rcad-modeling)
- Implemented in `libs/rcad-modeling/src/builder.rs`.
- Provides user-facing construction helpers for analytic geometry and primitive solids.
- API direction is intentionally aligned with OCCT constructor style:
  - Prefer direct public functions over fluent builder structs.
  - Keep validation at the modeling layer and return typed errors.
  - Return `Curve3`, `Surface3`, `PrimitiveSolid`, or `BRep` depending on the construction helper.

## 3. Visualization Pipeline (rcad-render)
- **Tessellation**: Converts analytic B-Rep surfaces to renderable mesh buffers on demand. When `Face.triangles` is pre-populated it is used as a cache; when absent the render pipeline tessellates from the analytic surface. Tessellation MUST NOT be triggered by modeling or export code.
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