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
- Analytic geometry coverage (Phase A + B):
  - **Curves (`Curve3`)**: `Line3`, `Circle3`, `Ellipse3`, `BSplineCurve3` (de Boor evaluation)
  - **Surfaces (`Surface3`)**: `Plane`, `CylindricalSurface`, `SphericalSurface`, `ConicalSurface`, `ToroidalSurface`, `BSplineSurface` (tensor-product de Boor), `LinearExtrusionSurface` (Phase K), `RevolutionSurface` (Phase K)
  - **2D Curves (`Curve2d`)**: `Line2d`, `Circle2d`, `BSplineCurve2` (Phase I — de Boor in 2D, for PCurves on B-spline surfaces), `Ellipse2d` (Phase J — 2D ellipse in parameter space)
  - **Evaluation traits**: `CurveEval` (`point_at`, `tangent_at`, `default_domain`) and `SurfaceEval` (`point_at`, `normal_at`, `default_domain`) — implemented for all analytic types; `Curve2dEval` (`point_at`) for all 2D variants
  - Primitive solids: `Box`, `Sphere`, `Cylinder`, `Cone`, `Torus`

### 2.2 Topological Structures
- Implemented in `libs/rcad-kernel/src/topology.rs`.
- `Vertex`: Represents a point in 3D space.
- `Edge`: Bounded portion of a curve (two vertex indices). Parameter range `[t1, t2]` stored in `GeomStore.edge_curve_range`.
- `WireEdge { idx: usize, forward: bool }`: An oriented edge reference inside a Wire.
- `Wire`: Ordered sequence of `WireEdge` entries with explicit orientation.
- `Face`: Bounded portion of a surface (outer wire + inner wires). The `triangles` field is **rendering metadata only** and must not influence modeling logic.
- `Shell`: Connected collection of faces.
- `Solid`: Bounded volume (one or more shells).

### 2.3 Topological Data Storage
- `BRep` stores topology arrays (`vertices`, `edges`, `solids`) plus geometric bindings (`geom: GeomStore`).
- `GeomStore` keeps curve/surface pools and mapping arrays from edges/faces to analytic geometry:
  - `curves: Vec<Curve3>` — analytic 3D curves
  - `edge_curve: Vec<Option<usize>>` — curve index per edge
  - `edge_curve_range: Vec<Option<[f64; 2]>>` — parameter range `[t1, t2]` per edge
  - `edge_degenerated: Vec<bool>` — degenerate edge flag (e.g., sphere pole)
  - `surfaces: Vec<Surface3>` — analytic 3D surfaces
  - `face_surface: Vec<Option<usize>>` — surface index per face
  - `curve2ds: Vec<Curve2d>` — 2D curves in parameter space (`Line2d`, `Circle2d`, `BSplineCurve2`, `Ellipse2d`)
  - `curve2d_range: Vec<Option<[f64; 2]>>` — parameter trim range per PCurve (Phase J); `None` = natural domain; parallel to `curve2ds`
  - `face_surface_range: Vec<Option<[f64; 4]>>` — per-face surface domain override `[u1, u2, v1, v2]` (Phase K); `None` = use `SurfaceEval::default_domain()`; parallel to `face_surface`
  - `edge_pcurves: Vec<Vec<PCurve>>` — per-edge PCurve bindings
  - `vertex_tolerance: Vec<f64>` — per-vertex tolerance (Phase I); falls back to `CONFUSION = 1e-7`
  - `edge_tolerance: Vec<f64>` — per-edge tolerance (Phase I); populated from STEP `UNCERTAINTY_MEASURE_WITH_UNIT` (Phase J)
  - `face_tolerance: Vec<f64>` — per-face tolerance (Phase I); populated from STEP `UNCERTAINTY_MEASURE_WITH_UNIT` (Phase J)
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
- `GeomStore.curve2ds: Vec<Curve2d>` — pool of 2D analytic curves (`Line2d`, `Circle2d`, `BSplineCurve2`, `Ellipse2d`)
- `GeomStore.curve2d_range: Vec<Option<[f64; 2]>>` — per-PCurve parameter trim range (Phase J); parallel to `curve2ds`; `None` = natural domain; `Some([t1, t2])` when originating from a STEP `TRIMMED_CURVE`
- `GeomStore.edge_pcurves: Vec<Vec<PCurve>>` — per-edge list of `PCurve { surface_idx, curve2d_idx }`
- Seam edges on closed surfaces (sphere, cylinder, torus) have **two** PCurves — one for each boundary side

**STEP mapping:**
```
EDGE_CURVE → SURFACE_CURVE(#3d_curve, (#pcurve1, #pcurve2)) → EDGE_CURVE
PCURVE('', #surface, DEFINITIONAL_REPRESENTATION(...#2d_curve...))
```
- `Curve2d::Line` → `LINE` (2D) entity in STEP
- `Curve2d::Circle` → `CIRCLE` with `AXIS2_PLACEMENT_2D` entity
- `Curve2d::Ellipse` → `ELLIPSE` with `AXIS2_PLACEMENT_2D` entity (Phase J)
- `Curve2d::BSpline` → `B_SPLINE_CURVE_WITH_KNOTS` with 2D control points (Phase J)

PCurves are required for full OCCT/CAE interoperability. Edges without PCurves fall back to 3D-curve-only STEP export, which is valid but loses parametric surface information.

## 2.5 Modeling Entry Layer (rcad-modeling)
- Implemented in `libs/rcad-modeling/src/builder/`.
- Provides user-facing construction helpers for analytic geometry and primitive solids.
- API direction is intentionally aligned with OCCT constructor style:
  - Prefer direct public functions over fluent builder structs.
  - Keep validation at the modeling layer and return typed errors (`BuildError`).
  - Return `Curve3`, `Surface3`, `PrimitiveSolid`, or `BRep` depending on the construction helper.
- Free-form BRep construction (Phase B):
  - `make_edge(brep, curve, t1, t2, v0, v1)` — add an edge with explicit curve and parameter range
  - `make_wire(edges: Vec<WireEdge>)` — construct an oriented Wire
  - `make_face(brep, surface, outer, inner_wires)` — add a face with analytic surface
  - `make_solid(brep, shells)` — add a solid
- Sweep operations (Phase B):
  - `extrude(profile, direction, distance)` — linear prism
  - `revolve(profile, axis_origin, axis_dir, angle)` — solid of revolution
- Multi-profile operations (Phase E):
  - `loft(profiles: &[Vec<DVec3>])` — connect N cross-section polygons with ruled lateral faces and planar caps
  - `sweep_pipe(profile_2d: &[DVec2], spine: &[DVec3])` — sweep 2D profile along a 3D polyline spine using Frenet-like frames; delegates to `loft`
  - `sweep_pipe_variable(profiles: &[Vec<DVec2>], spine: &[DVec3])` — variable-section sweep: a different 2D profile at each spine station (Phase L); analogous to OCCT `BRepOffsetAPI_MakePipeShell`
- Edge modification operations (Phase F):
  - `chamfer_edge(brep, edge_idx, dist)` — flat bevel; replaces edge with planar quad + 2 closing triangles; returns new BRep
  - `fillet_edge(brep, edge_idx, radius)` — cylindrical blend; setback = `radius / tan(β/2)` from exterior dihedral angle; returns new BRep
  - `fillet_edges(brep, edges: &[(usize, f64)])` — batch fillet API: applies `fillet_edge` for each entry, sorted by index descending (Phase L); safe for non-adjacent edges

## 2.7 Topology Query Layer (rcad-kernel / topo_query.rs)
- Analogous to OCCT `TopExp_Explorer` and `TopExp::MapShapesAndAncestors`.
- All functions operate on `solids[0].shells[0]`; safe on empty BRep.
- `edge_adjacent_faces(brep, edge_idx) -> Vec<usize>` — faces sharing an edge
- `face_edges(brep, face_idx) -> Vec<usize>` — edges of a face's outer wire
- `vertex_adjacent_edges(brep, vertex_idx) -> Vec<usize>` — edges incident on a vertex
- `face_count / edge_count / vertex_count` — shape size queries

## 2.8 Curvature Analysis (rcad-kernel / curvature.rs)
- Analogous to OCCT `GeomLProp_SLProps`.
- `principal_curvatures(surface, u, v) -> (k1, k2)`:
  - **Plane**: (0, 0)
  - **Cylinder(r)**: (1/r, 0)
  - **Sphere(r)**: (1/r, 1/r)
  - **Cone(α, v)**: (sin(α)/r_at, 0) where r_at = v·sin(α)
  - **Torus(R, r, v)**: (1/r, cos(v)/(R+r·cos(v)))
  - **BSpline**: numerical finite-difference via fundamental forms (I and II)
- `gaussian_curvature(surface, u, v) -> f64` — K = k1·k2
- `mean_curvature(surface, u, v) -> f64` — H = (k1+k2)/2
- `Color { r, g, b }` — sRGB color with preset constants (RED, GREEN, BLUE, …)
- `FaceColor { face_index, color }` — per-face color override
- `StepColor { solid_color, face_colors }` — color assignments for a BRep; used by `StepWriter::write_string_colored`

## 2.9 Analysis and Algorithms (rcad-algorithms)
- **Boolean operations** (`builder`): union, intersection, difference on convex BReps
- **Shape validity** (`brep_check`): `check(brep) -> CheckResult` — reports degenerate/invalid topology
- **Global properties** (`rcad-kernel/properties`): `surface_area`, `volume`, `centroid`, `inertia_tensor`
- **Section** (`section`): `section_polylines(brep, plane)` — cross-section line set
- **HLR** (`hlr`): `hlr(brep, camera, samples) -> HlrResult` — hidden-line removal via ray-triangle occlusion; `hlr_to_svg(result, scale, margin)` — SVG rendering

## 2.10 Curve Arc Length (rcad-kernel / arc_length.rs)
- Analogous to OCCT `GCPnts_AbscissaPoint` / `CPnts_AbscissaPoint::Length`.
- `arc_length(curve: &Curve3, t1: f64, t2: f64) -> f64` — signed arc length over `[t1, t2]`
  - `Line3`: exact — `t2 − t1` (direction is always unit)
  - `Circle3`: exact — `r · (t2 − t1)`
  - `Ellipse3`, `BSplineCurve3`: 16-point Gauss-Legendre quadrature of `|dP/dt|` (finite-difference speed)
- Returns signed value; caller takes `.abs()` for unsigned length.

## 2.11 Moment of Inertia Tensor (rcad-kernel / properties.rs)
- Analogous to OCCT `BRepGProp_VolumeProperties`.
- `inertia_tensor(brep: &BRep) -> InertiaTensor` — symmetric 3×3 tensor about the world origin
- `InertiaTensor { ixx, iyy, izz, ixy, ixz, iyz }` with `to_matrix() -> [[f64;3];3]`
- Uses divergence-theorem tetrahedral integration, consistent with `volume` / `centroid`.
- Diagonal terms `Ixx = ∫(y²+z²)dV`, etc.; off-diagonal `Ixy = −∫xy dV`, etc.
- Assumes unit density; multiply by material density for physical inertia.


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
- **Mapping**: Converts common entities (point/direction/line/circle/ellipse/B-spline/surface + topology entities) into internal `BRep`.
- **B-Spline curve support** (Phase D): `B_SPLINE_CURVE_WITH_KNOTS` parsed into `Curve3::BSpline`; written with compressed knot vector (`multiplicities + values`).
- **B-Spline surface support** (Phase E): `B_SPLINE_SURFACE_WITH_KNOTS` parsed into `Surface3::BSpline`; the STEP `[v][u]` control grid is transposed to `BSplineSurface.control_points[u][v]`; UV-grid triangulation via `SurfaceEval::point_at` for rendering.
- **Color export** (Phase D): `StepWriter::write_string_colored(brep, &StepColor)` emits the full `COLOUR_RGB → STYLED_ITEM` chain per STEP AP214.
- **Assembly export** (Phase D): `write_assembly(name, &[AssemblyComponent])` produces a multi-BRep STEP file with `PRODUCT` / `NEXT_ASSEMBLY_USAGE_OCCURRENCE` hierarchy; each component can carry a translation and color.
- **Curve2d export** (Phase J): `Curve2d::Ellipse` → `ELLIPSE` + `AXIS2_PLACEMENT_2D`; `Curve2d::BSpline` → `B_SPLINE_CURVE_WITH_KNOTS` with 2D control points.
- **Tolerance import** (Phase J): `UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(val), ...)` → `GeomStore.{vertex,edge,face}_tolerance` filled with `val`; falls back to `CONFUSION = 1e-7` when absent.
- **BSpline surface export** (Phase K): `Surface3::BSpline` → `B_SPLINE_SURFACE_WITH_KNOTS` with full control-point grid and knot vectors; kernel [u][v] grid transposed to STEP [v][u] order (was falling back to PLANE).
- **Swept surface import** (Phase K): `SURFACE_OF_LINEAR_EXTRUSION` → `Surface3::LinearExtrusion`; `SURFACE_OF_REVOLUTION` → `Surface3::Revolution`. Profile curve resolved via existing `resolve_curve`; direction/axis resolved via `direction_from_ref` / `placement_from_ref`.
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