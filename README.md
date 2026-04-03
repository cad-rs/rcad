# RCAD

A generic CAD engine written in pure Rust, targeting feature parity with Open CASCADE Technology (OCCT). Compiles to both native and WebAssembly.

## Current Status (April 2026)

Twelve development phases (A–L) of the OCCT parity roadmap are complete.

**Phase A — Geometry/topology foundations**
- `CurveEval` / `SurfaceEval` traits with `point_at`, `tangent_at`, `normal_at` for all analytic types.
- Edge parameter ranges `[t1, t2]` in `GeomStore.edge_curve_range`.
- `WireEdge { idx, forward }` orientation flags on every Wire edge.
- `BSplineCurve3` / `BSplineSurface` data types with de Boor evaluation.
- `BRep::bounding_box()` for axis-aligned bounds.

**Phase B — Modeling API**
- `make_edge` / `make_wire` / `make_face` / `make_solid` (analogous to `BRepBuilderAPI_Make*`).
- `extrude(profile, direction, distance)` — linear prism.
- `revolve(profile, axis, angle)` — solid of revolution.

**Phase C — Algorithms and analysis**
- `surface_area`, `volume`, `centroid` global properties.
- `BRepCheck` shape validity checker.
- `section(brep, plane)` → cross-section polylines.

**Phase D — Data exchange and visualization**
- STEP colored export: per-face and solid-level (`COLOUR_RGB` → `STYLED_ITEM` chain).
- STEP assembly: multi-BRep with `PRODUCT` / `NEXT_ASSEMBLY_USAGE_OCCURRENCE` hierarchy.
- `B_SPLINE_CURVE_WITH_KNOTS` read/write in STEP.
- HLR (Hidden-Line Removal): ray-triangle occlusion testing, SVG output via `hlr_to_svg`.

**Phase E — Loft, pipe sweep, B-Spline surface STEP read**
- `loft(profiles: &[Vec<DVec3>])` — connect N cross-section polygons into a closed solid.
- `sweep_pipe(profile_2d: &[DVec2], spine: &[DVec3])` — sweep 2D profile along 3D polyline via Frenet-like frames; delegates to `loft`.
- `B_SPLINE_SURFACE_WITH_KNOTS` STEP read: parses 2D control-point grid + knot vectors into `Surface3::BSpline`; UV-grid triangulation for rendering.

**Phase F — Chamfer and fillet**
- `chamfer_edge(brep, edge_idx, dist)` — flat bevel on a convex BRep edge; replaces the edge with a planar quad face and two triangular closing faces.
- `fillet_edge(brep, edge_idx, radius)` — cylindrical blend on a convex BRep edge; setback computed from the exterior dihedral angle (`radius / tan(β/2)`).
- Both operations rebuild the full BRep (non-destructive), limited to manifold edges shared by exactly two planar faces.

**Phase G — Topology query API and curvature analysis**
- `edge_adjacent_faces(brep, edge_idx)` / `face_edges(brep, face_idx)` / `vertex_adjacent_edges(brep, vertex_idx)` — public topology traversal (analogous to OCCT `TopExp_Explorer`).
- `face_count` / `edge_count` / `vertex_count` — shape size queries.
- `principal_curvatures(surface, u, v)` → `(k1, k2)` — analytic for Plane/Cylinder/Sphere/Cone/Torus; numerical finite-difference for BSpline.
- `gaussian_curvature` / `mean_curvature` — derived from principal curvatures (analogous to OCCT `GeomLProp_SLProps`).

**Phase H — Arc length and moment of inertia tensor**
- `arc_length(curve, t1, t2)` — exact for `Line3` (`|t2−t1|`) and `Circle3` (`r·|t2−t1|`); 16-point Gauss-Legendre quadrature for `Ellipse3` and `BSplineCurve3`. Returns signed value; `.abs()` for unsigned length. Analogous to OCCT `GCPnts_AbscissaPoint`.
- `inertia_tensor(brep)` → `InertiaTensor { ixx, iyy, izz, ixy, ixz, iyz }` — symmetric 3×3 moment of inertia about the world origin, computed via divergence-theorem tetrahedral integration over BRep triangles. Analogous to OCCT `BRepGProp_VolumeProperties`.

**Phase I — 2D B-Spline PCurve and per-entity tolerance system**
- `BSplineCurve2` — non-uniform rational B-spline in 2D parameter space; added as `Curve2d::BSpline` variant. Evaluated via de Boor's algorithm (2D analog of `de_boor`). Analogous to OCCT `Geom2d_BSplineCurve`.
- `Curve2dEval` dispatch updated for all three variants (`Line2d`, `Circle2d`, `BSplineCurve2`).
- Per-entity tolerance: `GeomStore.vertex_tolerance` / `edge_tolerance` / `face_tolerance` — parallel `Vec<f64>` arrays. Query helpers: `vertex_tolerance(brep, idx)`, `edge_tolerance`, `face_tolerance`, `model_tolerance`. Global constants: `CONFUSION = 1e-7`, `ANGULAR = 1e-12`, `APPROXIMATION = 1e-4`. Analogous to OCCT `Precision` class + `BRep_Tool::Tolerance`.

**Phase J — Ellipse2d PCurve, curve2d_range, STEP Curve2d I/O, STEP tolerance import**
- `Ellipse2d` — 2D ellipse parametric curve in parameter space; added as `Curve2d::Ellipse` variant. Parametric form `center + major_dir·a·cos(t) + minor_dir·b·sin(t)`, domain `[0, 2π]`. Analogous to OCCT `Geom2d_Ellipse`.
- `GeomStore.curve2d_range: Vec<Option<[f64; 2]>>` — per-PCurve parameter trim range, parallel to `curve2ds`. Stores `[t1, t2]` from STEP `TRIMMED_CURVE`; `None` = use natural domain. Analogous to `edge_curve_range` for 3D curves.
- STEP Curve2d export: `Curve2d::Ellipse` → `ELLIPSE` entity; `Curve2d::BSpline` → `B_SPLINE_CURVE_WITH_KNOTS` with 2D control points. Both use existing writer helpers.
- STEP tolerance import: `UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(val), ...)` → fills `GeomStore.vertex_tolerance`, `edge_tolerance`, `face_tolerance` with the file-specified value; falls back to `CONFUSION` when absent.

**Phase K — Swept surfaces, face domain, BSpline surface STEP export**
- `LinearExtrusionSurface` — `S(u,v) = profile.point_at(u) + v·direction`; normal = `tangent(u) × direction`. Analogous to OCCT `Geom_SurfaceOfLinearExtrusion`. STEP import: `SURFACE_OF_LINEAR_EXTRUSION` → `Surface3::LinearExtrusion`.
- `RevolutionSurface` — `S(u,v) = rotate(profile.point_at(v), axis_origin, axis_dir, angle=u)`; u ∈ [0, 2π], v from profile. Normal via finite-difference. Analogous to OCCT `Geom_SurfaceOfRevolution`.
- `GeomStore.face_surface_range: Vec<Option<[f64; 4]>>` — per-face surface parameter domain override `[u1, u2, v1, v2]`, parallel to `face_surface`. `face_domain(brep, idx)` returns the override when set, else `SurfaceEval::default_domain()`. Analogous to OCCT `BRep_Face::UVBounds()`.
- `BSplineSurface` STEP export: `write_surface` now emits `B_SPLINE_SURFACE_WITH_KNOTS` with full control-point grid and knot vectors (was falling back to PLANE). Control points transposed from kernel [u][v] to STEP [v][u] order.

**Phase L — Variable-section sweep and multi-edge fillet API**
- `sweep_pipe_variable(profiles: &[Vec<DVec2>], spine: &[DVec3])` — variable-section pipe sweep: a different 2D profile at each spine station, transformed via the same Frenet-like frame as `sweep_pipe`; delegates to `loft`. Analogous to OCCT `BRepOffsetAPI_MakePipeShell` with multiple sections.
- `fillet_edges(brep, edges: &[(usize, f64)])` — fillet multiple edges in a single call; sorts by index descending before applying so earlier fillets don't shift later indices (safe for non-adjacent edges). Analogous to adding multiple edges to `BRepFilletAPI_MakeFillet` before `Build()`.

**Phase M — All remaining P2 items**
- **SameParameter / SameRange edge flags**: `GeomStore.edge_same_parameter: Vec<bool>` + `edge_same_range: Vec<bool>`; STEP reader extracts the `same_parameter` field from `SURFACE_CURVE`; helper functions `edge_same_parameter(brep, idx)` / `edge_same_range(brep, idx)` default to `true` for RCAD-generated primitives.
- **Bezier curves and surfaces**: `BezierCurve3 / BezierSurface / BezierCurve2` added to `Curve3` / `Surface3` / `Curve2d` enums; evaluation via de Casteljau algorithm in homogeneous coordinates (supports rational weights). Analogous to OCCT `Geom_BezierCurve` / `Geom_BezierSurface`.
- **Offset curve and surface**: `OffsetCurve3 { basis, offset_distance, offset_dir }` — lateral offset `P + d·(tangent × dir).normalize()`; `OffsetSurface { basis, offset_distance }` — normal offset `P + d·normal`. Added as `Curve3::Offset` / `Surface3::Offset`. Analogous to OCCT `Geom_OffsetCurve` / `Geom_OffsetSurface`.
- **Boolean operation history**: `BooleanHistory { face_origins: Vec<FaceOrigin> }` maps each result face to `FaceOrigin::FromA(idx)`, `FromB(idx)`, or `Generated`. Convenience functions `union_with_history / intersection_with_history / difference_with_history` return `(BRep, BooleanHistory)`. Analogous to OCCT `BRepAlgoAPI_BuilderShape::Modified/Generated/Deleted`.
- **Corner blending**: `corner_blend(brep, vertex_idx, radius)` — blends a 3-valence convex corner by setting back each incident edge by `radius` and inserting a planar triangular closing patch. Eliminates gaps left at corners after `fillet_edges`. Analogous to OCCT `BRepFilletAPI_MakeFillet` corner resolution.

**Phase N — Curve fitting, point projection, analytic boolean intersections**
- **B-spline curve fitting** (`rcad_kernel::interpolate_points / approximate_points`): exact interpolation via collocation matrix + Gaussian elimination with partial pivoting; least-squares approximation via normal equations `(AᵀA)x = Aᵀb` with pinned endpoints. Chord-length parameterization; clamped cubic knot vectors. Analogous to OCCT `GeomAPI_Interpolate` / `GeomAPI_PointsToBSpline`.
- **Closest-point projection** (`rcad_kernel::closest_point_on_curve / closest_point_on_surface`): analytic closed-form projection for Plane, Sphere, Cylinder, Cone, Torus; Newton-Raphson refinement for all curves and parametric surfaces. Handles infinite-domain curves (Line). Analogous to OCCT `GeomAPI_ProjectPointOnCurve` / `GeomAPI_ProjectPointOnSurf`.
- **Analytic Plane×Sphere and Plane×Cylinder intersections** in the boolean PaveFiller: FF pass now dispatches these surface-type pairs to `inttools::plane_sphere` / `inttools::plane_cylinder` before falling back to marching; improved surface sampling for Plane and Cone geometries.

**Core infrastructure (ongoing)**
- STEP import: LINE / CIRCLE / ELLIPSE / B-SPLINE curves; PLANE / CYL / SPHERE / CONE / TORUS surfaces; PCurve / SURFACE_CURVE chains; GEOMETRIC_CURVE_SET.
- `rcad-kernel` separates analytic geometry (`geom`) and connectivity topology (`topology`).
- `rcad-render` centralizes picking, selection state, and highlight rendering.
- `rcad-scene` centralizes creation command state machines (Box/Sphere flow, preview, confirm/cancel/undo).
- `creator-egui` and `creator-iced` are interaction-aligned:
    - Left drag: rotate
    - Middle drag: pan
    - Mouse wheel: zoom
    - Click: select face/edge (with additive multi-select)

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
