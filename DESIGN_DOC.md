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
  - **Curves (`Curve3`)**: `Line3`, `Circle3`, `Ellipse3`, `BSplineCurve3` (de Boor evaluation), `BezierCurve3` (Phase M — de Casteljau), `OffsetCurve3` (Phase M — lateral offset)
  - **Surfaces (`Surface3`)**: `Plane`, `CylindricalSurface`, `SphericalSurface`, `ConicalSurface`, `ToroidalSurface`, `BSplineSurface` (tensor-product de Boor), `LinearExtrusionSurface` (Phase K), `RevolutionSurface` (Phase K), `BezierSurface` (Phase M — de Casteljau), `OffsetSurface` (Phase M — normal offset)
  - **2D Curves (`Curve2d`)**: `Line2d`, `Circle2d`, `BSplineCurve2` (Phase I — de Boor in 2D, for PCurves on B-spline surfaces), `Ellipse2d` (Phase J — 2D ellipse in parameter space), `BezierCurve2` (Phase M — de Casteljau 2D)
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
  - `edge_same_parameter: Vec<bool>` — per-edge SameParameter flag (Phase M); default `true` for RCAD-generated primitives; STEP reader populates from `SURFACE_CURVE` 4th field
  - `edge_same_range: Vec<bool>` — per-edge SameRange flag (Phase M); default `true`
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
  - `corner_blend(brep, vertex_idx, radius)` — blend a 3-valence convex corner (Phase M): sets back each incident edge by `radius` and inserts a planar triangular closing patch; eliminates gaps at corners after `fillet_edges`

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
- **Boolean history** (`history`): `BooleanHistory { face_origins: Vec<FaceOrigin> }` — maps each result face to `FaceOrigin::FromA(idx)`, `FromB(idx)`, or `Generated`; returned by `union/intersection/difference_with_history()` (Phase M); analogous to OCCT `BRepAlgoAPI_BuilderShape`
- **Shape validity** (`brep_check`): `check(brep) -> CheckResult` — reports degenerate/invalid topology
- **Global properties** (`rcad-kernel/properties`): `surface_area`, `volume`, `centroid`, `inertia_tensor`
- **Section** (`section`): `section_polylines(brep, plane)` — cross-section line set
- **HLR** (`hlr`): `hlr(brep, camera, samples) -> HlrResult` — hidden-line removal via ray-triangle occlusion; `hlr_to_svg(result, scale, margin)` — SVG rendering

## 2.10 Curve Arc Length (rcad-kernel / arc_length.rs)
- Analogous to OCCT `GCPnts_AbscissaPoint` / `CPnts_AbscissaPoint::Length`.
- `arc_length(curve: &Curve3, t1: f64, t2: f64) -> f64` — signed arc length over `[t1, t2]`
  - `Line3`: exact — `t2 − t1` (direction is always unit)
  - `Circle3`: exact — `r · (t2 − t1)`
  - `Ellipse3`, `BSplineCurve3`, `BezierCurve3`, `OffsetCurve3`: 16-point Gauss-Legendre quadrature of `|dP/dt|` (finite-difference speed)
- Returns signed value; caller takes `.abs()` for unsigned length.

## 2.11 Moment of Inertia Tensor (rcad-kernel / properties.rs)
- Analogous to OCCT `BRepGProp_VolumeProperties`.
- `inertia_tensor(brep: &BRep) -> InertiaTensor` — symmetric 3×3 tensor about the world origin
- `InertiaTensor { ixx, iyy, izz, ixy, ixz, iyz }` with `to_matrix() -> [[f64;3];3]`
- Uses divergence-theorem tetrahedral integration, consistent with `volume` / `centroid`.
- Diagonal terms `Ixx = ∫(y²+z²)dV`, etc.; off-diagonal `Ixy = −∫xy dV`, etc.
- Assumes unit density; multiply by material density for physical inertia.

## 2.12 Curve Fitting (rcad-kernel / fit.rs)
- Analogous to OCCT `GeomAPI_Interpolate` and `GeomAPI_PointsToBSpline`.
- `interpolate_points(pts: &[DVec3]) -> Result<BSplineCurve3, FitError>` — exact interpolation through all points
  - Chord-length parameterization: `t[i] = Σ|chord_i| / total`, normalized to [0, 1]
  - Clamped cubic knot vector: interior knots via Piegl & Tiller §9.3 averaging formula
  - Collocation matrix solved by Gaussian elimination with partial pivoting
  - Degree: `min(3, n-1)` so 2 points → linear, 3 → quadratic, ≥4 → cubic
- `approximate_points(pts: &[DVec3], n_ctrl: usize) -> Result<BSplineCurve3, FitError>` — least-squares B-spline with `n_ctrl` control points
  - Normal equations `(AᵀA)x = Aᵀb` solved per coordinate component
  - Endpoints pinned to first/last data points
  - Falls back to `interpolate_points` when `n_ctrl >= pts.len()`
- `FitError::TooFewPoints` / `FitError::DegeneratePoints`

## 2.13 Closest-Point Projection (rcad-kernel / projection.rs)
- Analogous to OCCT `GeomAPI_ProjectPointOnCurve` and `GeomAPI_ProjectPointOnSurf`.
- `closest_point_on_curve(curve, query, n_samples) -> CurveProjection { point, param, distance }`
  - Analytic for Line (infinite domain: centers sampling around closest point analytically)
  - Newton-Raphson with finite-difference tangent + second-order curvature correction
- `closest_point_on_surface(surface, query, n_samples) -> SurfaceProjection { point, params, distance }`
  - **Analytic**: Plane (dot product), Sphere (radial normalize), Cylinder (collapse axis + normalize), Cone (minimize 1D along generator), Torus (major-ring then tube projection)
  - **Numerical fallback**: uniform grid sampling + Gauss-Newton Newton refinement (2×2 system) for BSpline, Bezier, Offset, LinearExtrusion, Revolution surfaces

## 2.14 Shape Distance (rcad-kernel / distance.rs)
- Analogous to OCCT `BRepExtrema_DistShapeShape`.
- `min_distance(a: &BRep, b: &BRep) -> ShapeDistance { distance, point_on_a, point_on_b }`
  - Samples each face: 4×4 (u,v) grid + wire vertex positions
  - Projects each sample onto all analytic surfaces of the other BRep via `closest_point_on_surface`
  - Symmetric A→B and B→A passes; returns global minimum
- `point_to_shape_distance(query: DVec3, brep: &BRep) -> ShapeDistance`
  - Projects query onto every analytic face surface; returns closest result

## 2.15 Shell Sewing (rcad-modeling / sewing.rs)
- Analogous to OCCT `BRepOffsetAPI_Sewing`.
- `sew_shells(breps: &[BRep], tolerance: f64) -> SewingResult { brep, stitched_pairs, free_edges }`
  - **Step 1**: Concatenates all vertices/edges/faces from every input BRep, reindexing.
  - **Step 2**: Union-find vertex merge: pairs within `tolerance` are merged to one representative.
  - **Step 3**: Edge deduplication: edges sharing both (merged) endpoint vertices are stitched.
  - **Step 4**: Compacts vertex/edge arrays; assembles single shell with all faces.
  - **Step 5**: Reports free edges (only 1 incident face) for open-boundary diagnosis.
  - GeomStore surfaces and face_surface mappings are concatenated (not de-duplicated).

## 2.16 Analytic Section Curves (rcad-algorithms / section.rs)
- Analogous to OCCT `BRepAlgoAPI_Section` returning proper edge geometry.
- `section_curves(brep: &BRep, plane: &Plane) -> Vec<SectionCurve>`
  - `SectionCurve::Analytic(Curve3)` — exact result for analytic surfaces:
    - `Surface3::Plane` → `inttools::plane_plane` → `Curve3::Line`
    - `Surface3::Sphere` → `inttools::plane_sphere` → `Curve3::Circle`
    - `Surface3::Cylinder` → `inttools::plane_cylinder` → `Curve3::Circle / Ellipse / Line`
    - `Surface3::Cone` → `inttools::plane_cone` → `Curve3::Circle / Ellipse / Line`
  - `SectionCurve::Polyline(Vec<DVec3>)` — triangle-mesh fallback for Torus, BSpline, Bezier, Offset
- Existing `section()` and `section_polylines()` unchanged (backward compatible).

## 2.17 BRep Transform (rcad-kernel / lib.rs)
- Analogous to OCCT `BRepBuilderAPI_Transform` / `TopLoc_Location`.
- `BRep::apply_transform(mat: DAffine3)` — modifies in place:
  - All `vertices[i].point` via `mat.transform_point3`
  - All `Curve3` variants: Line origin/direction, Circle/Ellipse center/normal/major_dir, BSpline/Bezier control points, Offset basis curve (recursive)
  - All `Surface3` variants: Plane/Cylinder/Sphere/Cone/Torus origins and axes, BSpline/Bezier control-point grids, LinearExtrusion direction, Revolution axis; directions normalized after transform
  - Face normals (stored as `Plane.normal` in surfaces)
- `BRep::transformed(mat: DAffine3) -> BRep` — clone + apply (non-destructive).

## 2.18 Curve-Curve Extrema (rcad-kernel / extrema.rs)
- Analogous to OCCT `GeomAPI_ExtremaCurveCurve`.
- `extrema_curve_curve(c1: &Curve3, c2: &Curve3, n_samples: usize) -> CurveCurveExtrema`
  - **Coarse grid**: n×n grid over both curve domains; local minima are collected as seeds.
  - **Boundary seeds**: all four corner (s,t) combinations included to catch boundary minima.
  - **Newton-Raphson refinement**: finite-difference gradient `[2(C1-C2)·C1', -2(C1-C2)·C2']`; Gauss-Newton diagonal Hessian; backtracking line search (≤8 halvings).
  - **Deduplication**: pairs within `DEDUP_TOL = 1e-4` in parameter space are merged.
  - **Output**: `CurveCurveExtrema { pairs: Vec<ExtremaPair> }` sorted distance ascending; `ExtremaPair { param1, param2, point1, point2, distance }`. `min_distance()` convenience method.
  - Line domain clamped to `[-1e6, 1e6]` for infinite-line sampling.

## 2.19 Surface-Surface Intersection (rcad-algorithms / inttools/intss.rs)
- Analogous to OCCT `GeomAPI_IntSS`.
- `intersect_surfaces(s1: &Surface3, s2: &Surface3) -> SurfaceSurfaceIntersection`
- `SurfaceSurfaceIntersection { curves: Vec<SurfaceCurve> }` with `is_empty()` convenience.
- `SurfaceCurve` enum: `Circle(Circle3)`, `Ellipse(Ellipse3)`, `Line(Line3)`, `Point(DVec3)`, `Polyline(Vec<DVec3>)`.
- **Analytic dispatch**:
  - Plane × Plane: reuses `plane_plane_intersection` → Line or None.
  - Plane × Sphere: `plane_sphere_intersection` → Circle, Point, or None.
  - Plane × Cylinder: `plane_cylinder_intersection` → Circle, Ellipse, or None (perpendicular or oblique cut).
  - Sphere × Sphere: radical plane formula — `a = (d²+r1²-r2²)/(2d)` gives axial distance from s1 center to intersection circle; reports Circle or Point.
  - Cylinder × Cylinder (parallel axes): law-of-cosines angle → two generator lines (intersecting) or one tangent line.
- **Numeric fallback** (`numeric_intss`): 48×48 grid on s1 surface, 32×32 implicit-value cache on s2; sign-change detection → `SurfaceCurve::Polyline`.

## 2.20 Rectangular Trimmed Surface (rcad-kernel / geom.rs)
- Analogous to OCCT `Geom_RectangularTrimmedSurface`.
- `TrimmedSurface { basis: Box<Surface3>, trim: [f64; 4] }` with `TrimmedSurface::new(basis, u1, u2, v1, v2)`.
- Added as `Surface3::Trimmed(TrimmedSurface)` variant.
- `SurfaceEval` impl: `point_at`/`normal_at` delegate to basis; `default_domain()` returns `self.trim`.
- `apply_transform`: transforms basis geometry only — trim domain is in parameter space, not world space.
- STEP import: `RECTANGULAR_TRIMMED_SURFACE(name, #basis, u1, u2, v1, v2, .T., .T.)` → `Surface3::Trimmed(...)`.
- STEP export: writer strips to basis surface (trim bounds are encoded in face wire topology).

 Converts analytic B-Rep surfaces to renderable mesh buffers on demand. When `Face.triangles` is pre-populated it is used as a cache; when absent the render pipeline tessellates from the analytic surface. Tessellation MUST NOT be triggered by modeling or export code.
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
- **Color import** (Phase P): `StepReader::parse_string_with_color / read_file_with_color` → `(BRep, Option<StepColor>)`. Parses the full `STYLED_ITEM → PRESENTATION_STYLE_ASSIGNMENT → SURFACE_STYLE_USAGE → SURFACE_SIDE_STYLE → SURFACE_STYLE_FILL_AREA → FILL_AREA_STYLE → FILL_AREA_STYLE_COLOUR → COLOUR_RGB` chain. Maps STEP `ADVANCED_FACE` entity id → flat face index via `face_id_map` built during BRep assembly. Backward-compatible: existing `parse_string` / `read_file` unchanged.
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