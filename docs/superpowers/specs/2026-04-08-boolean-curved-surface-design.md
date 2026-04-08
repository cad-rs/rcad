# Boolean Operations: Curved Surface Support — Phase 1 Design

**Date:** 2026-04-08
**Status:** Approved
**Scope:** Parametric face splitting infrastructure + Plane×Curved and Curved×Curved Boolean operations with FEA/CAE precision
**Phase:** 1 of N — Foundations (parametric trim, IntSS PCurve, face splitting, test coverage)

## Overview

Upgrade rcad's Boolean operations from approximate boundary-based curved face splitting to precise parameter-space 2D clipping. This is the foundational phase that all subsequent curved Boolean improvements build on.

**Current state:** Boolean ops work exactly for planar bodies (box-box). Curved surfaces have an approximate `split_curved_face()` that splits the 3D boundary polygon without creating trimmed parametric sub-domains. IntSS returns only 3D curves, no PCurves.

**Target state:** IntSS returns PCurves alongside 3D curves. Face splitting operates in (u,v) parameter space for curved faces, producing geometrically precise sub-faces. All common surface pairs (Plane×Sphere, Plane×Cylinder, Sphere×Sphere, Cylinder×Cylinder, etc.) pass Boolean tests with FEA/CAE precision.

## Architecture

### Data Flow

```
IntSS: Surface A × Surface B
  → 3D curve (Circle3, Ellipse3, Line3, Polyline)
  → PCurve on A (Curve2d in A's parameter space)
  → PCurve on B (Curve2d in B's parameter space)

DS Loading:
  → For each curved face: compute uv_boundary from edge PCurves or projection

Face splitting:
  → Face boundary → uv_boundary (2D polygon)
  → Intersection curve → PCurve → 2D trim curve
  → 2D polygon clipping in parameter space → sub-face regions
  → Each sub-face → original Surface3 + new Wire boundary (3D remapped)
```

## Component 1: IntSS PCurve Output

### Change to IntSS result type

Extend the IntSS output to carry PCurves alongside 3D intersection curves.

**New struct** (in `rcad-algorithms/src/inttools/intss.rs`):

```rust
pub struct SurfaceIntersectionResult {
    pub curve_3d: SurfaceCurve,           // Existing 3D intersection curve
    pub pcurve_on_a: Option<Curve2d>,     // Projection in surface A's (u,v) domain
    pub pcurve_on_b: Option<Curve2d>,     // Projection in surface B's (u,v) domain
}
```

`SurfaceSurfaceIntersection.curves` changes from `Vec<SurfaceCurve>` to `Vec<SurfaceIntersectionResult>`.

### Analytic PCurve derivation

For analytic intersection pairs, PCurves are derived exactly (no numerical projection):

| Pair | 3D Curve | PCurve on A | PCurve on B |
|------|----------|-------------|-------------|
| Plane × Sphere | Circle3 | Line2d (projected to plane uv) | Circle2d (latitude circle on sphere (θ,φ) domain) |
| Plane × Cylinder | Ellipse3/Circle3 | Line2d or Ellipse2d in plane uv | Line2d in cylinder (θ,h) domain |
| Plane × Cone | Circle3/Ellipse3 | Line2d or Ellipse2d in plane uv | Curve in cone (θ,z) domain |
| Sphere × Sphere | Circle3 | Circle2d (latitude on sphere A) | Circle2d (latitude on sphere B) |
| Sphere × Cylinder (axial) | Circle3 | Circle2d on sphere | Line2d on cylinder (h=const) |
| Cylinder × Cylinder (parallel) | Line3 | Line2d on cyl A | Line2d on cyl B |

### Marching pair PCurves

For pairs solved by numerical marching (non-axial Sphere×Cylinder, general Cylinder×Cylinder, Cylinder×Cone, etc.):

1. The 3D polyline points are projected onto each surface's parameter domain using `closest_point_on_surface` (returns (u,v))
2. The resulting (u,v) point sequences are fit into `BSplineCurve2` using the existing `interpolate_points` infrastructure (adapted for 2D)

## Component 2: DSFace Parameter Domain Boundary

### Extension to DSFace

Add an optional (u,v) boundary to DSFace in `rcad-algorithms/src/bopds/ds.rs`:

```rust
pub struct DSFace {
    // Existing fields (unchanged)
    pub surface: Surface3,
    pub boundary_verts: Vec<usize>,
    pub boundary_edges: Vec<usize>,
    pub normal: DVec3,
    pub face_info: FaceInfo,

    // New field
    pub uv_boundary: Option<Vec<DVec2>>,  // Parameter-space boundary polygon
}
```

### UV boundary computation

Computed during `DS::new()` when loading faces:

- **Plane faces:** `uv_boundary = None` — continue using existing 2D projection logic in `split_planar_face()`
- **Primitive BRep faces (from make_sphere_brep, etc.):** Known parameterization, compute directly. Example: sphere face → `[(0,0), (2π,0), (2π,π), (0,π)]` for full sphere
- **Faces with PCurves (STEP import):** Extract Curve2d from `GeomStore.edge_pcurves`, sample into (u,v) points
- **Fallback (no PCurves):** Project 3D boundary vertices onto surface parameter domain via `closest_point_on_surface` returning (u,v)

## Component 3: Parametric Face Splitting

### Replacement of split_curved_face()

Replace the current approximate `split_curved_face()` in `builder.rs` with `split_curved_face_parametric()`.

**Algorithm:**

```
split_curved_face_parametric(face_idx):
  1. Retrieve DSFace.uv_boundary → 2D polygon P
  2. Collect all IntersectionCurve PCurves for this face → 2D trim curves
  3. Sample each trim curve into a 2D polyline:
     - Analytic PCurves (Line2d, Circle2d): small fixed sample count
     - BSpline2 PCurves: adaptive based on curvature
  4. Split P using trim curves into sub-regions
  5. For each sub-region → SubFace:
     - surface: original Surface3 (wire limits the active region)
     - boundary: new Wire (2D sub-region points remapped to 3D via surface.point_at(u,v))
     - sample_point: interior point of sub-region, mapped to 3D for classification
```

### 2D splitting strategy

**Simple case (single trim curve bisecting polygon):**
1. Find 2 intersection points between trim curve and polygon boundary
2. Split polygon at intersection points
3. Half A + trim curve forward = sub-face A
4. Half B + trim curve reversed = sub-face B

**Complex case (multiple trim curves):**
- Merge trim curves and boundary into a planar graph
- Extract face regions via monotone subdivision
- Each region becomes a sub-face

### 3D remapping

Sub-region (u,v) boundary points → `surface.point_at(u, v)` → 3D points → new Edges/Wires.

Trim curve edges carry both:
- 3D curve (from IntersectionCurve.curve_3d)
- PCurve (from IntSS output, stored directly)

## Component 4: Test Coverage

### Test matrix

| Test Case | Boolean Ops | Key Verification |
|-----------|-------------|------------------|
| Box × Sphere (overlapping) | Union, Intersection, Difference | Planar + spherical sub-faces; volume conservation |
| Box × Cylinder (through-hole) | Difference | Classic hole-in-box; 6 planar + 1 cylindrical face |
| Sphere × Sphere (overlapping) | Union, Intersection, Difference | Circle3 intersection; two spherical caps |
| Cylinder × Cylinder (orthogonal) | Intersection | Steinmetz solid; marching intersection |
| Sphere × Cylinder (axial) | Difference | Analytic pair; two latitude circles |
| Sphere × Cylinder (offset) | Intersection | Marching pair; BSpline2 PCurves |

### Verification criteria (FEA/CAE level)

1. **Topological integrity:** `brep_check()` passes — all faces closed, edges have two endpoints, wires continuous
2. **Geometric precision:** Sub-face boundary point deviation from original surface < `CONFUSION` (1e-7)
3. **Volume conservation:** `volume(A) + volume(B) - volume(A∩B) ≈ volume(A∪B)`, error < 1%
4. **STEP round-trip:** Export to STEP → reimport → `brep_check()` passes
5. **Classification correctness:** Each sub-face sample point classification consistent with Boolean semantics

### Test structure

Each test case:
1. Construct input BReps (primitives with known geometry)
2. Execute boolean_op
3. brep_check on result
4. Volume sanity check (compare against analytical or known value)
5. STEP export verification

## Files Affected

### Modified files

| File | Change |
|------|--------|
| `libs/rcad-algorithms/src/inttools/intss.rs` | New `SurfaceIntersectionResult` struct; all analytic pair functions return PCurves; marching wrapper computes PCurves from polyline projection |
| `libs/rcad-algorithms/src/bopds/ds.rs` | Add `uv_boundary: Option<Vec<DVec2>>` to DSFace; compute in `DS::new()` |
| `libs/rcad-algorithms/src/builder.rs` | Replace `split_curved_face()` with `split_curved_face_parametric()`; update `split_face()` dispatch |
| `libs/rcad-algorithms/src/pave_filler.rs` | Update FF pass to propagate PCurves from IntSS into IntersectionCurve |
| `libs/rcad-algorithms/src/bopds/ds.rs` (IntersectionCurve) | Add `pcurve_on_a: Option<Curve2d>`, `pcurve_on_b: Option<Curve2d>` fields |
| `libs/rcad-algorithms/src/lib.rs` | Add curved Boolean test cases |

### Possibly new files

| File | Purpose |
|------|---------|
| `libs/rcad-algorithms/src/inttools/pcurve_derive.rs` | Analytic PCurve derivation functions per surface pair |
| `libs/rcad-algorithms/src/split_parametric.rs` | Parametric face splitting algorithm (if builder.rs grows too large) |

### Unchanged

| Component | Reason |
|-----------|--------|
| `rcad-kernel` topology types | Face/Wire/Edge structs unchanged; PCurve infrastructure already exists |
| `rcad-kernel` Curve2d enum | Already has all needed 2D types (Line2d, Circle2d, Ellipse2d, BSpline2, Bezier2) |
| `classify.rs` | Point classification already supports Sphere/Cylinder/Cone; Torus is Phase 2 |
| `rcad-step`, `rcad-render`, `rcad-scene` | No changes needed for Phase 1 |

## Out of Scope (Phase 2+)

- Torus Boolean support (classification + IntSS pairs)
- BSpline surface Boolean support
- Edge-Face intersection for non-Line/Circle curves
- Vertex-Face classification for curved faces
- Coplanar curved face handling
- Performance optimization (parallel marching, adaptive grid)

## Decisions

| Decision | Rationale |
|----------|-----------|
| PCurves in IntSS output | Enables parameter-space splitting; avoids lossy 3D-to-2D re-projection |
| Analytic PCurve derivation per pair | Exact, no numerical error; matches OCCT approach |
| BSpline2 for marching PCurves | Smooth parametric representation; fits existing Curve2d infrastructure |
| Optional uv_boundary (None for planes) | Zero overhead for planar faces; existing logic unchanged |
| Original Surface3, not TrimmedSurface, for sub-faces | Wire boundary already limits the active region; TrimmedSurface adds unnecessary indirection |
| Monotone subdivision for complex cases | Handles arbitrary curve topology; proven algorithm for planar graphs |
