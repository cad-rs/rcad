# Boolean Curved Surface Support (Phase 1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Upgrade Boolean operations from approximate boundary-based curved face splitting to precise parameter-space 2D clipping, enabling FEA/CAE-grade Boolean results for all common analytic surface pairs.

**Architecture:** Extend IntSS to return PCurves alongside 3D curves. Add UV boundary computation to DSFace. Replace `split_curved_face()` with parameter-space 2D polygon clipping. Add comprehensive tests for curved Boolean operations.

**Tech Stack:** Rust, glam (DVec2/DVec3), rcad-kernel (Curve2d, SurfaceEval, projection), rcad-algorithms (BOPDS, IntSS, PaveFiller, BooleanBuilder)

**Spec:** `docs/superpowers/specs/2026-04-08-boolean-curved-surface-design.md`

---

## File Structure

| File | Responsibility | Action |
|------|---------------|--------|
| `libs/rcad-algorithms/src/inttools/intss.rs` | IntSS dispatch + result types | Modify: new `SurfaceIntersectionResult` struct, update return types |
| `libs/rcad-algorithms/src/inttools/pcurve_derive.rs` | Analytic PCurve derivation per surface pair | Create |
| `libs/rcad-algorithms/src/inttools/mod.rs` | Module registration | Modify: add `pub mod pcurve_derive` |
| `libs/rcad-algorithms/src/bopds/ds.rs` | DS data structures | Modify: add `uv_boundary` to DSFace, `pcurve_on_a/b` to IntersectionCurve |
| `libs/rcad-algorithms/src/pave_filler.rs` | FF pass propagates PCurves | Modify: use new IntSS result, store PCurves in IntersectionCurve |
| `libs/rcad-algorithms/src/builder.rs` | Face splitting | Modify: replace `split_curved_face()` with `split_curved_face_parametric()` |
| `libs/rcad-kernel/src/fit.rs` | 2D point interpolation | Modify: add `interpolate_points_2d()` |
| `libs/rcad-algorithms/src/lib.rs` | Tests | Modify: add curved Boolean tests |

---

### Task 1: Add `interpolate_points_2d` to rcad-kernel

Marching PCurves need 2D B-spline fitting. The existing `interpolate_points` works on `DVec3`. We need a `DVec2` variant.

**Files:**
- Modify: `libs/rcad-kernel/src/fit.rs`
- Test: `libs/rcad-kernel/src/fit.rs` (inline tests)

- [ ] **Step 1: Write the failing test**

Add to the existing `#[cfg(test)] mod tests` block in `libs/rcad-kernel/src/fit.rs`:

```rust
#[test]
fn interpolate_2d_line() {
    let pts = vec![
        DVec2::new(0.0, 0.0),
        DVec2::new(1.0, 1.0),
        DVec2::new(2.0, 2.0),
    ];
    let bsp = interpolate_points_2d(&pts).unwrap();
    // Check endpoints
    let p0 = bsp.point_at(0.0);
    let p2 = bsp.point_at(1.0);
    assert!((p0 - DVec2::new(0.0, 0.0)).length() < 1e-10);
    assert!((p2 - DVec2::new(2.0, 2.0)).length() < 1e-10);
    // Check midpoint
    let pm = bsp.point_at(0.5);
    assert!((pm - DVec2::new(1.0, 1.0)).length() < 1e-6);
}

#[test]
fn interpolate_2d_circle_arc() {
    use std::f64::consts::PI;
    // Quarter circle in 2D
    let n = 9;
    let pts: Vec<DVec2> = (0..=n)
        .map(|i| {
            let t = PI / 2.0 * i as f64 / n as f64;
            DVec2::new(t.cos(), t.sin())
        })
        .collect();
    let bsp = interpolate_points_2d(&pts).unwrap();
    // Check that interpolated midpoint is close to circle
    let pm = bsp.point_at(0.5);
    let r = pm.length();
    assert!((r - 1.0).abs() < 0.01, "midpoint radius {r} should be ~1.0");
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p rcad-kernel interpolate_2d -- --nocapture`
Expected: FAIL — `interpolate_points_2d` not found.

- [ ] **Step 3: Implement `interpolate_points_2d`**

Add to `libs/rcad-kernel/src/fit.rs`, after the existing `interpolate_points` function. The logic mirrors `interpolate_points` but operates on `DVec2` and returns `BSplineCurve2`:

```rust
use glam::DVec2;
use crate::geom::BSplineCurve2;

/// Interpolate a B-spline curve through the given 2D points.
///
/// Uses chord-length parameterization and cubic degree with clamped knots.
/// Analogous to `interpolate_points` but for 2D parameter-space curves (PCurves).
pub fn interpolate_points_2d(points: &[DVec2]) -> Result<BSplineCurve2, FitError> {
    let n = points.len();
    if n < 2 {
        return Err(FitError::TooFewPoints);
    }
    if n == 2 {
        // Degenerate to linear
        return Ok(BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: points.to_vec(),
            weights: vec![1.0; 2],
        });
    }

    let degree = 3.min(n - 1);

    // Chord-length parameterization
    let mut params = vec![0.0f64; n];
    let mut total = 0.0;
    for i in 1..n {
        total += (points[i] - points[i - 1]).length();
        params[i] = total;
    }
    if total < 1e-15 {
        return Err(FitError::DegeneratePoints);
    }
    for p in &mut params {
        *p /= total;
    }

    // Clamped knot vector
    let m = n + degree + 1;
    let mut knots = vec![0.0; m];
    for k in &mut knots[..=degree] {
        *k = 0.0;
    }
    for k in &mut knots[m - degree - 1..] {
        *k = 1.0;
    }
    let denom = (n - degree) as f64;
    for j in 1..n - degree {
        let mut sum = 0.0;
        for i in j..j + degree {
            sum += params[i];
        }
        knots[j + degree] = sum / degree as f64;
    }

    // Build and solve the interpolation linear system
    // N(i,degree)(params[j]) * ctrl[i] = points[j]
    let mut mat = vec![vec![0.0; n]; n];
    for (row, &t) in params.iter().enumerate() {
        let basis = all_basis_fns_2d(t, degree, &knots, n);
        for (col, &b) in basis.iter().enumerate() {
            mat[row][col] = b;
        }
    }

    // Solve for x and y separately using Gaussian elimination
    let mut rhs_x: Vec<f64> = points.iter().map(|p| p.x).collect();
    let mut rhs_y: Vec<f64> = points.iter().map(|p| p.y).collect();

    let mut aug_x = mat.clone();
    let mut aug_y = mat;

    // Forward elimination with partial pivoting
    for col in 0..n {
        let (pivot_row, _) = aug_x[col..]
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a[col].abs().partial_cmp(&b[col].abs()).unwrap())
            .unwrap();
        let pivot_row = pivot_row + col;
        aug_x.swap(col, pivot_row);
        aug_y.swap(col, pivot_row);
        rhs_x.swap(col, pivot_row);
        rhs_y.swap(col, pivot_row);

        let diag = aug_x[col][col];
        if diag.abs() < 1e-15 {
            return Err(FitError::DegeneratePoints);
        }
        for row in col + 1..n {
            let factor = aug_x[row][col] / diag;
            for k in col..n {
                let val = aug_x[col][k];
                aug_x[row][k] -= factor * val;
                let val_y = aug_y[col][k];
                aug_y[row][k] -= factor * val_y;
            }
            rhs_x[row] -= factor * rhs_x[col];
            rhs_y[row] -= factor * rhs_y[col];
        }
    }

    // Back substitution
    let mut ctrl_x = vec![0.0; n];
    let mut ctrl_y = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sx = rhs_x[i];
        let mut sy = rhs_y[i];
        for j in i + 1..n {
            sx -= aug_x[i][j] * ctrl_x[j];
            sy -= aug_y[i][j] * ctrl_y[j];
        }
        ctrl_x[i] = sx / aug_x[i][i];
        ctrl_y[i] = sy / aug_y[i][i];
    }

    let control_points: Vec<DVec2> = ctrl_x
        .iter()
        .zip(ctrl_y.iter())
        .map(|(&x, &y)| DVec2::new(x, y))
        .collect();

    Ok(BSplineCurve2 {
        degree,
        knots,
        control_points,
        weights: vec![1.0; n],
    })
}

/// Evaluate all n basis functions of given degree at parameter t.
fn all_basis_fns_2d(t: f64, degree: usize, knots: &[f64], n: usize) -> Vec<f64> {
    let mut basis = vec![0.0; n];
    // de Boor-Cox recurrence
    let m = knots.len();
    let mut n0 = vec![0.0; m - 1];
    for i in 0..m - 1 {
        n0[i] = if knots[i] <= t && t < knots[i + 1] {
            1.0
        } else {
            0.0
        };
    }
    // Handle right endpoint
    if (t - knots[m - 1]).abs() < 1e-15 {
        if m >= 2 {
            n0[m - 2] = 1.0;
        }
    }

    let mut prev = n0;
    for d in 1..=degree {
        let mut curr = vec![0.0; prev.len()];
        for i in 0..prev.len() - d {
            let denom1 = knots[i + d] - knots[i];
            let denom2 = knots[i + d + 1] - knots[i + 1];
            let left = if denom1.abs() > 1e-15 {
                (t - knots[i]) / denom1 * prev[i]
            } else {
                0.0
            };
            let right = if denom2.abs() > 1e-15 {
                (knots[i + d + 1] - t) / denom2 * prev[i + 1]
            } else {
                0.0
            };
            curr[i] = left + right;
        }
        prev = curr;
    }
    for (i, b) in basis.iter_mut().enumerate() {
        *b = prev[i];
    }
    basis
}
```

Also add the `Curve2dEval` implementation for `BSplineCurve2` `point_at` if not already present — check `geom.rs` for existing impl. (It should already exist since `Curve2dEval` is implemented for `Curve2d` which includes `BSpline(BSplineCurve2)`.)

- [ ] **Step 4: Add necessary imports at top of fit.rs**

Add `use glam::DVec2;` and `use crate::geom::BSplineCurve2;` to the imports in `fit.rs`. Also add the `Curve2dEval` trait import in the test module.

- [ ] **Step 5: Run tests to verify they pass**

Run: `cargo test -p rcad-kernel interpolate_2d -- --nocapture`
Expected: both tests PASS.

- [ ] **Step 6: Run full test suite**

Run: `cargo test --workspace`
Expected: all tests pass.

- [ ] **Step 7: Commit**

```bash
git add libs/rcad-kernel/src/fit.rs
git commit -m "feat(kernel): add interpolate_points_2d for 2D B-spline fitting"
```

---

### Task 2: Extend IntersectionCurve and SurfaceSurfaceIntersection with PCurves

Add PCurve fields to the data structures that carry intersection results through the pipeline.

**Files:**
- Modify: `libs/rcad-algorithms/src/bopds/ds.rs`
- Modify: `libs/rcad-algorithms/src/inttools/intss.rs`

- [ ] **Step 1: Add PCurve fields to IntersectionCurve**

In `libs/rcad-algorithms/src/bopds/ds.rs`, add to `IntersectionCurve`:

```rust
use rcad_kernel::geom::Curve2d;

/// An intersection curve from F-F intersection, bounded by vertices.
#[derive(Debug, Clone)]
pub struct IntersectionCurve {
    pub curve: Curve3,
    /// Sampled points from numerical marching (non-empty for marched curves).
    /// When non-empty this takes priority over `curve` for face splitting.
    pub polyline: Vec<DVec3>,
    pub start_vertex: usize,
    pub end_vertex: usize,
    pub t_range: [f64; 2],
    /// PCurve on the first face's surface (parameter-space representation).
    pub pcurve_on_a: Option<Curve2d>,
    /// PCurve on the second face's surface (parameter-space representation).
    pub pcurve_on_b: Option<Curve2d>,
}
```

- [ ] **Step 2: Add SurfaceIntersectionResult to intss.rs**

In `libs/rcad-algorithms/src/inttools/intss.rs`, add the new result type and update `SurfaceSurfaceIntersection`:

```rust
use rcad_kernel::geom::Curve2d;

/// A single intersection component with optional PCurves on both surfaces.
#[derive(Debug, Clone)]
pub struct SurfaceIntersectionResult {
    /// The 3D intersection curve.
    pub curve_3d: SurfaceCurve,
    /// Image of the intersection curve in surface A's (u,v) parameter domain.
    pub pcurve_on_a: Option<Curve2d>,
    /// Image of the intersection curve in surface B's (u,v) parameter domain.
    pub pcurve_on_b: Option<Curve2d>,
}

/// All intersection curves / components found between two surfaces.
#[derive(Debug, Clone, Default)]
pub struct SurfaceSurfaceIntersection {
    pub curves: Vec<SurfaceIntersectionResult>,
}
```

- [ ] **Step 3: Update all existing IntSS functions to return `SurfaceIntersectionResult`**

Every place that pushes to `out.curves` needs updating. For now, wrap with `pcurve_on_a: None, pcurve_on_b: None` — PCurve derivation is Task 3.

For example, in `plane_x_sphere`:

```rust
fn plane_x_sphere(p: &Plane, s: &SphericalSurface) -> SurfaceSurfaceIntersection {
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_plane_sphere(p, s) {
        PlaneSphereResult::Circle(c) => {
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(c),
                pcurve_on_a: None,
                pcurve_on_b: None,
            });
        }
        PlaneSphereResult::TangentPoint(pt) => {
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Point(pt),
                pcurve_on_a: None,
                pcurve_on_b: None,
            });
        }
        PlaneSphereResult::NoIntersection => {}
    }
    out
}
```

Apply the same pattern to: `plane_x_plane`, `plane_x_cylinder`, `plane_x_cone`, `sphere_x_sphere`, `sphere_x_cylinder`, `sphere_x_cone`, `cylinder_x_cylinder`, `numeric_intss`.

- [ ] **Step 4: Update all callers of SurfaceSurfaceIntersection.curves**

In `pave_filler.rs`, update `intersect_ff_by_marching` and anywhere that accesses `.curves[i]` to use `.curves[i].curve_3d` instead. Also update the IntersectionCurve construction to pass `pcurve_on_a: None, pcurve_on_b: None`.

In `intss.rs` tests, update assertions: `r.curves[0].curve_3d` instead of `r.curves[0]`.

- [ ] **Step 5: Fix all compilation errors**

Run: `cargo build --workspace`
Expected: compiles cleanly. Fix any remaining references to old field names.

- [ ] **Step 6: Run tests**

Run: `cargo test --workspace`
Expected: all tests pass (behavior unchanged, just struct shapes changed).

- [ ] **Step 7: Commit**

```bash
git add libs/rcad-algorithms/src/bopds/ds.rs libs/rcad-algorithms/src/inttools/intss.rs libs/rcad-algorithms/src/pave_filler.rs
git commit -m "refactor(algorithms): extend IntersectionCurve and IntSS with PCurve fields"
```

---

### Task 3: Implement analytic PCurve derivation

Derive exact PCurves for each analytic IntSS pair. This is the mathematical core of the improvement.

**Files:**
- Create: `libs/rcad-algorithms/src/inttools/pcurve_derive.rs`
- Modify: `libs/rcad-algorithms/src/inttools/mod.rs`
- Modify: `libs/rcad-algorithms/src/inttools/intss.rs`

- [ ] **Step 1: Create pcurve_derive.rs with module structure and first function**

Create `libs/rcad-algorithms/src/inttools/pcurve_derive.rs`:

```rust
//! Analytic PCurve derivation for surface-surface intersection pairs.
//!
//! Given a 3D intersection curve between two surfaces, compute the corresponding
//! 2D curve (PCurve) in each surface's (u,v) parameter domain.

use std::f64::consts::{PI, TAU};
use glam::{DVec2, DVec3};
use rcad_kernel::geom::*;
use rcad_kernel::projection::closest_point_on_surface;

/// Derive the PCurve of a Circle3 on a Plane's (u,v) domain.
///
/// A circle intersecting a plane is projected to a circle or ellipse
/// in the plane's local (u,v) coordinates.
pub fn circle_pcurve_on_plane(circle: &Circle3, plane: &Plane) -> Curve2d {
    // The plane's local basis: u_axis and v_axis
    let u_axis = any_perpendicular(plane.normal);
    let v_axis = plane.normal.cross(u_axis).normalize();

    // The circle center in plane (u,v)
    let d = circle.center - plane.origin;
    let cu = d.dot(u_axis);
    let cv = d.dot(v_axis);

    // The circle lies in the plane (since it's Plane×Sphere/Cylinder intersection).
    // Its image in (u,v) is a circle with the same radius.
    // Circle parametrization: center + r*(cos(t)*e1 + sin(t)*e2)
    // where e1, e2 are orthonormal in the circle's plane.
    // Since the circle lies in the plane, e1 and e2 project to vectors of length 1.

    // Find the circle's local axes
    let e1 = any_perpendicular(circle.normal);
    let e2 = circle.normal.cross(e1).normalize();

    // Project e1, e2 into plane (u,v) space
    let e1u = e1.dot(u_axis);
    let e1v = e1.dot(v_axis);
    let e2u = e2.dot(u_axis);
    let e2v = e2.dot(v_axis);

    // The 2D parametric curve is:
    // P(t) = (cu, cv) + r * (cos(t)*(e1u,e1v) + sin(t)*(e2u,e2v))
    // This is an ellipse in general form. Check if it's a circle.
    let r = circle.radius;
    let a_sq = e1u * e1u + e1v * e1v;
    let b_sq = e2u * e2u + e2v * e2v;
    let dot = e1u * e2u + e1v * e2v;

    if (a_sq - b_sq).abs() < 1e-10 && dot.abs() < 1e-10 {
        // It's a circle in 2D
        Curve2d::Circle(Circle2d {
            center: DVec2::new(cu, cv),
            radius: r * a_sq.sqrt(),
        })
    } else {
        // General ellipse — use BSpline approximation via sampling
        let n = 33;
        let pts: Vec<DVec2> = (0..n)
            .map(|i| {
                let t = TAU * i as f64 / (n - 1) as f64;
                DVec2::new(
                    cu + r * (t.cos() * e1u + t.sin() * e2u),
                    cv + r * (t.cos() * e1v + t.sin() * e2v),
                )
            })
            .collect();
        match rcad_kernel::fit::interpolate_points_2d(&pts) {
            Ok(bsp) => Curve2d::BSpline(bsp),
            Err(_) => fallback_pcurve_by_projection(
                &Curve3::Circle(circle.clone()),
                &[0.0, TAU],
                &Surface3::Plane(*plane),
            ),
        }
    }
}

/// Derive the PCurve of a Circle3 on a SphericalSurface's (θ,φ) domain.
///
/// For a sphere parameterized as (θ,φ) where θ ∈ [0,2π] (longitude) and
/// φ ∈ [0,π] (colatitude from north pole along `axis`), a circle of
/// intersection at a given latitude maps to a horizontal line φ = φ₀.
pub fn circle_pcurve_on_sphere(circle: &Circle3, sphere: &SphericalSurface) -> Curve2d {
    // The circle's center-to-sphere-center vector projected onto sphere axis
    let d = circle.center - sphere.center;
    let along_axis = d.dot(sphere.axis);

    // Colatitude: φ = acos(along_axis / radius)
    let cos_phi = (along_axis / sphere.radius).clamp(-1.0, 1.0);
    let phi = cos_phi.acos();

    // The circle maps to a horizontal line at v = phi in the sphere's (θ,φ) domain.
    // θ ranges over [0, 2π].
    Curve2d::Line(Line2d {
        origin: DVec2::new(0.0, phi),
        direction: DVec2::new(1.0, 0.0),
    })
}

/// Derive the PCurve of a Circle3 on a CylindricalSurface's (θ,h) domain.
///
/// A circle perpendicular to the cylinder axis at height h maps to
/// a horizontal line h = h₀ in the cylinder's (θ,h) domain.
pub fn circle_pcurve_on_cylinder(circle: &Circle3, cyl: &CylindricalSurface) -> Curve2d {
    let d = circle.center - cyl.origin;
    let h = d.dot(cyl.axis);

    // The circle at height h maps to a horizontal line in (θ,h) space.
    Curve2d::Line(Line2d {
        origin: DVec2::new(0.0, h),
        direction: DVec2::new(1.0, 0.0),
    })
}

/// Derive the PCurve of an Ellipse3 on a Plane's (u,v) domain.
pub fn ellipse_pcurve_on_plane(ellipse: &Ellipse3, plane: &Plane) -> Curve2d {
    let u_axis = any_perpendicular(plane.normal);
    let v_axis = plane.normal.cross(u_axis).normalize();

    let d = ellipse.center - plane.origin;
    let cu = d.dot(u_axis);
    let cv = d.dot(v_axis);

    // Project ellipse major axis into plane (u,v)
    let mu = ellipse.major_dir.dot(u_axis);
    let mv = ellipse.major_dir.dot(v_axis);
    let major_dir_2d = DVec2::new(mu, mv);
    let len = major_dir_2d.length();

    if len < 1e-12 {
        // Degenerate
        return fallback_pcurve_by_projection(
            &Curve3::Ellipse(*ellipse),
            &[0.0, TAU],
            &Surface3::Plane(*plane),
        );
    }

    Curve2d::Ellipse(Ellipse2d {
        center: DVec2::new(cu, cv),
        major_dir: major_dir_2d / len,
        major_radius: ellipse.semi_major * len,
        minor_radius: ellipse.semi_minor,
    })
}

/// Derive the PCurve of a Line3 on a CylindricalSurface's (θ,h) domain.
///
/// A line parallel to the cylinder axis at angle θ₀ maps to a vertical
/// line in (θ,h) space.
pub fn line_pcurve_on_cylinder(line: &Line3, cyl: &CylindricalSurface) -> Curve2d {
    // Project line origin onto the cylinder's cross-section
    let d = line.origin - cyl.origin;
    let d_perp = d - cyl.axis * d.dot(cyl.axis);

    let x_dir = any_perpendicular(cyl.axis);
    let y_dir = cyl.axis.cross(x_dir).normalize();

    let theta = d_perp.dot(y_dir).atan2(d_perp.dot(x_dir));
    let h = d.dot(cyl.axis);

    // Line parallel to axis → vertical line in (θ,h) space
    Curve2d::Line(Line2d {
        origin: DVec2::new(theta, h),
        direction: DVec2::new(0.0, 1.0),
    })
}

/// Derive the PCurve of a Line3 on a Plane's (u,v) domain.
pub fn line_pcurve_on_plane(line: &Line3, plane: &Plane) -> Curve2d {
    let u_axis = any_perpendicular(plane.normal);
    let v_axis = plane.normal.cross(u_axis).normalize();

    let d = line.origin - plane.origin;
    let ou = d.dot(u_axis);
    let ov = d.dot(v_axis);
    let du = line.direction.dot(u_axis);
    let dv = line.direction.dot(v_axis);

    Curve2d::Line(Line2d {
        origin: DVec2::new(ou, ov),
        direction: DVec2::new(du, dv),
    })
}

/// Fallback: compute PCurve by projecting 3D curve sample points onto a surface.
///
/// Used when analytic derivation is not available (e.g., marching results,
/// or complex surface pairs).
pub fn fallback_pcurve_by_projection(
    curve: &Curve3,
    t_range: &[f64; 2],
    surface: &Surface3,
) -> Curve2d {
    let n = 33;
    let pts: Vec<DVec2> = (0..n)
        .map(|i| {
            let t = t_range[0] + (t_range[1] - t_range[0]) * i as f64 / (n - 1) as f64;
            let p3d = curve.point_at(t);
            let proj = closest_point_on_surface(p3d, surface);
            DVec2::new(proj.params.0, proj.params.1)
        })
        .collect();

    match rcad_kernel::fit::interpolate_points_2d(&pts) {
        Ok(bsp) => Curve2d::BSpline(bsp),
        Err(_) => {
            // Last resort: linear 2D curve between endpoints
            let start = pts[0];
            let end = *pts.last().unwrap();
            Curve2d::Line(Line2d {
                origin: start,
                direction: (end - start).normalize_or_zero(),
            })
        }
    }
}

/// Derive PCurve from a 3D polyline by projecting each point onto the surface.
///
/// Returns a BSpline2 fitted through the projected (u,v) points.
pub fn polyline_pcurve_by_projection(polyline: &[DVec3], surface: &Surface3) -> Option<Curve2d> {
    if polyline.len() < 2 {
        return None;
    }
    let uv_pts: Vec<DVec2> = polyline
        .iter()
        .map(|&p| {
            let proj = closest_point_on_surface(p, surface);
            DVec2::new(proj.params.0, proj.params.1)
        })
        .collect();

    rcad_kernel::fit::interpolate_points_2d(&uv_pts)
        .ok()
        .map(Curve2d::BSpline)
}
```

- [ ] **Step 2: Register the module**

In `libs/rcad-algorithms/src/inttools/mod.rs`, add:

```rust
pub mod pcurve_derive;
```

- [ ] **Step 3: Run build to verify compilation**

Run: `cargo build -p rcad-algorithms`
Expected: compiles cleanly.

- [ ] **Step 4: Write tests for PCurve derivation**

Add a test file section at the bottom of `pcurve_derive.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::Curve2dEval;

    #[test]
    fn circle_on_plane_is_circle() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let circle = Circle3 {
            center: DVec3::new(1.0, 2.0, 0.0),
            normal: DVec3::Z,
            radius: 0.5,
        };
        let pc = circle_pcurve_on_plane(&circle, &plane);
        // Circle in XY plane → circle in (u,v)
        if let Curve2d::Circle(c2d) = &pc {
            assert!((c2d.radius - 0.5).abs() < 1e-6);
        } else {
            panic!("expected Circle2d, got {:?}", pc);
        }
    }

    #[test]
    fn circle_on_sphere_is_latitude() {
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
        };
        // Circle at z = 1.0 on sphere of radius 2
        let circle = Circle3 {
            center: DVec3::new(0.0, 0.0, 1.0),
            normal: DVec3::Z,
            radius: (4.0_f64 - 1.0).sqrt(), // sqrt(3)
        };
        let pc = circle_pcurve_on_sphere(&circle, &sphere);
        if let Curve2d::Line(l) = &pc {
            // φ = acos(1/2) = π/3
            let expected_phi = (0.5_f64).acos();
            assert!((l.origin.y - expected_phi).abs() < 1e-6, "phi should be π/3");
            assert!(l.direction.x.abs() > 0.9, "direction should be along θ");
        } else {
            panic!("expected Line2d for latitude, got {:?}", pc);
        }
    }

    #[test]
    fn circle_on_cylinder_is_h_line() {
        let cyl = CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        };
        let circle = Circle3 {
            center: DVec3::new(0.0, 0.0, 3.0),
            normal: DVec3::Z,
            radius: 1.0,
        };
        let pc = circle_pcurve_on_cylinder(&circle, &cyl);
        if let Curve2d::Line(l) = &pc {
            assert!((l.origin.y - 3.0).abs() < 1e-6, "h should be 3.0");
            assert!(l.direction.x.abs() > 0.9, "direction should be along θ");
        } else {
            panic!("expected Line2d, got {:?}", pc);
        }
    }

    #[test]
    fn fallback_projection_produces_bspline() {
        let circle = Curve3::Circle(Circle3 {
            center: DVec3::ZERO,
            normal: DVec3::Z,
            radius: 1.0,
        });
        let sphere = Surface3::Sphere(SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
        });
        let pc = fallback_pcurve_by_projection(&circle, &[0.0, TAU], &sphere);
        assert!(matches!(pc, Curve2d::BSpline(_)), "fallback should produce BSpline");
    }
}
```

- [ ] **Step 5: Run tests**

Run: `cargo test -p rcad-algorithms pcurve_derive -- --nocapture`
Expected: all 4 tests pass.

- [ ] **Step 6: Commit**

```bash
git add libs/rcad-algorithms/src/inttools/pcurve_derive.rs libs/rcad-algorithms/src/inttools/mod.rs
git commit -m "feat(algorithms): add analytic PCurve derivation for IntSS surface pairs"
```

---

### Task 4: Wire PCurves into IntSS and PaveFiller

Connect the PCurve derivation to the actual IntSS dispatch and propagate through PaveFiller into IntersectionCurve.

**Files:**
- Modify: `libs/rcad-algorithms/src/inttools/intss.rs`
- Modify: `libs/rcad-algorithms/src/pave_filler.rs`

- [ ] **Step 1: Update plane_x_sphere to derive PCurves**

In `intss.rs`, update `plane_x_sphere`:

```rust
fn plane_x_sphere(p: &Plane, s: &SphericalSurface) -> SurfaceSurfaceIntersection {
    let mut out = SurfaceSurfaceIntersection::default();
    match intersect_plane_sphere(p, s) {
        PlaneSphereResult::Circle(c) => {
            use crate::inttools::pcurve_derive::*;
            let pca = circle_pcurve_on_plane(&c, p);
            let pcb = circle_pcurve_on_sphere(&c, s);
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Circle(c),
                pcurve_on_a: Some(pca),
                pcurve_on_b: Some(pcb),
            });
        }
        PlaneSphereResult::TangentPoint(pt) => {
            out.curves.push(SurfaceIntersectionResult {
                curve_3d: SurfaceCurve::Point(pt),
                pcurve_on_a: None,
                pcurve_on_b: None,
            });
        }
        PlaneSphereResult::NoIntersection => {}
    }
    out
}
```

- [ ] **Step 2: Update plane_x_cylinder to derive PCurves**

Similar pattern — for `Circle` result use `circle_pcurve_on_plane` + `circle_pcurve_on_cylinder`; for `Ellipse` result use `ellipse_pcurve_on_plane` + `fallback_pcurve_by_projection`; for `TwoLines` use `line_pcurve_on_plane` + `line_pcurve_on_cylinder`.

- [ ] **Step 3: Update sphere_x_sphere to derive PCurves**

Use `circle_pcurve_on_sphere` for both surfaces.

- [ ] **Step 4: Update remaining analytic functions**

For `plane_x_plane` (Line result): `line_pcurve_on_plane` for both.
For `plane_x_cone`: `circle_pcurve_on_plane` + `fallback_pcurve_by_projection` (cone PCurve is complex).
For `sphere_x_cylinder`, `sphere_x_cone`, `cylinder_x_cylinder`: use available analytic functions or `fallback_pcurve_by_projection`.

- [ ] **Step 5: Update numeric_intss to derive PCurves from polylines**

```rust
fn numeric_intss(s1: &Surface3, s2: &Surface3) -> SurfaceSurfaceIntersection {
    // ... existing code to produce intersection_pts ...
    let mut out = SurfaceSurfaceIntersection::default();
    if intersection_pts.len() >= 2 {
        use crate::inttools::pcurve_derive::polyline_pcurve_by_projection;
        let pca = polyline_pcurve_by_projection(&intersection_pts, s1);
        let pcb = polyline_pcurve_by_projection(&intersection_pts, s2);
        out.curves.push(SurfaceIntersectionResult {
            curve_3d: SurfaceCurve::Polyline(intersection_pts),
            pcurve_on_a: pca,
            pcurve_on_b: pcb,
        });
    }
    out
}
```

- [ ] **Step 6: Update PaveFiller to propagate PCurves into IntersectionCurve**

In `pave_filler.rs`, wherever `IntersectionCurve` is constructed, pass through the PCurves. The key is that the PaveFiller calls `intersect_surfaces` indirectly via the face-face functions. Update `intersect_plane_sphere_faces`, `intersect_plane_cylinder_faces`, and `intersect_ff_by_marching` to:

1. Call `intersect_surfaces(&s1, &s2)` which now returns `SurfaceIntersectionResult` with PCurves
2. Store `pcurve_on_a` and `pcurve_on_b` from the result into the `IntersectionCurve`

For `intersect_plane_sphere_faces`, change the circle handling to use IntSS:

```rust
PlaneSphereResult::Circle(circle) => {
    let intss_result = crate::inttools::intss::intersect_surfaces(
        &Surface3::Plane(*plane),
        &Surface3::Sphere(*sphere),
    );
    if let Some(r) = intss_result.curves.first() {
        // ... create IntersectionCurve with r.pcurve_on_a, r.pcurve_on_b
    }
}
```

Alternatively, for simpler refactoring, call the pcurve_derive functions directly in the PaveFiller face-face methods and pass them into IntersectionCurve construction.

- [ ] **Step 7: Run tests**

Run: `cargo test --workspace`
Expected: all existing tests pass. PCurves are now flowing through the pipeline.

- [ ] **Step 8: Commit**

```bash
git add libs/rcad-algorithms/src/inttools/intss.rs libs/rcad-algorithms/src/pave_filler.rs
git commit -m "feat(algorithms): wire PCurve derivation into IntSS and PaveFiller pipeline"
```

---

### Task 5: Add UV boundary computation to DSFace

Compute parameter-space boundary polygons for curved faces during DS loading.

**Files:**
- Modify: `libs/rcad-algorithms/src/bopds/ds.rs`

- [ ] **Step 1: Add uv_boundary field to DSFace**

```rust
use glam::DVec2;

pub struct DSFace {
    // ... existing fields ...
    /// Parameter-space (u,v) boundary polygon. None for planar faces.
    pub uv_boundary: Option<Vec<DVec2>>,
}
```

Update the DSFace construction in `load_brep` to initialize `uv_boundary: None`.

- [ ] **Step 2: Implement UV boundary computation**

Add a method to DS:

```rust
impl DS {
    /// Compute UV boundary for all curved faces by projecting 3D boundary
    /// vertices onto the face surface's parameter domain.
    pub fn compute_uv_boundaries(&mut self) {
        for fi in 0..self.faces.len() {
            if matches!(self.faces[fi].surface, Surface3::Plane(_)) {
                continue; // Planar faces use existing 2D projection
            }

            let surface = self.faces[fi].surface.clone();
            let boundary_pts: Vec<DVec3> = self.faces[fi]
                .boundary_verts
                .iter()
                .map(|&vi| self.vertices[vi].point)
                .collect();

            if boundary_pts.is_empty() {
                continue;
            }

            let uv_pts: Vec<DVec2> = boundary_pts
                .iter()
                .map(|&p| {
                    let proj = rcad_kernel::projection::closest_point_on_surface(p, &surface);
                    DVec2::new(proj.params.0, proj.params.1)
                })
                .collect();

            self.faces[fi].uv_boundary = Some(uv_pts);
        }
    }
}
```

- [ ] **Step 3: Call compute_uv_boundaries in DS::new**

At the end of `DS::new`, after `load_brep(b, ...)`:

```rust
ds.compute_uv_boundaries();
ds
```

- [ ] **Step 4: Write test**

Add to the existing tests in `ds.rs`:

```rust
#[test]
fn ds_sphere_has_uv_boundary() {
    use rcad_modeling::make_sphere_brep;

    let a = make_sphere_brep(DVec3::ZERO, 1.0).unwrap();
    let b = make_sphere_brep(DVec3::new(1.0, 0.0, 0.0), 1.0).unwrap();
    let ds = DS::new(&a, &b);

    // Sphere faces should have uv_boundary computed
    let sphere_faces: Vec<_> = ds
        .faces
        .iter()
        .filter(|f| matches!(f.surface, Surface3::Sphere(_)))
        .collect();
    assert!(!sphere_faces.is_empty(), "should have sphere faces");
    for f in &sphere_faces {
        assert!(f.uv_boundary.is_some(), "sphere face should have uv_boundary");
        let uv = f.uv_boundary.as_ref().unwrap();
        assert!(uv.len() >= 3, "uv boundary should have at least 3 points");
    }
}
```

- [ ] **Step 5: Run tests**

Run: `cargo test -p rcad-algorithms ds_sphere_has_uv -- --nocapture`
Expected: PASS.

Run: `cargo test --workspace`
Expected: all tests pass.

- [ ] **Step 6: Commit**

```bash
git add libs/rcad-algorithms/src/bopds/ds.rs
git commit -m "feat(algorithms): compute UV boundary for curved DSFaces during DS loading"
```

---

### Task 6: Implement parametric face splitting

Replace the approximate `split_curved_face()` with parameter-space 2D polygon clipping.

**Files:**
- Modify: `libs/rcad-algorithms/src/builder.rs`

- [ ] **Step 1: Add the new split_curved_face_parametric method**

Add to `BooleanBuilder` impl in `builder.rs`:

```rust
use glam::DVec2;
use rcad_kernel::geom::Curve2dEval;

/// Split a curved face using parameter-space (u,v) 2D polygon clipping.
///
/// Algorithm:
/// 1. Get the face's UV boundary polygon
/// 2. For each intersection curve, get its PCurve on this face
/// 3. Sample the PCurve into a 2D polyline
/// 4. Split the UV polygon by the 2D trim polyline
/// 5. Map each sub-region back to 3D
fn split_curved_face_parametric(&self, face_idx: usize) -> Vec<SubFace> {
    let face = &self.ds.faces[face_idx];
    let surface = face.surface.clone();
    let normal = face.normal;

    // Get UV boundary
    let uv_boundary = match &face.uv_boundary {
        Some(uv) if uv.len() >= 3 => uv.clone(),
        _ => {
            // Fallback: no UV boundary, return whole face
            let boundary = face
                .boundary_verts
                .iter()
                .map(|&vi| self.ds.vertices[vi].point)
                .collect();
            return vec![SubFace { boundary, surface, normal }];
        }
    };

    // Collect PCurve trim lines for this face
    let mut trim_polylines_2d: Vec<Vec<DVec2>> = Vec::new();

    for &ci in &face.face_info.curves_in {
        let ic = &self.ds.intersection_curves[ci];

        // Determine which PCurve belongs to this face.
        // The IntersectionCurve stores pcurve_on_a (for face f1) and pcurve_on_b (for face f2).
        // We need to figure out if this face_idx was f1 or f2 in the interference.
        let pcurve = self.find_pcurve_for_face(ci, face_idx);

        if let Some(pc) = pcurve {
            // Sample the PCurve into 2D points
            let domain = pc.default_domain();
            let n = 32;
            let pts: Vec<DVec2> = (0..=n)
                .map(|i| {
                    let t = domain[0] + (domain[1] - domain[0]) * i as f64 / n as f64;
                    pc.point_at(t)
                })
                .collect();
            trim_polylines_2d.push(pts);
        } else {
            // No PCurve available — fallback to projecting the 3D polyline/curve
            let pts_3d = if !ic.polyline.is_empty() {
                ic.polyline.clone()
            } else {
                (0..=16)
                    .map(|i| {
                        let t = ic.t_range[0]
                            + (ic.t_range[1] - ic.t_range[0]) * i as f64 / 16.0;
                        use rcad_kernel::CurveEval;
                        ic.curve.point_at(t)
                    })
                    .collect()
            };
            let pts_2d: Vec<DVec2> = pts_3d
                .iter()
                .map(|&p| {
                    let proj =
                        rcad_kernel::projection::closest_point_on_surface(p, &surface);
                    DVec2::new(proj.params.0, proj.params.1)
                })
                .collect();
            if pts_2d.len() >= 2 {
                trim_polylines_2d.push(pts_2d);
            }
        }
    }

    if trim_polylines_2d.is_empty() {
        let boundary = face
            .boundary_verts
            .iter()
            .map(|&vi| self.ds.vertices[vi].point)
            .collect();
        return vec![SubFace { boundary, surface, normal }];
    }

    // Split the UV polygon by each trim polyline
    let mut uv_regions = vec![uv_boundary];

    for trim in &trim_polylines_2d {
        let mut next_regions = Vec::new();
        for region in &uv_regions {
            let split = split_uv_polygon_by_trim(region, trim);
            next_regions.extend(split);
        }
        uv_regions = next_regions;
    }

    // Map each UV sub-region back to 3D
    uv_regions
        .into_iter()
        .filter(|r| r.len() >= 3)
        .map(|uv_poly| {
            let boundary: Vec<DVec3> = uv_poly
                .iter()
                .map(|&uv| surface.point_at(uv.x, uv.y))
                .collect();
            SubFace {
                boundary,
                surface: surface.clone(),
                normal,
            }
        })
        .collect()
}

/// Find the PCurve that corresponds to a given face in an intersection curve.
fn find_pcurve_for_face(&self, curve_idx: usize, face_idx: usize) -> Option<&Curve2d> {
    // Search interferences to determine if this face was f1 or f2
    for intf in &self.ds.interferences {
        if let Interference::FaceFace { f1, f2, curves, .. } = intf {
            if curves.contains(&curve_idx) {
                let ic = &self.ds.intersection_curves[curve_idx];
                if *f1 == face_idx {
                    return ic.pcurve_on_a.as_ref();
                } else if *f2 == face_idx {
                    return ic.pcurve_on_b.as_ref();
                }
            }
        }
    }
    None
}
```

- [ ] **Step 2: Add the 2D polygon splitting helper**

Add as a free function in `builder.rs`:

```rust
/// Split a 2D polygon by a 2D trim polyline.
///
/// The trim polyline should cross the polygon boundary at 2 points,
/// splitting it into 2 sub-regions.
fn split_uv_polygon_by_trim(polygon: &[DVec2], trim: &[DVec2]) -> Vec<Vec<DVec2>> {
    let n = polygon.len();
    if n < 3 || trim.len() < 2 {
        return vec![polygon.to_vec()];
    }

    let trim_start = trim[0];
    let trim_end = *trim.last().unwrap();

    // Find the polygon edge indices closest to trim start and end
    let find_closest_edge = |target: DVec2| -> Option<(usize, DVec2)> {
        let mut best_dist = f64::INFINITY;
        let mut best_idx = 0;
        let mut best_pt = target;

        for i in 0..n {
            let j = (i + 1) % n;
            let a = polygon[i];
            let b = polygon[j];
            // Project target onto segment a-b
            let ab = b - a;
            let ab_len_sq = ab.length_squared();
            if ab_len_sq < 1e-20 {
                continue;
            }
            let t = ((target - a).dot(ab) / ab_len_sq).clamp(0.0, 1.0);
            let proj = a + ab * t;
            let dist = (target - proj).length();
            if dist < best_dist {
                best_dist = dist;
                best_idx = i;
                best_pt = proj;
            }
        }
        if best_dist < f64::INFINITY {
            Some((best_idx, best_pt))
        } else {
            None
        }
    };

    let (idx_start, pt_start) = match find_closest_edge(trim_start) {
        Some(v) => v,
        None => return vec![polygon.to_vec()],
    };
    let (idx_end, pt_end) = match find_closest_edge(trim_end) {
        Some(v) => v,
        None => return vec![polygon.to_vec()],
    };

    if idx_start == idx_end {
        return vec![polygon.to_vec()]; // Degenerate: both on same edge
    }

    // Ensure consistent ordering
    let (ia, ib, pa, pb, trim_fwd) = if idx_start <= idx_end {
        (idx_start, idx_end, pt_start, pt_end, true)
    } else {
        (idx_end, idx_start, pt_end, pt_start, false)
    };

    // Sub-region A: polygon[0..=ia] + pa + trim + pb + polygon[ib+1..]
    let mut sub_a: Vec<DVec2> = polygon[..=ia].to_vec();
    sub_a.push(pa);
    if trim_fwd {
        sub_a.extend_from_slice(&trim[1..trim.len() - 1]);
    } else {
        for &p in trim[1..trim.len() - 1].iter().rev() {
            sub_a.push(p);
        }
    }
    sub_a.push(pb);
    if ib + 1 < n {
        sub_a.extend_from_slice(&polygon[ib + 1..]);
    }

    // Sub-region B: pa + polygon[ia+1..=ib] + pb + reverse_trim
    let mut sub_b: Vec<DVec2> = vec![pa];
    sub_b.extend_from_slice(&polygon[ia + 1..=ib]);
    sub_b.push(pb);
    if trim_fwd {
        for &p in trim[1..trim.len() - 1].iter().rev() {
            sub_b.push(p);
        }
    } else {
        sub_b.extend_from_slice(&trim[1..trim.len() - 1]);
    }

    // Dedup near-coincident points
    let dedup_2d = |v: Vec<DVec2>| -> Vec<DVec2> {
        let mut r: Vec<DVec2> = Vec::new();
        for p in v {
            if r.is_empty() || (p - *r.last().unwrap()).length() > 1e-10 {
                r.push(p);
            }
        }
        if r.len() > 1 && (r[0] - *r.last().unwrap()).length() < 1e-10 {
            r.pop();
        }
        r
    };

    let sub_a = dedup_2d(sub_a);
    let sub_b = dedup_2d(sub_b);

    let mut result = Vec::new();
    if sub_a.len() >= 3 {
        result.push(sub_a);
    }
    if sub_b.len() >= 3 {
        result.push(sub_b);
    }
    if result.is_empty() {
        vec![polygon.to_vec()]
    } else {
        result
    }
}
```

- [ ] **Step 3: Update split_face dispatch to use the new method**

In `split_face()`, change the curved surface branch:

```rust
Surface3::Cylinder(_)
| Surface3::Sphere(_)
| Surface3::Cone(_)
| Surface3::Torus(_) => self.split_curved_face_parametric(face_idx),
```

Keep the old `split_curved_face` method as `_split_curved_face_legacy` for reference but don't call it. Or remove it entirely.

- [ ] **Step 4: Run existing tests**

Run: `cargo test --workspace`
Expected: all existing box-box tests still pass (planar path unchanged).

- [ ] **Step 5: Commit**

```bash
git add libs/rcad-algorithms/src/builder.rs
git commit -m "feat(algorithms): replace split_curved_face with parameter-space 2D clipping"
```

---

### Task 7: Add curved Boolean operation tests

Comprehensive tests for all target surface pair combinations.

**Files:**
- Modify: `libs/rcad-algorithms/src/lib.rs`

- [ ] **Step 1: Add Box × Sphere tests**

Add to the test module in `libs/rcad-algorithms/src/lib.rs`:

```rust
#[test]
fn boolean_box_sphere_intersection() {
    use rcad_modeling::{make_box_brep, make_sphere_brep};

    let a = make_box_brep(
        DVec3::ZERO,
        DVec3::X,
        DVec3::Y,
        2.0, 2.0, 2.0,
    ).unwrap();
    let b = make_sphere_brep(DVec3::new(1.0, 1.0, 1.0), 1.5).unwrap();

    let result = boolean_op(BooleanOpType::Intersection, &a, &b);
    assert!(result.is_ok(), "box-sphere intersection should succeed");
    let brep = result.unwrap();
    assert!(!brep.solids[0].shells[0].faces.is_empty(), "result should have faces");

    // Volume sanity: intersection must be smaller than both inputs
    let v_result = rcad_kernel::properties::volume(&brep);
    let v_box = rcad_kernel::properties::volume(&a);
    let v_sphere = rcad_kernel::properties::volume(&b);
    assert!(v_result > 0.0, "result volume should be positive");
    assert!(v_result < v_box, "intersection should be smaller than box");
    assert!(v_result < v_sphere, "intersection should be smaller than sphere");
}

#[test]
fn boolean_box_sphere_difference() {
    use rcad_modeling::{make_box_brep, make_sphere_brep};

    let a = make_box_brep(
        DVec3::ZERO,
        DVec3::X,
        DVec3::Y,
        2.0, 2.0, 2.0,
    ).unwrap();
    let b = make_sphere_brep(DVec3::new(1.0, 1.0, 1.0), 0.5).unwrap();

    let result = boolean_op(BooleanOpType::Difference, &a, &b);
    assert!(result.is_ok(), "box-sphere difference should succeed");
    let brep = result.unwrap();

    let v_result = rcad_kernel::properties::volume(&brep);
    let v_box = rcad_kernel::properties::volume(&a);
    let v_sphere = rcad_kernel::properties::volume(&b);
    // Box minus sphere should be approximately v_box - v_sphere
    let expected = v_box - v_sphere;
    let error = ((v_result - expected) / expected).abs();
    assert!(error < 0.15, "volume error {error:.2} should be < 15%");
}

#[test]
fn boolean_box_sphere_union() {
    use rcad_modeling::{make_box_brep, make_sphere_brep};

    let a = make_box_brep(
        DVec3::ZERO,
        DVec3::X,
        DVec3::Y,
        2.0, 2.0, 2.0,
    ).unwrap();
    let b = make_sphere_brep(DVec3::new(1.0, 1.0, 2.0), 1.0).unwrap();

    let result = boolean_op(BooleanOpType::Union, &a, &b);
    assert!(result.is_ok(), "box-sphere union should succeed");
    let brep = result.unwrap();
    let v_result = rcad_kernel::properties::volume(&brep);
    let v_box = rcad_kernel::properties::volume(&a);
    assert!(v_result > v_box, "union should be larger than box");
}
```

- [ ] **Step 2: Add Sphere × Sphere tests**

```rust
#[test]
fn boolean_sphere_sphere_intersection() {
    use rcad_modeling::make_sphere_brep;

    let a = make_sphere_brep(DVec3::ZERO, 1.0).unwrap();
    let b = make_sphere_brep(DVec3::new(1.0, 0.0, 0.0), 1.0).unwrap();

    let result = boolean_op(BooleanOpType::Intersection, &a, &b);
    assert!(result.is_ok(), "sphere-sphere intersection should succeed");
    let brep = result.unwrap();
    let v = rcad_kernel::properties::volume(&brep);
    // Lens-shaped intersection of two unit spheres separated by 1.0
    // Analytical: V = (5π/12) ≈ 1.309
    assert!(v > 0.5 && v < 2.0, "volume {v} should be in reasonable range");
}

#[test]
fn boolean_sphere_sphere_difference() {
    use rcad_modeling::make_sphere_brep;

    let a = make_sphere_brep(DVec3::ZERO, 2.0).unwrap();
    let b = make_sphere_brep(DVec3::new(1.0, 0.0, 0.0), 1.0).unwrap();

    let result = boolean_op(BooleanOpType::Difference, &a, &b);
    assert!(result.is_ok(), "sphere-sphere difference should succeed");
    let brep = result.unwrap();
    let v = rcad_kernel::properties::volume(&brep);
    let v_a = rcad_kernel::properties::volume(&a);
    assert!(v > 0.0 && v < v_a, "difference volume should be positive and less than A");
}
```

- [ ] **Step 3: Add Box × Cylinder test (through-hole)**

```rust
#[test]
fn boolean_box_cylinder_hole() {
    use rcad_modeling::{make_box_brep, make_cylinder_brep};

    // Box centered at origin
    let a = make_box_brep(
        DVec3::new(-1.0, -1.0, -1.0),
        DVec3::X,
        DVec3::Y,
        2.0, 2.0, 2.0,
    ).unwrap();
    // Cylinder along Z axis, smaller than box
    let b = make_cylinder_brep(
        DVec3::new(0.0, 0.0, -2.0),
        DVec3::Z,
        DVec3::X,
        0.5,
        4.0,
    ).unwrap();

    let result = boolean_op(BooleanOpType::Difference, &a, &b);
    assert!(result.is_ok(), "box-cylinder difference should succeed");
    let brep = result.unwrap();
    let v = rcad_kernel::properties::volume(&brep);
    let v_box = 8.0; // 2x2x2
    let v_cyl_hole = std::f64::consts::PI * 0.25 * 2.0; // π*r²*h through the box height
    let expected = v_box - v_cyl_hole;
    let error = ((v - expected) / expected).abs();
    assert!(error < 0.15, "volume error {error:.2} should be < 15%");
}
```

- [ ] **Step 4: Add Cylinder × Cylinder test (Steinmetz)**

```rust
#[test]
fn boolean_cylinder_cylinder_intersection() {
    use rcad_modeling::make_cylinder_brep;

    // Two perpendicular cylinders
    let a = make_cylinder_brep(
        DVec3::new(0.0, 0.0, -2.0),
        DVec3::Z,
        DVec3::X,
        1.0,
        4.0,
    ).unwrap();
    let b = make_cylinder_brep(
        DVec3::new(-2.0, 0.0, 0.0),
        DVec3::X,
        DVec3::Y,
        1.0,
        4.0,
    ).unwrap();

    let result = boolean_op(BooleanOpType::Intersection, &a, &b);
    assert!(result.is_ok(), "cylinder-cylinder intersection should succeed");
    let brep = result.unwrap();
    let v = rcad_kernel::properties::volume(&brep);
    // Steinmetz solid volume = 16r³/3 ≈ 5.333 for r=1
    assert!(v > 2.0 && v < 8.0, "Steinmetz volume {v} should be reasonable");
}
```

- [ ] **Step 5: Run all new tests**

Run: `cargo test -p rcad-algorithms boolean_box_sphere boolean_sphere_sphere boolean_box_cylinder boolean_cylinder_cylinder -- --nocapture`
Expected: all new tests PASS.

- [ ] **Step 6: Run full test suite to verify no regressions**

Run: `cargo test --workspace`
Expected: all tests pass (old + new).

- [ ] **Step 7: Run clippy and fmt**

```bash
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
```

Expected: clean.

- [ ] **Step 8: Commit**

```bash
git add libs/rcad-algorithms/src/lib.rs
git commit -m "test(algorithms): add curved Boolean operation tests for box-sphere, sphere-sphere, box-cylinder, cylinder-cylinder"
```

---

### Task 8: Volume conservation validation tests

Add tests that verify the Boolean identity: `V(A) + V(B) - V(A∩B) ≈ V(A∪B)`.

**Files:**
- Modify: `libs/rcad-algorithms/src/lib.rs`

- [ ] **Step 1: Write volume conservation test**

```rust
#[test]
fn boolean_volume_conservation_box_sphere() {
    use rcad_modeling::{make_box_brep, make_sphere_brep};

    let a = make_box_brep(
        DVec3::ZERO,
        DVec3::X,
        DVec3::Y,
        2.0, 2.0, 2.0,
    ).unwrap();
    let b = make_sphere_brep(DVec3::new(1.0, 1.0, 1.0), 1.0).unwrap();

    let v_a = rcad_kernel::properties::volume(&a);
    let v_b = rcad_kernel::properties::volume(&b);

    let union = boolean_op(BooleanOpType::Union, &a, &b).unwrap();
    let inter = boolean_op(BooleanOpType::Intersection, &a, &b).unwrap();

    let v_union = rcad_kernel::properties::volume(&union);
    let v_inter = rcad_kernel::properties::volume(&inter);

    // V(A∪B) = V(A) + V(B) - V(A∩B)
    let expected_union = v_a + v_b - v_inter;
    let error = ((v_union - expected_union) / expected_union).abs();
    assert!(
        error < 0.05,
        "Volume conservation failed: V(A∪B)={v_union:.4}, V(A)+V(B)-V(A∩B)={expected_union:.4}, error={error:.4}"
    );
}

#[test]
fn boolean_volume_conservation_spheres() {
    use rcad_modeling::make_sphere_brep;

    let a = make_sphere_brep(DVec3::ZERO, 1.5).unwrap();
    let b = make_sphere_brep(DVec3::new(1.0, 0.0, 0.0), 1.0).unwrap();

    let v_a = rcad_kernel::properties::volume(&a);
    let v_b = rcad_kernel::properties::volume(&b);

    let union = boolean_op(BooleanOpType::Union, &a, &b).unwrap();
    let inter = boolean_op(BooleanOpType::Intersection, &a, &b).unwrap();

    let v_union = rcad_kernel::properties::volume(&union);
    let v_inter = rcad_kernel::properties::volume(&inter);

    let expected_union = v_a + v_b - v_inter;
    let error = ((v_union - expected_union) / expected_union).abs();
    assert!(
        error < 0.05,
        "Volume conservation failed: V(A∪B)={v_union:.4}, expected={expected_union:.4}, error={error:.4}"
    );
}
```

- [ ] **Step 2: Run tests**

Run: `cargo test -p rcad-algorithms boolean_volume_conservation -- --nocapture`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add libs/rcad-algorithms/src/lib.rs
git commit -m "test(algorithms): add volume conservation validation for curved Boolean ops"
```
