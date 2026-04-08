# Benchmark System Design

**Date:** 2026-04-08
**Status:** Approved
**Scope:** Add criterion.rs benchmarks for core algorithm paths across 3 crates

## Overview

Add performance benchmarks using criterion.rs to establish baselines and detect regressions for critical code paths: Boolean operations, surface-surface intersection, properties computation, point projection, and STEP I/O.

## Organization

Each benchmarked crate has its own `benches/` directory with criterion as a workspace dev-dependency. This follows standard Rust convention and allows `cargo bench -p <crate>` for targeted runs.

## Benchmarks

### rcad-algorithms — `libs/rcad-algorithms/benches/boolean.rs`

| Benchmark | Input | Measures |
|-----------|-------|----------|
| `boolean_union_boxes` | Two overlapping 1×1×1 boxes | Full Boolean union pipeline (DS + PaveFiller + Builder) |
| `boolean_diff_box_sphere` | Box(2×2×2) - Sphere(r=0.8) | Curved Boolean with parametric face splitting |
| `intss_plane_sphere` | Plane × Sphere (r=1) | Analytic IntSS dispatch + PCurve derivation |

### rcad-kernel — `libs/rcad-kernel/benches/properties.rs`

| Benchmark | Input | Measures |
|-----------|-------|----------|
| `volume_sphere` | Sphere BRep (r=1) | Volume via divergence theorem |
| `closest_point_on_sphere` | 100 random points → Sphere(r=1) | Analytic point projection throughput |

### rcad-step — `libs/rcad-step/benches/step_io.rs`

| Benchmark | Input | Measures |
|-----------|-------|----------|
| `step_roundtrip_box` | Box BRep → STEP string → parse back | Serialization + deserialization |

## Configuration

- `criterion = { version = "0.5" }` added to `[workspace.dev-dependencies]` in root Cargo.toml
- Each crate's Cargo.toml gets `criterion.workspace = true` under `[dev-dependencies]` and a `[[bench]]` section with `harness = false`
- CI: benchmarks NOT added to CI workflow (too slow for PR checks). Can be added later as a separate optional job.

## Decisions

| Decision | Rationale |
|----------|-----------|
| criterion over built-in bench | Stable Rust, statistical analysis, HTML reports |
| Per-crate benches/ | Standard convention, allows targeted `cargo bench -p` |
| 6 benchmarks total | Covers critical paths without excessive noise |
| No CI integration yet | Benchmark variance on shared runners makes it unreliable for regression detection |
