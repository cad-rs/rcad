# CI/CD Pipeline Design for rcad

**Date:** 2026-04-08
**Status:** Approved
**Scope:** GitHub Actions CI pipeline — first improvement item in the rcad quality roadmap

## Overview

Add a single GitHub Actions workflow (`.github/workflows/ci.yml`) that runs on every push to `main` and every pull request targeting `main`. The workflow contains 4 parallel jobs covering linting, testing, WASM build verification, and example compilation.

## Trigger

```yaml
on:
  push:
    branches: [main]
  pull_request:
    branches: [main]
```

## Jobs

### 1. lint (ubuntu-latest)

Purpose: Enforce code formatting and static analysis.

Steps:
1. `actions/checkout@v4`
2. `dtolnay/rust-toolchain@stable` with components: `rustfmt, clippy`
3. `Swatinem/rust-cache@v2`
4. `cargo fmt --all --check`
5. `cargo clippy --workspace --all-targets -- -D warnings`

Clippy treats all warnings as errors (`-D warnings`). No custom `.clippy.toml` — use default rules.

### 2. test (matrix: ubuntu-latest, macos-latest, windows-latest)

Purpose: Run the full test suite across platforms.

Steps:
1. `actions/checkout@v4`
2. `dtolnay/rust-toolchain@stable`
3. `Swatinem/rust-cache@v2`
4. `cargo test --workspace`

Three platforms run in parallel to catch platform-specific issues (path separators, floating-point behavior, etc.).

### 3. wasm (ubuntu-latest)

Purpose: Verify that both creator apps compile to WebAssembly.

Steps:
1. `actions/checkout@v4`
2. `dtolnay/rust-toolchain@stable` with targets: `wasm32-unknown-unknown`
3. `Swatinem/rust-cache@v2`
4. Install trunk via `jetli/trunk-action@v0.1.0` (pre-built binary, avoids slow `cargo install`)
5. `trunk build` in `apps/creator-egui`
6. `trunk build` in `apps/creator-iced`

Note: wasm-bindgen is pinned to `=0.2.100` in the workspace. Trunk version must be compatible with this pin.

### 4. examples (ubuntu-latest)

Purpose: Ensure all 21+ examples compile successfully.

Steps:
1. `actions/checkout@v4`
2. `dtolnay/rust-toolchain@stable`
3. `Swatinem/rust-cache@v2`
4. `cargo build --examples`

Compile-only — examples are not executed in CI (some produce STEP output files).

## Caching Strategy

All jobs use `Swatinem/rust-cache@v2`, which keys on `Cargo.lock` by default. The WASM job caches separately due to different compilation target (`wasm32-unknown-unknown`).

## Decisions

| Decision | Rationale |
|----------|-----------|
| Single workflow file | Sufficient for this project size; avoids config duplication across multiple files |
| Parallel jobs (not serial) | Faster feedback; failures in one area don't block diagnosis of others |
| `-D warnings` on clippy | Enforces clean code; can be relaxed per-lint if needed |
| `trunk build` not `cargo build --target wasm32` | Trunk is the actual build tool used; validates the real build path including index.html |
| `jetli/trunk-action` over `cargo install trunk` | Pre-built binary saves ~2-3 min compilation time |
| Examples: build-only | Running examples would generate files and add complexity with no testing value |

## Future Extensions

These are explicitly out of scope for this iteration but noted for later:
- **Benchmarks job** — add when criterion benchmarks are introduced
- **Coverage job** — add with tarpaulin or llvm-cov when coverage targets are set
- **Release workflow** — add when crate publishing is needed
