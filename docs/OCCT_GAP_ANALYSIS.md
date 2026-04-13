# RCAD vs OCCT Gap Analysis

Date: 2026-04-12

## Purpose

This document turns the current RCAD capability inventory into an OCCT-oriented gap analysis and an execution roadmap. The goal is not to chase feature-count parity blindly, but to identify the areas where RCAD still differs most from OCCT in practical engineering use.

The key conclusion is:

RCAD already covers a large portion of the geometric foundation layer, but still trails OCCT most in kernel maturity layers:

- robust boolean architecture
- healing and tolerance management
- richer topology and document model
- industrial-grade data exchange depth
- advanced post-processing and result simplification

## Executive Summary

### Areas where RCAD is already strong

- Analytic and spline geometry coverage is broad.
- Core B-Rep entities are in place.
- Primitive creation, sweeps, loft, fillet, chamfer, thicken, draft, mirror, and array operations exist.
- STEP read/write, assembly IO, color IO, OBJ IO, and an IGES mesh bridge exist.
- HLR, projections, curvature, extrema, sectioning, and shape properties are implemented.

### Areas where RCAD is still behind OCCT

- Boolean operations are not yet at OCCT's robustness and breadth, especially on curved solids and mixed-dimension workflows.
- Healing and tolerance propagation are still basic compared with ShapeFix / ShapeAnalysis / ShapeUpgrade style workflows.
- Topology containers are narrower than OCCT's full shape graph, especially around Compound and CompSolid workflows.
- Post-operation simplification and same-domain unification are still limited.
- Exchange support is good for STEP AP203/AP214, but not yet at OCCT's full document and interoperability depth.

## Capability Matrix

| Domain | RCAD today | Relative to OCCT | Gap level | Recommended direction |
|---|---|---|---|---|
| Geometry kernel | Strong analytic + spline coverage | Close on base geometry | Low | Maintain correctness and performance |
| Core B-Rep model | Vertex/Edge/Wire/Face/Shell/Solid | Behind on richer shape container graph | Medium | Add Compound, CompSolid, better non-manifold support |
| Primitive and sweep modeling | Good breadth | Competitive on basic features | Low | Improve edge cases and history consistency |
| Fillet/chamfer/thicken/draft | Present, but narrower robustness envelope | Behind on maturity | Medium | Strengthen corner cases, variable laws, self-intersection handling |
| Boolean framework | Functional, especially planar; curved partial | Significantly behind | High | Build OCCT-style split-first boolean framework and advanced options |
| Healing and validation | Basic repair + check available | Significantly behind | High | Expand into full healing/tolerance pipeline |
| Post-op simplification | Limited | Behind | High | Add unify-same-domain, small-edge cleanup, tolerance correction |
| Topology history | Face-oriented history exists | Behind | Medium | Extend to Deleted/Modified/Generated across vertices, edges, faces, solids |
| STEP exchange | Good AP203/AP214 support | Solid but not full parity | Medium | Add AP242 depth, richer metadata, stronger healing on import |
| IGES exchange | Mesh bridge only | Far behind | High | Add analytic/B-Rep IGES support or de-prioritize explicitly |
| Assembly/document model | Assembly tree + colors exist | Behind XCAF/OCAF depth | Medium | Add names, layers, materials, persistent document metadata |
| Meshing and visualization | Good rendering baseline | Behind industrial meshing controls | Medium | Add tunable meshing and stronger HLR accuracy |

## Highest-Priority Gaps

## 1. Robust Boolean Architecture

This is the largest practical gap.

RCAD already has a DS / PaveFiller / Builder shape, but OCCT's advantage is not just that it supports Fuse/Common/Cut. Its real strength is the broader boolean platform:

- General Fuse as the common split engine
- Splitter as a first-class operation
- Cells Builder for reusable split parts
- MakerVolume for volume construction from split shells/faces
- advanced options such as fuzzy tolerance, gluing, non-destructive mode, inverted-solid checks, and OBB acceleration
- result simplification and same-domain unification
- rich history information

### What RCAD should add

- General-fuse style split-first core separate from final boolean classification.
- Splitter operation over arbitrary objects/tools groups.
- Cells-builder style reusable split-cell graph for repeated boolean expressions.
- Volume-maker workflow for constructing solids from split faces/shells.
- Advanced options layer:
  - fuzzy tolerance
  - gluing mode
  - non-destructive mode
  - check-inverted toggle
  - OBB pruning in addition to current acceleration
- Result simplification stage:
  - unify same-domain faces
  - unify tangent edges/faces where valid
  - remove internal faces after fuse
  - small-edge cleanup and tolerance correction
- Full topological history for vertices, edges, faces, and solids.

### Why this matters

Without this layer, RCAD can demonstrate many boolean examples but will continue to lag on real imported models, near-coincident geometry, repeated boolean editing, and downstream naming/history use.

## 2. Healing and Tolerance System

RCAD already stores per-vertex, per-edge, and per-face tolerances. That is the right foundation, but the surrounding workflow is still much thinner than OCCT.

### Current RCAD baseline

- close-vertex merge
- degenerate-face removal
- face normal recomputation
- wire orientation repair
- basic validity checks

### What is still missing

- SameParameter repair for edge 3D curve vs PCurve consistency
- SameRange repair where imported ranges disagree
- edge-on-surface consistency checks
- seam-edge and degenerated-edge repair workflows
- shell closure and solid closure validation
- self-interference detection beyond basic topology checks
- tolerance refinement and propagation after split, sew, and boolean operations
- import-heal-export workflows for external CAD data
- stronger analyzers for imported bad geometry

### Why this matters

OCCT survives many poor inputs not because the geometry is simpler, but because its analyzer/healer stack is deep. If RCAD aims to replace OCCT in CAE preprocessing or CAD interoperability, this area is mandatory.

## 3. Richer Topology and Document Model

RCAD's shape model is good for solid modeling, but narrower than OCCT's full shape/document ecosystem.

### Gaps

- Compound support
- CompSolid support
- more explicit open-shell and mixed-dimension container workflows
- stronger non-manifold topology handling
- persistent naming infrastructure beyond lightweight history
- document metadata comparable to XCAF/OCAF usage patterns:
  - names
  - layers
  - materials
  - attributes
  - instance references
  - richer assembly metadata

### Why this matters

Once workflows go beyond single modeled solids and into imported assemblies, CAE preparation, reusable feature editing, or multi-object operations, this gap becomes visible quickly.

## 4. Exchange Depth and Interoperability

STEP support is already one of RCAD's better areas, but OCCT still leads in depth, metadata coverage, and recovery quality.

### Recommended improvements

- strengthen STEP AP242 support
- preserve more metadata through import/export
- add import healing as part of default read pipelines
- improve validation of exported PCurves and tolerances
- decide strategically whether to expand IGES beyond the current mesh bridge

If IGES B-Rep support is not a priority, document that explicitly and focus on STEP AP242 instead.

## 5. Meshing, HLR, and Result Quality

RCAD has working rendering and HLR, but OCCT still has deeper production controls around tessellation and engineering drawing output.

### Recommended improvements

- configurable meshing deflection and angular tolerances
- incremental or cached tessellation for edited models
- post-boolean remeshing quality improvements
- stronger analytic edge preservation in HLR output
- more exact silhouette/hidden-line handling on curved surfaces

## Roadmap

## P0: Make RCAD Credible as an OCCT-Class Boolean Kernel

Target outcome: imported and modeled solids go through boolean workflows with much better robustness and cleaner results.

### Deliverables

- Introduce a reusable general-fuse core in `rcad-algorithms`.
- Add Splitter API for object/tool split workflows.
- Add fuzzy tolerance option to boolean, section, and split operations.
- Add non-destructive processing mode.
- Add result simplification pass:
  - same-domain face unification
  - internal face removal for fuse
  - small-edge cleanup
- Expand history tracking to edges and solids.
- Close current curved-boolean gaps called out in production-readiness notes.

### Acceptance criteria

- Curved-solid boolean tests are no longer marked partial in project capability docs.
- Near-coincident boolean regression cases become deterministic with fuzzy mode.
- Boolean results expose Deleted/Modified/Generated mappings for all topological dimensions that RCAD supports.

## P1: Build Industrial Healing and Import Robustness

Target outcome: imported CAD data can be analyzed, repaired, and pushed into modeling/boolean workflows more safely.

### Deliverables

- Add ShapeFix-style healing module set.
- Add SameParameter and SameRange repair utilities.
- Add shell closure, manifoldness, and face-on-surface analyzers.
- Add import pipeline hooks:
  - parse
  - analyze
  - optional heal
  - report
- Add stronger tolerance propagation rules after sew/split/boolean.

### Acceptance criteria

- STEP import of imperfect but recoverable models succeeds through an analyze/heal pipeline.
- Healing reports expose concrete diagnostics instead of only pass/fail.
- `repair` evolves from a convenience helper into a structured healing workflow.

## P2: Expand Topology and Document Model

Target outcome: RCAD can support more OCCT-like assembly, compound, and document-centric workflows.

### Deliverables

- Add Compound and CompSolid topology containers.
- Introduce richer document metadata model for names, layers, materials, and attributes.
- Add persistent naming strategy tied to operation history.
- Expand assembly graph semantics for external references and instance attributes.

### Acceptance criteria

- Multi-body and mixed container workflows can be represented without flattening into ad hoc solids.
- Export/import preserves more semantic metadata.
- Application code can query stable object identity beyond transient face indices.

## Module-Level Task Breakdown

## `libs/rcad-algorithms`

### Highest value work

- General fuse core
- Splitter API
- fuzzy / glue / non-destructive boolean options
- result simplification and same-domain unification
- better curved face-face analytic intersection coverage
- stronger post-op PCurve rebuild
- richer history graph

## `libs/rcad-kernel`

### Highest value work

- Compound / CompSolid support
- richer validity analysis
- tolerance propagation rules
- persistent naming hooks
- stronger topology graph services

## `libs/rcad-modeling`

### Highest value work

- stabilize advanced fillet/chamfer/thicken edge cases
- preserve analytic geometry and history through more feature edits
- integrate healing after risky operations when needed

## `libs/rcad-step`

### Highest value work

- AP242-oriented improvements
- richer metadata preservation
- import diagnostics and optional heal pipeline integration
- stronger validation of generated PCurves / surface bindings / tolerances

## `libs/rcad-scene`

### Highest value work

- consume richer topology history and persistent naming
- prepare object/document layer for semantic metadata instead of only geometry payloads

## `libs/rcad-render`

### Highest value work

- tunable meshing controls
- improved HLR quality on curved geometry
- better cache invalidation for edited models

## Suggested Milestone Order

1. Strengthen boolean core before adding more feature breadth.
2. Add healing/tolerance workflows immediately after boolean core stabilization.
3. Expand topology/document model once history and healing are credible.
4. Deepen exchange support after the kernel can recover and preserve more semantics.

## Recommended Non-Goals for Now

These are valid future areas, but should not outrank the kernel-maturity work above.

- chasing many new geometry primitive types
- expanding UI features first
- adding niche interchange formats before AP242/healing maturity
- premature optimization of rendering before boolean/healing/document layers improve

## Bottom Line

RCAD is no longer missing the obvious basic geometry layer. The next stage is to become harder to break.

To move materially closer to OCCT, the project should prioritize:

1. boolean robustness and reusable split architecture
2. healing and tolerance management
3. richer topology and document semantics

If these three areas improve, RCAD will start to feel less like a promising geometric engine and more like a credible industrial CAD kernel platform.