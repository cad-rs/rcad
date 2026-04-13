# RCAD vs OCCT Gap Analysis

Date: 2026-04-13 (revised against OCCT 7.9.3 stable / 8.0.0-rc5, incorporating P3–P7 partial completions)

## Purpose

This document turns the current RCAD capability inventory into an OCCT-oriented gap analysis and an execution roadmap. The goal is not to chase feature-count parity blindly, but to identify the areas where RCAD still differs most from OCCT in practical engineering use.

The key conclusion is:

RCAD has made substantial progress since the last revision. P3–P7 partial work has added small-edge cleanup, SameParameter/SameRange diagnosis+repair, shell and wire diagnostics, draft prism + revolution features, STEP material/layer/general-property metadata plumbing, and new helix/involute/spiral curve evaluators. However, OCCT 8.0.0 (currently at rc5, release imminent) has itself leapt forward significantly -- it introduces an entirely new graph-based topology representation (BRepGraph, 49 000+ lines), a Gordon surface transfinite interpolation framework, a fully redesigned geometry evaluation architecture, the TKHelix toolkit, and defeaturing + connected-shape APIs. The production gap is narrowing but OCCT's target is moving.

RCAD still trails OCCT most in:

- boolean robustness and architectural breadth (largest gap, unchanged)
- topology representation depth: OCCT 8.0.0 adds BRepGraph as a second, richer topology layer
- healing and tolerance management (partially improved in P3–P4)
- industrial-grade data exchange depth
- geometry evaluation breadth: helicoids, spirals, ellipsoids, Gordon surfaces

## OCCT Version Reference

| Metric | Value |
|---|---|
| Latest stable | **7.9.3** (December 2025) |
| Master pre-release | **8.0.0-rc5** (April 2026, ~60–70 days from final release) |
| 8.0.0 total changes vs 7.9.0 | 460+ improvements and bug fixes across all rc stages |
| New toolkits in 8.0.0 | **TKHelix** (helix geometry), **BRepGraph** API (graph topology), **GeomBndLib** (geometry-aware bboxes) |
| Architectural changes | Unified TKDE data exchange, `TopoDS_TShape` contiguous-array overhaul, `BRepGraph` graph BRep, devirtualized EvalD\* geometry, C++17 minimum, `std::exception` inheritance, thread-local error handlers |
| New algorithm classes | `BRepAlgoAPI_Defeaturing`, `BOPAlgo_MakeConnected`, `GeomFill_Gordon`, `ExtremaPC`, `PointSetLib_Props/Equation` |
| Handle API | Return-by-value APIs, deprecated out-parameter overloads, removed `TColGeom`/`TColGeom2d` packages |
| Toolkit count | ~60 toolkits across 7 modules (3 new in 8.0.0) |

## Executive Summary

### Areas where RCAD is already strong (updated after P1.4–P6)

- Broad analytic and B-spline geometry coverage.
- Core B-Rep types: Vertex, Edge, Wire, Face, Shell, Solid.
- Primitive creation, extrude, revolve, sweep multi-section, loft, variable-radius fillet, chamfer, thicken, draft, mirror, and array.
- STEP AP214/AP242 write, assembly IO, color/material/layer attributes, OBJ IO, IGES mesh bridge.
- HLR with dense silhouette sampling on curved surfaces (P2.2).
- Analytic extrema/distance fast paths for sphere-sphere, plane-sphere, parallel planes (P2.3).
- Curvature, sectioning, shape properties, projections.
- mesh_dirty incremental mesh caching (P2.1).
- **P3**: Small-edge / degenerate-edge cleanup integrated into `simplify_brep_post_ops`.
- **P4**: SameParameter and SameRange diagnosis+repair; shell manifold / open-edge analyzer (`analyze_shell_topology`); wire gap/self-intersection report (`analyze_wire_issues`).
- **P5**: Prism, draft prism, and revolution feature operations.
- **P6**: STEP material/layer extraction, GENERAL_PROPERTY extraction and export API.
- **P7 (partial)**: Circular helix, circle involute, Archimedean spiral, logarithmic spiral evaluators.

### Areas where RCAD is still behind OCCT

- **Boolean robustness** on curved solids and near-coincident geometry remains the largest practical gap. Splitter, CellsBuilder, fuzzy tolerance, defeaturing, and MakeConnected are all absent.
- **BRepGraph**: OCCT 8.0.0 ships an entirely new graph-based BRep API (49 000+ lines) with history tracking, mutation guards, deduplication, and validation. RCAD has no equivalent second-tier topology layer.
- **Healing pipeline**: ShapeUpgrade_UnifySameDomain (same-domain face merging) and ShapeProcess (batch repair chain) are still missing; these dominate post-boolean result quality.
- **Exchange**: STEP AP242 write is present but AP242 complete read, GDT/DimTol, kinematics, and FEA entities are not.
- **Geometry evaluation breadth**: OCCT 8.0.0 adds helicoid, spiral, ellipsoid, and parametric curve evaluation classes not present in RCAD.
- **Gordon surface**: N×M transfinite interpolation of curve networks (`GeomFill_Gordon`) is new in OCCT 8.0.0; RCAD supports only 4-boundary Coons patches.
- **Topology containers**: No Compound / CompSolid; no non-manifold topology.
- **Feature library**: Prism, draft prism, revolution feature, and cylindrical hole are present; rib/slot and SplitShape remain.

## Capability Matrix

| Domain | RCAD today | OCCT reference | Gap level | Recommended direction |
|---|---|---|---|---|
| Geometry kernel | Strong analytic + spline coverage | TKG2d / TKG3d / TKGeomBase | Low | Maintain correctness; helix curves now in OCCT TKHelix |
| Geometry evaluation breadth | Basic D0/D1 on standard types | GeomEval/Geom2dEval (helicoids, spirals, ellipsoids, parametric) | Medium | Add helicoid, spiral, ellipsoid eval classes |
| Gordon / N×M surface fill | Coons 4-boundary only | GeomFill_Gordon (8.0) | Medium | Implement transfinite N×M interpolation |
| Point cloud analysis | Missing | PointSetLib (8.0) -- PCA, inertia, dimensionality | Low | Low priority; add if needed for import QC |
| Core B-Rep model | Vertex/Edge/Wire/Face/Shell/Solid | TKBRep | Medium | Add Compound, CompSolid, non-manifold support |
| Graph-based topology layer | None | BRepGraph (8.0, 49k lines) -- history, mutation, dedup, validate | High | Long-term: design graph topology API for history + editing |
| Primitive and sweep modeling | Good breadth (extrude/revolve/loft/sweep) | TKPrim + TKOffset | Low | Edge cases; N-sided fill, normal projection |
| Fillet / chamfer | Variable-radius fillet, chamfer | TKFillet | Medium | Angle-mode chamfer, 2-D fillet API, corner cases |
| Thicken / draft / offset | Present | TKOffset BRepOffset | Medium | Shell offset (MakeOffsetShape), evolved surface |
| Feature library | Prism + draft prism + revolution + cylindrical hole | TKFeat (boss/pocket/rib/hole) | Medium | Rib/slot, SplitShape |
| Boolean framework | Fuse/Cut/Common/Section + imprint | TKBO | High | CellsBuilder, Splitter, fuzzy, glue, defeaturing |
| Post-op simplification | Small-edge cleanup (P3) | ShapeUpgrade_UnifySameDomain, BOPAlgo cleanup | High | Same-domain unification (highest impact), internal-face removal |
| Healing and validation | SameParameter/SameRange + shell + wire diagnostics (P4 partial) | TKShHealing (10 packages) | High | ShapeUpgrade_UnifySameDomain, ShapeProcess, tolerance rules |
| Topology history | Face-level history | BOPAlgo history, BRepGraph_History (8.0), OCAF naming | Medium | Extend to edges/vertices; persistent naming; BRepGraph history |
| STEP exchange: write | AP214 + AP242 + material/layer + GENERAL_PROPERTY | TKDESTEP STEPCAFControl | Medium | GDT write, property_definition relations, PCurve validation |
| STEP exchange: read | Basic import | TKDESTEP + STEPCAFControl_Reader | High | Full AP242 read, GDT, kinematics, FEA entities |
| IGES exchange | Mesh bridge only | TKDEIGES | High | Add analytic/B-Rep IGES or document as non-goal |
| Assembly / document model | Colors + shape tree + material + layer (P6 partial) | TKXCAF (XCAFDoc_*) | Medium | DimTol, GDT annotations, notes, persistent attributes |
| Meshing and visualization | mesh_dirty caching (P2.1), HLR dense silhouettes (P2.2) | TKMesh + TKHLR | Medium | Tunable deflection/angular tolerances, incremental remesh |
| Thread safety | Single-threaded | BRepCheck thread-safe (8.0), thread-local error handlers (8.0) | Low | Not urgent unless parallel workflows are added |

## Highest-Priority Gaps

## 0. BRepGraph: Graph-Based Topology Layer (NEW in OCCT 8.0.0 -- gap: Medium-High)

OCCT 8.0.0-rc5 introduces `BRepGraph`, an entirely new graph-based representation of topology and BRep geometry as an alternative to the traditional `TopoDS_Shape` linked structure. This is 49 000+ lines of new code with 20+ GTest files.

Key capabilities of BRepGraph that RCAD does not have:

| BRepGraph feature | Purpose |
|---|---|
| `BRepGraph_NodeId`-typed incidence tables | Graph traversal without pointer chasing |
| `BRepGraph_History` | Persistent shape history (which new faces came from which old faces) |
| `BRepGraph_MutGuard` | Safe mutation of topology with invariant checking |
| `BRepGraph_Deduplicate` | Remove coincident geometry across copies |
| `BRepGraph_Validate` | Full topology validity checking |
| `BRepGraph_Compact` | Compact sparse topology after edits |
| `BRepGraph_Builder` | Construct graphs programmatically (no TopoDS needed) |
| `BRepGraph_Tool` | Geometry access analogous to `BRep_Tool`, but over graph nodes |

**Impact on RCAD**: The `BRepGraph` layer is what OCCT will use for persistent naming, richer history (edges, vertices, solids), and safer boolean post-processing. It is the foundation for future defeaturing, rib, and history-based re-feature workflows. RCAD's topology is currently a flat Rust struct with no equivalent graph API -- this will become a growing architectural gap as OCCT users expect history-aware operations.

**Recommended direction**: Design a lightweight graph topology wrapper (`rcad-kernel` or new `rcad-graph` crate) that maps existing `BRep` topology to a history-capable graph. This can start as read-only traversal and grow to support mutation.

## 1. Robust Boolean Architecture (gap: High)

This remains the largest practical gap.

OCCT's TKBO provides a tiered boolean platform beyond simple Fuse/Cut/Common:

| OCCT class | Purpose | RCAD equivalent |
|---|---|---|
| BOPAlgo_PaveFiller | Interference computation core | Present (pave_filler.rs) |
| BOPAlgo_Builder | Shape assembly from pave data | Partial (builder.rs) |
| BRepAlgoAPI_Fuse/Cut/Common/Section | Standard boolean API | Present |
| BRepAlgoAPI_Splitter | Split objects by tools | Missing |
| BOPAlgo_CellsBuilder | Reusable split-cell graph | Missing |
| BOPAlgo_MakerVolume | Solid from split faces/shells | Missing |
| BRepAlgoAPI_Defeaturing | Remove interior features (8.0) | Missing |
| BOPAlgo_MakeConnected | Connect disconnected geometry (8.0) | Missing |
| BOPAlgo_CheckerSI | Self-intersection checker | Partial (brep_check.rs) |
| Fuzzy tolerance option | Near-coincident robustness | Missing |
| Gluing option | Shared-face fast path | Missing |
| Result simplification | Same-domain unify after boolean | Missing |

### What RCAD should add (ordered by impact)

1. General-fuse split-first core independent of final boolean classification.
2. Splitter API (split objects by arbitrary tool shapes).
3. Fuzzy tolerance option for near-coincident geometry.
4. Result simplification pass (same-domain face merging, internal face removal after fuse, small-edge cleanup).
5. CellsBuilder for reusable split-cell expressions (needed for CAE partitioning).
6. Defeaturing (remove small interior pockets/bosses from imported parts).

## 2. Healing and Tolerance System (gap: High, partially improved)

OCCT's TKShHealing comprises 10 packages. After P3–P4, RCAD has improved from "one sub-package equivalent" to approximately two.

| OCCT package | Key capability | RCAD status |
|---|---|---|
| ShapeFix_Face | Repair degenerated/invalid faces | Partial (brep_repair.rs) |
| ShapeFix_Wire | Reorder wires, close gaps, remove degenerate edges | Partial |
| ShapeFix_Edge | SameParameter, SameRange, degenerated edges | **SameParameter: Present (P4)**; **SameRange: Present (baseline scan+repair)** |
| ShapeFix_Shell | Repair shell orientation, manifoldness | **Shell analyzer: Present (P4)**; fix: Missing |
| ShapeFix_Solid | Solid closure, shell orientation | Missing |
| ShapeAnalysis_Surface | UV consistency, surface bounds analysis | Missing |
| ShapeAnalysis_Wire | Wire gap, self-intersection, area | **Partial (gap/self-intersection report present)** |
| ShapeUpgrade_UnifySameDomain | Merge co-planar/co-cylindrical faces | **Missing (critical for boolean results)** |
| ShapeCustom | BSpline restriction, convert to indirect | Missing |
| ShapeProcess | Batch pipeline with operator chain | Missing |

### What RCAD should add next (ordered by impact)

- **ShapeUpgrade_UnifySameDomain equivalent**: merges adjacent co-planar or co-cylindrical faces -- this alone cleans up most boolean result artifacts. Highest single-item impact.
- SameRange repair deepening: extend beyond baseline range alignment to richer edge/surface consistency.
- Wire analysis deepening: add area/orientation/quality metrics and automatic fix strategies.
- Import healing pipeline: promote current JSON diagnostics wiring to full staged analyze→diagnose→heal orchestration.
- Tolerance propagation rules after boolean/split/sew operations.

## 3. Feature Library (gap: Medium, partially implemented)

OCCT's TKFeat provides parametric feature operations built on top of boolean operations. P5 added prism (boss/pocket) and cylindrical hole.

| OCCT class | Feature type | RCAD status |
|---|---|---|
| BRepFeat_MakePrism | Boss/pocket (blind/through/up-to) | **Present (P5) -- polygon profile only** |
| BRepFeat_MakeDPrism | Draft prism | **Present (P5+)** |
| BRepFeat_MakeRevol | Revolution feature | **Present (P5+)** |
| BRepFeat_MakeCylindricalHole | Cylindrical hole | **Present (P5)** |
| BRepFeat_MakeLinearForm | Rib/slot (linear) | Missing |
| BRepFeat_MakeRevolutionForm | Rib/slot (revolved) | Missing |
| BRepFeat_Gluer | Glue shapes at interface | Missing |
| BRepFeat_SplitShape | Split face by wire | Missing |

Prism, draft prism, revolution, and hole now cover most common mechanical operations. Remaining items (rib/slot, SplitShape) are still important for sheet-metal and structural workflows.

## 4. Document Model and Exchange Depth (gap: Medium-High, partially improved)

OCCT's XDE/XCAF layer provides a structured document model far deeper than RCAD's current support. P6 added material and layer parsing.

| OCCT attribute | Purpose | RCAD status |
|---|---|---|
| XCAFDoc_Color / ColorTool | Per-face/part colors | Present (step IO) |
| XCAFDoc_ShapeTool | Assembly tree | Present |
| XCAFDoc_Material / MaterialTool | Material assignment + density | **Present (P6)** |
| XCAFDoc_Layer / LayerTool | Layer/group membership | **Present (P6)** |
| XCAFDoc_DimTol / GeomTolerance | GDT annotations | Missing |
| XCAFDimTolObjects | Dimensional tolerance objects | Missing |
| XCAFNoteObjects | Notes and markup | Missing |
| XCAFView | View definitions | Missing |
| StepKinematics | Joint/mechanism metadata | Missing |

For STEP AP242 round-trip:

| AP242 area | RCAD status |
|---|---|
| Geometry + topology write | Present |
| Color/style write | Present |
| Assembly write | Present |
| Material/layer write | **Present (P6)** |
| AP242 GDT write | Missing |
| AP242 read (import) | Basic |
| Kinematics read/write | Missing |
| FEA entity read | Missing |

**Notable in OCCT 8.0.0**: STEP general properties export (`property_definition` entities for arbitrary string metadata) and stream-based DE_Wrapper read/write are new. RCAD now has baseline `GENERAL_PROPERTY` + `PROPERTY_DEFINITION` read/write linkage, but still lacks deeper AP242 relationship coverage and DE_Wrapper-style stream APIs.

## 5. Post-Operation Simplification (gap: High, small-edge cleanup now present)

OCCT's ShapeUpgrade_UnifySameDomain and BOPAlgo_Builder cleanup passes automatically clean up boolean results. P3 added small-edge cleanup; same-domain face merging is still missing.

RCAD currently produces boolean results that may contain:
- ~~tiny edges from near-coincident intersections~~ (small-edge cleanup added in P3)
- many small adjacent co-planar faces where one merged face would suffice **(still missing)**
- dangling internal faces after fuse operations **(still missing)**
- mismatched tolerances at operation boundaries **(still missing)**

The most impactful remaining item is **same-domain face unification** (merging co-planar or co-cylindrical adjacent faces). This is the single change that most visibly improves boolean output quality for downstream meshing and exchange.

## 6. Geometry Evaluation Breadth (NEW gap from OCCT 8.0.0)

OCCT 8.0.0 introduces `GeomEval` / `Geom2dEval` evaluation classes that extend the geometry hierarchy with new curve and surface types:

| New OCCT type | Description | RCAD status |
|---|---|---|
| `GeomEval_CircularHelixCurve` | Circular helix curve (TKHelix) | **Present (P7 partial)** |
| `GeomEval_SineWaveCurve` / `Geom2dEval_SineWaveCurve` | Sine wave along a baseline | Missing |
| `Geom2dEval_ArchimedeanSpiralCurve` | Archimedean spiral | **Present (P7 partial)** |
| `Geom2dEval_LogarithmicSpiralCurve` | Logarithmic spiral | **Present (P7 partial)** |
| `Geom2dEval_CircleInvoluteCurve` | Circle involute (gear tooth profile) | **Present (P7 partial)** |
| `GeomEval_TBezierSurface` / `AHTBezierSurface` | Parametric generalized Bezier surfaces | Missing |
| `GeomFill_Gordon` | N×M transfinite surface from curve network | Missing |
| `ExtremaPC` | Point-to-curve extrema with per-type dispatch | Partial (basic projection) |

The gear-tooth involute curve is particularly important for mechanical design. Gordon surface and helix are the others most likely to be demanded by RCAD users in production scenarios.

## Roadmap

### ✅ Completed (P3–P6)

| Deliverable | Status |
|---|---|
| Small-edge / degenerate-edge cleanup (`remove_small_edges`, `SimplifyOptions`) | **Done** |
| SameParameter diagnosis (`diagnose_same_parameter`) | **Done** |
| SameParameter repair (`fix_same_parameter_with_scan`) | **Done** |
| Shell manifold / open-edge analyzer (`analyze_shell_topology`, `ShellTopologyReport`) | **Done** |
| Polygon prism boss/pocket feature (`make_prism`, `build_polygon_prism`) | **Done** |
| STEP material extraction (`StepMaterial`, density-aware) | **Done** |
| STEP layer extraction (`StepLayer`) | **Done** |

### P3 Remaining: Boolean Robustness and Result Simplification (highest priority)

Target: imported and modeled solids through boolean workflows with much better robustness and cleaner results.

| Deliverable | Effort | Status |
|---|---|---|
| **Same-domain face unification** | Medium | **Not started (highest single impact)** |
| Internal-face removal after fuse | Small | Not started |
| Splitter API (object/tool split) | Medium | Not started |
| Fuzzy tolerance option | Small | Not started |
| CellsBuilder (split-cell graph) | Large | **Done (baseline expression evaluator)** |
| Defeaturing pass | Large | Not started |
| Full history to edges and solids | Medium | Not started |

### P4 Remaining: Industrial Healing Pipeline

Target: imported CAD data can be analyzed, repaired, and pushed into modeling/boolean workflows safely.

| Deliverable | Effort | Status |
|---|---|---|
| **ShapeUpgrade_UnifySameDomain equivalent** | Medium | **Not started (critical)** |
| SameRange repair | Medium | **Done (baseline scan+repair)** |
| Face-on-surface consistency checker | Medium | **Done (baseline diagnosis API)** |
| Wire gap / self-intersection analyzer | Medium | **Done (wire report API)** |
| Import analyze/heal/report pipeline | Medium | **Partial (healing JSON + wire stats)** |
| Tolerance propagation after boolean/sew | Large | Not started |

### P5 Remaining: Feature Library

Target: full range of parametric features.

| Deliverable | Effort | Status |
|---|---|---|
| Draft prism | Medium | **Done** |
| Revolution feature | Medium | **Done** |
| SplitShape (face by wire) | Small | Not started |
| MakeLinearForm (rib/slot) | Large | Not started |

### P6 Remaining: Document Model and AP242 Depth

Target: XCAF-comparable document model; AP242 GDT and kinematics round-trip.

| Deliverable | Effort | Status |
|---|---|---|
| GDT / DimTol write | Medium | Not started |
| AP242 read (full import) | Large | Not started |
| Kinematics read | Large | Not started |
| Persistent naming in history | Large | Not started |
| STEP general property export (string metadata) | Small | **Done (baseline `GENERAL_PROPERTY` + `PROPERTY_DEFINITION`)** |

### P7 (NEW): Graph Topology and Evaluation Breadth

These are new gaps opened by OCCT 8.0.0's architectural leaps.

| Deliverable | Effort | Notes |
|---|---|---|
| Graph-based topology wrapper (`BRepGraph` equivalent) | Large | Foundation for persistent naming + richer history |
| Circular helix curve (for spring/coil modeling) | Small | **Done (kernel + arc-length support)** |
| Circle involute curve (gear tooth profiles) | Small | **Done (kernel)** |
| Gordon surface (N×M transfinite fill) | Medium | Upgrades surface fill beyond 4-boundary Coons |
| Archimedean / logarithmic spiral curves | Small | **Done (kernel)** |

## Production Readiness Gap Assessment

Answering the question: **how far is RCAD from being production-grade?**

### Definition used here

"Production-grade" = a downstream application (CAE preprocessor, manufacturing toolpath planner, or PDM system) can rely on RCAD for:
1. Importing real-world STEP files without manual intervention.
2. Running boolean and feature operations on those imports reliably.
3. Exporting results with correct metadata.
4. Supporting history-based editing and re-feature.

### Gap summary by category (updated after P3–P6)

| Category | Current state | Production bar | Delta | Change since last revision |
|---|---|---|---|---|
| Geometry kernel | 90% | 95% | Minor gaps (helix, involute, some surface types) | Unchanged |
| B-Rep model | 75% | 90% | Compound/CompSolid, non-manifold | Unchanged |
| Graph topology layer | 0% | 50% | BRepGraph equivalent needed for history/naming | **New gap (OCCT 8.0.0)** |
| Sweep / loft / extrude | 80% | 90% | N-sided fill (Gordon), medial axis | Unchanged |
| Fillet / chamfer | 70% | 85% | Corner cases, angle-mode chamfer | Unchanged |
| Boolean core | 55% | 85% | Fuzzy, splitter, cleanup, defeaturing | Unchanged |
| Healing pipeline | **35%** | 80% | SameRange, ShapeUpgrade, wire analyzers | **+15% (P3–P4)** |
| Feature library | **35%** | 60% | SplitShape, rib/slot | **+35% (P5+)** |
| Document / XCAF | **40%** | 75% | GDT, DimTol, persistent naming | **+10% (P6)** |
| STEP exchange depth | 60% | 80% | AP242 full read, GDT round-trip, property_definition chain | **+10% (P6)** |
| Meshing controls | 40% | 75% | Tunable deflection, incremental update | Unchanged |
| Geometry eval breadth | 78% | 85% | Gordon fill and advanced Eval surfaces; OCCT 8.0.0 raised the bar | **Improved (P7 partial)** |

### Estimated work remaining (revised)

Assuming 1 developer working full-time on kernel work:

| Phase | Focus | Calendar estimate | Status |
|---|---|---|---|
| P3 remaining | Same-domain unification, splitter, fuzzy | 2 months | In progress |
| P4 remaining | ShapeUpgrade, face consistency, tolerance propagation | 1-2 months | In progress |
| P5 remaining | SplitShape, rib/slot | 1-2 months | In progress |
| P6 remaining | GDT write, AP242 full read | 2-3 months | Partial |
| P7 | Graph topology + Gordon surface | 2-3 months | In progress |
| Hardening, edge cases, test coverage | Ongoing | +2 months across all phases | Ongoing |

**Total to reach credible production baseline: approximately 8–12 months** (revised down from 9–13 due to delivered P4/P5/P7 partial items; still bounded by same-domain unification, AP242 read depth, and Gordon/BRepGraph scope).

### What would most accelerate the timeline

1. **Same-domain face unification** (P3 remaining) -- immediately makes all existing boolean results cleaner and is the single highest-ROI item remaining.
2. **Fuzzy boolean tolerance** (P3) -- removes the most common class of near-coincident import failures.
3. **SameRange + wire gap analysis** (P4) -- together with the SameParameter work already done, completes the healing pipeline for import.
4. **Deferring BRepGraph / graph topology** (P7) until a second developer is available -- it is architecturally important but not blocking immediate production use.

## Module-Level Task Breakdown

### libs/rcad-algorithms

| Task | Priority | Status |
|---|---|---|
| General fuse split-first core | P3 - High | Not started |
| Splitter API | P3 - High | Not started |
| Fuzzy / glue boolean options | P3 - High | Not started |
| **Result simplification: same-domain unification** | **P3 - High** | **Not started (highest impact)** |
| Result simplification: internal face removal | P3 - Medium | Not started |
| CellsBuilder (split-cell graph) | P3 - Medium | **✅ Done (baseline)** |
| Defeaturing pass | P3 - Medium | Not started |
| Richer history graph (edges, solids) | P3 - Medium | Not started |
| ~~SameParameter / SameRange repair~~ | ~~P4 - High~~ | **✅ SameParameter done (P4)** |
| SameRange repair | P4 - High | **✅ Done (scan+repair)** |
| ShapeUpgrade_UnifySameDomain equivalent | P4 - High | Not started |
| Face-on-surface consistency checker | P4 - Medium | **✅ Done (diagnose_face_surface_consistency)** |
| ~~Shell / manifoldness analyzer~~ | ~~P4 - Medium~~ | **✅ Done (P4)** |
| Wire gap / self-intersection analyzer | P4 - Medium | **✅ Done (analyze_wire_issues)** |
| Import analyze/heal pipeline | P4 - Medium | **Partial (JSON diagnostics wired)** |
| Tolerance propagation rules | P4 - Medium | Not started |
| ~~Small-edge cleanup~~ | ~~P3 - Medium~~ | **✅ Done (P3)** |
| ~~Feature prism / cylindrical hole~~ | ~~P5 - Medium~~ | **✅ Done (P5)** |
| Draft prism feature | P5 - Medium | **✅ Done** |
| Revolution feature | P5 - Medium | **✅ Done** |
| Rib / slot feature | P5 - Low | Not started |

### libs/rcad-kernel

| Task | Priority | Status |
|---|---|---|
| Compound / CompSolid topology | P6 - Medium | Not started |
| Non-manifold topology support | P3 - Low | Not started |
| Persistent naming hooks | P6 - Medium | Not started |
| Richer validity analysis | P4 - Medium | Not started |
| Tolerance propagation rules | P4 - High | Not started |
| **Graph topology wrapper (BRepGraph equivalent)** | **P7 - High** | **Not started** |

### libs/rcad-modeling

| Task | Priority | Status |
|---|---|---|
| **N-sided surface fill -- Gordon N×M transfinite** | **P7 - Medium** | **Not started** |
| **Circle involute curve (gear tooth)** | **P7 - Medium** | **✅ Done** |
| **Circular helix curve** | **P7 - Low** | **✅ Done** |
| Normal projection of wire onto surface | P3 - Low | Not started |
| Stabilize advanced fillet corner cases | P3 - Medium | Not started |
| Angle-mode chamfer | P3 - Low | Not started |
| 2-D fillet/chamfer API | Convenience | Not started |

### libs/rcad-step

| Task | Priority | Status |
|---|---|---|
| AP242 full read (import) | P6 - High | Not started |
| GDT / DimTol read and write | P6 - Medium | Not started |
| ~~Material / layer read and write~~ | ~~P6 - Low~~ | **✅ Done (P6)** |
| STEP general property export (arbitrary metadata) | P6 - Low | **Partial (`GENERAL_PROPERTY` + `PROPERTY_DEFINITION` baseline; deeper chains pending)** |
| Stream-based read/write (DE_Wrapper style) | P6 - Low | Not started |
| Kinematics read | P6 - Low | Not started |
| Import healing pipeline integration | P4 - Medium | **Partial (healing report JSON + wire stats)** |
| Stronger PCurve / tolerance validation on export | P4 - Medium | Not started |

### libs/rcad-render

| Task | Priority | Status |
|---|---|---|
| Tunable meshing deflection / angular tolerances | P3 - Medium | Not started |
| Incremental cache invalidation for edited models | P2.1 followup - Medium | Not started |
| Better HLR on curved geometry | P2.2 followup - Low | Not started |

## Recommended Non-Goals for Now

- **Helix geometry (TKHelix)**: OCCT 8.0.0 ships a full helix toolkit. It is still a low-priority RCAD addition unless spiral/spring modeling is required. Revisit in P7.
- VRML / PLY / glTF output -- low engineering value; focus on STEP AP242 depth first.
- Full FEA entity round-trip via AP209 -- too specialized.
- Comprehensive IGES B-Rep support -- STEP AP242 is superior for the same use cases.
- More geometry primitive types -- the current set is sufficient.
- Premature rendering optimization -- wait until boolean/healing/document layers improve.
- **BRepGraph full implementation** -- important long-term, but implementing 49 000+ lines of equivalent code is a multi-quarter project best done after P3–P6 remaining items are closed.

## Bottom Line

After P3–P6, RCAD has meaningfully advanced its healing, feature, and metadata coverage. The situation as of April 2026:

- The geometric foundation is solid and competitive.
- The feature toolbox covers the most common operations (extrude, revolve, sweep, fillet, boolean, prism, cylindrical hole).
- The weak points remain **boolean robustness and result simplification**, **healing depth**, and **AP242 exchange completeness**.
- OCCT 8.0.0 (imminent) raises the bar further: BRepGraph, Gordon surfaces, evaluation geometry breadth, and the TKHelix toolkit create new gaps that did not exist 6 months ago.

The three areas that dominate the remaining production gap, in order:

1. **Boolean robustness and result simplification** -- same-domain face merge is the single most impactful open item.
2. **Healing depth** -- SameRange repair and wire analysis complete the most critical path.
3. **BRepGraph / graph topology** -- not blocking today, but will be required for history-based re-feature and persistent naming as applications grow.

At the current trajectory, RCAD can become a credible industrial CAD kernel with **8–12 months of focused kernel work**, assuming the priority order above and one dedicated developer. If BRepGraph and Gordon/eval-breadth work are deferred to a later phase, the immediate production baseline can be reached in roughly **5–7 months** (same-domain unification + healing pipeline completion + AP242 read depth).

