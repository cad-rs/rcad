# RCAD vs OCCT Gap Analysis

Date: 2026-04-13 (revised against OCCT 7.9.3 stable / 8.0.0-rc5, incorporating P3–P7 completions)

## Purpose

This document turns the current RCAD capability inventory into an OCCT-oriented gap analysis and an execution roadmap. The goal is not to chase feature-count parity blindly, but to identify the areas where RCAD still differs most from OCCT in practical engineering use.

The key conclusion is:

RCAD has made substantial progress since the last revision. P3–P7 work has added same-domain face unification Phase 1 (plane/cylinder/cone/torus/sphere) and Phase 2 (topological + geometric double-validation with UV-region/edge-continuity checks), tolerance propagation, shell and wire diagnostics, SplitShape and rib/slot baselines, AP242 metadata chain parsing (PROPERTY_DEFINITION_REPRESENTATION + DIMENSIONAL_LOCATION/SIZE + GEOMETRIC_TOLERANCE), a baseline read-only BRepGraph traversal API, and Gordon surface transfinite interpolation support. However, OCCT 8.0.0 (currently at rc5, release imminent) has itself leapt forward significantly -- it introduces an entirely new graph-based topology representation (BRepGraph, 49 000+ lines), a production-grade Gordon framework, a fully redesigned geometry evaluation architecture, the TKHelix toolkit, and defeaturing + connected-shape APIs. The production gap is narrowing but OCCT's target is moving.

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
- **P6**: STEP material/layer extraction, GENERAL_PROPERTY extraction and export API, AP242 metadata-chain expansion (DATUM / DATUM_SYSTEM / kinematic pair entities) with read/write baselines.
- **P7 (partial)**: Circular helix, circle involute, Archimedean/logarithmic spiral, 2D and 3D sine-wave evaluators.

### Areas where RCAD is still behind OCCT

- **Boolean robustness** on curved solids and near-coincident geometry remains the largest practical gap. Splitter/CellsBuilder/fuzzy tolerance are now baseline-present (including adaptive retry escalation), and MakeConnected-style cleanup now has iterative growth-aware baseline passes with tolerance-cap safety plus scoped mode with semantic seed strategies (short-edge / near-duplicate / tolerance-tagged / multi-PCurve / topology-seam-candidates / hybrid), plus issue-driven pre-make-connected and healing-stall fallback integration in the boolean pipeline; defeaturing and richer connectivity rebuilding semantics are still absent.
- **BRepGraph depth**: OCCT 8.0.0 ships an extensive graph-based BRep API (49 000+ lines) with history tracking, mutation guards, deduplication, and validation. RCAD now has a baseline topology graph with history events, a full RAII `BRepGraphMutGuard` (checkpoint/commit/rollback/validate), `validate_invariants`, and compact/dedup primitives. The remaining gap is rich persistent-history naming.
- **Healing pipeline depth**: same-domain unification now has Phase 1 baseline (plane/cylinder/cone/torus/sphere) and Phase 2 double-validation (topological edge-continuity + geometric UV-region checks); healing orchestration now includes issue-driven pre-make-connected, SameRange/SameParameter parametric consistency prepass with iterative reconciliation, plus stalled fallback, but broader ShapeFix/ShapeProcess coverage remains missing.
- **Exchange**: STEP AP242 write is present and AP242 read now covers an expanded metadata baseline (GDT/DimTol + DATUM / DATUM_SYSTEM + kinematic pair entities), but full semantic AP242 import and FEA entities are still not covered.
- **Geometry evaluation breadth**: OCCT 8.0.0 adds helicoid, spiral, ellipsoid, and parametric curve evaluation classes not present in RCAD.
- **Gordon surface robustness**: N×M transfinite interpolation baseline exists in RCAD, but OCCT's GeomFill_Gordon remains richer and more production-tuned.
- **Topology containers**: No Compound / CompSolid; no non-manifold topology.
- **Feature library**: Prism, draft prism, revolution feature, cylindrical hole, rib/slot baseline, and SplitShape baseline are present; advanced constraints and robustness remain.

## Capability Matrix

| Domain | RCAD today | OCCT reference | Gap level | Recommended direction |
|---|---|---|---|---|
| Geometry kernel | Strong analytic + spline coverage | TKG2d / TKG3d / TKGeomBase | Low | Maintain correctness; helix curves now in OCCT TKHelix |
| Geometry evaluation breadth | Basic D0/D1 on standard types | GeomEval/Geom2dEval (helicoids, spirals, ellipsoids, parametric) | Medium | Add helicoid, spiral, ellipsoid eval classes |
| Gordon / N×M surface fill | Baseline transfinite Gordon support | GeomFill_Gordon (8.0) | Medium | Improve robustness and continuity controls |
| Point cloud analysis | Missing | PointSetLib (8.0) -- PCA, inertia, dimensionality | Low | Low priority; add if needed for import QC |
| Core B-Rep model | Vertex/Edge/Wire/Face/Shell/Solid | TKBRep | Medium | Add Compound, CompSolid, non-manifold support |
| Graph-based topology layer | Baseline traversal + history events + **RAII mutation guard (`BRepGraphMutGuard`)** + checkpoint/rollback + `validate_invariants` + validate/compact/dedup primitives | BRepGraph (8.0, 49k lines) -- history, mutation, dedup, validate | Medium | Add richer persistent-history semantics and persistent naming |
| Primitive and sweep modeling | Good breadth (extrude/revolve/loft/sweep) | TKPrim + TKOffset | Low | Edge cases; N-sided fill, normal projection |
| Fillet / chamfer | Variable-radius fillet, chamfer | TKFillet | Medium | Angle-mode chamfer, 2-D fillet API, corner cases |
| Thicken / draft / offset | Present | TKOffset BRepOffset | Medium | Shell offset (MakeOffsetShape), evolved surface |
| Feature library | Prism + draft prism + revolution + cylindrical hole + rib/slot + SplitShape baseline | TKFeat (boss/pocket/rib/hole) | Medium | Expand feature constraints and robustness |
| Boolean framework | Fuse/Cut/Common/Section + imprint + splitter/cells + MakerVolume baseline + adaptive fuzzy retry + iterative/growth-aware/capped make-connected cleanup (global+scoped, semantic seeds incl. tolerance-tagged + multi-PCurve + topology seam candidates, with history-informed seed-edge preference plus low-coverage heuristic augmentation, and seed source/count metadata + stable edge labels reporting) | TKBO | High | glue/deeper scoped connectivity rebuilding semantics, deeper failure recovery |
| Post-op simplification | Small-edge cleanup + same-domain unification baseline | ShapeUpgrade_UnifySameDomain, BOPAlgo cleanup | Medium-High | Internal-face removal and richer same-domain criteria |
| Healing and validation | SameParameter/SameRange + shell + wire diagnostics (P4 partial) + staged healing with issue-driven pre-make-connected, iterative parametric consistency pass, and make-connected-on-stall fallback | TKShHealing (10 packages) | High | ShapeUpgrade_UnifySameDomain, ShapeProcess, tolerance rules |
| Topology history | Face-level history | BOPAlgo history, BRepGraph_History (8.0), OCAF naming | Medium | Extend to edges/vertices; persistent naming; BRepGraph history |
| STEP exchange: write | AP214 + AP242 + material/layer + GENERAL_PROPERTY | TKDESTEP STEPCAFControl | Medium | GDT write, property_definition relations, PCurve validation |
| STEP exchange: read | Basic import + expanded AP242 metadata-chain baseline (PDR/DimLoc/DimSize/GeomTol/Datum/DatumSystem/KinematicPair) | TKDESTEP + STEPCAFControl_Reader | High | Full semantic AP242 read and FEA entities |
| IGES exchange | Mesh bridge only | TKDEIGES | High | Add analytic/B-Rep IGES or document as non-goal |
| Assembly / document model | Colors + shape tree + material + layer (P6 partial) | TKXCAF (XCAFDoc_*) | Medium | DimTol, GDT annotations, notes, persistent attributes |
| Meshing and visualization | mesh_dirty caching (P2.1), HLR dense silhouettes (P2.2) | TKMesh + TKHLR | Medium | Tunable deflection/angular tolerances, incremental remesh |
| Thread safety | Single-threaded | BRepCheck thread-safe (8.0), thread-local error handlers (8.0) | Low | Not urgent unless parallel workflows are added |

## Highest-Priority Gaps

## 0. BRepGraph: Graph-Based Topology Layer (NEW in OCCT 8.0.0 -- gap: Medium-High)

OCCT 8.0.0-rc5 introduces `BRepGraph`, an entirely new graph-based representation of topology and BRep geometry as an alternative to the traditional `TopoDS_Shape` linked structure. This is 49 000+ lines of new code with 20+ GTest files.

Key capabilities exposed by OCCT `BRepGraph`, with RCAD status:

| BRepGraph feature | Purpose | RCAD status |
|---|---|---|
| `BRepGraph_NodeId`-typed incidence tables | Graph traversal without pointer chasing | **Baseline present** (index-based adjacency tables) |
| `BRepGraph_History` | Persistent shape history (which new faces came from which old faces) | Partial |
| `BRepGraph_MutGuard` | Safe mutation of topology with invariant checking | **Done** |
| `BRepGraph_Deduplicate` | Remove coincident geometry across copies | **Baseline present** |
| `BRepGraph_Validate` | Full topology validity checking | **Baseline present** (`validate_invariants`) |
| `BRepGraph_Compact` | Compact sparse topology after edits | **Baseline present** |
| `BRepGraph_Builder` | Construct graphs programmatically (no TopoDS needed) | **Done** |
| `BRepGraph_Tool` | Geometry access analogous to `BRep_Tool`, but over graph nodes | **Done** |

**Impact on RCAD**: The `BRepGraph` layer is what OCCT will use for persistent naming, richer history (edges, vertices, solids), and safer boolean post-processing. It is the foundation for future defeaturing, rib, and history-based re-feature workflows. RCAD now has a baseline topology graph traversal API with history events, a full RAII `BRepGraphMutGuard` (commit/rollback/checkpoint/validate_invariants), programmatic graph construction (`BRepGraphBuilder`), graph-node geometry access (`BRepGraphTool`), and compact/dedup primitives. The remaining gap is richer persistent-history semantics and stable persistent naming.

**Recommended direction**: Design a lightweight graph topology wrapper (`rcad-kernel` or new `rcad-graph` crate) that maps existing `BRep` topology to a history-capable graph. This can start as read-only traversal and grow to support mutation.

## 1. Robust Boolean Architecture (gap: High)

This remains the largest practical gap.

OCCT's TKBO provides a tiered boolean platform beyond simple Fuse/Cut/Common:

| OCCT class | Purpose | RCAD equivalent |
|---|---|---|
| BOPAlgo_PaveFiller | Interference computation core | Present (pave_filler.rs) |
| BOPAlgo_Builder | Shape assembly from pave data | Partial (builder.rs) |
| BRepAlgoAPI_Fuse/Cut/Common/Section | Standard boolean API | Present |
| BRepAlgoAPI_Splitter | Split objects by tools | Present (baseline split-first API) |
| BOPAlgo_CellsBuilder | Reusable split-cell graph | Present (baseline expression evaluator) |
| BOPAlgo_MakerVolume | Solid from split faces/shells | **Done (baseline `MakerVolume` API: region-mask / explicit-index / `CellExpr` assembly over reusable split cells, with history variant)** |
| BRepAlgoAPI_Defeaturing | Remove interior features (8.0) | **Done (baseline `defeature_brep` + cylindrical feature detection/fill-remove + small-face identification)** |
| BOPAlgo_MakeConnected | Connect disconnected geometry (8.0) | Partial (iterative + growth-aware + tolerance-capped baseline merge/small-edge cleanup; scoped mode with semantic short/near-dup/tolerance-tagged/multi-PCurve/topology-seam/hybrid seeds) |
| BOPAlgo_CheckerSI | Self-intersection checker | Partial (brep_check.rs) |
| Fuzzy tolerance option | Near-coincident robustness | Present (baseline + robust retry ladder) |
| Gluing option | Shared-face fast path | **Partial (`BooleanOptions.use_glue` + shared-face skip/merge baseline in filler/builder)** |
| Result simplification | Same-domain unify + internal-face cleanup after boolean | Done (same-domain Phase 1+2; internal-face removal Phase 1+2 baseline) |

### What RCAD should add (ordered by impact)

1. Deepen make-connected cleanup from iterative baseline (merge/small-edge) to scoped connectivity rebuilding.
2. Add gluing/shared-face fast path for robust and faster near-contact operations.
3. Deepen fuzzy strategy with failure-class-aware escalation tuning (baseline adaptive policy now present).
4. Result simplification pass (same-domain face merging, internal face removal after fuse, small-edge cleanup).
5. Deepen defeaturing robustness (feature-class breadth, topology healing after suppression).

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
| ShapeUpgrade_UnifySameDomain | Merge co-planar/co-cylindrical faces | **Done (Phase 1: plane/cylinder/cone/torus/sphere; Phase 2: topological guards)** |
| ShapeCustom | BSpline restriction, convert to indirect | Missing |
| ShapeProcess | Batch pipeline with operator chain | Missing |

### What RCAD should add next (ordered by impact)

- Deepen same-domain unification beyond current baseline (more surface classes, safer topology guards, stronger merge diagnostics).
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
| BRepFeat_MakeLinearForm | Rib/slot (linear) | **Present (baseline)** |
| BRepFeat_MakeRevolutionForm | Rib/slot (revolved) | **Present (baseline)** |
| BRepFeat_Gluer | Glue shapes at interface | Missing |
| BRepFeat_SplitShape | Split face by wire | **Present (baseline)** |

Prism, draft prism, revolution, hole, rib/slot, and SplitShape now cover a broader baseline feature set. Remaining gap is mainly robustness and advanced feature variants.

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
| StepKinematics | Joint/mechanism metadata | Partial (baseline AP242 kinematic-pair metadata read/write support) |

For STEP AP242 round-trip:

| AP242 area | RCAD status |
|---|---|
| Geometry + topology write | Present |
| Color/style write | Present |
| Assembly write | Present |
| Material/layer write | **Present (P6)** |
| AP242 GDT write | Missing |
| AP242 read (import) | Basic + expanded metadata-chain baseline (property-definition representation + dimensional location/size + geometric tolerance entities + DATUM / DATUM_SYSTEM / kinematic pair entities) |
| Kinematics read/write | Partial (metadata-level read/write baseline for pair entities) |
| FEA entity read | Missing |

**Notable in OCCT 8.0.0**: STEP general properties export (`property_definition` entities for arbitrary string metadata) and stream-based DE_Wrapper read/write are new. RCAD now has baseline `GENERAL_PROPERTY` + `PROPERTY_DEFINITION` read/write linkage and `Read`/`Write` stream APIs, but still lacks deeper AP242 relationship coverage.

## 5. Post-Operation Simplification (gap: High, small-edge cleanup now present)

OCCT's ShapeUpgrade_UnifySameDomain and BOPAlgo_Builder cleanup passes automatically clean up boolean results. RCAD now has baseline same-domain face unification plus small-edge cleanup.

RCAD currently produces boolean results that may contain:
- ~~tiny edges from near-coincident intersections~~ (small-edge cleanup added in P3)
- many small adjacent same-domain faces where one merged face would suffice **(partially improved by baseline unification)**
- dangling internal faces after fuse operations **(still missing)**
- mismatched tolerances at operation boundaries **(still missing)**

The most impactful remaining simplification items are deeper same-domain criteria, internal-face removal after fuse, and stronger tolerance reconciliation across operation boundaries.

## 6. Geometry Evaluation Breadth (NEW gap from OCCT 8.0.0)

OCCT 8.0.0 introduces `GeomEval` / `Geom2dEval` evaluation classes that extend the geometry hierarchy with new curve and surface types:

| New OCCT type | Description | RCAD status |
|---|---|---|
| `GeomEval_CircularHelixCurve` | Circular helix curve (TKHelix) | **Present (P7 partial)** |
| `GeomEval_SineWaveCurve` / `Geom2dEval_SineWaveCurve` | Sine wave along a baseline | **Present (P7: 2D + 3D evaluators)** |
| `Geom2dEval_ArchimedeanSpiralCurve` | Archimedean spiral | **Present (P7 partial)** |
| `Geom2dEval_LogarithmicSpiralCurve` | Logarithmic spiral | **Present (P7 partial)** |
| `Geom2dEval_CircleInvoluteCurve` | Circle involute (gear tooth profile) | **Present (P7 partial)** |
| `GeomEval_TBezierSurface` / `AHTBezierSurface` | Parametric generalized Bezier surfaces | Missing |
| `GeomFill_Gordon` | N×M transfinite surface from curve network | **Present (baseline)** |
| `ExtremaPC` | Point-to-curve extrema with per-type dispatch | **Partial → improved** (analytic O(1) dispatch for Line/Circle/Ellipse; Newton-Raphson fallback for all other curve types; 11 unit tests) |

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
| **Same-domain face unification** | Medium | **Done (baseline plane/cylinder/cone/torus/sphere)** |
| Internal-face removal after fuse | Small | Done (Phase 1 threshold-based cleanup; Phase 2 true-duplicate detection) |
| Splitter API (object/tool split) | Medium | **Done (baseline split-first API)** |
| Fuzzy tolerance option | Small | **Done (boolean/split options + adaptive retry policy baseline)** |
| CellsBuilder (split-cell graph) | Large | **Done (baseline expression evaluator)** |
| MakeConnected baseline pass | Medium | **Done (iterative + growth-aware + tolerance-capped merge-near vertices + small-edge cleanup; global+scoped modes with semantic short/near-dup/tolerance-tagged/multi-PCurve/topology-seam/hybrid seeds, plus history-informed seed-edge preference with low-coverage heuristic augmentation and scoped seed source/count metadata + stable edge-label reporting)** |
| Defeaturing pass | Large | Done (baseline: cylindrical feature detection + boolean fill/remove and small-face identification) |
| Full history to edges and solids | Medium | **Done (DS a_vertex_count/a_edge_count boundary tracking; annotate_history_from_ds position-matches result vertices/edges to DS origin ranges; populates BooleanHistory.vertex_origins and edge_origins in both sequential and parallel build paths; VertexOrigin: FromA/FromB/Intersection; EdgeOrigin: FromA/FromB/Generated/SplitFromA/SplitFromB)** |

### P4 Remaining: Industrial Healing Pipeline

Target: imported CAD data can be analyzed, repaired, and pushed into modeling/boolean workflows safely.

| Deliverable | Effort | Status |
|---|---|---|
| **ShapeUpgrade_UnifySameDomain equivalent** | Medium | **Done (Phase 1 analytic same-domain merge; Phase 2 topological + UV guards)** |
| SameRange repair | Medium | **Done (baseline scan+repair)** |
| Face-on-surface consistency checker | Medium | **Done (baseline diagnosis API)** |
| Wire gap / self-intersection analyzer | Medium | **Done (wire report API)** |
| Import analyze/heal/report pipeline | Medium | **Done (staged analyze/repair/final report + JSON integration)** |
| Tolerance propagation after boolean/sew | Large | **Done (baseline propagation API)** |

### P5 Remaining: Feature Library

Target: full range of parametric features.

| Deliverable | Effort | Status |
|---|---|---|
| Draft prism | Medium | **Done** |
| Revolution feature | Medium | **Done** |
| SplitShape (face by wire) | Small | **Done (baseline)** |
| MakeLinearForm (rib/slot) | Large | **Done (baseline)** |

### P6 Remaining: Document Model and AP242 Depth

Target: XCAF-comparable document model; AP242 GDT and kinematics round-trip.

| Deliverable | Effort | Status |
|---|---|---|
| GDT / DimTol write | Medium | **Done (baseline AP242 metadata entities write path, including datum-reference entities)** |
| AP242 read (full import) | Large | **Done (baseline metadata-chain scope: PDR/DimLoc/DimSize/GeomTol/Datum/DatumSystem/KinematicPair extraction)** |
| Kinematics read | Large | **Done (baseline metadata extraction for kinematic pair entities)** |
| Persistent naming in history | Large | Not started |
| STEP general property export (string metadata) | Small | **Done (baseline `GENERAL_PROPERTY` + `PROPERTY_DEFINITION`)** |

### P7 (NEW): Graph Topology and Evaluation Breadth

These are new gaps opened by OCCT 8.0.0's architectural leaps.

| Deliverable | Effort | Notes |
|---|---|---|
| Graph-based topology wrapper (`BRepGraph` equivalent) | Large | **Done (O(1) adjacency, DFS/BFS traversal, dirty tracking in `brep_graph.rs`)** |
| Circular helix curve (for spring/coil modeling) | Small | **Done (kernel + arc-length support)** |
| Circle involute curve (gear tooth profiles) | Small | **Done (kernel)** |
| Gordon surface (N×M transfinite fill) | Medium | **Done (baseline evaluator)** |
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
| Geometry kernel | 92% | 95% | Minor gaps (advanced surface families and robustness) | Improved |
| B-Rep model | 75% | 90% | Compound/CompSolid, non-manifold | **+3% (non-manifold repair hints)** |
| Graph topology layer | **42%** | 50% | RAII MutGuard + checkpoint/rollback/validate + Builder/Tool done; persistent-history naming still pending | **+17% (BRepGraphMutGuard + BRepGraphCheckpoint + BRepGraphBuilder + BRepGraphTool + validate_invariants + ExtremaPC analytic dispatch)** |
| Sweep / loft / extrude | **83%** | 90% | medial axis | **+3% (wire surface projection)** |
| Fillet / chamfer | **80%** | 85% | Further corner-case depth | **+10% (angle-mode chamfer + safe wrappers + 2-D API)** |
| Boolean core | 65% | 85% | splitter/fuzzy baseline done; defeaturing and stronger cleanup remain | Improved |
| Healing pipeline | **50%** | 80% | staged reports + tolerance propagation baseline added; deep ShapeProcess semantics pending | Improved |
| Feature library | **70%** | 80% | Advanced constraints and robustness hardening | **+35% (P5+)** |
| Document / XCAF | **43%** | 75% | GDT, DimTol, full persistent naming propagation | **+13% (P6 + propagate_through_remap / identity_map / iter)** |
| STEP exchange depth | 68% | 80% | AP242 metadata-chain read + metadata-entity write baseline; full AP242 read/write depth pending | Improved |
| Meshing controls | 40% | 75% | Tunable deflection, incremental update | Unchanged |
| Geometry eval breadth | 82% | 85% | Gordon baseline + sine-wave evaluator baseline; advanced families still pending | Improved |

### Estimated work remaining (revised)

Assuming 1 developer working full-time on kernel work:

| Phase | Focus | Calendar estimate | Status |
|---|---|---|---|
| P3 remaining | defeaturing robustness + deeper make-connected/glue hardening | 1-2 months | In progress |
| P4 remaining | ShapeProcess-grade healing semantics and tolerance policy hardening | 1-2 months | In progress |
| P5 remaining | feature constraints/variants and robustness hardening | 1-2 months | In progress |
| P6 remaining | AP242 semantic read depth + broader XCAF-style document semantics | 2-3 months | Partial |
| P7 remaining | persistent history/naming semantics (mutation guard + Builder/Tool now **done**) | 1-2 months | In progress |
| Hardening, edge cases, test coverage | Ongoing | +2 months across all phases | Ongoing |

**Total to reach credible production baseline: approximately 7–11 months** (revised down as P3-P7 baseline items are now largely landed; still bounded by AP242 semantic depth, healing-pipeline depth, and persistent-history semantics).

### What would most accelerate the timeline

1. **Defeaturing robustness hardening** (P3) -- highest remaining ROI in boolean cleanup quality.
2. **Deeper make-connected + glue policies** (P3) -- biggest reduction in near-coincident failure fallout.
3. **AP242 semantic read depth** (P6) -- largest interoperability blocker for production import workflows.
4. **BRepGraph persistent-history semantics** (P7) -- foundation for stable history-driven editing and naming.

## Module-Level Task Breakdown

### libs/rcad-algorithms

| Task | Priority | Status |
|---|---|---|
| General fuse split-first core | P3 - High | Done (baseline split-first general_fuse API + per-object splitter/fuse reporting) |
| Splitter API | P3 - High | Done (baseline split_brep and grouped object/tool variants) |
| Fuzzy / glue boolean options | P3 - High | Done (BooleanOptions glue path + fuzzy analytic coverage baseline) |
| **Result simplification: same-domain unification** | **P3 - High** | **Done (Phase 1: plane/cylinder/cone/torus/sphere; Phase 2: topological+UV validation)** |
| Result simplification: internal face removal | P3 - Medium | Done (Phase 1: threshold-based + same-domain checks; Phase 2: topological true-duplicate detection) |
| CellsBuilder (split-cell graph) | P3 - Medium | **✅ Done (baseline)** |
| MakerVolume (solid from split cells/shells) | P3 - Medium | **Done (`MakerVolume` + `make_solid_from_region` + history variant over reusable cells / `CellExpr`)** |
| Defeaturing pass | P3 - Medium | Done (baseline: cylindrical feature detection + boolean fill/remove and small-face identification) |
| Richer history graph (edges, solids) | P3 - Medium | Done (edge/vertex + aggregated shell/solid origins with persistent labels) |
| ~~SameParameter / SameRange repair~~ | ~~P4 - High~~ | **✅ SameParameter done (P4)** |
| SameRange repair | P4 - High | **✅ Done (scan+repair)** |
| ShapeUpgrade_UnifySameDomain equivalent | P4 - High | Done (Phase 1 baseline: plane/cylinder/cone/torus/sphere; Phase 2 topological guards) |
| Face-on-surface consistency checker | P4 - Medium | **✅ Done (diagnose_face_surface_consistency)** |
| ~~Shell / manifoldness analyzer~~ | ~~P4 - Medium~~ | **✅ Done (P4)** |
| Wire gap / self-intersection analyzer | P4 - Medium | **✅ Done (analyze_wire_issues)** |
| Import analyze/heal pipeline | P4 - Medium | **Done (staged healing report + JSON integration)** |
| Tolerance propagation rules | P4 - Medium | **Done (baseline bottom-up/top-down propagation API)** |
| ~~Small-edge cleanup~~ | ~~P3 - Medium~~ | **✅ Done (P3)** |
| ~~Feature prism / cylindrical hole~~ | ~~P5 - Medium~~ | **✅ Done (P5)** |
| Draft prism feature | P5 - Medium | **✅ Done** |
| Revolution feature | P5 - Medium | **✅ Done** |
| Rib / slot feature | P5 - Low | Done (baseline: make_linear_rib + make_revolution_rib via prism/revolve boolean) |

### libs/rcad-kernel

| Task | Priority | Status |
|---|---|---|
| Compound / CompSolid topology | P6 - Medium | Done (baseline: `Compound` + `CompSolid` structs in topology.rs; add/iter API) |
| Non-manifold topology support | P3 - Low | **Done (`RepairHint` enum + `ManifoldRepairHints` struct; `BRepGraph::edge_valence`, `vertex_degree`, `repair_hints` — classifies StitchablePair / UnmatchedBoundaryEdge / OrphanEdge / MultiManifoldEdge / NonManifoldVertex; 7 tests)** |
| Persistent naming hooks | P6 - Medium | **Done (`PersistentNamingHooks::propagate_through_remap` + `propagate_face_remap` + `identity_map` + `iter`; empty-map = passthrough semantics; 5 propagation tests)** |
| Richer validity analysis | P4 - Medium | Done (baseline: Euler characteristic + genus + orientation consistency + RicherValidityReport) |
| Tolerance propagation rules | P4 - High | Done (kernel write API: `set/update_vertex/edge/face_tolerance`, `finalize_tolerance_hierarchy`, `resize_tolerance_arrays`) |
| **Graph topology wrapper (BRepGraph equivalent)** | **P7 - High** | **Done (O(1) adjacency + DFS/BFS + dirty tracking; `BRepGraph::from_brep`, iterators exported; `BRepGraphMutGuard` RAII + `BRepGraphCheckpoint` + `validate_invariants`; `BRepGraphBuilder` + `BRepGraphTool`; 35 graph tests total)** |

### libs/rcad-modeling

| Task | Priority | Status |
|---|---|---|
| **N-sided surface fill -- Gordon N×M transfinite** | **P7 - Medium** | **Done (baseline Gordon surface evaluator in rcad-kernel/rcad-modeling)** |
| **Circle involute curve (gear tooth)** | **P7 - Medium** | **✅ Done** |
| **Circular helix curve** | **P7 - Low** | **✅ Done** |
| Normal projection of wire onto surface | P3 - Low | **Done (`project_wire_onto_surface` in `wire_ops.rs`: projects Wire vertices via `closest_point_on_surface`, reconnects with Line3 edges; tests cover identity projection on plane + z-elevation drop)** |
| Stabilize advanced fillet corner cases | P3 - Medium | **Done (`fillet_edge_safe` + `chamfer_edge_safe`: auto-clamp radius/dist to `0.49 × min(edge_len, shortest_adj_edge)`, `SafeFilletResult` carries `was_clamped` flag; 6 tests)** |
| Angle-mode chamfer | P3 - Low | **Done (`chamfer_edge_angle` + `chamfer_edge_angle_with_history`; setback formula sb1 = dist·sin(α)/sin(β−α); 9 unit tests covering asymmetric setbacks, 45° symmetry, error cases)** |
| 2-D fillet/chamfer API | Convenience | **Done (`fillet_wire_2d` + `chamfer_wire_2d` in `wire_ops.rs`: shared `round_corners_2d` core; setback = param/tan(θ/2) with per-corner fallthrough when edge too short; 8 polygon tests)** |

### libs/rcad-step

| Task | Priority | Status |
|---|---|---|
| AP242 full read (import) | P6 - High | **Done (baseline metadata-chain import scope includes PDR/DimLoc/DimSize/GeomTol/Datum/DatumSystem/KinematicPair extraction)** |
| GDT / DimTol read and write | P6 - Medium | **Done (baseline DIMENSIONAL_LOCATION / DIMENSIONAL_SIZE / GEOMETRIC_TOLERANCE + DATUM + GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE read/write metadata support)** |
| ~~Material / layer read and write~~ | ~~P6 - Low~~ | **✅ Done (P6)** |
| STEP general property export (arbitrary metadata) | P6 - Low | **Done (baseline `GENERAL_PROPERTY` + `PROPERTY_DEFINITION` linkage)** |
| Stream-based read/write (DE_Wrapper style) | P6 - Low | **Done (`StepReader::parse_reader*` + `StepWriter::write_to*` stream APIs over `Read`/`Write`, with stream round-trip tests)** |
| Kinematics read | P6 - Low | **Done (baseline metadata extraction for AP242 pair entities + AP242 metadata writer emit path)** |
| Import healing pipeline integration | P4 - Medium | **Done (staged healing report + JSON wiring)** |
| Stronger PCurve / tolerance validation on export | P4 - Medium | **Done (`validate_export_readiness` in rcad-step: PCurve index-bounds + cardinality + missing-PCurve + tolerance-floor checks; `ExportReadinessReport` with `summary()`)** |

### libs/rcad-render

| Task | Priority | Status |
|---|---|---|
| Tunable meshing deflection / angular tolerances | P3 - Medium | **Done (`TessellationOptions = TessellationParams` re-export + `Tessellator::tessellate_with_options` in rcad-render; wires `mesh_brep` with chord/angle params to GPU mesh builder)** |
| Incremental cache invalidation for edited models | P2.1 followup - Medium | **Done (`EditedModelDelta` + `Tessellator::invalidate_cache_for_edits` + `tessellate_after_edits`; adjacency-driven face invalidation via `BRepGraph` for vertex/edge/face edits)** |
| Better HLR on curved geometry | P2.2 followup - Low | Not started |

## Recommended Non-Goals for Now

- **Helix geometry (TKHelix)**: baseline evaluators are now present in RCAD; a full dedicated toolkit equivalent remains a later-phase item.
- VRML / PLY / glTF output -- low engineering value; focus on STEP AP242 depth first.
- Full FEA entity round-trip via AP209 -- too specialized.
- Comprehensive IGES B-Rep support -- STEP AP242 is superior for the same use cases.
- More geometry primitive types -- the current set is sufficient.
- Premature rendering optimization -- wait until boolean/healing/document layers improve.
- **BRepGraph full implementation** -- important long-term, but implementing 49 000+ lines of equivalent code is a multi-quarter project best done after P3–P6 remaining items are closed.

## Bottom Line

After P3–P7 baseline implementation, RCAD has meaningfully advanced its healing, feature, graph traversal, and metadata coverage. The situation as of April 2026:

- The geometric foundation is solid and competitive.
- The feature toolbox covers the most common operations (extrude, revolve, sweep, fillet, boolean, prism, cylindrical hole, split-face, baseline rib/slot).
- The weak points remain **boolean robustness hardening**, **healing depth**, and **full AP242 exchange completeness**.
- OCCT 8.0.0 (imminent) still raises the bar further: deep BRepGraph mutation/history, production Gordon tooling, and broader evaluation geometry remain significant moving targets.

The three areas that dominate the remaining production gap, in order:

1. **Boolean robustness and result simplification hardening** -- splitter/fuzzy/same-domain baselines exist; robustness and defeaturing remain the biggest gap.
2. **Healing depth** -- staged pipeline exists, but ShapeProcess-like coverage and richer repair semantics remain.
3. **BRepGraph / graph topology depth** -- RAII mutation guard (`BRepGraphMutGuard`), checkpoint/rollback, `validate_invariants`, `BRepGraphBuilder`, and `BRepGraphTool` are now in place; persistent-history naming remains.

At the current trajectory, RCAD can become a credible industrial CAD kernel with **8–12 months of focused kernel work**, assuming the priority order above and one dedicated developer. If BRepGraph and Gordon/eval-breadth work are deferred to a later phase, the immediate production baseline can be reached in roughly **5–7 months** (same-domain unification + healing pipeline completion + AP242 read depth).

