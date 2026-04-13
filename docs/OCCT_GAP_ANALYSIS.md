# RCAD vs OCCT Gap Analysis

Date: 2026-04-13 (revised against OCCT 7.9.3 / 8.0.0-rc5)

## Purpose

This document turns the current RCAD capability inventory into an OCCT-oriented gap analysis and an execution roadmap. The goal is not to chase feature-count parity blindly, but to identify the areas where RCAD still differs most from OCCT in practical engineering use.

The key conclusion is:

RCAD already covers a large portion of the geometric foundation layer, and recent work (P1.4-P2.3) has materially improved variable fillet, STEP AP242 output, HLR silhouette quality, incremental mesh caching, and analytic distance. However, RCAD still trails OCCT most in:

- boolean robustness and architectural breadth
- healing and tolerance management
- richer topology and document model
- industrial-grade data exchange depth
- advanced post-processing and result simplification

## OCCT Version Reference

| Metric | Value |
|---|---|
| Latest stable | **7.9.3** (December 2025) |
| Master (pre-release) | **8.0.0-rc5** (April 2026) |
| New in 8.0.0 | TKHelix toolkit, TKDE unified DE framework, BRepGraph graph API, BRepAlgoAPI_Defeaturing, Handle API modernization, Windows ARM64 support |
| Toolkit count | ~57 toolkits across 7 modules |

## Executive Summary

### Areas where RCAD is already strong (updated after P1.4-P2.3)

- Broad analytic and B-spline geometry coverage.
- Core B-Rep types: Vertex, Edge, Wire, Face, Shell, Solid.
- Primitive creation, extrude, revolve, sweep multi-section, loft, variable-radius fillet, chamfer, thicken, draft, mirror, and array.
- STEP AP214/AP242 write, assembly IO, color attributes, OBJ IO, IGES mesh bridge.
- HLR with dense silhouette sampling on curved surfaces (P2.2).
- Analytic extrema/distance fast paths for sphere-sphere, plane-sphere, parallel planes (P2.3).
- Curvature, sectioning, shape properties, projections.
- mesh_dirty incremental mesh caching (P2.1).

### Areas where RCAD is still behind OCCT

- Boolean robustness on curved solids and near-coincident geometry is the largest practical gap.
- OCCT 8.0.0 adds BRepAlgoAPI_Defeaturing (remove interior features) and BOPAlgo_MakeConnected that RCAD has no equivalent for.
- Healing and tolerance pipelines are basic compared to OCCT's 10-package TKShHealing stack.
- Post-operation simplification (same-domain unification, small-face/edge cleanup) is missing.
- Exchange: STEP AP242 write is present but AP242 read, kinematics, FEA, and GDT round-trips are not.
- No Compound / CompSolid topology containers.
- No feature library (boss, pocket, rib, hole -- OCCT's TKFeat).
- No N-sided surface filling (MakeFilling) or medial-axis extraction.
- No persistent naming or document metadata beyond colors and assembly placement.

## Capability Matrix

| Domain | RCAD today | OCCT reference | Gap level | Recommended direction |
|---|---|---|---|---|
| Geometry kernel | Strong analytic + spline coverage | TKG2d / TKG3d / TKGeomBase | Low | Maintain correctness; add helix curves if needed |
| Core B-Rep model | Vertex/Edge/Wire/Face/Shell/Solid | TKBRep | Medium | Add Compound, CompSolid, non-manifold support |
| Primitive and sweep modeling | Good breadth (extrude/revolve/loft/sweep) | TKPrim + TKOffset | Low | Improve edge cases; add N-sided fill, normal projection |
| Fillet / chamfer | Variable-radius fillet, chamfer | TKFillet | Medium | Angle-mode chamfer, 2-D fillet API, more corner cases |
| Thicken / draft / offset | Present | TKOffset BRepOffset | Medium | Shell offset (MakeOffsetShape), evolved surface |
| Feature library | None | TKFeat (boss/pocket/rib/hole) | High | Implement: prism/revolution feature, rib, hole |
| Boolean framework | Fuse/Cut/Common/Section + imprint | TKBO | High | CellsBuilder, Splitter, fuzzy, glue, defeaturing |
| Post-op simplification | None | ShapeUpgrade_UnifySameDomain, BOPAlgo cleanup | High | Same-domain unification, internal-face removal |
| Healing and validation | Basic repair + check | TKShHealing (10 packages) | High | ShapeFix-style pipeline, SameParameter/SameRange repair |
| Topology history | Face-level history | BOPAlgo history, OCAF naming | Medium | Extend to edges/vertices/solids; add persistent naming |
| STEP exchange: write | AP214 + AP242 basic | TKDESTEP STEPCAFControl | Medium | Richer AP242 metadata, better PCurve/tolerance validation |
| STEP exchange: read | Basic import | TKDESTEP + STEPCAFControl_Reader | High | Full AP242 read, GDT, kinematics, FEA entities |
| IGES exchange | Mesh bridge only | TKDEIGES | High | Add analytic/B-Rep IGES or document as non-goal |
| Assembly / document model | Colors + shape tree | TKXCAF (XCAFDoc_*) | Medium | Materials, layers, DimTol, notes, persistent attributes |
| Meshing and visualization | mesh_dirty caching (P2.1), HLR dense silhouettes (P2.2) | TKMesh + TKHLR | Medium | Tunable deflection/angular tolerances, incremental remesh |

## Highest-Priority Gaps

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
| BRepAlgoAPI_Defeaturing | Remove interior features (new 8.0) | Missing |
| BOPAlgo_MakeConnected | Connect disconnected geometry (new 8.0) | Missing |
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

## 2. Healing and Tolerance System (gap: High)

OCCT's TKShHealing comprises 10 packages. RCAD's healing module (brep_repair.rs + healing.rs) covers roughly the work of one sub-package.

| OCCT package | Key capability | RCAD status |
|---|---|---|
| ShapeFix_Face | Repair degenerated/invalid faces | Partial in brep_repair.rs |
| ShapeFix_Wire | Reorder wires, close gaps, remove degenerate edges | Partial |
| ShapeFix_Edge | SameParameter, SameRange, degenerated edges | Missing |
| ShapeFix_Shell | Repair shell orientation, manifoldness | Missing |
| ShapeFix_Solid | Solid closure, shell orientation | Missing |
| ShapeAnalysis_Surface | UV consistency, surface bounds analysis | Missing |
| ShapeAnalysis_Wire | Wire gap, self-intersection, area | Missing |
| ShapeUpgrade_UnifySameDomain | Merge co-planar/co-cylindrical faces | Missing (critical for boolean results) |
| ShapeCustom | BSpline restriction, convert to indirect | Missing |
| ShapeProcess | Batch pipeline with operator chain | Missing |

### What RCAD should add

- SameParameter repair: verify and fix edge 3D curve vs PCurve consistency.
- SameRange repair: fix mismatched curve parameter ranges across edges.
- ShapeUpgrade_UnifySameDomain equivalent: merges adjacent co-planar or co-cylindrical faces -- this alone cleans up most boolean result artifacts.
- Import healing pipeline: analyze -> diagnose -> heal -> report on imported STEP data.
- Tolerance propagation rules after boolean/split/sew operations.

## 3. Feature Library (gap: High)

OCCT's TKFeat provides parametric feature operations built on top of boolean operations. RCAD has no equivalent.

| OCCT class | Feature type | RCAD status |
|---|---|---|
| BRepFeat_MakePrism | Boss/pocket (blind/through/up-to) | Missing |
| BRepFeat_MakeDPrism | Draft prism | Missing |
| BRepFeat_MakeRevol | Revolution feature | Missing |
| BRepFeat_MakeCylindricalHole | Cylindrical hole | Missing |
| BRepFeat_MakeLinearForm | Rib/slot (linear) | Missing |
| BRepFeat_MakeRevolutionForm | Rib/slot (revolved) | Missing |
| BRepFeat_Gluer | Glue shapes at interface | Missing |
| BRepFeat_SplitShape | Split face by wire | Missing |

Feature operations are heavily used in mechanical design workflows. Without them, applications must compose boolean operations manually, which is fragile and history-unfriendly.

## 4. Document Model and Exchange Depth (gap: Medium-High)

OCCT's XDE/XCAF layer provides a structured document model far deeper than RCAD's current color+assembly support.

| OCCT attribute | Purpose | RCAD status |
|---|---|---|
| XCAFDoc_Color / ColorTool | Per-face/part colors | Present (step IO) |
| XCAFDoc_ShapeTool | Assembly tree | Present |
| XCAFDoc_Material / MaterialTool | Material assignment | Missing |
| XCAFDoc_Layer / LayerTool | Layer/group membership | Missing |
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
| AP242 GDT write | Missing |
| AP242 read (import) | Basic |
| Kinematics read/write | Missing |
| FEA entity read | Missing |

## 5. Post-Operation Simplification (gap: High)

OCCT's ShapeUpgrade_UnifySameDomain and BOPAlgo_Builder cleanup passes automatically clean up boolean results.

RCAD currently produces boolean results that may contain:
- many small adjacent co-planar faces where one merged face would suffice
- dangling internal faces after fuse operations
- tiny edges from near-coincident intersections
- mismatched tolerances at operation boundaries

Producing clean results without this layer means downstream consumers (meshing, rendering, exchange) see unnecessarily complex topology.

## Roadmap

### P3: Boolean Robustness and Result Simplification (highest priority)

Target: imported and modeled solids through boolean workflows with much better robustness and cleaner results.

Estimated effort: 2-3 months of focused work.

| Deliverable | Effort |
|---|---|
| Splitter API (object/tool split) | Medium |
| Fuzzy tolerance option | Small |
| Same-domain face unification | Medium |
| Internal-face removal after fuse | Small |
| Small-edge / small-face cleanup | Medium |
| CellsBuilder (split-cell graph) | Large |
| Defeaturing pass | Large |
| Full history to edges and solids | Medium |

### P4: Industrial Healing Pipeline

Target: imported CAD data can be analyzed, repaired, and pushed into modeling/boolean workflows safely.

Estimated effort: 3-4 months.

| Deliverable | Effort |
|---|---|
| SameParameter and SameRange repair | Medium |
| ShapeUpgrade_UnifySameDomain equivalent | Medium |
| Shell closure and manifoldness analyzer | Small |
| Face-on-surface consistency checker | Medium |
| Import analyze/heal/report pipeline | Medium |
| Tolerance propagation after boolean/sew | Large |

### P5: Feature Library

Target: parametric boss/pocket/rib/hole features built on top of boolean + history operations.

Estimated effort: 2-3 months.

| Deliverable | Effort |
|---|---|
| MakePrism (blind/through/up-to) | Medium |
| MakeDPrism (draft prism) | Medium |
| MakeCylindricalHole | Small |
| MakeLinearForm (rib/slot) | Large |
| SplitShape | Small |

### P6: Document Model and AP242 Depth

Target: XCAF-comparable document model; AP242 GDT and kinematics round-trip.

Estimated effort: 2-3 months.

| Deliverable | Effort |
|---|---|
| Material attribute | Small |
| Layer attribute | Small |
| GDT / DimTol write | Medium |
| AP242 read (full) | Large |
| Persistent naming in history | Large |

## Production Readiness Gap Assessment

Answering the question: **how far is RCAD from being production-grade?**

### Definition used here

"Production-grade" = a downstream application (CAE preprocessor, manufacturing toolpath planner, or PDM system) can rely on RCAD for:
1. Importing real-world STEP files without manual intervention.
2. Running boolean and feature operations on those imports reliably.
3. Exporting results with correct metadata.
4. Supporting history-based editing and re-feature.

### Gap summary by category

| Category | Current state | Production bar | Delta |
|---|---|---|---|
| Geometry kernel | 90% | 95% | Minor gaps (helix, some surface types) |
| B-Rep model | 75% | 90% | Compound/CompSolid, non-manifold |
| Sweep / loft / extrude | 80% | 90% | N-sided fill, medial axis |
| Fillet / chamfer | 70% | 85% | Corner cases, angle-mode chamfer |
| Boolean core | 55% | 85% | Fuzzy, splitter, cleanup, defeaturing |
| Healing pipeline | 20% | 80% | SameParameter, ShapeUpgrade, analyzers |
| Feature library | 0% | 60% | Still useful but not blocking for all apps |
| Document / XCAF | 30% | 75% | Materials, layers, GDT, persistent naming |
| STEP exchange depth | 50% | 80% | AP242 read, GDT round-trip |
| Meshing controls | 40% | 75% | Tunable deflection, incremental update |

### Estimated work remaining

Assuming 1 developer working full-time on kernel work:

| Phase | Focus | Calendar estimate |
|---|---|---|
| P3 | Boolean robustness + result simplification | 2-3 months |
| P4 | Healing pipeline + import robustness | 3-4 months |
| P5 | Feature library (boss/pocket/rib/hole) | 2-3 months |
| P6 | Document model + AP242 round-trip depth | 2-3 months |
| Hardening, edge cases, test coverage | Ongoing | +2 months across all phases |

**Total to reach credible production baseline: approximately 11-15 months.**

This estimate assumes no fundamental algorithmic gaps in the current geometry kernel -- the main work is architecture breadth and robustness, not from-scratch algorithm research.

### What would most accelerate the timeline

1. Prioritizing **same-domain unification** (P3) early -- it immediately makes all existing boolean results cleaner and makes P4/P5 easier.
2. Implementing **fuzzy boolean tolerance** (P3) -- removes the most common class of near-coincident import failures.
3. Starting a **healing framework scaffold** (P4) even before full SameParameter support -- an analyze/diagnose/report pipeline is usable before all fixes exist.
4. Deferring the full **feature library** (P5) if the target application does not require history-based parametric editing.

## Module-Level Task Breakdown

### libs/rcad-algorithms

| Task | Priority |
|---|---|
| General fuse split-first core | P3 - High |
| Splitter API | P3 - High |
| Fuzzy / glue boolean options | P3 - High |
| Result simplification: same-domain unification | P3 - High |
| Result simplification: internal face removal | P3 - Medium |
| CellsBuilder (split-cell graph) | P3 - Medium |
| Defeaturing pass | P3 - Medium |
| Richer history graph (edges, solids) | P3 - Medium |
| SameParameter / SameRange repair | P4 - High |
| ShapeUpgrade_UnifySameDomain equivalent | P4 - High |
| Shell / manifoldness / surface analyzers | P4 - Medium |
| Import analyze/heal pipeline | P4 - Medium |
| Tolerance propagation rules | P4 - Medium |
| Feature prism / revolution | P5 - Medium |
| Cylindrical hole feature | P5 - Low |
| Rib / slot feature | P5 - Low |

### libs/rcad-kernel

| Task | Priority |
|---|---|
| Compound / CompSolid topology | P6 - Medium |
| Non-manifold topology support | P3 - Low |
| Persistent naming hooks | P6 - Medium |
| Richer validity analysis | P4 - Medium |
| Tolerance propagation rules | P4 - High |

### libs/rcad-modeling

| Task | Priority |
|---|---|
| N-sided surface fill (MakeFilling equivalent) | P3 - Low |
| Normal projection of wire onto surface | P3 - Low |
| Stabilize advanced fillet corner cases | P3 - Medium |
| Angle-mode chamfer | P3 - Low |
| 2-D fillet/chamfer API | Convenience |

### libs/rcad-step

| Task | Priority |
|---|---|
| AP242 full read (import) | P6 - High |
| GDT / DimTol read and write | P6 - Medium |
| Material / layer read and write | P6 - Low |
| Kinematics read | P6 - Low |
| Import healing pipeline integration | P4 - Medium |
| Stronger PCurve / tolerance validation on export | P4 - Medium |

### libs/rcad-render

| Task | Priority |
|---|---|
| Tunable meshing deflection / angular tolerances | P3 - Medium |
| Incremental cache invalidation for edited models | P2.1 followup - Medium |
| Better HLR on curved geometry | P2.2 followup - Low |

## Recommended Non-Goals for Now

- Helix geometry (TKHelix) -- only needed for spring/coil modeling; low priority.
- VRML / PLY / glTF output -- low engineering value; focus on STEP AP242 depth first.
- Full FEA entity round-trip via AP209 -- too specialized.
- Comprehensive IGES B-Rep support -- STEP AP242 is superior for the same use cases.
- More geometry primitive types -- the current set is sufficient.
- Premature rendering optimization -- wait until boolean/healing/document layers improve.

## Bottom Line

RCAD is no longer missing the basic geometric layer. The situation after P1.4-P2.3 is:

- The foundation is solid.
- The toolbox of operations is competitive at the feature level.
- The weak points are **depth and robustness**, not breadth.

The three areas that dominate the remaining production gap, in order:

1. **Boolean robustness and result simplification** -- affects every non-trivial workflow.
2. **Healing and import robustness** -- affects every production import scenario.
3. **Document model and AP242 depth** -- affects every downstream consumer of RCAD output.

At the current trajectory, RCAD can become a credible industrial CAD kernel with **11-15 months of focused kernel work**, assuming the priority order above and one dedicated developer. If the target application scope is narrower (e.g., new-model-only, no import), the required work drops to roughly **P3 + partial P4 = 5-7 months**.
