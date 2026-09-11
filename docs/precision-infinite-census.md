# Precision::Infinite coordinated switch — census record (E3-T)

Status: **LANDED** (kernel producers + consumers by agent A6; rcad-algo/rcad-brep
producers + Group A/B consumers by agent A7; three anchor misreads found by the
debug gate run and fixed 1:1 — see §6). This document records the census the
switch was planned from, the landing, and the residual-noise classification.

## 1. The three OCCT infinite-value layers (do not conflate)

| Layer | Value | Source |
|---|---|---|
| `Precision::Infinite()` | **2e100** | Precision.hxx L371 |
| `Precision::IsInfinite(R)` / `IsPositiveInfinite` / `IsNegativeInfinite` | `\|R\| >= 1e100` / `R >= 1e100` / `R <= -1e100` | Precision.hxx L350-367 |
| Package-local constants | Bnd `THE_BND_PRECISION_INFINITE = 1e100` (Bnd_Box.cxx L27, Bnd_Box2d.cxx L27); `RealFirst()/RealLast()` = ±f64::MAX | per package |

rcad infrastructure: `rcad-kernel/src/core/precision.rs` — `INFINITE_VALUE = 2e100`,
`is_infinite_value` (`|R| >= 1e100`), `is_positive_infinite_value`, `is_negative_infinite_value`,
`REAL_FIRST/REAL_LAST = ±f64::MAX`, `BND_PRECISION_INFINITE = 1e100`.

## 2. Producers switched to ±INFINITE_VALUE (OCCT anchors cited)

Kernel (geom/eval.rs): Line3 (Geom_Line.cxx L137/144), Hyperbola3 (L94/101),
Parabola3 (L114/128 — was a `[-1e4,1e4]` stand-in), Plane (L181-184),
CylindricalSurface (L162-163), **ConicalSurface (L212-213 — V1 corrected
0.0 → −2e100)**, LinearExtrusionSurface (L142-145); Curve2dEval trait default
(geom/mod.rs); Line2d (Geom2d_Line.cxx L144/151), Parabola2d, Hyperbola2d.

rcad-algo: brep_fill_axe.rs:41 (BRepFill.cxx L726/L844),
proj_lib_h_comp_projected_curve_b.rs:368-376 (ProjLib_CompProjectedCurve.cxx
L331-341), brep_blend_walking.rs:869/874 + brep_blend_walking_b.rs:1337 +
chfi3d.rs:2898-2899 + chfi3d_filbuilder_c3.rs:972-973 (ChFiDS_ElSpine.cxx
L42-43), chfi2d_builder_0.rs:49 (Geom_Line), chfi3d_builder_0.rs:638-645
(plane/cyl/cone natural bounds; cone V1 = −2e100),
brep_fill_sweep_b.rs:178 (Geom2d_Line.cxx L142-151), bop/ds/mod.rs
face_actual_uv_bounds fallbacks, bean_face_intersector.rs ExtremaExtCS u/v init
(Extrema_GenExtCS.cxx Initialize chain).

RealLast trio (bop/int_tools/face_make_curve.rs:1115/1128/1141): OCCT
GeomInt_LineConstructor.cxx L940-976 uses exact `RealLast()` equality →
`REAL_LAST`, never INFINITE_VALUE.

## 3. Consumers flipped to precision predicates

Group A (exact IEEE `==`/`!=` compares, all zero-residual): tool_rehost.rs
Plane() (GeomAdaptor_SurfaceOfLinearExtrusion.cxx L341-352), hlr/contap/
surface_adaptor.rs converters (output STAYS ±f64::MAX per Contap_HContTool.cxx
L136-162), brep_offset_inter2d_b.rs ×8, brep_offset_make_offset.rs (dead
`to_occt_infinite` deleted), fclass2d.rs PerformInfinitePoint
(IntTools_FClass2d.cxx L627 `== RealLast()`), face_make_curve.rs:1141.

Group B (IEEE is_finite/is_infinite predicates on domain values): the ElCLib
guard quartet (pave_filler.rs, shape_build/edge.rs, chfi2d_builder.rs,
chfi2d_ana_fillet_algo.rs — ElCLib.cxx L121-128), bop context.rs
(IntTools_Context.cxx L818/L875), ds/mod.rs, builder.rs surface_is_closed ×7 +
domain-center fallbacks ×3 (GeomLib.cxx L2717-2813), face_explorer.rs
(BRepClass_FaceExplorer.cxx L151-166), make_blocks.rs, brep_sweep_builder.rs
(BRep_Builder.cxx L1224-1227), tool_rehost.rs, chfi3d.rs, walking_b.rs
(BRepBlend_Walking.cxx L2198/L2259), approx_curve_on_surface.rs,
extrema_gen_ext_cs.rs clamp, int_quad_quad.rs bounds branch, w_line_tool.rs,
brep_offset_inter2d.rs, make_offset_e.rs (BRepOffset_MakeOffset.cxx
L5475/L5492), shape_analysis/mod.rs ×7, brep_check/e1.rs, solid_explorer.rs
(SolidExplorer.cxx L454-467), g_inter.rs, bnd_box_tree.rs, geomfill/fixed.rs
(anchor verified: GeomFill_Fixed.cxx L104-105 = Precision::Infinite),
rcad-brep/adaptor.rs predicate form.

## 4. Deliberate exceptions (verified against anchors, NOT flipped)

- int_quad_quad.rs `is_first_open/is_last_open`: OCCT IntAna_Curve.cxx L251-261
  returns the STORED `firstbounded/lastbounded` flags, not infinity tests —
  rcad keeps `!param.is_finite()` (infinite bound ⟺ open side), the faithful
  encoding for bounded curves.
- int_quad_quad.rs TrigFunctionRoots infeasible-domain test: OCCT
  math_TrigonometricFunctionRoots.cxx L96-108 compares against
  `RealFirst()/RealLast()` literals (±f64::MAX), not the 1e100 predicate.
- geom_lprop REAL_LAST/REAL_FIRST uses, bnd_lib, distance sentinels, BBox
  accumulators, 3D-point finiteness checks: IEEE by design (noise class).

## 5. Residual-noise classes (left IEEE, verified harmless)

BBox accumulator seeds (bnd_lib ×32, gprop, algo_tools, tolerance, brep_repair),
best-distance sentinels (extrema.rs, gprop/tri.rs, bvh, core/color.rs,
math/curvature.rs), approx_int tol sentinels, standard_epsilon branches,
3D/UV point finiteness checks, point_line.rs (`1/INF=0`).

## 6. Debug-gate corrections (2026-09-11, after the first full-gate run)

The switch exposed three anchor misreads in its own consumer flips, each fixed
1:1 (26 gate failures → green):

1. `extrema_locate_ext_pc.rs` (Perform's not-found interval arm): OCCT
   Extrema_GLocateExtPC.hxx L129/L153-160/L166-168 keeps loop-carried
   `myintuinf/myintusup` for the final try — rcad re-indexed the degraded
   `inter` and panicked (index out of bounds), starving paving of one WLine in
   cone/sphere pairs (pf_stage 9-vs-10 counts, 15 builder_stage failures).
2. `int_quad_quad.rs`: the A7 first-flip of `is_first_open/is_last_open` to
   infinity predicates was polarity-inverted vs IntAna_Curve.cxx L251-261
   (stored bounded flags); restored to the `!is_finite()` flag encoding, and
   the TrigFunctionRoots infeasible-domain test to RealFirst/RealLast literals
   per math_TrigonometricFunctionRoots.cxx L96-108.
3. `pave_filler_make_blocks.rs` PutClosingPaveOnCurve: De Morgan inversion of
   OCCT BOPAlgo_PaveFiller_6.cxx L3513-3516 `if (!aIC.HasBounds()) return;` —
   bounded curves (every circle) were being rejected; restored
   `if is_infinite_value(t0) || is_infinite_value(t1) { return; }`.

## 7. Completeness greps (must stay at the classified-noise baseline)

```
grep -rnE "[!=]= (f64::|std::f64::)?(NEG_)?INFINITY" libs/rcad-kernel/src libs/rcad-algo/src libs/rcad-brep/src --include="*.rs" | grep -v to_bits
grep -rn "is_infinite()\|is_finite()" libs/rcad-algo/src --include="*.rs" | grep -vE "\.x\.|\.y\.|\.z\."
```
Post-landing residuals are §5 noise only.
