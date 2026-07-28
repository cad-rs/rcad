# TKMath 模块重构交接

## 一句话提示词（在新 session 中直接使用）

```
继续 rcad-kernel TKMath 重构：创建 math/int_ana/（从 rcad-algorithms/src/int_ana.rs 迁移解析求交算法）和 math/extrema/（从 base/geom_api/project.rs 提取 Extrema_ExtPC/Extrema_ExtPS 的低层极值逻辑），使模块对齐 OCCT TKMath 的 IntAna 和 Extrema 包。
```

## 已完成

### `rcad-kernel/src/math/` — OCCT TKMath（FoundationClasses）

| 模块 | OCCT 包 | 说明 |
|------|---------|------|
| `bnd/` | Bnd | BndBox struct + curve/surface bbox |
| `bspl` | BSplCLib + BSplSLib | de_boor, bspline_tangent, find_knot_span |
| `bvh/` | BVH | Bvh struct (SAH, ray_cast, candidate_pairs) |
| `convert/` | Convert | sphere/cylinder/plane → BSpline |
| `cs_lib` | CSLib | Surface normal computation (from derivatives, singular-case rescue) |
| `el/` | ElCLib + ElSLib | elementary curve/surface analytic eval |
| `gcpnts/` | GCPnts | arc_length (AbscissaPoint) |
| `gprop/` | GProp | surface_area, volume, inertia, plate |
| `math_poly` | MathPoly | solve_linear~quartic, laguerre_roots |
| `plib` | PLib | de_casteljau, eval_polynomial, power→Bezier |
| `poly/` | Poly | Triangulation, Connect (mesh data) |
| `top_loc` | TopLoc | TopLoc struct with chaining |
| `root/` | MathRoot | newton_raphson, bisection, trig_roots |
| `opt/` | MathOpt | BFGS, FRPR, NewtonMin, Powell, etc. |
| `lin/` | MathLin | SVD, Gauss, Crout, eigenvalues |
| `integ/` | MathInteg | simpson, gaussian_quadrature |
| `sys/` | MathSys | newton_2d, newton_3d |
| `curvature` | LProp | principal_curvatures, gaussian/mean |

### `rcad-kernel/src/base/` — OCCT TKGeomBase（ModelingData）

| 模块 | OCCT 包 | 说明 |
|------|---------|------|
| `int_ana` | IntAna | Analytic line/plane/cylinder/sphere/cone/torus intersection |
| `extrema` | Extrema | Point-curve, point-surface, and curve-curve extrema |
| `geom_api/` | GeomAPI | Projection, Interpolate, Extrema (GeomAPI-level wrappers) |
| `extend/` | — | curve/surface extension |

注：`math/` 中的 `int_ana` 和 `extrema` 现在是转发到 `base/` 的桩模块，保持编译兼容。

## 未完成（待迁移）

### 1. PointSetLib（低优先级）

OCCT 较新加的包，点云 PCA 分析。rcad 目前无对应实现，可后期再加。

## 当前 git 状态

在 `rcad` 子模块的 `main` 分支，最新提交在 `ddf8d28`。所有 parent repo 的 submodule 指针已更新。
