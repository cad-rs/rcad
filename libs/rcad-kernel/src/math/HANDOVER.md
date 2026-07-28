# TKMath 模块重构交接

## 一句话提示词（在新 session 中直接使用）

```
继续 rcad-kernel TKMath 重构：创建 math/int_ana/（从 rcad-algorithms/src/int_ana.rs 迁移解析求交算法）和 math/extrema/（从 base/geom_api/project.rs 提取 Extrema_ExtPC/Extrema_ExtPS 的低层极值逻辑），使模块对齐 OCCT TKMath 的 IntAna 和 Extrema 包。
```

## 已完成

`rcad-kernel/src/math/` 已有以下 OCCT TKMath 模块：

| 模块 | OCCT 包 | 说明 |
|------|---------|------|
| `bnd/` | Bnd | BndBox struct + curve/surface bbox |
| `bspl` | BSplCLib + BSplSLib | de_boor, bspline_tangent, find_knot_span |
| `bvh/` | BVH | Bvh struct (SAH, ray_cast, candidate_pairs) |
| `convert/` | Convert | sphere/cylinder/plane → BSpline |
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

`base/geom_api/` 覆盖 GeomAPI（TKGeomBase）。

## 未完成（待迁移）

### 1. IntAna — 解析求交

OCCT TKMath 的 `IntAna` 包提供：
- `IntAna_QuadQuad` — 二次曲面-二次曲面求交
- `IntAna_QuadQuadGeo` — 二次曲面-二次曲面几何解
- `IntAna_IntLinTorus` — 直线-圆环面求交
- `IntAna_IntLinCylinder` / `IntAna_IntLinSphere` 等

rcad 当前实现在：`rcad-algorithms/src/int_ana.rs`
引用路径：`crate::int_ana::intersect_line_*`

迁移方式：
1. 创建 `rcad-kernel/src/math/int_ana.rs`
2. 将 `rcad-algorithms/src/int_ana.rs` 中的函数移到 kernel
3. 在 `rcad-algorithms/src/int_ana.rs` 中创建转发 `pub use rcad_kernel::math::int_ana::*;`
4. 更新所有 `crate::int_ana::` → `rcad_kernel::math::int_ana::`（或通过转发保持兼容）

### 2. Extrema — 低层极值算法

OCCT TKMath 的 `Extrema` 包提供：
- `Extrema_ExtCC` — 曲线-曲线极值（当前在 `base/geom_api/extrema.rs`）
- `Extrema_ExtPC` — 点-曲线极值（嵌入在 `base/geom_api/project.rs` 的 `closest_point_on_curve`）
- `Extrema_ExtPS` — 点-曲面极值（嵌入在 `base/geom_api/project.rs` 的 `numeric_surface_projection`）
- `Extrema_ExtPElC` — 点到解析曲线（Line/Circle/Ellipse 的解析路径在 `project.rs` 里）

当前这些逻辑都嵌入在 `base/geom_api/` 的 GeomAPI 封装中，没有独立的低层模块。

迁移方式：
1. 创建 `rcad-kernel/src/math/extrema.rs`
2. 从 `base/geom_api/project.rs` 提取 `Extrema_ExtPElC` 的解析逻辑（Line/Circle/Ellipse 的解析最近点）
3. 从 `base/geom_api/project.rs` 提取 `Extrema_ExtPS` 的数值网格+Newton 细化逻辑
4. `base/geom_api/extrema.rs` 的 `Extrema_ExtCC` 保持不变或提取底层
5. 在 `base/geom_api/project.rs` 中改为调用 `math/extrema.rs`

### 3. PointSetLib（低优先级）

OCCT 较新加的包，点云 PCA 分析。rcad 目前无对应实现，可后期再加。

## 当前 git 状态

在 `rcad` 子模块的 `main` 分支，最新提交在 `ddf8d28`。所有 parent repo 的 submodule 指针已更新。
