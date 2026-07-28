# rcad-kernel TKMath/TKGeomBase 模块重构交接

## 一句话提示词（在新 session 中直接使用）

```
继续 rcad-kernel TKGeomBase 补全：创建 base/gc/（从 geom/ 提取 GC_MakeCircle/MakeLine/MakePlane 等构造算法）和 base/intana2d/（2D 解析求交），使模块对齐 OCCT TKGeomBase 的 GC 和 IntAna2d 包。
```

## TKMath 状态（FoundationClasses — 已全部对齐）

`rcad-kernel/src/math/` 所有 OCCT TKMath 包全覆盖：

| 模块 | OCCT 包 | 说明 |
|------|---------|------|
| `bnd/` | Bnd | BndBox struct + curve/surface bbox |
| `bspl` | BSplCLib + BSplSLib | de_boor, bspline_tangent, find_knot_span |
| `bvh/` | BVH | Bvh struct (SAH, ray_cast, candidate_pairs) |
| `convert/` | Convert | sphere/cylinder/plane → BSpline |
| `cs_lib` | CSLib | Surface normal computation (derivatives, singular-case rescue) |
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

## TKGeomBase 状态（ModelingData — 补全中）

`rcad-kernel/src/base/`：

| 模块 | OCCT 包 | 说明 | 优先级 |
|------|---------|------|--------|
| `int_ana` | IntAna | 3D 解析求交（线/平面/柱/球/锥/环） | ✅ 已对齐 |
| `extrema` | Extrema | 点-曲/点-面/曲-曲极值 | ✅ 已对齐 |
| `geom_api/` | GeomAPI | 投影/插值/Extrema 高层封装 | ✅ 已对齐 |
| `extend/` | — | 曲线/曲面延伸 | ⚠️ 局部 |
| `geom/` | — | 几何类型 + SurfaceEval | ✅ 基础 |
| — | **GC** | **几何构造（MakeCircle/MakeLine/MakePlane 等）** | **🔴 待建** |
| — | **IntAna2d** | **2D 解析求交** | **🔴 待建** |
| — | **ProjLib/GeomProjLib** | **曲线-曲面投影** | **🟡 部分** |
| — | **Approx/AppParCurves** | **拟合/逼近** | **🟡 部分** |
| — | **GeomLib** | **几何工具函数** | **🟡 部分** |
| — | GCPnts | 曲线上取点 | ✅ math/gcpnts/ |
| — | GeomConvert | 几何转换 | ✅ math/convert/ |
| — | GProp | 全局属性 | ✅ math/gprop/ |
| — | LProp | 局部属性 | ✅ math/curvature |
| — | BndLib/GeomBndLib | 包围盒 | ✅ math/bnd/ |
| — | CPnts, AdvApp2Var, AppCont, Hermit, FEmTool, GeomTools, GC2d, gce | 其他 | ⬜ 低优先级 |

## 下阶段启动建议

### 推荐：GC（几何构造）

OCCT `GC` 包提供 `GC_MakeCircle`、`GC_MakeLine`、`GC_MakePlane`、`GC_MakeSegment` 等。目前 rcad 的这些功能散落在 `geom/` 模块或 Example 代码中。迁移方式：

1. 创建 `base/gc.rs`（或 `base/gc/` 目录）
2. 从 OCCT `GC_MakeXxx.hxx` 翻译：MakeCircle（圆心+法向+半径/三点/两点）、MakeLine（点+方向/两点）、MakePlane（点+法向/三点）等
3. 搜刮 rcad 现有代码中重复的构造逻辑，统一到 GC 模块

### 次选：IntAna2d

OCCT `IntAna2d` 提供 2D 解析求交（线-线、线-圆、圆-圆、线-椭圆、椭圆-椭圆等）。`rcad-algorithms/src/inttools/` 中已有一些 2D 求交逻辑可提取。

### 更高优先级？

如果布尔对齐方向需要更精确的曲线-曲面投影，优先做 **ProjLib**（`base/proj_lib.rs`），从 `base/geom_api/project.rs` 的 `make_pcurve_on_surface` 提取。

## 当前 git 状态

在 `rcad` 子模块的 `main` 分支，最新提交在 `60e5e6fb`。Parent repo submodule 指针未更新。
