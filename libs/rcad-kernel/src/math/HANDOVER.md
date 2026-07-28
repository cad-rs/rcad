# rcad-kernel TKMath/TKGeomBase 模块重构交接

## 一句话提示词（在新 session 中直接使用）

```
继续 rcad 布尔管线形式对齐：从 `BOPAlgo_Builder::PerformInternal1` 顶层开始，逐层递归对齐 OCCT TKBO 源码到 rcad `build_with_history`。当前已对齐到 `fill_images_edges` / `fill_images_faces` 层。使用 OCCT 8.0 源码在 `C:\Users\lilu\works\OCCT\src\ModelingAlgorithms\TKBO\`。编译：`cargo check -p rcad-algorithms`。
```

## TKMath 状态（FoundationClasses — 全部对齐）

`rcad-kernel/src/math/` 所有 OCCT TKMath 包全覆盖：

| 模块 | OCCT 包 |
|------|---------|
| `bnd/` | Bnd (BndBox 数据) |
| `bspl` | BSplCLib + BSplSLib |
| `bvh/` | BVH |
| `cs_lib` | CSLib |
| `el/` | ElCLib + ElSLib |
| `math_poly` | MathPoly |
| `plib` | PLib |
| `poly/` | Poly |
| `top_loc` | TopLoc |
| `root/` | MathRoot |
| `opt/` | MathOpt |
| `lin/` | MathLin |
| `integ/` | MathInteg |
| `sys/` | MathSys |
| `curvature` | LProp 曲率 |
| `arc_length` | 弧长（委托 base::gcpnts） |

## TKGeomBase 状态（ModelingData — 100% 全部对齐）

`rcad-kernel/src/base/` 所有 OCCT TKGeomBase 包全覆盖：

| 模块 | OCCT 包 | 说明 |
|------|---------|------|
| `gc/` | GC + GCE2d + gce | MakeCircle/MakeLine/MakePlane 等构造算法 |
| `int_ana.rs` | IntAna | 3D 解析求交（线/平面/柱/球/锥/环） |
| `int_ana2d/` | IntAna2d | 2D 解析求交（线/圆/二次曲线） |
| `geom_lib/` | GeomLib | Tool, IsPlanarSurface, CheckCurveOnSurface |
| `geom_lprop/` | GeomLProp | CLProps(曲线局部), SLProps(曲面局部) |
| `proj_lib/` | ProjLib | 投影器（平面/柱/球/锥/环） |
| `geom_proj_lib/` | GeomProjLib | 高层投影 + 采样回退 |
| `geom_api/` | GeomAPI | project, interpolate, extrema, IntCS, IntSS |
| `convert/` | GeomConvert | 解析→BSpline + 反向(BSpline→解析/Bezier) |
| `geom2d_convert/` | Geom2dConvert | 2D BSpline↔Bezier + 组合/逼近 |
| `bnd_lib/` | BndLib | curve_bounding_box + surface_bounding_box |
| `geom_bnd_lib/` | GeomBndLib | AddCurve/AddSurface 到 BndBox |
| `gcpnts/` | GCPnts | 弧长取点 |
| `cpnts/` | CPnts | 均匀取点 |
| `gprop/` | GProp | 表面积/体积/惯性矩/板 |
| `extrema.rs` | Extrema | ExtPC/ExtPS/ExtCC + GenLocateExtPS |
| `lprop/` | LProp | 曲率极值/拐点分析 |
| `geom_tools/` | GeomTools | Dump 曲线/曲面/2D曲线 |
| `hermit/` | Hermit | Hermite 插值曲线 |
| `approx/` | Approx+AppCont+AppDef+AppParCurves | 曲线/曲面逼近 + 多线并行逼近 |

## 对齐原则

布尔管线全面对齐期间，所有底层模块**优先采用 OCCT 形式对齐**。布尔管线全部通过后，可逐步转换为 Rust 惯用 API。详见 `AGENTS.md`「后续重构指导」。

## 当前 git

在 `rcad` 子模块的 `main` 分支。最新提交包含 TKGeomBase 全部 23 个包的完对齐。
