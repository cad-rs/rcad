# RCAD 功能文档

> 版本：2026-04 · Phase P3 完成（含扫掠/圆角历史追踪）

---

## 目录

1. [概述](#1-概述)
2. [几何内核](#2-几何内核)
3. [拓扑结构](#3-拓扑结构)
4. [建模 API](#4-建模-api)
5. [算法层](#5-算法层)
6. [数据交换](#6-数据交换)
7. [渲染与可视化](#7-渲染与可视化)
8. [精度与容差](#8-精度与容差)
9. [与 OCCT 的功能对照](#9-与-occt-的功能对照)
10. [草图约束](#10-草图约束rcad-constraints)

---

## 1. 概述

RCAD 是一个用 Rust 实现的参数化实体建模内核，目标是在工程软件、CAE 前处理和电磁仿真前处理等场景中提供可嵌入、高安全性的几何计算能力。

**技术选型原则**

| 方面 | 选择 |
|------|------|
| 语言 | Rust（内存安全、无 GC、C ABI 友好）|
| 数学库 | `glam`（SIMD 友好的向量/矩阵，f64）|
| 序列化 | `serde` + JSON / bincode |
| 渲染 | `wgpu`（跨平台 GPU，WebGPU API）|
| 包管理 | Cargo workspace（独立编译、清晰边界）|

**Workspace 结构**

```
rcad-kernel        几何类型、拓扑、BRep、变换、投影、曲率、形状属性
rcad-modeling      建模 API：扫掠、圆角、偏移、历史树
rcad-algorithms    布尔运算、面印记、截面、HLR、IntSS、BRep 检查
rcad-step          STEP AP203/AP214 读写、装配体 IO
rcad-render        wgpu 渲染、拾取、HLR 描边
rcad-scene         创建命令状态机、装配体管理
rcad-constraints   2D/3D 草图约束求解、草图 → BRep
apps/creator-egui  egui 桌面应用
apps/creator-iced  iced 桌面应用
```

---

## 2. 几何内核

### 2.1 三维曲线（`Curve3`）

```
Line3        — 无限直线（origin + direction）
Circle3      — 圆（center + normal + radius）
Ellipse3     — 椭圆（center + normal + major_dir + semi_a + semi_b）
Hyperbola3   — 双曲线（center + normal + major_dir + semi_major + semi_minor）
Parabola3    — 抛物线（vertex + normal + axis_dir + focal_param）
BSplineCurve3— 非均匀有理 B-Spline（de Boor 算法，任意阶、任意节点向量）
BezierCurve3 — Bezier 曲线（de Casteljau）
OffsetCurve3 — 偏移曲线（基曲线 + 偏移距离 + 偏移平面法向）
```

所有曲线实现 `CurveEval` trait：

| 方法 | 说明 |
|------|------|
| `point_at(t) → DVec3` | 参数求值 |
| `tangent_at(t) → DVec3` | 切向量（有限差分）|
| `default_domain() → [f64; 2]` | 参数域 `[t0, t1]` |

曲线参数域通过 `GeomStore.edge_curve_range` 截断，等价于 OCCT `Geom_TrimmedCurve`。

---

### 2.2 三维曲面（`Surface3`）

```
Plane                  — 无限平面（origin + normal）
CylindricalSurface     — 圆柱面（origin + axis + radius；u=方位角，v=高度）
SphericalSurface       — 球面（center + axis + radius；u=经度，v=纬度余角）
ConicalSurface         — 锥面（apex + axis + half_angle）
ToroidalSurface        — 环面（center + axis + R + r）
BSplineSurface         — 张量积 B-Spline 曲面（双向 de Boor）
BezierSurface          — 张量积 Bezier 曲面（de Casteljau）
LinearExtrusionSurface — 线性拉伸面（轮廓曲线 + 方向向量）
RevolutionSurface      — 旋转面（轮廓曲线 + 轴）
OffsetSurface          — 偏移面（基面 + 法向偏移距离）
TrimmedSurface         — 矩形裁剪面（基面 + [u1,u2,v1,v2] 参数框）
```

所有曲面实现 `SurfaceEval` trait：

| 方法 | 说明 |
|------|------|
| `point_at(u, v) → DVec3` | 参数求值 |
| `normal_at(u, v) → DVec3` | 单位法向量 |
| `default_domain() → [f64; 4]` | 参数域 `[u0,u1,v0,v1]` |

---

### 2.3 二维曲线（PCurve，`Curve2d`）

用于将三维曲线在曲面参数域内的表示（STEP 必须）：

```
Line2d   — 二维直线
Circle2d — 二维圆
Ellipse2d— 二维椭圆
BSplineCurve2 — 二维 B-Spline
BezierCurve2  — 二维 Bezier
```

每条边可通过 `GeomStore.edge_pcurves` 存储若干 PCurve（每个相邻面一条）。

---

### 2.4 曲率与局部属性

| 功能 | API |
|------|-----|
| 主曲率 κ₁, κ₂ | `principal_curvatures(surface, u, v)` |
| 平均曲率 H | `mean_curvature(surface, u, v)` |
| 高斯曲率 K | `gaussian_curvature(surface, u, v)` |
| 曲线曲率向量 | `curvature_vector(curve, t)` |
| 密切圆半径 | `osculating_radius(curve, t)` |
| 弧长 | `arc_length(curve, t0, t1, n)` |

---

### 2.5 投影与极值

| 功能 | API |
|------|-----|
| 点投影到曲线 | `closest_point_on_curve(curve, query, n)` → `CurveProjection` |
| 点投影到曲面 | `closest_point_on_surface(surface, query, n)` → `SurfaceProjection` |
| 曲线-曲线极值 | `extrema_curve_curve(c1, c2, n)` → `CurveCurveExtrema { pairs }` |

`CurveCurveExtrema.pairs` 包含所有局部最小距离对 `(param1, param2, point1, point2, distance)`，按距离升序排列。

---

### 2.6 曲面-曲面交线（`intersect_surfaces`）

```rust
intersect_surfaces(s1: &Surface3, s2: &Surface3) -> SurfaceSurfaceIntersection
intersect_surfaces_with_density(s1, s2, density: usize) -> SurfaceSurfaceIntersection
```

`density` 控制数值 marching 网格分辨率（默认 48），对高曲率曲面可适当增大。

解析对：

| 组合 | 结果 |
|------|------|
| Plane × Plane | Line 或 None |
| Plane × Sphere | Circle、Point 或 None |
| Plane × Cylinder | Circle 或 None |
| Sphere × Sphere | Circle、Point 或 None |
| Cylinder × Cylinder（平行轴）| 1–2 条 Line |

其余组合退化为 48×48 网格 + marching 数值折线。

### 2.7 NURBS 互操作（`nurbs_convert`）

| 功能 | API | 说明 |
|------|-----|------|
| 任意曲线 → BSpline | `curve_to_bspline(curve, n_samples)` | 解析精确或采样插值 |
| 任意曲面 → BSpline | `surface_to_bspline(surface, n_u, n_v)` | 解析精确或双线性采样 |
| Line3 → BSpline | `line_to_bspline(line)` | 度-1，2 控制点 |
| Circle3 → BSpline | `circle_to_bspline(circle)` | 度-2，9 控制点，精确 |
| Ellipse3 → BSpline | `ellipse_to_bspline(ellipse)` | 同上，按半轴缩放 |
| Plane → BSpline | `plane_to_bspline(plane)` | 度-(1,1)，4 控制点 |
| Cylinder → BSpline | `cylinder_to_bspline(cyl)` | 度-(2,1)，精确 |
| Sphere → BSpline | `sphere_to_bspline(sphere)` | 度-(2,2)，5×9 控制网格，精确 |
| Bezier → BSpline | `bezier_curve_to_bspline` / `bezier_surface_to_bspline` | 插入端点节点 |

`BSplineSurface::point_at` 使用全有理张量积求值（`de_boor_homo` + 有理 de Boor），正确传递权重。

---

### 2.8 曲线/曲面裁剪与延伸（`extend`）

| 操作 | API | 类比 |
|------|-----|------|
| 曲线裁剪 | `trim_curve(curve, t0, t1)` | `Geom_TrimmedCurve` |
| 曲线延伸到点 | `extend_curve_to_point(curve, end, target)` | `GeomAPI_ExtendCurveToPoint` |
| 曲线延伸指定长度 | `extend_curve_by_length(curve, end, length)` | — |
| 曲面裁剪 | `trim_surface(surface, u0, u1, v0, v1)` | `Geom_RectangularTrimmedSurface` |
| BSpline 曲面延伸 | `extend_bspline_surface(surface, boundary, dist)` | `GeomAPI_ExtendSurfaceToShape` |

`trim_curve` 通过 Boehm 节点插入（重数 = degree+1）精确裁剪，结果参数域归一化到 `[0, 1]`。

---

### 2.9 形状属性（`properties`）

类比 OCCT `BRepGProp`。优先使用已有三角面片；无三角化时退化为 64×64 UV 网格采样（曲面体体积误差 < 0.1%），最终回退到外线框顶点扇形三角化。

| 功能 | API | 说明 |
|------|-----|------|
| 表面积 | `surface_area(brep: &BRep) → f64` | 所有面三角面积之和 |
| 体积 | `volume(brep: &BRep) → f64` | 有符号四面体法 |
| 质心 | `centroid(brep: &BRep) → DVec3` | 体积加权质心 |
| 惯性张量 | `inertia_tensor(brep: &BRep) → InertiaTensor` | `{ ixx, iyy, izz, ixy, ixz, iyz }` + `to_matrix()` |

---

## 3. 拓扑结构

### BRep 数据模型

```
BRep
 ├── vertices: Vec<Vertex>          3D 点
 ├── edges: Vec<Edge>               顶点对 (start, end)
 ├── solids: Vec<Solid>
 │     └── shells: Vec<Shell>
 │           └── faces: Vec<Face>
 │                 ├── outer_wire: Wire
 │                 │     └── edges: Vec<WireEdge { idx, forward }>
 │                 ├── inner_wires: Vec<Wire>  (孔洞)
 │                 ├── normal: DVec3
 │                 └── triangles: Vec<[usize;3]>  (渲染缓存)
 └── geom: GeomStore
       ├── curves: Vec<Curve3>
       ├── surfaces: Vec<Surface3>
       ├── curve2ds: Vec<Curve2d>
       ├── edge_curve: Vec<Option<usize>>
       ├── face_surface: Vec<Option<usize>>
       ├── edge_pcurves: Vec<Vec<PCurve>>
       ├── edge_curve_range: Vec<Option<[f64;2]>>
       ├── face_surface_range: Vec<Option<[f64;4]>>
       ├── edge_degenerated: Vec<bool>
       ├── vertex/edge/face_tolerance: Vec<f64>
       └── edge_same_parameter: Vec<bool>
```

`WireEdge.forward` 表示边的遍历方向（等价于 OCCT `FORWARD/REVERSED`），使同一条几何边可被相邻两个面以相反方向引用，无需复制。

---

### 变换

```rust
brep.apply_transform(mat: DAffine3)         // 原地变换
brep.transformed(mat: DAffine3) -> BRep     // 返回新实例，原 BRep 不变
```

变换覆盖所有顶点坐标及曲线/曲面的解析几何参数（origin、axis、控制点等）。`TrimmedSurface` 的参数域不受影响（参数空间坐标）。

---

### 拓扑查询（`topo_query`）

类比 OCCT `TopExp_Explorer` 和 `TopExp::MapShapesAndAncestors`：

| 功能 | API |
|------|-----|
| 边的邻接面 | `edge_adjacent_faces(brep, edge_idx) → Vec<usize>` |
| 面的所有边 | `face_edges(brep, face_idx) → Vec<usize>` |
| 顶点的邻接边 | `vertex_adjacent_edges(brep, vertex_idx) → Vec<usize>` |
| 面数 | `face_count(brep) → usize` |
| 边数 | `edge_count(brep) → usize` |
| 顶点数 | `vertex_count(brep) → usize` |
| 退化边判断 | `is_degenerate_edge(brep, edge_idx) → bool` |

所有函数对空 BRep 安全（返回 0 或空 Vec）。

---

### 外观属性（`appearance`）

类比 OCCT `XCAFDoc_ColorTool`，颜色与几何拓扑解耦存储：

```rust
pub struct Color { pub r: f64, pub g: f64, pub b: f64 }

pub struct StepColor {
    pub solid_color: Option<Color>,   // 实体级默认色
    pub face_colors: Vec<FaceColor>,  // 面级覆盖
}
```

内置预设色：`RED / GREEN / BLUE / YELLOW / CYAN / MAGENTA / WHITE / GRAY / SILVER / GOLD / ORANGE / BLACK`

`StepColor::color_for_face(face_index)` 优先返回面级覆盖，回退到实体级默认色。

---

## 4. 建模 API

### 4.1 基本体（`PrimitiveSolid`）

```rust
BRep::from_primitive(PrimitiveSolid::Box { width, height, depth })
BRep::from_primitive(PrimitiveSolid::Sphere { radius })
BRep::from_primitive(PrimitiveSolid::Cylinder { radius, height })
BRep::from_primitive(PrimitiveSolid::Cone { base_radius, height })
BRep::from_primitive(PrimitiveSolid::Torus { major_radius, minor_radius })
```

所有基本体的 `GeomStore` 在构造时完整填充，无需额外调用。

---

### 4.2 扫掠类

| 操作 | API | 说明 |
|------|-----|------|
| 拉伸 | `extrude(profile, dir, dist)` | 轮廓沿方向线性拉伸 |
| 拉伸（历史） | `extrude_with_history(...)` | 返回 `SweepHistory`（bottom_cap / top_cap / lateral_faces / profile_edge_to_lateral）|
| 旋转 | `revolve(profile, axis, angle)` | 轮廓绕轴旋转 |
| 旋转（历史） | `revolve_with_history(...)` | 返回 `SweepHistory`（全旋转时 cap 为空）|
| 管道扫掠 | `sweep_pipe(profile, spine)` | 轮廓沿脊线路径扫掠 |
| 放样 | `loft(profiles)` | 多截面放样 |
| 放样（历史） | `loft_with_history(...)` | 返回 `LoftHistory`（bottom_cap / top_cap / lateral_faces）|
| 变截面扫掠 | `sweep_pipe_variable(profiles, spine)` | 截面沿路径线性插值缩放 |

---

### 4.3 倒角与圆角

| 操作 | API |
|------|-----|
| 等距圆角 | `fillet(brep, edge_indices, radius)` |
| 等距圆角（历史） | `fillet_edge_with_history(brep, edge_idx, radius)` → `(BRep, FilletHistory)` |
| 等距倒角 | `chamfer(brep, edge_indices, dist)` |
| 等距倒角（历史） | `chamfer_edge_with_history(brep, edge_idx, dist)` → `(BRep, FilletHistory)` |
| 批量圆角 | `fillet_edges(brep, [(edge_idx, radius), ...])` |
| 批量圆角（历史） | `fillet_edges_with_history(...)` → `(BRep, MultiFilletHistory)` |
| 顶点混合圆角 | `corner_blend(brep, vertex_idx, radius)` |
| 顶点混合（历史） | `corner_blend_with_history(...)` → `(BRep, CornerBlendHistory)` |

`FilletHistory` 包含 `modified_faces`（被修剪的原始面）、`fillet_face`（新生成的圆角/倒角面）、`closing_faces`（端点处的三角封闭面）。

---

### 4.4 偏移与薄壁

| 操作 | API |
|------|-----|
| 面法向偏移 | `offset_surface(face, dist)` → `OffsetSurface` |
| 壳体偏移 | `offset_shell(brep, dist)` |

---

### 4.5 操作历史

每次布尔运算可附带历史追踪：

```rust
let (result, history) = boolean_op_with_history(op, &a, &b)?;
// history.faces[i] = FaceOrigin::FromA(src) | FromB(src) | Generated
```

`FaceOrigin` 可用于将结果面映射回原始输入面，支持参数化重建。

扫掠操作同样支持历史追踪：

```rust
// 拉伸历史
let (brep, hist) = extrude_with_history(&profile, face, dir, dist)?;
// hist.bottom_cap    → 底面索引
// hist.top_cap       → 顶面索引
// hist.lateral_faces → 侧面索引（按轮廓边顺序）
// hist.profile_edge_to_lateral → 轮廓边 → 侧面映射

// 放样历史
let (brep, hist) = loft_with_history(&profiles)?;
// hist.bottom_cap    → 底面
// hist.top_cap       → 顶面
// hist.lateral_faces → 侧面（按 section × edge 顺序）

// 圆角/倒角历史
let (brep, hist) = fillet_edge_with_history(&brep, edge_idx, radius)?;
// hist.modified_faces → 被修剪的原始面
// hist.fillet_face    → 新生成的圆角面
// hist.closing_faces  → 端点三角封闭面
```

---

## 5. 算法层

### 5.1 B-Rep 修复（`repair`）

```rust
repair(brep: &BRep, tolerance: f64) -> (BRep, RepairReport)
```

| 操作 | API | 说明 | OCCT 等价 |
|------|-----|------|-----------|
| 合并近邻顶点 | `merge_close_vertices(brep, tol)` | 并查集合并，重映射边/线框索引 | `BRepOffsetAPI_Sewing` / `ShapeFix_Wire` |
| 删除退化面 | `remove_degenerate_faces(brep)` | <3 边或 Newell 面积≈0 | `ShapeFix_Shape` |
| 重算面法向 | `recompute_face_normals(brep)` | Newell 方法，返回修改计数 | `BRepLib` 法向修复 |
| 修复线框方向 | `fix_wire_orientation(brep, tol)` | 翻转断链处的边方向 | `ShapeFix_Wire::FixClosed` |
| 全部修复 | `repair(brep, tol)` | 四步合一，返回 `RepairReport` | `ShapeFix_Shape::Perform` |

---

### 5.2 布尔运算

```rust
boolean_op(BooleanOpType::{Union|Intersection|Difference}, &a, &b) -> Result<BRep>
boolean_op_with_history(op, &a, &b) -> Result<(BRep, BooleanHistory)>
boolean_op_par(op, &a, &b) -> Result<(BRep, BooleanHistory)>   // Rayon 并行版本
```

**实现架构（OCCT BOPAlgo 风格）：**

```
DS::new(a, b)          构建顶点/边/面工作数据集
PaveFiller::perform()  六趟求交：VV → VE → EE → VF → EF → FF
  ├── 解析交：Plane×Plane/Sphere/Cylinder，Sphere×Sphere，Cylinder×Cylinder
  └── 数值交：marching（折线存储，面 AABB 边界约束）
BooleanBuilder::build() 光线法分类 → 子面分割 → 组装结果 BRep
```

**当前能力：**
- 纯平面体（长方体等）：完整 Union/Intersection/Difference，多场景测试通过
- 曲面体（球/柱/锥）：FF pass 解析 + marching，面分类支持曲面，子面分割仍在改善中

---

### 5.3 面印记（`imprint_brep`）

```rust
imprint_brep(target: &BRep, tool: &BRep) -> ImprintResult {
    brep: BRep,                         // target 各面被 tool 边界分割后的新 BRep
    seam_edges: Vec<(usize, usize)>,    // (target face idx, tool face idx) 共享边界对
}
```

将 `tool` 的几何边界印记到 `target` 的面上，不做布尔分类、保留所有 target 面。是 FEM/FDTD 共形网格生成的前置步骤。

---

### 5.4 间隙/重叠检测（`detect_gaps_overlaps`）

```rust
detect_gaps_overlaps(a: &BRep, b: &BRep, tolerance: f64) -> GapOverlapReport {
    gaps: Vec<Gap>,                     // 0 < d ≤ tolerance 的面对
    overlaps: Vec<Overlap>,             // d < 0（估计穿透深度）的面对
    shared_faces: Vec<(usize, usize)>,  // d ≈ 0 且法向反平行的共面对
}
```

流程：面 AABB 预筛选 → 采样 ≤5 点 → `closest_point_on_surface` 测距 → 分类。

---

### 5.5 截面

```rust
section_polylines(brep, plane) -> Vec<Vec<DVec3>>   // 折线截面（始终可用）
section_curves(brep, plane) -> Vec<SectionCurve>    // 解析截面
```

`SectionCurve::Analytic(Curve3)` 对 Plane/Sphere/Cylinder/Cone 面返回精确圆/椭圆/直线；其余退化为 `Polyline`。

---

### 5.6 形状距离

```rust
min_distance(a: &BRep, b: &BRep) -> ShapeDistance {
    distance: f64,
    point_on_a: DVec3,
    point_on_b: DVec3,
}
```

暴力面对循环，每个面 4×4 采样 + 线框顶点，双向 A→B、B→A 取最小值。

---

### 5.7 壳缝合（`sew_shells`）

```rust
sew_shells(breps: &[BRep], tolerance: f64) -> SewingResult {
    brep: BRep,
    stitched_pairs: Vec<(usize, usize)>,
    free_edges: Vec<usize>,
}
```

并查集顶点合并 + 边去重 + 单壳组装。

---

### 5.8 曲线拟合与插值

| 操作 | API |
|------|-----|
| 点列插值 | `interpolate_points(pts, degree)` → `BSplineCurve3` |
| 最小二乘近似 | `approximate_points(pts, degree, n_ctrl)` → `BSplineCurve3` |

---

### 5.9 消隐线渲染（HLR）

```rust
hlr(brep: &BRep, camera: &HlrCamera, samples: usize) -> HlrResult
hlr_to_svg(result: &HlrResult, scale: f64, margin: f64) -> String
```

投影所有边，分类为可见/遮挡，输出 SVG 线框图。

**相机预设：**

```rust
HlrCamera::isometric(distance)   // 等轴测（+X+Y+Z 方向）
HlrCamera::front(distance)       // 正视图（沿 +Y 看，up = +Z）
HlrCamera::top(distance)         // 俯视图（沿 -Z 看）
HlrCamera::right(distance)       // 右视图（沿 -X 看，up = +Z）
```

**解析轮廓线（silhouette）：**

| 曲面类型 | 轮廓线生成方式 |
|---------|--------------|
| `CylindricalSurface` | 两条平行于轴的直线（解析，视方向垂直于轴时有效）|
| `SphericalSurface` | 垂直于视方向的大圆（64 段密集折线）|
| `ConicalSurface` | 从顶点出发的两条母线（解析）|
| `ToroidalSurface` | 两条轮廓曲线（每条 64 段，按 `normal·view=0` 数值求解 v 角）|

**`HlrResult` 输出：**

```rust
result.visible()   // 可见线段迭代器
result.hidden()    // 遮挡线段迭代器
// 每个 HlrSegment 含 start/end (DVec2)、visible (bool)、curve_hint (Option<CurveHint>)
```

`CurveHint::Circle` 使 SVG 导出器输出 `<path>` 弧线而非折线。

---

### 5.10 曲面-曲面交线（IntSS）

见 [2.6 节](#26-曲面-曲面交线intersect_surfaces)。

---

### 5.11 BRep 检查（`brep_check`）

类比 OCCT `BRepCheck_Analyzer`，只读检查，不修改 BRep：

```rust
check(brep: &BRep) -> CheckResult { issues: Vec<CheckIssue> }
```

| 检查项 | `CheckIssue` 变体 | 说明 |
|--------|-----------------|------|
| 线框闭合性 | `OpenWire { solid, shell, face, wire_pos }` | 相邻边端点不连续 |
| 面法向有效性 | `ZeroNormal { solid, shell, face }` | 法向量为零向量 |
| 退化面 | `DegenerateFace { solid, shell, face }` | 外线框边数 < 3 |
| 边索引越界 | `InvalidEdgeIndex { solid, shell, face, edge_idx }` | WireEdge 引用越界 |
| 顶点索引越界 | `InvalidVertexIndex { edge, vertex_idx }` | 边引用顶点越界 |

`CheckResult::is_valid()` 当 `issues` 为空时返回 `true`。

---

### 5.12 曲线-曲面交点（`curve_surface`）

解析求交，类比 OCCT `IntCurvesFace_ShapeIntersector`：

```rust
pub struct CurveSurfaceHit {
    pub point: DVec3,
    pub curve_param: f64,   // 交点在曲线上的参数值
}
```

| 组合 | API |
|------|-----|
| Line × Cylinder | `intersect_line_cylinder(line, t_range, cyl) → Vec<CurveSurfaceHit>` |
| Line × Sphere | `intersect_line_sphere(line, t_range, sphere) → Vec<CurveSurfaceHit>` |
| Line × Cone | `intersect_line_cone(line, t_range, cone) → Vec<CurveSurfaceHit>` |

每个函数返回参数域 `t_range` 内所有交点（0–2 个），按 `curve_param` 升序排列。

---

## 6. 数据交换

### STEP 读取（AP203 / AP214）

支持实体类型：

| 几何 | 拓扑 |
|------|------|
| `LINE` → `Line3` | `ADVANCED_FACE` |
| `CIRCLE` → `Circle3` | `FACE_OUTER_BOUND` / `FACE_BOUND` |
| `ELLIPSE` → `Ellipse3` | `EDGE_CURVE` / `ORIENTED_EDGE` |
| `HYPERBOLA` → `Hyperbola3` | `VERTEX_POINT` |
| `PARABOLA` → `Parabola3` | `SHELL_BASED_SURFACE_MODEL` |
| `OFFSET_CURVE_3D` → `OffsetCurve3` | `CLOSED_SHELL` |
| `B_SPLINE_CURVE_WITH_KNOTS` → `BSplineCurve3` | `MANIFOLD_SOLID_BREP` |
| `PLANE` → `Plane` | |
| `CYLINDRICAL_SURFACE` → `CylindricalSurface` | |
| `SPHERICAL_SURFACE` → `SphericalSurface` | |
| `CONICAL_SURFACE` → `ConicalSurface` | |
| `TOROIDAL_SURFACE` → `ToroidalSurface` | |
| `B_SPLINE_SURFACE_WITH_KNOTS` → `BSplineSurface` | |
| `SURFACE_OF_LINEAR_EXTRUSION` → `LinearExtrusionSurface` | |
| `SURFACE_OF_REVOLUTION` → `RevolutionSurface` | |
| `RECTANGULAR_TRIMMED_SURFACE` → `TrimmedSurface` | |
| `SURFACE_CURVE` / PCurve 链 | `GEOMETRIC_CURVE_SET` |
| `UNCERTAINTY_MEASURE_WITH_UNIT` → 容差字段 | |

**颜色导入：**

```rust
StepReader::parse_string_with_color(s) -> Result<(BRep, Option<StepColor>)>
StepReader::read_file_with_color(path) -> Result<(BRep, Option<StepColor>)>
```

解析链：`STYLED_ITEM → PSA → SSU → SSS → SSFA → FAS → FASC → COLOUR_RGB`

---

### STEP 写出

```rust
StepWriter::write_string(&brep, ExportSelection) -> String
StepWriter::write_file(&brep, path, ExportSelection)
StepWriter::write_string_colored(&brep, &color, ExportSelection) -> String
```

支持导出全部已知几何类型，`TrimmedSurface` 写基面（裁剪由面线框拓扑隐含）。

---

### OBJ 网格读写

```rust
ObjReader::parse_string(obj) -> Result<BRep, ObjError>
ObjReader::read_file(path) -> Result<BRep, ObjError>

ObjWriter::write_string(&brep) -> String
ObjWriter::write_file(&brep, path) -> Result<usize>
write_obj(&brep, &mut writer) -> Result<usize>
```

范围：三角网格级别交换（`v` / `f`），`f` 多边形按扇形三角化，支持正/负索引。

---

### IGES 网格读写（Type 106）

```rust
IgesReader::parse_string(iges) -> Result<BRep, IgesError>
IgesReader::read_file(path) -> Result<BRep, IgesError>

IgesWriter::write_string(&brep) -> String
IgesWriter::write_file(&brep, path) -> Result<usize>
```

范围：通过 IGES Type 106（copious-data polyline）桥接三角面片，定位于网格互操作，不包含解析曲面/拓扑语义。

---

### 装配体 IO（`rcad-step assembly`）

类比 OCCT `XCAFDoc_ShapeTool`，支持平坦装配体和嵌套树形装配体的 STEP 往返。

**平坦装配体（`AssemblyComponent`）：**

```rust
// 写出
let comp_a = AssemblyComponent::new("part_a", brep_a)
    .with_translation(DVec3::new(10.0, 0.0, 0.0));
let step_str = write_assembly("my_asm", &[comp_a, comp_b]);

// 读回（每个组件几何独立隔离）
let components: Vec<AssemblyComponent> = read_assembly(&step_str)?;
// components[i].name, components[i].brep
```

**嵌套树形装配体（`AssemblyNode`）：**

```rust
// 构建树
let leaf_a = AssemblyNode::leaf("part_a", brep_a);
let leaf_b = AssemblyNode::leaf("part_b", brep_b);
let sub    = AssemblyNode::branch("sub_asm", vec![leaf_a, leaf_b]);
let root   = AssemblyNode::branch("root_asm", vec![sub, leaf_c]);

// 写出（递归生成 NAUO 层级）
let step_str = write_assembly_tree("root_asm", &root);

// 读回（保留树结构）
let tree: AssemblyNode = read_assembly_tree(&step_str)?;
// tree.name, tree.children, tree.brep (叶节点有 Some(BRep)，分支节点为 None)
```

`AssemblyNode` 支持 `.with_transform(DAffine3)` / `.with_translation(DVec3)` / `.with_color(Color)`。

STEP 结构：每个组件对应一组 `PRODUCT_DEFINITION` + `SHAPE_DEFINITION_REPRESENTATION`，层级关系通过 `NEXT_ASSEMBLY_USAGE_OCCURRENCE` 表达。读取时通过 BFS 从各组件的 `SHAPE_REPRESENTATION` 出发收集可达实体，实现几何隔离。

---

## 7. 渲染与可视化

### 渲染管线（wgpu）

| 功能 | 说明 |
|------|-----|
| 实体着色 | 三角面片渲染，Blinn-Phong 光照（可配置光照方向） |
| 线框 | 边的可见性着色（选中/悬停/普通）|
| 显示模式 | SolidWithEdges / Solid / Wireframe / Transparent 四种模式 |
| 坐标系可视化 | XYZ 轴箭头（红/绿/蓝），带锥形箭头和线段轴杆 |
| 背景网格 | XZ 平面网格，主线/次线区分，可开关 |
| 拾取 | 鼠标点击 → face/edge 索引（光线投射拾取）|
| 多选 | 累积选择模式（additive_select）|
| 相机 | 透视投影，左键旋转，中键平移，滚轮缩放 |
| 可配置光照 | 光照方向通过 uniform 传入，支持 headlight 模式（光跟随相机）|
| Per-object 颜色 | 动态设置模型颜色（set_model_color）|
| 截图导出 | 离屏渲染到 RGBA 纹理 → PNG 文件（screenshot_to_file）|
| HLR 描边 | `hlr_to_svg` 生成工程图风格 SVG |

---

## 8. 精度与容差

| 参数 | 值 | 说明 |
|------|----|------|
| `CONFUSION` (`rcad-kernel`) | `1e-7` | 点重合判断（OCCT `Precision::Confusion()`）|
| `ANGULAR` (`rcad-kernel`) | `1e-12` | 角度精度（OCCT `Precision::Angular()`）|
| `APPROXIMATION` (`rcad-kernel`) | `1e-4` | 曲面近似容差（OCCT `Precision::Approximation()`）|
| `TOLERANCE_ABS` (`rcad-algorithms`) | `1e-7` | 算法层点重合判断，与 `CONFUSION` 对齐 |
| `TOLERANCE_ANG` (`rcad-algorithms`) | `1e-9` | 算法层平行向量判断；比 `ANGULAR` 宽松，容许交叉计算中的浮点积累误差 |
| `vertex_tolerance` | per-vertex | 从 STEP `UNCERTAINTY_MEASURE_WITH_UNIT` 读取 |
| `edge_tolerance` | per-edge | 同上 |
| `face_tolerance` | per-face | 同上 |
| `edge_same_parameter` | bool | PCurve 与 3D 曲线的参数一致性标记 |
| `edge_same_range` | bool | 参数域一致性标记 |

---

## 9. 与 OCCT 的功能对照

### 几何类型

| 类别 | OCCT | RCAD | 状态 |
|------|------|------|------|
| 直线 | `Geom_Line` | `Line3` | ✅ |
| 圆 | `Geom_Circle` | `Circle3` | ✅ |
| 椭圆 | `Geom_Ellipse` | `Ellipse3` | ✅ |
| B-Spline 曲线 | `Geom_BSplineCurve` | `BSplineCurve3` | ✅ |
| Bezier 曲线 | `Geom_BezierCurve` | `BezierCurve3` | ✅ |
| 裁剪曲线 | `Geom_TrimmedCurve` | `edge_curve_range` | ✅ 参数域截断 |
| 偏移曲线 | `Geom_OffsetCurve` | `OffsetCurve3` | ✅ |
| 双曲线 | `Geom_Hyperbola` | `Hyperbola3` | ✅ |
| 抛物线 | `Geom_Parabola` | `Parabola3` | ✅ |
| 平面 | `Geom_Plane` | `Plane` | ✅ |
| 圆柱面 | `Geom_CylindricalSurface` | `CylindricalSurface` | ✅ |
| 球面 | `Geom_SphericalSurface` | `SphericalSurface` | ✅ |
| 锥面 | `Geom_ConicalSurface` | `ConicalSurface` | ✅ |
| 环面 | `Geom_ToroidalSurface` | `ToroidalSurface` | ✅ |
| B-Spline 曲面 | `Geom_BSplineSurface` | `BSplineSurface` | ✅ |
| Bezier 曲面 | `Geom_BezierSurface` | `BezierSurface` | ✅ |
| 裁剪曲面 | `Geom_RectangularTrimmedSurface` | `TrimmedSurface` | ✅ |
| 拉伸面 | `Geom_SurfaceOfLinearExtrusion` | `LinearExtrusionSurface` | ✅ |
| 旋转面 | `Geom_SurfaceOfRevolution` | `RevolutionSurface` | ✅ |
| 偏移面 | `Geom_OffsetSurface` | `OffsetSurface` | ✅ |

### 建模算法

| 功能 | OCCT | RCAD | 状态 |
|------|------|------|------|
| 基本体 | `BRepPrimAPI_Make*` | `BRep::from_primitive` | ✅ |
| 拉伸 | `BRepPrimAPI_MakePrism` | `extrude` | ✅ |
| 拉伸历史 | `BRepPrimAPI_MakePrism::Generated()` | `extrude_with_history` | ✅ |
| 旋转 | `BRepPrimAPI_MakeRevol` | `revolve` | ✅ |
| 旋转历史 | `BRepPrimAPI_MakeRevol::Generated()` | `revolve_with_history` | ✅ |
| 管道扫掠 | `BRepOffsetAPI_MakePipe` | `sweep_pipe` | ✅ |
| 放样 | `BRepOffsetAPI_ThruSections` | `loft` | ✅ |
| 放样历史 | — | `loft_with_history` | ✅ |
| 变截面扫掠 | `BRepOffsetAPI_MakePipeShell` | `sweep_pipe_variable` | ✅ |
| 等距圆角 | `BRepFilletAPI_MakeFillet` | `fillet` | ✅ |
| 圆角历史 | `BRepFilletAPI_MakeFillet::Modified()` | `fillet_edge_with_history` | ✅ |
| 等距倒角 | `BRepFilletAPI_MakeChamfer` | `chamfer` | ✅ |
| 倒角历史 | `BRepFilletAPI_MakeChamfer::Modified()` | `chamfer_edge_with_history` | ✅ |
| 顶点混合圆角 | `ChFi3d_Builder` | `corner_blend` | ✅ |
| 布尔 Union | `BRepAlgoAPI_Fuse` | `boolean_op(Union, ...)` | ✅ 平面体；曲面体部分 |
| 布尔 Intersection | `BRepAlgoAPI_Common` | `boolean_op(Intersection, ...)` | ✅ 平面体；曲面体部分 |
| 布尔 Difference | `BRepAlgoAPI_Cut` | `boolean_op(Difference, ...)` | ✅ 平面体；曲面体部分 |
| 面印记 | `BRepAlgoAPI_Splitter` | `imprint_brep` | ✅ 平面体完整 |
| 壳缝合 | `BRepOffsetAPI_Sewing` | `sew_shells` | ✅ |
| 变换 | `BRepBuilderAPI_Transform` | `apply_transform / transformed` | ✅ |

### 分析与查询

| 功能 | OCCT | RCAD | 状态 |
|------|------|------|------|
| 点投影到曲线 | `GeomAPI_ProjectPointOnCurve` | `closest_point_on_curve` | ✅ |
| 点投影到曲面 | `GeomAPI_ProjectPointOnSurf` | `closest_point_on_surface` | ✅ |
| 曲线-曲线极值 | `GeomAPI_ExtremaCurveCurve` | `extrema_curve_curve` | ✅ |
| 曲面-曲面交线 | `GeomAPI_IntSS` | `intersect_surfaces` | ✅ 解析+数值 |
| 形状最小距离 | `BRepExtrema_DistShapeShape` | `min_distance` | ✅ |
| 曲率 | `GeomLProp_SLProps` | `principal/mean/gaussian_curvature` | ✅ |
| 弧长 | `GCPnts_AbscissaPoint` | `arc_length` | ✅ |
| 截面曲线 | `BRepAlgoAPI_Section` | `section_curves` | ✅ 解析 |
| 间隙/重叠检测 | — | `detect_gaps_overlaps` | ✅ |
| 颜色属性（读） | `XCAFDoc_ColorTool` | `parse_string_with_color` | ✅ |
| 颜色属性（写） | `XCAFDoc_ColorTool` | `write_string_colored` | ✅ |
| 体积/表面积/质心/惯性 | `BRepGProp` | `volume / surface_area / centroid / inertia_tensor` | ✅ |
| 拓扑查询 | `TopExp_Explorer` | `edge_adjacent_faces / face_edges / vertex_adjacent_edges` | ✅ |
| BRep 有效性检查 | `BRepCheck_Analyzer` | `check` | ✅ |
| 曲线-曲面交点 | `IntCurvesFace_ShapeIntersector` | `intersect_line_cylinder/sphere/cone` | ✅ 解析 |

### 数据交换

| 格式 | OCCT | RCAD | 状态 |
|------|------|------|------|
| STEP AP203/AP214 读 | `STEPControl_Reader` | `StepReader` | ✅ |
| STEP AP203/AP214 写 | `STEPControl_Writer` | `StepWriter` | ✅ |
| IGES（网格桥接） | `IGESControl_Reader/Writer` | `IgesReader` / `IgesWriter` | ✅ Type 106 |
| OBJ（网格） | `RWObj` | `ObjReader` / `ObjWriter` | ✅ |
| GDS-II | — | — | ❌ |

### 目前不支持的 OCCT 功能

| 功能 | OCCT 类 | 备注 | 状态 |
|------|---------|------|------|
| 体网格 (TetGen) | — | 独立集成 | ❌ |

---

## 10. 草图约束（`rcad-constraints`）

2D 参数化草图约束求解器，类比 OCCT `GCS` / FreeCAD Sketcher。

### 实体与草图 API

```rust
Sketch::new() -> Self
Sketch::add_point(x, y) -> EntityId
Sketch::add_line(x1, y1, x2, y2) -> EntityId
Sketch::add_circle(cx, cy, r) -> EntityId
Sketch::add_arc(cx, cy, r, start_angle, end_angle) -> EntityId
Sketch::fix_param(param_idx)       // 锁定单个参数
Sketch::fix_entity(id)             // 锁定实体全部参数
Sketch::add_constraint(c: Constraint)
Sketch::dof() -> i64               // 当前自由度数（0 = 完全约束）
Sketch::solve() -> SolveResult     // Newton-Raphson 非线性迭代求解
```

访问求解后的几何：

```rust
Sketch::point(id) -> PointRef
Sketch::line_start(id) / line_end(id) -> PointRef
Sketch::center(id) -> PointRef
Sketch::point_coords(p: PointRef) -> DVec2
Sketch::entity_params(id) -> &[f64]
```

### 约束类型

**点约束（2 方程）：**

| 约束 | 构造 | 说明 |
|------|------|------|
| `Fixed` | `Constraint::fix_point(point, x, y)` | 固定点坐标 |
| `Coincident` | `Constraint::coincident(p1, p2)` | 两点重合 |

**距离/长度（1 方程）：**

| 约束 | 构造 | 说明 |
|------|------|------|
| `PointDistance` | `Constraint::point_distance(p1, p2, d)` | 两点欧氏距离 |
| `LineLength` | `LineLength { line, length }` | 线段长度 |

**线方向（1 方程）：**

| 约束 | 说明 |
|------|------|
| `Horizontal(line)` | y1 == y2 |
| `Vertical(line)` | x1 == x2 |

**线-线关系（1 方程）：**

| 约束 | 说明 |
|------|------|
| `Parallel(l1, l2)` | 方向向量叉积 = 0 |
| `Perpendicular(l1, l2)` | 方向向量点积 = 0 |
| `EqualLength(l1, l2)` | 长度平方差 = 0 |
| `Angle { l1, l2, angle_rad }` | 两线夹角 |

**圆/弧约束（1 方程）：**

| 约束 | 说明 |
|------|------|
| `Radius { circle, radius }` | 固定半径 |
| `EqualRadius(c1, c2)` | 两圆/弧半径相等 |
| `PointOnCircle { point, circle }` | 点在圆上 |
| `Tangent { circle, line }` | 圆与直线相切 |
| `CircleCircleTangent { c1, c2, external }` | 两圆外切（`external=true`）或内切 |
| `ArcArcTangent { a1, a2, external }` | 两弧外切或内切（与圆-圆相切同残差）|

**点-线关系（1 方程）：**

| 约束 | 说明 |
|------|------|
| `PointOnLine { point, line }` | 点在直线（无限延伸）上 |

**对称（2 方程）：**

| 约束 | 说明 |
|------|------|
| `Symmetric { p1, p2, line }` | p1、p2 关于 line 对称（中点在线上 + 连线垂直于线）|

### 求解器

Newton-Raphson 迭代，数值 Jacobian（有限差分），固定参数通过掩码跳过。

```rust
pub struct SolveResult {
    pub converged: bool,
    pub residual: f64,   // 最终残差 L2 范数
    pub iterations: u32,
}
```

### 草图 → BRep

```rust
// 线框 BRep（所有实体 → 边，无面）
Sketch::to_wire_brep() -> BRep

// 实体 BRep（闭合线段多边形 → 拉伸实体）
Sketch::to_solid_brep(height: f64) -> Option<BRep>
```

`to_solid_brep` 自动将草图中的 `Line` 实体链接成闭合多边形，构建平面轮廓面，沿 `+Z` 方向拉伸 `height`，返回含 6 个面（底/顶/4 侧）的实体 BRep（矩形轮廓）。多边形顶点按端点连接自动排序，不要求输入顺序。

---

## 11. 3D 草图约束（`rcad-constraints::space3d`）

3D 参数化空间约束求解器，类比 2D 草图约束但操作于三维空间。

### 实体类型

| 实体 | 参数 | 说明 |
|------|------|------|
| `SpacePoint` | 3 (x, y, z) | 空间点 |
| `SpaceLine` | 6 (x0, y0, z0, x1, y1, z1) | 空间线段 |
| `Plane` | 4 (nx, ny, nz, d) | 平面（法向 + 距离）|
| `Sphere` | 4 (cx, cy, cz, r) | 球面 |

### 约束类型

| 约束 | 方程数 | 说明 |
|------|--------|------|
| `Fixed { point, x, y, z }` | 3 | 固定点坐标 |
| `Coincident(p1, p2)` | 3 | 两点重合 |
| `PointDistance { p1, p2, distance }` | 1 | 两点距离 |
| `PointOnLine { point, line }` | 2 | 点在直线上 |
| `PointOnPlane { point, plane }` | 1 | 点在平面上 |
| `LineParallelLine(l1, l2)` | 2 | 两线平行 |
| `LinePerpendicularLine(l1, l2)` | 1 | 两线垂直 |
| `LineLength { line, length }` | 1 | 线段长度 |
| `PlaneNormal { plane }` | 2 | 法向单位化 |
| `PlaneParallel(p1, p2)` | 2 | 两平面平行 |
| `PlaneAngle { p1, p2, angle }` | 1 | 两平面夹角 |
| `SphereRadius { sphere, radius }` | 1 | 球面半径 |
| `SphereTangent { s1, s2, external }` | 1 | 两球相切 |
| `SphereSphereAngle { s1, s2, angle }` | 1 | 两球心连线方向约束 |

### 求解器 API

```rust
SpaceSketch::new() -> Self
SpaceSketch::add_point(x, y, z) -> SpaceEntityId
SpaceSketch::add_line(x0, y0, z0, x1, y1, z1) -> SpaceEntityId
SpaceSketch::add_plane(nx, ny, nz, d) -> SpaceEntityId
SpaceSketch::add_sphere(cx, cy, cz, r) -> SpaceEntityId
SpaceSketch::add_constraint(c: SpaceConstraint)
SpaceSketch::fix_param(param_idx)
SpaceSketch::fix_entity(id)
SpaceSketch::solve() -> SpaceSolveResult
```

求解器使用 Newton-Raphson 迭代 + Tikhonov 正则化（处理欠约束系统），数值 Jacobian（中心有限差分）。

---

*文档更新于 2026-04-11（P3 更新：§5.9 HLR 解析轮廓线新增锥面/环面；§6 装配体 IO 重写，新增嵌套树形装配体 `AssemblyNode` / `write_assembly_tree` / `read_assembly_tree`；§10 草图约束全面扩充，新增 `ArcArcTangent`、`Symmetric`、`to_solid_brep` 拉伸，约束类型从 3 种扩展到 16 种；§11 新增 3D 空间约束求解器；§4.2-4.3 扫掠/圆角操作新增 `*_with_history` 变体，支持面来源追踪；P1: `mesh_brep` FEM 质量网格生成等价于 `BRepMesh_IncrementalMesh`；P1: 消除 `intss.rs` 中 NaN panic；已移除"目前不支持的 OCCT 功能"表中已完成项）*
