# RCAD 功能文档

> 版本：2026-04 · Phase T 完成

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
rcad-kernel        几何类型、拓扑、BRep、变换、投影、曲率
rcad-modeling      建模 API：扫掠、圆角、偏移、历史树
rcad-algorithms    布尔运算、面印记、截面、HLR、IntSS
rcad-step          STEP AP203/AP214 读写
rcad-render        wgpu 渲染、拾取、HLR 描边
rcad-scene         创建命令状态机（Box/Sphere 流程）
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
```

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
| 旋转 | `revolve(profile, axis, angle)` | 轮廓绕轴旋转 |
| 管道扫掠 | `sweep_pipe(profile, spine)` | 轮廓沿脊线路径扫掠 |
| 放样 | `loft(profiles)` | 多截面放样 |
| 变截面扫掠 | `sweep_variable(profiles, spine)` | 截面沿路径线性插值缩放 |

---

### 4.3 倒角与圆角

| 操作 | API |
|------|-----|
| 等距圆角 | `fillet(brep, edge_indices, radius)` |
| 等距倒角 | `chamfer(brep, edge_indices, dist)` |
| 批量圆角 | `fillet_edges(brep, [(edge_idx, radius), ...])` |

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
hlr(brep: &BRep, camera: &HlrCamera) -> HlrResult
hlr_to_svg(result: &HlrResult, width, height) -> String
```

投影所有边，分类为可见/遮挡，输出 SVG 线框图。

---

### 5.10 曲面-曲面交线（IntSS）

见 [2.6 节](#26-曲面-曲面交线intersect_surfaces)。

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
| 旋转 | `BRepPrimAPI_MakeRevol` | `revolve` | ✅ |
| 管道扫掠 | `BRepOffsetAPI_MakePipe` | `sweep_pipe` | ✅ |
| 放样 | `BRepOffsetAPI_ThruSections` | `loft` | ✅ |
| 变截面扫掠 | `BRepOffsetAPI_MakePipeShell` | `sweep_variable` | ✅ |
| 等距圆角 | `BRepFilletAPI_MakeFillet` | `fillet` | ✅ |
| 等距倒角 | `BRepFilletAPI_MakeChamfer` | `chamfer` | ✅ |
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

### 数据交换

| 格式 | OCCT | RCAD | 状态 |
|------|------|------|------|
| STEP AP203/AP214 读 | `STEPControl_Reader` | `StepReader` | ✅ |
| STEP AP203/AP214 写 | `STEPControl_Writer` | `StepWriter` | ✅ |
| IGES（网格桥接） | `IGESControl_Reader/Writer` | `IgesReader` / `IgesWriter` | ✅ Type 106 |
| OBJ（网格） | `RWObj` | `ObjReader` / `ObjWriter` | ✅ |
| GDS-II | — | — | ❌ |

### 目前不支持的 OCCT 功能（低优先或尚未规划）

| 功能 | OCCT 类 | 备注 | 状态 |
|------|---------|------|------|
| B-Rep 修复/清理 | `ShapeFix_*` | `merge_close_vertices` / `repair` 已实现 | ✅ |
| NURBS 互操作 | `GeomConvert` | `curve_to_bspline` / `surface_to_bspline` 已实现 | ✅ |
| 曲线裁剪/延伸 | `GeomAPI_ExtendCurveToPoint` | `trim_curve` / `extend_curve_*` 已实现 | ✅ |
| 曲面裁剪/延伸 | `Geom_RectangularTrimmedSurface` | `trim_surface` / `extend_bspline_surface` 已实现 | ✅ |
| 装配体/实例化 | `XCAFDoc_ShapeTool` | 无场景图 | ❌ |
| 参数化约束求解 | `GCS`, Sketcher | 独立模块 | ❌ |
| FEM 网格生成 | `BRepMesh_IncrementalMesh` | 仅渲染用三角化 | ❌ |
| 体网格 (TetGen) | — | 独立集成 | ❌ |

---

*文档更新于 2026-04-09（UV 缝合修复、rcad-scene/rcad-modeling/rcad-kernel 集成测试新增、Clippy 清理完成）*
