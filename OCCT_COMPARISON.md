# RCAD vs OCCT 对比文档

> 本文档系统比较 RCAD 与 Open CASCADE Technology (OCCT) 的数据结构、功能实现与接口设计，
> 帮助团队明确差距、规划后续开发方向。
>
> **文档状态：** 基于 RCAD 当前代码（2026-04）生成，随代码演进应同步更新。
> **Phase A–I 全部完成（2026-04-03）。**

---

## 目录

1. [整体架构对比](#1-整体架构对比)
2. [几何层对比](#2-几何层对比)
   - 2.1 三维曲线
   - 2.2 三维曲面
   - 2.3 二维曲线（PCurve）
   - 2.4 几何评估与局部属性
3. [拓扑层对比](#3-拓扑层对比)
   - 3.1 拓扑实体
   - 3.2 拓扑属性与元数据
   - 3.3 拓扑遍历与查询
4. [建模 API 对比](#4-建模-api-对比)
   - 4.1 基本体（Primitives）
   - 4.2 曲线 / 曲面构造
   - 4.3 布尔运算
   - 4.4 扫掠 / 拉伸 / 旋转
   - 4.5 倒角 / 圆角
   - 4.6 薄壁 / 偏移
5. [数据交换对比](#5-数据交换对比)
6. [全局属性与分析](#6-全局属性与分析)
7. [精度与容差体系](#7-精度与容差体系)
8. [渲染与可视化](#8-渲染与可视化)
9. [综合差距汇总表](#9-综合差距汇总表)
10. [开发路线建议](#10-开发路线建议)

---

## 1. 整体架构对比

### OCCT 分层架构

```
FoundationClasses (数学、容差、内存)
    └── ModelingData (几何与拓扑数据结构)
            └── ModelingAlgorithms (布尔、圆角、扫掠…)
                    └── DataExchange (STEP、IGES、OBJ…)
                            └── Visualization (AIS、V3d、OpenGL)
```

| 层次 | OCCT 包 | RCAD 对应 |
|------|---------|-----------|
| 数学基础 | gp, TColgp | `glam::DVec3/DVec2/DMat3` (外部依赖) |
| 几何数据 | Geom, Geom2d, GeomAdaptor | `rcad-kernel/src/geom.rs` |
| 拓扑数据 | TopoDS, TopAbs, BRep | `rcad-kernel/src/lib.rs` + `topology.rs` |
| 建模入口 | BRepBuilderAPI, BRepPrimAPI | `rcad-modeling/src/builder/` |
| 算法 | BRepAlgoAPI, BRepFilletAPI, BRepOffsetAPI | `rcad-algorithms/` (仅布尔) |
| 数据交换 | STEPControl, IGESControl | `rcad-step/` (仅 STEP) |
| 可视化 | AIS, V3d, Graphic3d | `rcad-render/` (wgpu) |

**关键架构差异：**
- OCCT 为 C++ 单体库，内部紧耦合；RCAD 为 Cargo workspace，模块边界清晰，编译隔离。
- OCCT 使用 `Handle<T>`（引用计数智能指针）管理所有几何/拓扑对象；RCAD 目前以 `Vec` + index 代替，避免了指针图，但也缺少 OCCT 中任意子形状共享的能力。
- OCCT 拓扑对象携带 `TopLoc_Location`（变换矩阵）实现实例化共享；RCAD 直接在顶点坐标中存储变换后的位置，没有实例化层。

---

## 2. 几何层对比

### 2.1 三维曲线

#### OCCT `Geom_Curve` 继承树

```
Geom_Curve (abstract, parametric: t → Point3)
├── Geom_BoundedCurve
│   ├── Geom_TrimmedCurve      ← 裁剪曲线：任意 Curve + [t1, t2]
│   ├── Geom_BSplineCurve      ← B-Spline
│   └── Geom_BezierCurve       ← Bezier
├── Geom_Conic (abstract)
│   ├── Geom_Circle
│   ├── Geom_Ellipse
│   ├── Geom_Hyperbola
│   └── Geom_Parabola
├── Geom_Line
└── Geom_OffsetCurve           ← 偏移曲线（距离 d + 基准曲线）
```

每条曲线都有：
- `Value(t)` → Point3
- `D1(t)` → (Point3, Vec3) — 一阶导数
- `D2(t)` → (Point3, Vec3, Vec3) — 二阶导数
- `FirstParameter()` / `LastParameter()` — 参数域
- `IsClosed()` / `IsPeriodic()`

#### RCAD `Curve3` 现状

```rust
pub enum Curve3 {
    Line(Line3),
    Circle(Circle3),
    Ellipse(Ellipse3),
    BSpline(BSplineCurve3),  // Phase B 新增
}
```

| 曲线类型 | OCCT | RCAD | 差距说明 |
|---------|------|------|---------|
| 直线（无限） | `Geom_Line` | `Line3` | ✅ 对等 |
| 圆 | `Geom_Circle` | `Circle3` | ✅ 对等 |
| 椭圆 | `Geom_Ellipse` | `Ellipse3` | ✅ 对等 |
| **裁剪曲线** | `Geom_TrimmedCurve` | 以 `edge_curve_range [t1,t2]` 实现 | ✅ 参数范围已支持 |
| B-Spline | `Geom_BSplineCurve` | `BSplineCurve3`（de Boor）| ✅ Phase B |
| Bezier | `Geom_BezierCurve` | ❌ 缺失 | 中优先 |
| 双曲线 | `Geom_Hyperbola` | ❌ 缺失 | 低优先（工程场景罕见）|
| 抛物线 | `Geom_Parabola` | ❌ 缺失 | 低优先 |
| 偏移曲线 | `Geom_OffsetCurve` | ❌ 缺失 | 中优先（倒角需要）|
| 参数求值 | `Value(t)` → Point3 | `CurveEval::point_at(t)` | ✅ Phase A |
| 参数域 | `FirstParameter/LastParameter` | `CurveEval::default_domain()` + `edge_curve_range` | ✅ Phase A |

---

### 2.2 三维曲面

#### OCCT `Geom_Surface` 继承树

```
Geom_Surface (abstract, parametric: (u,v) → Point3)
├── Geom_BoundedSurface
│   ├── Geom_BSplineSurface
│   ├── Geom_BezierSurface
│   └── Geom_RectangularTrimmedSurface   ← 裁剪面
├── Geom_ElementarySurface
│   ├── Geom_Plane
│   ├── Geom_CylindricalSurface
│   ├── Geom_SphericalSurface
│   ├── Geom_ConicalSurface
│   └── Geom_ToroidalSurface
├── Geom_SweptSurface
│   ├── Geom_SurfaceOfLinearExtrusion    ← 线性扫掠面
│   └── Geom_SurfaceOfRevolution        ← 旋转面
└── Geom_OffsetSurface                  ← 偏移面
```

每个曲面都有：
- `Value(u, v)` → Point3
- `D1(u, v)` → (Point3, dU, dV)
- `UPeriod()` / `VPeriod()`
- `UIsoCurve(u)` / `VIsoCurve(v)` → `Geom_Curve`
- `Bounds(u1, u2, v1, v2)` — 参数域

#### RCAD `Surface3` 现状

```rust
pub enum Surface3 {
    Plane(Plane),
    Cylinder(CylindricalSurface),
    Sphere(SphericalSurface),
    Cone(ConicalSurface),
    Torus(ToroidalSurface),
    BSpline(BSplineSurface),  // Phase B 新增
}
```

| 曲面类型 | OCCT | RCAD | 差距说明 |
|---------|------|------|---------|
| 平面 | `Geom_Plane` | `Plane` | ✅ 对等 |
| 柱面 | `Geom_CylindricalSurface` | `CylindricalSurface` | ✅ 对等 |
| 球面 | `Geom_SphericalSurface` | `SphericalSurface` | ✅ 对等 |
| 锥面 | `Geom_ConicalSurface` | `ConicalSurface` | ✅ 对等 |
| 环面 | `Geom_ToroidalSurface` | `ToroidalSurface` | ✅ 对等 |
| **裁剪面** | `Geom_RectangularTrimmedSurface` | ❌ 缺失 | 中优先 |
| B-Spline 面 | `Geom_BSplineSurface` | `BSplineSurface`（张量积 de Boor）| ✅ Phase B |
| Bezier 面 | `Geom_BezierSurface` | ❌ 缺失 | 中优先 |
| 线性扫掠面 | `Geom_SurfaceOfLinearExtrusion` | `LinearExtrusionSurface`（Phase K）| ✅ Phase K |
| 旋转面 | `Geom_SurfaceOfRevolution` | `RevolutionSurface`（Phase K）| ✅ Phase K |
| 偏移面 | `Geom_OffsetSurface` | ❌ 缺失 | 低优先 |
| (u,v) 求值 | `Value(u,v)` → Point3 | `SurfaceEval::point_at(u,v)` | ✅ Phase A |
| 法向量 | `Normal(u,v)` | `SurfaceEval::normal_at(u,v)` | ✅ Phase A |
| 参数域查询 | `Bounds(u1,u2,v1,v2)` | `SurfaceEval::default_domain()` / `face_domain()` | ✅ Phase A + K |
| 面域覆盖 | `BRep_Face::UVBounds()` | `GeomStore.face_surface_range` + `face_domain()`（Phase K）| ✅ Phase K |

---

### 2.3 二维曲线（PCurve / Geom2d）

#### OCCT `Geom2d_Curve` 继承树

```
Geom2d_Curve
├── Geom2d_BoundedCurve
│   ├── Geom2d_TrimmedCurve
│   ├── Geom2d_BSplineCurve
│   └── Geom2d_BezierCurve
├── Geom2d_Conic
│   ├── Geom2d_Circle
│   ├── Geom2d_Ellipse
│   ├── Geom2d_Hyperbola
│   └── Geom2d_Parabola
├── Geom2d_Line
└── Geom2d_OffsetCurve
```

#### RCAD `Curve2d` 现状

```rust
pub enum Curve2d {
    Line(Line2d),
    Circle(Circle2d),
    Ellipse(Ellipse2d), // Phase J 新增
    BSpline(BSplineCurve2),  // Phase I 新增
}
```

| 类型 | OCCT | RCAD | 差距 |
|------|------|------|------|
| 2D 直线 | `Geom2d_Line` | `Line2d` | ✅ |
| 2D 圆 | `Geom2d_Circle` | `Circle2d` | ✅ |
| **2D B-Spline** | `Geom2d_BSplineCurve` | `BSplineCurve2`（de Boor 2D）| ✅ Phase I |
| **2D 椭圆** | `Geom2d_Ellipse` | `Ellipse2d`（Phase J）| ✅ Phase J |
| **2D 裁剪曲线** | `Geom2d_TrimmedCurve` | `curve2d_range: Vec<Option<[f64;2]>>`（Phase J）| ✅ Phase J（参数范围记录）|

---

### 2.4 几何评估与局部属性

OCCT 提供丰富的几何计算 API，RCAD Phase A 补全了核心求值接口：

| 功能 | OCCT 包/类 | RCAD 现状 |
|------|-----------|-----------|
| 曲线点求值 `C(t)` | `Geom_Curve::Value(t)` | ✅ `CurveEval::point_at(t)` |
| 曲线切向量 `C'(t)` | `Geom_Curve::D1(t)` | ✅ `CurveEval::tangent_at(t)` |
| 曲线曲率 | `GeomLProp_CLProps` | ❌ |
| 曲面点求值 `S(u,v)` | `Geom_Surface::Value(u,v)` | ✅ `SurfaceEval::point_at(u,v)` |
| 曲面法向量 | `GeomLProp_SLProps` | ✅ `SurfaceEval::normal_at(u,v)` |
| 曲线投影到曲面 | `GeomAPI_ProjectPointOnSurf` | ❌ |
| 曲线-曲线交点 | `GeomAPI_ExtremaCurveCurve` | ❌ |
| 曲面-曲面交线 | `GeomAPI_IntSS` | ❌ |
| 曲线插值 | `GeomAPI_Interpolate` | ❌ |
| 曲线近似（最小二乘）| `GeomAPI_PointsToBSpline` | ❌ |
---

## 3. 拓扑层对比

### 3.1 拓扑实体

#### OCCT `TopoDS` 实体层次（由简到繁）

```
TopoDS_Shape (基类，携带 Orientation + Location)
├── TopoDS_Vertex      — 点
├── TopoDS_Edge        — 有向边段（关联 Curve3 + 参数范围 + PCurves）
├── TopoDS_Wire        — 有序边环
├── TopoDS_Face        — 有界曲面片（关联 Surface3 + 边界 Wire）
├── TopoDS_Shell       — 面的集合
├── TopoDS_Solid       — 封闭 Shell 构成的体
├── TopoDS_CompSolid   — 多体组合
└── TopoDS_Compound    — 任意形状组合
```

每个 `TopoDS_Shape` 携带：
- `Orientation`：`FORWARD / REVERSED / INTERNAL / EXTERNAL`
- `TopLoc_Location`：变换矩阵（支持实例化共享）
- 子形状迭代：`TopExp_Explorer`

#### RCAD 拓扑现状

```rust
pub struct Vertex { pub point: DVec3 }
pub struct Edge   { pub start: usize, pub end: usize }
pub struct WireEdge { pub idx: usize, pub forward: bool }  // Phase A 新增
pub struct Wire   { pub edges: Vec<WireEdge> }             // 含方向
pub struct Face   { pub outer_wire: Wire, pub inner_wires: Vec<Wire>,
                    pub normal: DVec3, pub triangles: Vec<[usize;3]> }
pub struct Shell  { pub faces: Vec<Face> }
pub struct Solid  { pub shells: Vec<Shell> }
pub struct BRep   { pub vertices: Vec<Vertex>, pub edges: Vec<Edge>,
                    pub solids: Vec<Solid>, pub geom: GeomStore }
```

| 实体 | OCCT | RCAD | 差距说明 |
|------|------|------|---------|
| Vertex | `TopoDS_Vertex` + 容差 | `Vertex { point }` | 缺容差字段（P2）|
| Edge | `TopoDS_Edge` + 参数范围 + 方向 | `Edge { start, end }` + `edge_curve_range` | ✅ 参数范围已补全（Phase A）|
| Wire | `TopoDS_Wire` | `Wire { edges: Vec<WireEdge> }` | ✅ 含方向标志（Phase A）|
| Face | `TopoDS_Face` | `Face { outer_wire, inner_wires, normal }` | ✅ 基本对等 |
| Shell | `TopoDS_Shell` | `Shell { faces }` | ✅ 基本对等 |
| Solid | `TopoDS_Solid` | `Solid { shells }` | ✅ 基本对等 |
| **CompSolid** | `TopoDS_CompSolid` | ❌ | 低优先 |
| **Compound** | `TopoDS_Compound` | ❌ | 低优先 |
| **Location** | `TopLoc_Location` | ❌ | 中优先：实例化、装配 |
| 子形状迭代 | `TopExp_Explorer` | 需手动嵌套 | 中优先：可封装辅助函数 |

---

### 3.2 拓扑属性与元数据

#### OCCT `BRep_TEdge`（边的存储细节）

```
BRep_TEdge:
  Tolerance          f64
  SameParameter      bool    ← 3D曲线与PCurve参数化是否同步
  SameRange          bool    ← 不同PCurve是否使用相同参数范围
  Degenerated        bool    ← 退化边（如球极点）
  Curve3D            Geom_Curve + [t1, t2]
  CurvesOnSurface[]  (Geom2d_Curve, Geom_Surface, TopLoc_Location, [t1, t2])
  PolygonOnSurface[] (可视化网格缓存)
```

#### RCAD `Edge` + `GeomStore` 现状

```rust
Edge { start: usize, end: usize }

GeomStore {
    edge_curve:       Vec<Option<usize>>,      // 指向 Curve3
    edge_curve_range: Vec<Option<[f64; 2]>>,   // [t1, t2] — Phase A 新增
    edge_degenerated: Vec<bool>,               // Phase A 新增
    edge_pcurves:     Vec<Vec<PCurve>>,        // PCurve { surface_idx, curve2d_idx }
    vertex_tolerance: Vec<f64>,               // Phase I 新增 — per-vertex 容差
    edge_tolerance:   Vec<f64>,               // Phase I 新增 — per-edge 容差
    face_tolerance:   Vec<f64>,               // Phase I 新增 — per-face 容差
    // 仍缺: same_parameter, same_range
}
```

**关键字段状态：**

| OCCT 字段 | RCAD | 状态 |
|-----------|------|------|
| **Edge 参数范围 `[t1, t2]`** | `edge_curve_range` | ✅ Phase A |
| `Degenerated` | `edge_degenerated` | ✅ Phase A |
| **`Tolerance` (per edge/vertex/face)** | `vertex_tolerance / edge_tolerance / face_tolerance` + 查询函数 | ✅ Phase I |
| `SameParameter` | ❌ | P2 |
| `SameRange` | ❌ | P2 |

---

### 3.3 拓扑遍历与查询

| 功能 | OCCT | RCAD |
|------|------|------|
| 遍历所有边 | `TopExp_Explorer(shape, EDGE)` | `face_edges(brep, face_idx)` + 手动迭代 | ✅ Phase G |
| 查找边的相邻面 | `TopTools_IndexedDataMapOfShapeListOfShape` | `edge_adjacent_faces(brep, edge_idx)` | ✅ Phase G |
| 查找顶点共享的边 | `TopExp::MapShapesAndAncestors` | `vertex_adjacent_edges(brep, vertex_idx)` | ✅ Phase G |
| 形状比较 | `TopoDS_Shape::IsEqual/IsSame` | ❌ | 低优先 |
| 子形状计数 | `BRepTools::NbFaces` 等 | `face_count / edge_count / vertex_count` | ✅ Phase G |
| 形状有效性检查 | `BRepCheck_Analyzer` | ✅ `check(brep)` Phase C | ✅ |

---

## 4. 建模 API 对比

### 4.1 基本体（Primitives）

| 基本体 | OCCT API | RCAD API | 状态 |
|--------|----------|----------|------|
| 长方体 | `BRepPrimAPI_MakeBox` | `box_brep(origin, x_dir, y_dir, w, h, d)` | ✅ 对等 |
| 球 | `BRepPrimAPI_MakeSphere` | `sphere_brep(center, radius)` | ✅ 对等 |
| 柱 | `BRepPrimAPI_MakeCylinder` | `cylinder_brep(center, axis, ref, r, h)` | ✅ 对等 |
| 锥 | `BRepPrimAPI_MakeCone` | `cone_brep(center, axis, ref, r, h)` | ✅ 对等 |
| 环 | `BRepPrimAPI_MakeTorus` | `torus_brep(center, axis, ref, R, r)` | ✅ 对等 |
| 半空间 | `BRepPrimAPI_MakeHalfSpace` | ❌ | 低优先 |
| 楔体 | `BRepPrimAPI_MakeWedge` | ❌ | 低优先 |
| 棱柱 | `BRepPrimAPI_MakePrism` | ❌ (只有 box) | 中优先 |

### 4.2 曲线 / 曲面构造

| 功能 | OCCT API | RCAD API | 状态 |
|------|----------|----------|------|
| 直线 | `GC_MakeLine` | `line(origin, dir)` | ✅ |
| 圆 | `GC_MakeCircle` | `circle(center, normal, r)` | ✅ |
| 椭圆 | `GC_MakeEllipse` | `ellipse(...)` | ✅ |
| **裁剪弧** | `GC_MakeArcOfCircle` → `TrimmedCurve` | `edge_curve_range [t1,t2]` | ✅ Phase A |
| **插值曲线** | `GeomAPI_Interpolate` | ❌ | 中优先 |
| 平面 | `GC_MakePlane` | `plane(origin, normal)` | ✅ |
| 柱面 | 直接构造 | `cylindrical_surface(...)` | ✅ |
| **边构造** | `BRepBuilderAPI_MakeEdge` | `make_edge(brep, curve, t1, t2, v0, v1)` | ✅ Phase B |
| **线框构造** | `BRepBuilderAPI_MakeWire` | `make_wire(edges)` | ✅ Phase B |
| **面构造** | `BRepBuilderAPI_MakeFace` | `make_face(brep, surface, outer, inner_wires)` | ✅ Phase B |
| **体构造** | `BRepBuilderAPI_MakeSolid` | `make_solid(brep, shells)` | ✅ Phase B |

### 4.3 布尔运算

| 功能 | OCCT API | RCAD API | 状态 |
|------|----------|----------|------|
| 并集 | `BRepAlgoAPI_Fuse` | `boolean_op(Union, a, b)` | ✅ 基本实现 |
| 交集 | `BRepAlgoAPI_Common` | `boolean_op(Intersection, a, b)` | ✅ 基本实现 |
| 差集 | `BRepAlgoAPI_Cut` | `boolean_op(Difference, a, b)` | ✅ 基本实现 |
| 截面线 | `BRepAlgoAPI_Section` | ❌ | 中优先 |
| 多体布尔 | `BRepAlgoAPI_BooleanOperation` (n 个输入) | 仅支持 2 个输入 | 中优先 |
| 非流形结果 | 自动处理 | ❌ | 中优先 |
| 历史（面/边映射）| `BRepAlgoAPI::Modified/Generated` | ❌ | 中优先 |

### 4.4 扫掠 / 拉伸 / 旋转

| 功能 | OCCT API | RCAD | 状态 |
|------|----------|------|------|
| 线性拉伸 | `BRepPrimAPI_MakePrism` | `extrude(profile, direction, distance)` | ✅ Phase B |
| 旋转体 | `BRepPrimAPI_MakeRevol` | `revolve(profile, axis, angle)` | ✅ Phase B |
| 管道扫掠 | `BRepOffsetAPI_MakePipe` | `sweep_pipe(profile_2d, spine)` | ✅ Phase E |
| 变截面扫掠 | `BRepOffsetAPI_MakePipeShell` | ❌ | 中优先 |
| Loft（放样）| `BRepOffsetAPI_ThruSections` | `loft(profiles)` | ✅ Phase E |
| 加厚 | `BRepOffsetAPI_MakeThickSolid` | ❌ | 低优先 |

### 4.5 倒角 / 圆角

| 功能 | OCCT API | RCAD | 状态 |
|------|----------|------|------|
| 等半径圆角 | `BRepFilletAPI_MakeFillet` | `fillet_edge(brep, edge_idx, radius)` | ✅ Phase F（凸边，平面面）|
| 变半径圆角 | `BRepFilletAPI_MakeFillet` (变量版) | ❌ | 低优先 |
| 直倒角 | `BRepFilletAPI_MakeChamfer` | `chamfer_edge(brep, edge_idx, dist)` | ✅ Phase F（凸边，平面面）|

### 4.6 薄壁 / 偏移 / 缝合

| 功能 | OCCT API | RCAD | 状态 |
|------|----------|------|------|
| 偏移面 | `BRepOffset_MakeOffset` | ❌ | 低优先 |
| 抽壳 | `BRepOffsetAPI_MakeThickSolid` | ❌ | 低优先 |
| 开放壳缝合 | `BRepOffsetAPI_Sewing` | ❌ | 中优先 |
| 形状修复 | `ShapeFix_Shape` | ❌ | 中优先 |

---

## 5. 数据交换对比

### STEP 读写

| 实体类型 | OCCT 读 | OCCT 写 | RCAD 读 | RCAD 写 |
|---------|---------|---------|---------|---------|
| 解析几何曲线（LINE/CIRCLE/ELLIPSE）| ✅ | ✅ | ✅ | ✅ |
| B_SPLINE_CURVE_WITH_KNOTS | ✅ | ✅ | ✅ | ✅ Phase D |
| TRIMMED_CURVE | ✅ | ✅ | 解析参数范围 | ✅（以 `edge_curve_range` 导出）|
| 解析曲面（PLANE/CYL/SPHERE…）| ✅ | ✅ | ✅ | ✅ |
| B_SPLINE_SURFACE | ✅ | ✅ | ✅ Phase E | ❌ |
| PCURVE / SURFACE_CURVE | ✅ | ✅ | ✅ | ✅ |
| ADVANCED_FACE | ✅ | ✅ | ✅ | ✅ |
| MANIFOLD_SOLID_BREP | ✅ | ✅ | ✅ | ✅ |
| SHELL_BASED_SURFACE_MODEL | ✅ | ✅ | ✅ | ✅ |
| 装配体（PRODUCT / NAUO）| ✅ | ✅ | 部分 | ✅ Phase D |
| 颜色 / 材质 | ✅ | ✅ | ❌ | ✅ Phase D |
| PMI（标注尺寸）| ✅ | ✅ | ❌ | ❌ |
| 容差信息 | ✅ | ✅ | ❌ | ❌ |

### 其他格式

| 格式 | OCCT | RCAD |
|------|------|------|
| IGES | ✅ 完整 | ❌ |
| OBJ / STL | ✅ | ❌ (仅渲染内部) |
| GLTF | ✅ (7.x) | ❌ |
| DXF | 插件 | ❌ |
| BREP（原生格式）| ✅ | ❌ |

---

## 6. 全局属性与分析

Phase C 完成了核心属性计算，Phase G 补全了曲率分析和拓扑查询，Phase H 添加弧长与惯性张量：

| 功能 | OCCT 包/类 | RCAD |
|------|-----------|------|
| 面积计算 | `GProp_GProps` + `BRepGProp::SurfaceProperties` | ✅ `surface_area(brep)` Phase C |
| 体积计算 | `GProp_GProps` + `BRepGProp::VolumeProperties` | ✅ `volume(brep)` Phase C |
| 质心计算 | `GProp_GProps::CentreOfMass` | ✅ `centroid(brep)` Phase C |
| 惯性矩 | `GProp_GProps::MatrixOfInertia` | ✅ `inertia_tensor(brep)` Phase H |
| 包围盒 | `Bnd_Box` + `BRepBndLib` | ✅ `bounding_box()` Phase A |
| 曲线弧长 | `GCPnts_AbscissaPoint` | ✅ `arc_length(curve, t1, t2)` Phase H |
| 曲面曲率分析 | `BRepLProp_SLProps` / `GeomLProp_SLProps` | ✅ `principal_curvatures / gaussian_curvature / mean_curvature` Phase G |
| 最近点 | `BRepExtrema_DistShapeShape` | ❌ |
| 法向量场 | `BRep_Tool::Normal` | ✅ `SurfaceEval::normal_at` Phase A |
| 拓扑有效性 | `BRepCheck_Analyzer` | ✅ `check(brep)` Phase C |
| 形状连通性 | `TopTools_ConnectedIterator` | ✅ `edge_adjacent_faces / vertex_adjacent_edges` Phase G |

---

## 7. 精度与容差体系

### OCCT 的精度分层

```
Precision::Confusion()      = 1e-7  ← 点重合容差（默认）
Precision::Angular()        = 1e-12 ← 角度容差（弧度）
Precision::Intersection()   = 1e-7  ← 交运算容差
Precision::Approximation()  = 1e-4  ← 离散化近似容差

BRep_Vertex::Tolerance       ← 每个顶点独立容差
BRep_Edge::Tolerance         ← 每条边独立容差（≥ 顶点容差）
BRep_Face::Tolerance         ← 每个面独立容差
```

### RCAD 现状

```rust
// rcad-kernel/src/tolerance.rs — Phase I 新增
pub const CONFUSION: f64 = 1e-7;       // 点重合容差（对应 Precision::Confusion）
pub const ANGULAR: f64 = 1e-12;        // 角度容差
pub const APPROXIMATION: f64 = 1e-4;   // 离散化近似容差

// GeomStore 中的 per-entity 容差字段（Phase I 新增）
pub vertex_tolerance: Vec<f64>,        // 无存储值时回退至 CONFUSION
pub edge_tolerance: Vec<f64>,
pub face_tolerance: Vec<f64>,

// 查询函数：有值→返回存储值；值为 0 或越界→返回 CONFUSION
pub fn vertex_tolerance(brep: &BRep, vertex_idx: usize) -> f64 { ... }
pub fn edge_tolerance(brep: &BRep, edge_idx: usize) -> f64 { ... }
pub fn face_tolerance(brep: &BRep, face_flat_idx: usize) -> f64 { ... }
pub fn model_tolerance(brep: &BRep) -> f64 { ... }  // 返回所有实体容差的最大值
```

**差距影响：**
- 布尔运算中的退化情形（共面、共边、接触）判断不稳定（`SameParameter`/`SameRange` 仍缺失）
- STEP 导入后容差信息已从 `UNCERTAINTY_MEASURE_WITH_UNIT` 写入 GeomStore（Phase J 完成）
- 圆角算法对近零边无保护

---

## 8. 渲染与可视化

| 功能 | OCCT AIS/V3d | RCAD rcad-render |
|------|-------------|-----------------|
| 渲染后端 | OpenGL (旧) / Vulkan (7.x) | **wgpu**（跨平台 + WASM）✅ |
| 实时着色 | Phong/PBR | 基础 Phong ✅ |
| 面高亮 | AIS_Shape | SelectionState + face highlight ✅ |
| 边高亮 | AIS_Shape | 屏幕空间 edge overlay ✅ |
| 拾取 | `BRepIntCurveSurface` ray casting | 光线投射 + 屏幕空间边 ✅ |
| 装配体显示 | AIS_Shape + 子装配 | ❌ |
| 隐线消除 HLR | `HLRBRep` | ✅ `hlr()` + `hlr_to_svg()` Phase D |
| 动画 | `AIS_Animation` | ❌ |
| 文字标注 | `AIS_Text` | ❌ |
| 截面视图 | `Graphic3d_ClipPlane` | ✅ `section_polylines()` Phase C |
| 多视口 | `V3d_Viewer` | ❌ |
| WASM 支持 | ❌（OpenGL-ES 限制） | ✅（wgpu → WebGPU）|

> RCAD 的渲染层在跨平台（特别是 WASM）方面有明显优势；但高级可视化功能仍是空白。

---

## 9. 综合差距汇总表

以下优先级定义：
- 🔴 **P0（阻塞）**：缺少会导致当前功能崩溃或正确性问题
- 🟠 **P1（高优先）**：对 CAD 引擎核心功能不可缺少
- 🟡 **P2（中优先）**：提高 OCCT 兼容性和用户体验
- 🟢 **P3（低优先）**：高级特性，可延后

| 领域 | 缺失功能 | 优先级 | 说明 |
|------|---------|--------|------|
| **几何** | `Edge` 参数范围 `[t1, t2]` | ✅ 已完成 | Phase A |
| **几何** | 曲线求值 `C(t)` / 曲面求值 `S(u,v)` | ✅ 已完成 | Phase A |
| **几何** | `BSplineCurve3` / `BSplineSurface` | ✅ 已完成 | Phase B |
| **几何** | `Curve2d::Ellipse` + `TrimmedCurve2d` | ✅ 已完成 | Phase J（`Ellipse2d` + `curve2d_range`）|
| **几何** | `Curve2d::BSpline` (2D B-Spline PCurve) | ✅ 已完成 | Phase I（`BSplineCurve2`，de Boor 2D）|
| **几何** | 线性扫掠面 / 旋转面 | ✅ 已完成 | Phase K（`LinearExtrusionSurface` / `RevolutionSurface`）|
| **几何** | 面参数域覆盖 | ✅ 已完成 | Phase K（`face_surface_range` + `face_domain()`）|
| **几何** | Bezier 曲线 / 曲面 | 🟡 P2 | |
| **几何** | 偏移曲线 / 偏移面 | 🟡 P2 | 倒角需要 |
| **拓扑** | `WireEdge.forward` 方向标志 | ✅ 已完成 | Phase A |
| **拓扑** | `Degenerated` edge 标记 | ✅ 已完成 | Phase A |
| **拓扑** | per-vertex / per-edge / per-face 容差 | ✅ 已完成 | Phase I（`vertex/edge/face_tolerance`，回退 CONFUSION）|
| **拓扑** | 相邻面查询 API | ✅ 已完成 | Phase G（`edge_adjacent_faces`）|
| **拓扑** | `Compound` / `CompSolid` | 🟢 P3 | 装配体 |
| **建模** | `make_edge/wire/face/solid` | ✅ 已完成 | Phase B |
| **建模** | 线性拉伸 `extrude` | ✅ 已完成 | Phase B |
| **建模** | 旋转体 `revolve` | ✅ 已完成 | Phase B |
| **建模** | 圆角 `MakeFillet` | ✅ 已完成 | Phase F（`fillet_edge`，凸平面边）|
| **建模** | 管道扫掠 `MakePipe` | ✅ 已完成 | Phase E（`sweep_pipe`）|
| **建模** | 倒角 `MakeChamfer` | ✅ 已完成 | Phase F（`chamfer_edge`，凸平面边）|
| **建模** | Loft 放样 | ✅ 已完成 | Phase E（`loft`）|
| **建模** | 加厚 / 抽壳 | 🟢 P3 | |
| **布尔** | 截面线 `Section` | ✅ 已完成 | Phase C（`section_polylines`）|
| **布尔** | 形状历史映射 | 🟡 P2 | 特征树需要 |
| **数据交换** | B-Spline STEP 读写 | ✅ 已完成 | Phase D |
| **数据交换** | 颜色 / 材质 STEP | ✅ 已完成 | Phase D |
| **数据交换** | 装配体 STEP (NAUO) | ✅ 已完成 | Phase D |
| **数据交换** | IGES / OBJ / GLTF | 🟢 P3 | |
| **数据交换** | B_SPLINE_SURFACE STEP 读 | ✅ 已完成 | Phase E |
| **数据交换** | B_SPLINE_SURFACE STEP 写 | ✅ 已完成 | Phase K（`B_SPLINE_SURFACE_WITH_KNOTS` 导出）|
| **数据交换** | 扫掠面 STEP 读 | ✅ 已完成 | Phase K（`SURFACE_OF_LINEAR_EXTRUSION` / `SURFACE_OF_REVOLUTION`）|
| **分析** | 面积 / 体积 / 质心 | ✅ 已完成 | Phase C |
| **分析** | 包围盒 `Bnd_Box` | ✅ 已完成 | Phase A |
| **分析** | `BRepCheck` 有效性 | ✅ 已完成 | Phase C |
| **分析** | 曲率分析 | ✅ 已完成 | Phase G（`principal_curvatures`，解析+数值）|
| **分析** | 曲线弧长 | ✅ 已完成 | Phase H（`arc_length`，Line/Circle 解析，Ellipse/BSpline GL16）|
| **分析** | 惯性矩张量 | ✅ 已完成 | Phase H（`inertia_tensor`，散度定理三角积分）|
| **精度** | 全局 `Precision` 配置 + per-entity 容差 | ✅ 已完成 | Phase I（`CONFUSION/ANGULAR/APPROXIMATION` + `vertex/edge/face_tolerance` 查询）|
| **渲染** | 隐线消除 HLR | ✅ 已完成 | Phase D |
| **渲染** | 截面视图 | ✅ 已完成 | Phase C（section_polylines + SVG）|
| **渲染** | 多视口 | 🟢 P3 | |
| **渲染** | 动画 / 文字标注 | 🟢 P3 | |

---

## 10. 开发路线建议

基于上述差距分析，十一个阶段已全部完成。以下记录各阶段的实际产出，供后续规划参考。

### Phase A — 几何/拓扑基础加固 ✅ 已完成

1. `GeomStore.edge_curve_range: Vec<Option<[f64; 2]>>` — Edge 参数范围
2. `CurveEval` / `SurfaceEval` trait — `point_at`, `tangent_at`, `normal_at`, `default_domain`
3. `WireEdge { idx: usize, forward: bool }` — Wire 边方向
4. `GeomStore.edge_degenerated: Vec<bool>` — 退化边标记
5. `BRep::bounding_box()` — 包围盒

### Phase B — 建模能力扩展 ✅ 已完成

1. `make_edge / make_wire / make_face / make_solid` (`rcad-modeling/brep_builder`)
2. `extrude(profile, direction, distance)` — 线性拉伸
3. `revolve(profile, axis, angle)` — 旋转体
4. `BSplineCurve3 / BSplineSurface` — de Boor 求值；加入 `Curve3` / `Surface3` enum

### Phase C — 算法完善 ✅ 已完成

1. `surface_area(brep)`, `volume(brep)`, `centroid(brep)` (`rcad-kernel/properties`)
2. `check(brep)` — `BRepCheck` 形状有效性检查 (`rcad-algorithms/brep_check`)
3. `section(brep, plane)` / `section_polylines(brep, plane)` — 截面线 (`rcad-algorithms/section`)

### Phase D — 数据交换与高级功能 ✅ 已完成

1. `StepWriter::write_string_colored(brep, &StepColor)` — STEP 颜色导出（`COLOUR_RGB` → `STYLED_ITEM`）
2. `write_assembly(name, &[AssemblyComponent])` — STEP 装配体（`PRODUCT` + `NAUO`）
3. `B_SPLINE_CURVE_WITH_KNOTS` 读写（压缩节点向量 + de Boor 精确导出）
4. `hlr(brep, camera, samples)` + `hlr_to_svg(result, scale, margin)` — 隐线消除 + SVG 输出

### Phase E — 扫掠/放样与 B-Spline 面 STEP 读 ✅ 已完成

1. `B_SPLINE_SURFACE_WITH_KNOTS` STEP 读取 — 解析 2D 控制点网格 + 展开节点向量 → `Surface3::BSpline`；UV 网格三角化用于渲染
2. `loft(profiles: &[Vec<DVec3>])` — 多截面放样；所有截面顶点数须相等；生成封闭 BRep 实体
3. `sweep_pipe(profile_2d: &[DVec2], spine: &[DVec3])` — 管道扫掠；基于 Frenet-like 帧变换 2D 截面到 3D，委托给 `loft`

### Phase F — 倒角与圆角 ✅ 已完成

1. `chamfer_edge(brep, edge_idx, dist)` — 平面凸边倒角；`dist` 为每侧切入距离；返回新 BRep（6→9 面：6 原始 + 1 倒角四边形 + 2 封口三角形）
2. `fillet_edge(brep, edge_idx, radius)` — 圆角；退刀距离 `setback = radius / tan(β/2)`（β 为外侧二面角）；圆角面为 `CylindricalSurface`；返回新 BRep
3. 内部辅助：`find_adjacent_faces`（O(F) wire 扫描）、`setback_direction`（向内方向计算）、`copy_face_remapped`（顶点重映射重建）

### Phase G — 拓扑查询与曲率分析 ✅ 已完成

1. `edge_adjacent_faces(brep, edge_idx)` — 查找共享某条边的所有面（O(F) 扫描）
2. `face_edges(brep, face_idx)` — 获取面外边框的所有边索引
3. `vertex_adjacent_edges(brep, vertex_idx)` — 获取某顶点关联的所有边
4. `face_count / edge_count / vertex_count` — 形状尺寸查询
5. `principal_curvatures(surface, u, v)` — 主曲率 (k1, k2)：解析（Plane/Cylinder/Sphere/Cone/Torus）+ 数值有限差分（BSpline）
6. `gaussian_curvature / mean_curvature` — 由主曲率导出

### Phase H — 弧长与惯性矩张量 ✅ 已完成

1. `arc_length(curve, t1, t2)` — 有符号弧长；`Line3`/`Circle3` 解析精确，`Ellipse3`/`BSplineCurve3` 16 点 Gauss-Legendre 数值积分（有限差分速度场）
2. `InertiaTensor { ixx, iyy, izz, ixy, ixz, iyz }` — 对称 3×3 惯性张量；`to_matrix()` 返回行主矩阵
3. `inertia_tensor(brep)` — 散度定理四面体积分（与 `volume` / `centroid` 同模式）；单位密度，结果乘密度即得物理惯性矩

### Phase I — 2D B-Spline PCurve + 容差体系 ✅ 已完成

1. `BSplineCurve2` — 2D 参数空间非均匀有理 B-spline；de Boor 2D 算法（3 分量齐次坐标）；加入 `Curve2d::BSpline` 变体；类比 OCCT `Geom2d_BSplineCurve`
2. `Curve2dEval` dispatch 更新覆盖三个变体（`Line2d`、`Circle2d`、`BSplineCurve2`）
3. Per-entity 容差字段 `GeomStore.vertex_tolerance / edge_tolerance / face_tolerance` — 与拓扑数组平行的 `Vec<f64>`；缺失/零值回退至 `CONFUSION = 1e-7`
4. 精度常量 `CONFUSION = 1e-7` / `ANGULAR = 1e-12` / `APPROXIMATION = 1e-4`（类比 OCCT `Precision` 类）
5. 查询函数 `vertex_tolerance / edge_tolerance / face_tolerance / model_tolerance`（返回所有实体容差的最大值）

### Phase J — Ellipse2d PCurve + curve2d_range + STEP Curve2d I/O + STEP tolerance import ✅ 已完成

1. `Ellipse2d` — 2D 参数空间椭圆；加入 `Curve2d::Ellipse` 变体；参数化 `center + major_dir·a·cos(t) + minor_dir·b·sin(t)`，域 `[0, 2π]`；类比 OCCT `Geom2d_Ellipse`
2. `GeomStore.curve2d_range: Vec<Option<[f64; 2]>>` — per-PCurve 参数范围；`None` = 自然域；`Some([t1, t2])` 来自 STEP `TRIMMED_CURVE`；类比 `edge_curve_range` 用于 3D
3. STEP Curve2d 导出：`Curve2d::Ellipse` → `ELLIPSE` + `AXIS2_PLACEMENT_2D`；`Curve2d::BSpline` → `B_SPLINE_CURVE_WITH_KNOTS`（2D 控制点）
4. STEP tolerance 导入：`UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(val), ...)` → 填充 `GeomStore.{vertex,edge,face}_tolerance`；缺失时回退 `CONFUSION`

### 下一阶段候选（Phase K）

优先级较高的剩余工作：
- **变截面扫掠** `MakePipeShell` — P2
- **多边同时圆角（corner blending）** — P2

### Phase K — 扫掠面 + 面参数域 + BSplineSurface STEP 导出 ✅ 已完成

1. `LinearExtrusionSurface` — `S(u,v) = profile.point_at(u) + v·direction`；法向 = `tangent(u) × direction`；类比 OCCT `Geom_SurfaceOfLinearExtrusion`；STEP 导入：`SURFACE_OF_LINEAR_EXTRUSION` → `Surface3::LinearExtrusion`
2. `RevolutionSurface` — `S(u,v) = rotate(profile.point_at(v), axis_origin, axis_dir, angle=u)`；u ∈ [0, 2π]，v 来自 profile；法向数值差分；类比 OCCT `Geom_SurfaceOfRevolution`
3. `GeomStore.face_surface_range: Vec<Option<[f64; 4]>>` — 逐面曲面参数域覆盖 `[u1,u2,v1,v2]`；`face_domain()` 优先返回覆盖值，回退 `SurfaceEval::default_domain()`；类比 OCCT `BRep_Face::UVBounds()`
4. `BSplineSurface` STEP 导出：`write_surface` 现在输出 `B_SPLINE_SURFACE_WITH_KNOTS`（含完整控制点格和节点向量）；内核 [u][v] 网格转置为 STEP [v][u] 顺序（原先回退为 PLANE）

---

### 阶段时序（实际完成）

```
Phase A（几何基础）  ████████
Phase B（建模 API）  ░░░░░░░░████████
Phase C（算法）      ░░░░░░░░░░░░░░░░████████
Phase D（交换/高级） ░░░░░░░░░░░░░░░░░░░░░░░░████████
Phase E（扫掠/B面）  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░████████
Phase F（倒角/圆角） ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░████████
Phase G（拓扑/曲率） ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░████████
Phase H（弧长/惯性） ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░████████
Phase I（PCurve/容差）░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░████████
Phase J（Ellipse2d/STEP容差）░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░████████
Phase K（扫掠面/面域/BSpline导出）░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░████████
```

---

*文档更新于 2026-04-04，基于 RCAD Phase K 完成。*
