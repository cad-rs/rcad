# RCAD vs OCCT 对比文档

> 本文档系统比较 RCAD 与 Open CASCADE Technology (OCCT) 的数据结构、功能实现与接口设计，
> 帮助团队明确差距、规划后续开发方向。
>
> **文档状态：** 基于 RCAD 当前代码（2026-04）生成，随代码演进应同步更新。

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
    Line(Line3),     // origin + direction（无限长）
    Circle(Circle3), // center + normal + radius（完整圆）
    Ellipse(Ellipse3),
}
```

| 曲线类型 | OCCT | RCAD | 差距说明 |
|---------|------|------|---------|
| 直线（无限） | `Geom_Line` | `Line3` | ✅ 对等 |
| 圆 | `Geom_Circle` | `Circle3` | ✅ 对等 |
| 椭圆 | `Geom_Ellipse` | `Ellipse3` | ✅ 对等 |
| **裁剪曲线** | `Geom_TrimmedCurve` | ❌ 缺失 | **高优先**：边界边的参数域无法精确表达 |
| B-Spline | `Geom_BSplineCurve` | ❌ 缺失 | 必须支持才能导入通用 STEP |
| Bezier | `Geom_BezierCurve` | ❌ 缺失 | 中优先 |
| 双曲线 | `Geom_Hyperbola` | ❌ 缺失 | 低优先（工程场景罕见）|
| 抛物线 | `Geom_Parabola` | ❌ 缺失 | 低优先 |
| 偏移曲线 | `Geom_OffsetCurve` | ❌ 缺失 | 中优先（倒角需要）|
| 参数求值 | `Value(t)` → Point3 | ❌ 缺失 | **高优先**：布尔/圆角/偏移依赖 |
| 参数域 | `FirstParameter/LastParameter` | ❌ 缺失 | **高优先**：Edge 参数范围 |

**当前最大缺陷：`Curve3` 只存储几何形状，没有参数域信息。**  
OCCT 中一条边的实际曲线范围由 `BRep_Edge::Range()` 提供（对应 `TrimmedCurve` 的 `[t1, t2]`）；
RCAD 目前通过顶点位置隐式推断，无法精确表达弧长不到一整圈的圆弧边。

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
| B-Spline 面 | `Geom_BSplineSurface` | ❌ 缺失 | **高优先**（自由曲面） |
| Bezier 面 | `Geom_BezierSurface` | ❌ 缺失 | 中优先 |
| 线性扫掠面 | `Geom_SurfaceOfLinearExtrusion` | ❌ 缺失 | 中优先（拉伸算法需要）|
| 旋转面 | `Geom_SurfaceOfRevolution` | ❌ 缺失 | 中优先（旋转算法需要）|
| 偏移面 | `Geom_OffsetSurface` | ❌ 缺失 | 低优先 |
| (u,v) 求值 | `Value(u,v)` → Point3 | ❌ 缺失 | **高优先**：渲染、布尔需要 |
| 参数域查询 | `Bounds(u1,u2,v1,v2)` | ❌ 缺失 | **高优先** |

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
}
```

| 类型 | OCCT | RCAD | 差距 |
|------|------|------|------|
| 2D 直线 | `Geom2d_Line` | `Line2d` | ✅ |
| 2D 圆 | `Geom2d_Circle` | `Circle2d` | ✅ |
| **2D 椭圆** | `Geom2d_Ellipse` | ❌ | 中优先（椭圆边在椭球面的 PCurve）|
| 2D B-Spline | `Geom2d_BSplineCurve` | ❌ | 高优先（B-Spline 面上的边必须用）|
| 2D 裁剪曲线 | `Geom2d_TrimmedCurve` | ❌ | 高优先（PCurve 的参数域）|

---

### 2.4 几何评估与局部属性

OCCT 提供丰富的几何计算 API，RCAD 目前几乎全部缺失：

| 功能 | OCCT 包/类 | RCAD 现状 |
|------|-----------|-----------|
| 曲线点求值 `C(t)` | `Geom_Curve::Value(t)` | ❌ |
| 曲线切向量 `C'(t)` | `Geom_Curve::D1(t)` | ❌ |
| 曲线曲率 | `GeomLProp_CLProps` | ❌ |
| 曲面点求值 `S(u,v)` | `Geom_Surface::Value(u,v)` | ❌ |
| 曲面法向量 | `GeomLProp_SLProps` | ❌ |
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
pub struct Wire   { pub edges: Vec<usize> }
pub struct Face   { pub outer_wire: Wire, pub inner_wires: Vec<Wire>,
                    pub normal: DVec3, pub triangles: Vec<[usize;3]> }
pub struct Shell  { pub faces: Vec<Face> }
pub struct Solid  { pub shells: Vec<Shell> }
pub struct BRep   { pub vertices: Vec<Vertex>, pub edges: Vec<Edge>,
                    pub solids: Vec<Solid>, pub geom: GeomStore }
```

| 实体 | OCCT | RCAD | 差距说明 |
|------|------|------|---------|
| Vertex | `TopoDS_Vertex` + 容差 | `Vertex { point }` | RCAD 缺容差字段 |
| Edge | `TopoDS_Edge` + 参数范围 + 方向 | `Edge { start, end }` | **缺参数范围 `[t1,t2]`；缺方向标志** |
| Wire | `TopoDS_Wire` | `Wire { edges }` | 缺边方向（每条边可正/反向）|
| Face | `TopoDS_Face` | `Face { outer_wire, inner_wires, normal }` | ✅ 基本对等；缺容差 |
| Shell | `TopoDS_Shell` | `Shell { faces }` | ✅ 基本对等 |
| Solid | `TopoDS_Solid` | `Solid { shells }` | ✅ 基本对等 |
| **CompSolid** | `TopoDS_CompSolid` | ❌ | 低优先 |
| **Compound** | `TopoDS_Compound` | ❌ | 低优先 |
| **Orientation** | `TopAbs_Orientation` | ❌ (隐式) | **高优先**：布尔/圆角依赖 |
| **Location** | `TopLoc_Location` | ❌ | 中优先：实例化、装配 |
| 子形状迭代 | `TopExp_Explorer` | ❌ | 高优先：算法遍历 |

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
    edge_curve: Vec<Option<usize>>,      // 指向 Curve3，无参数范围
    edge_pcurves: Vec<Vec<PCurve>>,      // PCurve { surface_idx, curve2d_idx }
    // 缺: tolerance, same_parameter, same_range, degenerated
}
```

**关键缺失字段：**

| OCCT 字段 | RCAD | 影响 |
|-----------|------|------|
| `Tolerance` (per edge/vertex/face) | ❌ | 精度体系缺失，布尔操作容差判断困难 |
| `SameParameter` | ❌ | 导出/导入 STEP 精度损失 |
| `SameRange` | ❌ | PCurve 一致性检验缺失 |
| `Degenerated` | ❌ | 退化边（如极点）处理不安全 |
| **Edge 参数范围 `[t1, t2]`** | ❌ | **最高优先**：弧段无法精确表示 |

---

### 3.3 拓扑遍历与查询

| 功能 | OCCT | RCAD |
|------|------|------|
| 遍历所有边 | `TopExp_Explorer(shape, EDGE)` | 需手动嵌套遍历 `solid→shell→face→wire` |
| 查找边的相邻面 | `TopTools_IndexedDataMapOfShapeListOfShape` | ❌ 无内建 API |
| 查找顶点共享的边 | `TopExp::MapShapesAndAncestors` | ❌ |
| 形状比较 | `TopoDS_Shape::IsEqual/IsSame` | ❌ |
| 子形状计数 | `BRepTools::NbFaces` 等 | ❌ |
| 形状有效性检查 | `BRepCheck_Analyzer` | ❌ |

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
| **裁剪弧** | `GC_MakeArcOfCircle` → `TrimmedCurve` | ❌ | **高优先** |
| **插值曲线** | `GeomAPI_Interpolate` | ❌ | 高优先 |
| 平面 | `GC_MakePlane` | `plane(origin, normal)` | ✅ |
| 柱面 | 直接构造 | `cylindrical_surface(...)` | ✅ |
| **边构造** | `BRepBuilderAPI_MakeEdge` | ❌ 用户无法直接建边 | **高优先** |
| **线框构造** | `BRepBuilderAPI_MakeWire` | ❌ | 高优先 |
| **面构造** | `BRepBuilderAPI_MakeFace` | ❌ | 高优先 |
| **壳构造** | `BRepBuilderAPI_MakeShell` | ❌ | 中优先 |
| **体构造** | `BRepBuilderAPI_MakeSolid` | ❌ | 中优先 |

> **说明：** RCAD 目前只能通过 `BRep::from_primitive` 创建预定义形状；
> 用户无法从曲线/曲面自由组装拓扑形状，这是与 OCCT 最大的建模能力差距。

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
| 线性拉伸 | `BRepPrimAPI_MakePrism` | ❌ | **高优先** |
| 旋转体 | `BRepPrimAPI_MakeRevol` | ❌ | **高优先** |
| 管道扫掠 | `BRepOffsetAPI_MakePipe` | ❌ | 高优先 |
| 变截面扫掠 | `BRepOffsetAPI_MakePipeShell` | ❌ | 中优先 |
| Loft（放样）| `BRepOffsetAPI_ThruSections` | ❌ | 中优先 |
| 加厚 | `BRepOffsetAPI_MakeThickSolid` | ❌ | 低优先 |

### 4.5 倒角 / 圆角

| 功能 | OCCT API | RCAD | 状态 |
|------|----------|------|------|
| 等半径圆角 | `BRepFilletAPI_MakeFillet` | ❌ | 高优先 |
| 变半径圆角 | `BRepFilletAPI_MakeFillet` (变量版) | ❌ | 低优先 |
| 直倒角 | `BRepFilletAPI_MakeChamfer` | ❌ | 中优先 |

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
| B_SPLINE_CURVE_WITH_KNOTS | ✅ | ✅ | 解析但忽略 | ❌ |
| TRIMMED_CURVE | ✅ | ✅ | 解析但不建拓扑 | ❌ |
| 解析曲面（PLANE/CYL/SPHERE…）| ✅ | ✅ | ✅ | ✅ |
| B_SPLINE_SURFACE | ✅ | ✅ | ❌ | ❌ |
| PCURVE / SURFACE_CURVE | ✅ | ✅ | ✅ | ✅ |
| ADVANCED_FACE | ✅ | ✅ | ✅ | ✅ |
| MANIFOLD_SOLID_BREP | ✅ | ✅ | ✅ | ✅ |
| SHELL_BASED_SURFACE_MODEL | ✅ | ✅ | ✅ | ✅ |
| 装配体（PRODUCT hierarchy）| ✅ | ✅ | 部分 | 部分 |
| 颜色 / 材质 | ✅ | ✅ | ❌ | ❌ |
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

这是 RCAD 目前完全空白的一个领域，对 CAE 应用极为重要：

| 功能 | OCCT 包/类 | RCAD |
|------|-----------|------|
| 面积计算 | `GProp_GProps` + `BRepGProp::SurfaceProperties` | ❌ |
| 体积计算 | `GProp_GProps` + `BRepGProp::VolumeProperties` | ❌ |
| 质心计算 | `GProp_GProps::CentreOfMass` | 仅顶点平均 |
| 惯性矩 | `GProp_GProps::MatrixOfInertia` | ❌ |
| 包围盒 | `Bnd_Box` + `BRepBndLib` | ❌（需手动遍历顶点）|
| 曲线弧长 | `GCPnts_AbscissaPoint` | ❌ |
| 曲面曲率分析 | `BRepLProp_SLProps` | ❌ |
| 最近点 | `BRepExtrema_DistShapeShape` | ❌ |
| 法向量场 | `BRep_Tool::Normal` | 仅存储于 `Face.normal` |
| 拓扑有效性 | `BRepCheck_Analyzer` | ❌ |
| 形状连通性 | `TopTools_ConnectedIterator` | ❌ |

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
// rcad-kernel/src/lib.rs
// 仅有一个全局常量（在 BRep::is_analytic_valid 等检查中隐式用到）
// 无 per-vertex / per-edge / per-face 容差字段
// 无全局精度配置 API
```

**差距影响：**
- 布尔运算中的退化情形（共面、共边、接触）判断不稳定
- STEP 导入后容差信息丢失，无法做 ShapeFix
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
| 隐线消除 HLR | `HLRBRep` | ❌ |
| 动画 | `AIS_Animation` | ❌ |
| 文字标注 | `AIS_Text` | ❌ |
| 截面视图 | `Graphic3d_ClipPlane` | ❌ |
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
| **几何** | `Edge` 参数范围 `[t1, t2]` | 🔴 P0 | 弧段边无法精确表示 |
| **几何** | 曲线求值 `C(t)` / 曲面求值 `S(u,v)` | 🔴 P0 | 布尔、圆角、偏移依赖 |
| **几何** | `TrimmedCurve3` | 🟠 P1 | 任意弧段边构造 |
| **几何** | `BSplineCurve3` / `BSplineSurface` | 🟠 P1 | 自由曲面 STEP 导入必须 |
| **几何** | `Curve2d::Ellipse` + `TrimmedCurve2d` | 🟡 P2 | 椭圆体 PCurve |
| **几何** | `SurfaceOfRevolution` / `LinearExtrusion` | 🟡 P2 | 扫掠算法的输出曲面 |
| **拓扑** | `Orientation` 标志（Edge / Face 级别）| 🔴 P0 | Wire 方向依赖隐式推断，脆弱 |
| **拓扑** | `Degenerated` edge 标记 | 🟠 P1 | 球极点等退化边 |
| **拓扑** | per-vertex / per-edge 容差 | 🟠 P1 | 精度体系完整性 |
| **拓扑** | `TopExp_Explorer` 等遍历工具 | 🟠 P1 | 算法层需要 |
| **拓扑** | 相邻面查询 API | 🟡 P2 | 圆角需要 |
| **拓扑** | `Compound` / `CompSolid` | 🟢 P3 | 装配体 |
| **建模** | `MakeEdge(curve, t1, t2)` | 🟠 P1 | 用户建边入口 |
| **建模** | `MakeWire` / `MakeFace` / `MakeSolid` | 🟠 P1 | 自由形体构造 |
| **建模** | 线性拉伸 `MakePrism` | 🟠 P1 | 最常用建模操作 |
| **建模** | 旋转体 `MakeRevol` | 🟠 P1 | |
| **建模** | 圆角 `MakeFillet` | 🟠 P1 | |
| **建模** | 管道扫掠 `MakePipe` | 🟡 P2 | |
| **建模** | 倒角 `MakeChamfer` | 🟡 P2 | |
| **建模** | Loft 放样 | 🟡 P2 | |
| **建模** | 加厚 / 抽壳 | 🟢 P3 | |
| **布尔** | 截面线 `Section` | 🟡 P2 | |
| **布尔** | 形状历史映射 | 🟡 P2 | 特征树需要 |
| **数据交换** | B-Spline STEP 写出 | 🟠 P1 | |
| **数据交换** | 颜色 / 材质 STEP | 🟡 P2 | |
| **数据交换** | 装配体 STEP | 🟡 P2 | |
| **数据交换** | IGES / OBJ / GLTF | 🟢 P3 | |
| **分析** | 面积 / 体积 / 质心 | 🟠 P1 | CAE 前处理必须 |
| **分析** | 包围盒 `Bnd_Box` | 🟠 P1 | 布尔加速、拾取 |
| **分析** | 曲率分析 | 🟡 P2 | |
| **分析** | `BRepCheck_Analyzer` 有效性 | 🟡 P2 | |
| **精度** | 全局 `Precision` 配置 | 🟡 P2 | |
| **渲染** | 隐线消除 HLR | 🟢 P3 | |
| **渲染** | 截面视图 | 🟢 P3 | |

---

## 10. 开发路线建议

基于上述差距分析，建议分四个阶段推进：

### Phase A — 几何/拓扑基础加固（P0 + P1 核心）

**目标：** 使 BRep 模型在数学上完备，为后续算法提供正确基础。

1. **`Edge` 增加参数范围 `[t1, t2]`**
   - 在 `GeomStore` 或 `Edge` 结构体增加 `edge_curve_range: Vec<Option<[f64; 2]>>`
   - 更新 `create_*` 函数填写正确的参数范围
   - 使 STEP reader/writer 读写 `Edge_Curve` 的参数范围

2. **几何求值 trait `CurveEval` / `SurfaceEval`**
   ```rust
   pub trait CurveEval {
       fn point_at(&self, t: f64) -> DVec3;
       fn tangent_at(&self, t: f64) -> DVec3;
       fn domain(&self) -> [f64; 2];
   }
   pub trait SurfaceEval {
       fn point_at(&self, u: f64, v: f64) -> DVec3;
       fn normal_at(&self, u: f64, v: f64) -> DVec3;
       fn domain(&self) -> [f64; 4]; // [u1, u2, v1, v2]
   }
   ```

3. **`TrimmedCurve3` 类型**
   ```rust
   pub struct TrimmedCurve3 { pub basis: Box<Curve3>, pub t1: f64, pub t2: f64 }
   // 或者在 Curve3 enum 增加 Trimmed 变体
   pub enum Curve3 { Line(..), Circle(..), Ellipse(..), Trimmed(TrimmedCurve3) }
   ```

4. **`Orientation` 标志**
   - 在 `Wire` 的边引用上增加方向标志：`Vec<(usize, bool)>` （边 index + 是否正向）
   - 或专门的 `OrientedEdge { edge_idx: usize, forward: bool }` 类型（已在写出器中存在）

5. **`Degenerated` 标记**
   - `GeomStore.edge_degenerated: Vec<bool>`

6. **包围盒计算**
   ```rust
   // libs/rcad-kernel/src/lib.rs
   impl BRep {
       pub fn bounding_box(&self) -> [DVec3; 2]  // [min, max]
   }
   ```

---

### Phase B — 建模能力扩展（P1 建模操作）

**目标：** 用户可以从曲线/曲面自由组装形状，支持拉伸、旋转。

1. **`BRepBuilder` — 自由建模入口**（对应 OCCT `BRepBuilderAPI`）
   ```rust
   // libs/rcad-modeling/src/brep_builder.rs
   pub fn make_edge(brep: &mut BRep, curve: Curve3, t1: f64, t2: f64,
                    v0: usize, v1: usize) -> usize
   pub fn make_wire(edges: Vec<OrientedEdge>) -> Wire
   pub fn make_face(brep: &mut BRep, surface: Surface3, outer: Wire,
                    inner_wires: Vec<Wire>) -> usize
   pub fn make_solid(brep: &mut BRep, shells: Vec<Shell>) -> usize
   ```

2. **线性拉伸 `extrude(profile_face, direction, distance)`**
   - 输入：一个或多个闭合 Wire / Face
   - 输出：BRep Solid
   - 对应 OCCT `BRepPrimAPI_MakePrism`

3. **旋转体 `revolve(profile, axis, angle)`**
   - 输入：Wire / Face + 旋转轴 + 角度
   - 输出：BRep Solid
   - 对应 OCCT `BRepPrimAPI_MakeRevol`

4. **`BSplineCurve3` / `BSplineSurface` 基础支持**
   - 先只支持数据存储（control_points, knots, weights）
   - 再实现 `CurveEval` / `SurfaceEval`
   - 使 STEP reader 能正确填充而不是忽略

---

### Phase C — 算法完善（P1 算法 + P2 分析）

**目标：** 圆角、全局属性、形状分析。

1. **圆角 `fillet(brep, edge_idx, radius)`**

2. **面积 / 体积 / 质心计算**
   ```rust
   // libs/rcad-kernel/src/properties.rs
   pub fn surface_area(brep: &BRep) -> f64
   pub fn volume(brep: &BRep) -> f64
   pub fn centroid(brep: &BRep) -> DVec3
   ```

3. **`BRepCheck` — 形状有效性检查**

4. **截面线 `section(brep, plane)` → Wire**

---

### Phase D — 数据交换与高级功能（P2/P3）

**目标：** 更完整的 STEP 覆盖、IGES、装配体、可视化增强。

1. **B-Spline STEP 读写完整支持**
2. **STEP 颜色 / 材质**
3. **STEP 装配体（PRODUCT hierarchy）**
4. **HLR 隐线消除**
5. **多视口 / 截面视图**

---

### 阶段时序建议

```
Phase A（几何基础）  ████████░░░░░░░░░░░░░░░░░░░░
Phase B（建模 API）  ░░░░░░░░████████░░░░░░░░░░░░
Phase C（算法）      ░░░░░░░░░░░░░░░░████████░░░░
Phase D（交换/高级） ░░░░░░░░░░░░░░░░░░░░░░░░████
```

Phase A 和 Phase B 的前半段（`make_edge` / `make_wire` / `make_face`）可以并行推进，
因为它们是独立的数据结构扩展。Phase B 后半段（拉伸、旋转）依赖 Phase A 的几何求值 trait 完成。

---

*文档生成于 2026-04-03，基于 RCAD commit `7f603d5`。*
