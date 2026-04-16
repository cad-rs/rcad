# RCAD 开源/专业版拆分与 rmsh 集成设计

> 日期：2026-04-16
> 状态：实施完成

**实施进度**：
- [x] Phase 1: 印记功能拆分设计
- [x] Phase 2: rcad-pro 仓库创建
- [x] Phase 3: GDS/OAS/IO 模块迁移
- [x] Phase 4: rmsh 集成（依赖更新完成，编译通过）
- [x] Phase 5: 专业版功能开发
  - [x] CAE 格式导出（VTK、MSH、INP、CDB）
  - [x] 高级印记功能（多体装配印记、耦合接口检测）
  - [x] 仿真工作流模块（rcad-feflow）

---

## 1. 概述

### 1.1 背景

RCAD 项目需要拆分为开源版本和专业版本，同时支持与 rmsh（FEM 网格生成库）的集成。

### 1.2 目标

1. **开源版本**：提供完整的几何内核能力，支撑 rmsh 的网格生成需求
2. **专业版本**：提供多物理场仿真、CAE 网格格式等专业功能
3. **rmsh 集成**：rmsh 直接依赖 rcad 的 BRep 类型，无循环依赖

### 1.3 仓库架构

```
┌─────────────────────────────────────────────────────┐
│  rcad (开源仓库 - MIT/Apache 2.0)                   │
│                                                     │
│  完整几何内核 + 基础建模 + 数据交换                   │
└────────────────────────────┬────────────────────────┘
                             │ 直接依赖 BRep
                             ▼
┌─────────────────────────────────────────────────────┐
│  rmsh (独立开源库)                                   │
│                                                     │
│  FEM 网格生成（四面体/六面体/边界层）                  │
└─────────────────────────────────────────────────────┘
                             ▲
                             │ 依赖 rcad
┌────────────────────────────┴────────────────────────┐
│  rcad-pro (私有仓库 - 商业许可)                      │
│                                                     │
│  多物理场耦合 + CAE 网格格式 + 半导体设计格式          │
└─────────────────────────────────────────────────────┘
```

---

## 2. 仓库划分

### 2.1 rcad (开源版本)

| 模块 | 功能 | 说明 |
|------|------|------|
| `rcad-kernel` | 几何内核 | B-Rep 拓扑、解析几何、NURBS、曲率、投影 |
| `rcad-modeling` | 建模 API | 基本体、扫掠、圆角、倒角、偏移、薄壁 |
| `rcad-algorithms` | 算法层 | 布尔运算、几何修复、截面、HLR、基础印记 |
| `rcad-step` | STEP 读写 | AP203/AP214 完整支持、装配体 IO |
| `rcad-iges` | IGES 读写 | Type 106 网格桥接 |
| `rcad-render` | 渲染 | wgpu 渲染管线、拾取、显示模式 |
| `rcad-scene` | 场景 | 创建命令状态机、装配体管理 |
| `rcad-constraints` | 约束求解 | 2D/3D 草图约束、草图→BRep |

**开源版许可证**：MIT + Apache 2.0 双许可

### 2.2 rcad-pro (专业版本)

| 模块 | 功能 | 说明 |
|------|------|------|
| `rcad-coupling` | 多物理场耦合 | 高级印记、共享拓扑、耦合界面定义 |
| `rcad-cae-formats` | CAE 网格格式 | VTK、MSH (Gmsh)、INP (Abaqus)、CDB (ANSYS) |
| `rcad-gds` | GDS-II 读写 | 半导体版图格式（从开源版迁移） |
| `rcad-oas` | OASIS 读写 | 半导体版图格式（从开源版迁移） |
| `rcad-io-pro` | 专业 IO | VTK、MSH、INP、CDB 等格式 |
| `rcad-feflow` | 仿真工作流 | 边界条件定义、材料属性管理 |

**专业版许可证**：商业许可

### 2.3 rmsh (独立库)

| 功能 | 说明 |
|------|------|
| 表面网格生成 | 三角形/四边形网格 |
| 体网格生成 | 四面体/六面体网格 |
| 边界层网格 | CFD 边界层、结构化边界层 |
| 网格质量优化 | 平滑、质量改进 |

**依赖关系**：`rmsh → rcad-kernel`（仅依赖几何类型）

---

## 3. rmsh 集成设计

### 3.1 依赖关系

```toml
# rmsh/Cargo.toml
[dependencies]
rcad-kernel = { version = "0.x", path = "../rcad/libs/rcad-kernel" }
```

**关键原则**：rmsh 仅依赖 `rcad-kernel`，不依赖 `rcad-algorithms` 等高级模块。

### 3.2 rmsh 需要的几何能力

| 能力 | 来源 | 说明 |
|------|------|------|
| BRep 拓扑 | `rcad-kernel/topology` | Vertex/Edge/Wire/Face/Shell/Solid |
| 曲面类型 | `rcad-kernel/geom` | Plane/Cylinder/Sphere/Cone/Torus/NURBS |
| 曲线类型 | `rcad-kernel/geom` | Line/Circle/Ellipse/NURBS |
| 曲面求值 | `SurfaceEval::point_at` | 参数域→3D 点 |
| 曲线求值 | `CurveEval::point_at` | 参数→3D 点 |
| 曲面法向 | `SurfaceEval::normal_at` | 用于网格方向 |
| 曲率计算 | `rcad-kernel/curvature` | 用于网格自适应加密 |
| 参数域范围 | `default_domain()` | 曲面参数域边界 |
| 曲面三角化 | `rcad-kernel/triangulate` | 表面网格初始离散 |

### 3.3 rcad-kernel 为 rmsh 提供的 API

```rust
// rmsh 需要的核心接口

/// 遍历 BRep 的所有面
pub fn faces(brep: &BRep) -> impl Iterator<Item = (usize, &Face)>;

/// 获取面的曲面几何
pub fn face_surface(brep: &BRep, face_idx: usize) -> Option<&Surface3>;

/// 获取面的参数域范围
pub fn face_domain(brep: &BRep, face_idx: usize) -> [f64; 4];

/// 遍历面的边界环
pub fn face_wires(brep: &BRep, face_idx: usize) -> Vec<&Wire>;

/// 获取边的曲线几何
pub fn edge_curve(brep: &BRep, edge_idx: usize) -> Option<&Curve3>;

/// 曲面三角化（可选，用于初始离散）
pub fn tessellate_face(face: &Face, surface: &Surface3, tolerance: f64) -> Vec<[DVec3; 3]>;
```

### 3.4 避免循环依赖

```
❌ 错误的依赖：
   rcad-algorithms → rmsh (网格生成)
   rmsh → rcad-kernel (几何类型)
   rcad-algorithms → rcad-kernel
   → 循环：rmsh 既是依赖又是被依赖

✅ 正确的依赖：
   rmsh → rcad-kernel (仅几何类型)
   rcad-algorithms → rcad-kernel (独立)
   rcad-pro → rcad-algorithms + rmsh (应用层组合)
```

**设计原则**：
- rmsh 不依赖 `rcad-algorithms`（布尔运算等高级功能）
- rmsh 不依赖 `rcad-modeling`（建模 API）
- 网格生成算法完全在 rmsh 内部实现

---

## 4. 功能迁移计划

### 4.1 迁移到专业版

| 功能 | 源位置 | 目标位置 | 优先级 |
|------|--------|----------|--------|
| GDS-II 读写 | `rcad-gds` | `rcad-pro/rcad-gds` | 高 |
| OASIS 读写 | `rcad-oas` | `rcad-pro/rcad-oas` | 高 |
| 高级印记 | `rcad-algorithms/imprint` 部分 | `rcad-pro/rcad-coupling` | 中 |
| CAE 格式导出 | 无（新增） | `rcad-pro/rcad-cae-formats` | 中 |

### 4.2 保留在开源版

| 功能 | 位置 | 说明 |
|------|------|------|
| 基础印记 | `rcad-algorithms/imprint` | 简单面印记，rmsh 可用 |
| 去特征化 | `rcad-algorithms/defeature` | 小孔/圆角移除 |
| 几何修复 | `rcad-algorithms/brep_repair` | 缝合、法向修复 |
| 布尔运算 | `rcad-algorithms/builder` | Union/Intersection/Difference |

### 4.3 rcad-io 拆分

`rcad-io` 模块需要拆分为开源基础 IO 和专业格式支持：

**开源版（保留）**：
- `rcad-io` - 统一 IO 接口定义
- 基础格式：STEP、IGES、OBJ、STL
- 统一的 `Reader`/`Writer` trait 定义

**专业版（新增 `rcad-io-pro`）**：
- CAE 网格格式：VTK、MSH (Gmsh)、INP (Abaqus)、CDB (ANSYS)
- 继承开源版的 trait，扩展专业格式支持

```
rcad-io (开源)           - trait 定义 + 基础格式
    ↑
    └── rcad-io-pro (专业) - CAE 网格格式实现
```

### 4.4 印记功能拆分

**基础版（开源）**：
- `imprint_brep_basic(target, tool)` - 简单的两体印记
- 用于共形网格准备
- rmsh 可调用

**高级版（专业）**：
- `imprint_assembly(breps)` - 多体装配印记
- `imprint_with_tolerance(...)` - 容差控制
- `detect_coupling_interfaces(...)` - 耦合界面自动检测
- 用于多物理场耦合仿真

---

## 5. 目录结构

### 5.1 rcad (开源仓库)

```
rcad/
├── Cargo.toml
├── libs/
│   ├── rcad-kernel/
│   │   ├── src/
│   │   │   ├── geom.rs          # 曲线/曲面类型
│   │   │   ├── topology.rs      # B-Rep 拓扑
│   │   │   ├── triangulate.rs   # 曲面三角化
│   │   │   ├── curvature.rs     # 曲率计算
│   │   │   ├── projection.rs    # 投影计算
│   │   │   └── lib.rs
│   │   └── Cargo.toml
│   ├── rcad-modeling/
│   ├── rcad-algorithms/
│   ├── rcad-step/
│   ├── rcad-iges/
│   ├── rcad-render/
│   ├── rcad-scene/
│   ├── rcad-constraints/
│   └── rcad-io/
├── apps/
│   ├── creator-egui/
│   └── creator-iced/
└── examples/
```

### 5.2 rcad-pro (专业版仓库)

```
rcad-pro/
├── Cargo.toml
├── libs/
│   ├── rcad-coupling/
│   │   ├── src/
│   │   │   ├── imprint_advanced.rs
│   │   │   ├── shared_topology.rs
│   │   │   └── coupling_interface.rs
│   │   └── Cargo.toml
│   ├── rcad-cae-formats/
│   │   ├── src/
│   │   │   ├── vtk.rs
│   │   │   ├── gmsh_msh.rs
│   │   │   ├── abaqus_inp.rs
│   │   │   └── ansys_cdb.rs
│   │   └── Cargo.toml
│   ├── rcad-gds/
│   ├── rcad-oas/
│   └── rcad-feflow/
└── apps/
    └── cae-preprocessor/        # CAE 前处理应用
```

### 5.3 rmsh (独立仓库)

```
rmsh/
├── Cargo.toml
├── src/
│   ├── surface_mesh/
│   │   ├── triangle.rs
│   │   └── quad.rs
│   ├── volume_mesh/
│   │   ├── tetrahedron.rs
│   │   └── hexahedron.rs
│   ├── boundary_layer.rs
│   ├── quality.rs
│   └── lib.rs
└── tests/
```

---

## 6. 实施步骤

### Phase 1：印记功能拆分（优先，1 周）

1. 分析现有 `imprint.rs` 功能边界
2. 设计基础版与高级版 API
3. 在 `rcad-algorithms/imprint.rs` 中重构为基础功能
4. 创建 `rcad-pro/rcad-coupling/imprint_advanced.rs` 存放高级功能
5. 添加测试验证拆分正确性

### Phase 2：仓库准备（1 周）

1. 创建 `rcad-pro` 私有仓库
2. 设置 CI/CD（独立于开源版）
3. 确定专业版许可证条款

### Phase 3：模块迁移（2 周）

1. 将 `rcad-gds` 迁移到 `rcad-pro`
2. 将 `rcad-oas` 迁移到 `rcad-pro`
3. 拆分 `rcad-io` 为开源基础 + 专业格式
4. 更新 `rcad` 主仓库的 `Cargo.toml`

### Phase 4：rmsh 集成（3 周）

1. 确定 rmsh 需要的精确 API
2. 在 `rcad-kernel` 中添加/暴露必要接口
3. rmsh 添加对 `rcad-kernel` 的依赖
4. 集成测试

### Phase 5：专业版功能开发（持续）

1. 实现 VTK/MSH/INP 格式导出
2. 完善高级印记功能
3. 开发多物理场耦合支持

---

## 7. 许可证策略

| 仓库 | 许可证 | 说明 |
|------|--------|------|
| rcad | MIT + Apache 2.0 双许可 | 开源，可商用 |
| rmsh | MIT + Apache 2.0 双许可 | 开源，可商用 |
| rcad-pro | 商业许可 | 需购买授权 |

**开源版承诺**：
- 完整的几何内核能力
- 无功能阉割，适合 rmsh 集成
- 持续维护和更新

**专业版价值**：
- 多物理场仿真专用功能
- CAE 软件集成格式
- 半导体设计格式支持
- 技术支持和定制服务

---

## 8. 风险与缓解

| 风险 | 影响 | 缓解措施 |
|------|------|----------|
| 开源版功能过于完整，专业版价值不足 | 商业模式受影响 | 专业版聚焦专业领域，持续增加专有功能 |
| rmsh 与 rcad 版本不同步 | 集成问题 | 使用语义化版本，保证 API 稳定 |
| 印记功能拆分边界模糊 | 维护困难 | 明确定义基础/高级边界，文档化 |
| rcad-io 拆分影响现有用户 | 兼容性问题 | 保持开源版 API 不变，专业版作为扩展 |
| 专业版依赖开源版更新 | 同步成本 | 自动化依赖更新，定期同步 |

---

## 9. 后续工作

1. **印记拆分详细设计**：定义基础版与高级版的功能边界（优先）
2. **rmsh API 需求确认**：与 rmsh 开发者确认精确的几何接口需求
3. **rcad-io 拆分设计**：确定哪些格式留在开源，哪些放入专业版
4. **CAE 格式实现**：VTK/MSH/INP 格式的详细设计
5. **CI/CD 配置**：双仓库的持续集成流程

---

*文档版本：1.1*
*作者：Claude + 用户协作*
*最后更新：2026-04-16*
