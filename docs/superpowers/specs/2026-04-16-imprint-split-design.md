# 印记功能拆分设计

> 日期：2026-04-16
> 状态：设计完成，待实施

---

## 1. 概述

将 `rcad-algorithms/imprint.rs` 的功能拆分为开源基础版和专业高级版。

---

## 2. 当前功能分析

### 2.1 现有 API

| API | 行数 | 功能 | 复杂度 |
|-----|------|------|--------|
| `imprint_brep(target, tool)` | 83-173 | 两体印记 | 中 |
| `detect_gaps_overlaps(a, b, tol)` | 364-436 | 间隙/重叠检测 | 中 |
| `min_distance(a, b)` | 446-503 | BVH最小距离 | 中 |

### 2.2 内部依赖

```
imprint_brep
├── DS::new + PaveFiller::perform  (布尔基础设施)
├── split_face_by_curves
│   ├── split_planar_face_simple
│   └── split_curved_face
│       ├── split_curved_face_legacy
│       └── UV多边形分割函数
└── triangulate_polygon

detect_gaps_overlaps
├── collect_faces_with_surfaces
├── aabb (包围盒)
└── closest_point_on_surface (kernel)

min_distance
├── Bvh::build + query_aabb
└── closest_point_on_surface
```

---

## 3. 拆分方案

### 3.1 开源版保留

**文件**：`rcad-algorithms/src/imprint.rs`（保持不变）

**公开 API**：

```rust
/// 两体印记 - 将 tool 边界印记到 target 面上
/// 用于共形网格准备
pub fn imprint_brep(target: &BRep, tool: &BRep) -> ImprintResult;

/// 检测间隙和重叠
pub fn detect_gaps_overlaps(a: &BRep, b: &BRep, tolerance: f64) -> GapOverlapReport;

/// 计算最小距离
pub fn min_distance(a: &BRep, b: &BRep) -> f64;
```

**保留原因**：
- `imprint_brep`：rmsh 需要用于共形网格拓扑准备
- `detect_gaps_overlaps`：基础几何分析功能
- `min_distance`：基础几何分析功能

### 3.2 专业版新增

**文件**：`rcad-pro/libs/rcad-coupling/src/imprint_advanced.rs`

**新增 API**：

```rust
/// 多体装配印记
/// 对装配体中所有组件两两印记，生成共形拓扑
pub fn imprint_assembly(breps: &[BRep]) -> AssemblyImprintResult;

/// 带选项的印记
/// 支持容差控制、并行计算
pub fn imprint_with_options(
    target: &BRep,
    tool: &BRep,
    options: ImprintOptions,
) -> ImprintResult;

/// 耦合界面自动检测
/// 检测装配体中可能耦合的界面对
pub fn detect_coupling_interfaces(
    assembly: &[BRep],
    tolerance: f64,
) -> Vec<CouplingInterface>;

/// 共享拓扑分析
/// 分析印记后的BRep，提取共享边/顶点拓扑
pub fn shared_topology_analysis(brep: &BRep) -> SharedTopologyReport;

/// 印记选项
pub struct ImprintOptions {
    pub tolerance: f64,
    pub parallel: bool,
    pub keep_seam_history: bool,
}

/// 装配印记结果
pub struct AssemblyImprintResult {
    pub brep: BRep,                    // 合并后的印记结果
    pub component_mapping: Vec<usize>, // 结果面 → 源组件
    pub coupling_interfaces: Vec<CouplingInterface>,
}

/// 耦合界面
pub struct CouplingInterface {
    pub face_a: usize,
    pub face_b: usize,
    pub interface_type: InterfaceType,
    pub area: f64,
}

pub enum InterfaceType {
    Coincident,    // 共面
    Gap(f64),      // 间隙（带距离）
    Overlap(f64),  // 重叠（带深度）
}

/// 共享拓扑报告
pub struct SharedTopologyReport {
    pub shared_edges: Vec<(usize, usize)>,  // 共享边对
    pub shared_vertices: Vec<(usize, usize)>, // 共享顶点对
    pub manifold_edges: Vec<usize>,          // 流形边（2个面共享）
    pub non_manifold_edges: Vec<usize>,      // 非流形边（>2个面）
}
```

---

## 4. 实施步骤

### Step 1：验证开源版 API 稳定性

当前 `imprint.rs` 无需修改，验证测试通过即可。

### Step 2：创建 rcad-pro 仓库结构

```
rcad-pro/
├── Cargo.toml
└── libs/
    └── rcad-coupling/
        ├── Cargo.toml
        └── src/
            ├── lib.rs
            ├── imprint_advanced.rs
            ├── coupling_interface.rs
            └── shared_topology.rs
```

### Step 3：实现专业版功能

1. `imprint_assembly` - 遍历所有组件对，调用 `imprint_brep`
2. `imprint_with_options` - 封装 `imprint_brep` + 容差控制
3. `detect_coupling_interfaces` - 基于 `detect_gaps_overlaps` 扩展
4. `shared_topology_analysis` - 基于 `BRepGraph` 分析

### Step 4：添加集成测试

```rust
// rcad-pro/libs/rcad-coupling/tests/imprint_assembly.rs

#[test]
fn test_imprint_three_boxes() {
    let box1 = make_box(DVec3::ZERO, 1.0, 1.0, 1.0);
    let box2 = make_box(DVec3::new(0.9, 0.0, 0.0), 1.0, 1.0, 1.0);
    let box3 = make_box(DVec3::new(0.0, 0.9, 0.0), 1.0, 1.0, 1.0);
    
    let result = imprint_assembly(&[box1, box2, box3]);
    
    assert!(!result.coupling_interfaces.is_empty());
}
```

---

## 5. 依赖关系

```
rcad-algorithms (开源)
├── imprint_brep
├── detect_gaps_overlaps
└── min_distance

rcad-coupling (专业)
├── depends on: rcad-algorithms
├── imprint_assembly → 调用 imprint_brep
├── imprint_with_options → 调用 imprint_brep
├── detect_coupling_interfaces → 调用 detect_gaps_overlaps
└── shared_topology_analysis → 调用 BRepGraph
```

---

## 6. 文件变更清单

| 操作 | 文件 | 说明 |
|------|------|------|
| 保持 | `rcad-algorithms/src/imprint.rs` | 开源版不变 |
| 新增 | `rcad-pro/libs/rcad-coupling/Cargo.toml` | 专业版crate配置 |
| 新增 | `rcad-pro/libs/rcad-coupling/src/lib.rs` | 模块入口 |
| 新增 | `rcad-pro/libs/rcad-coupling/src/imprint_advanced.rs` | 高级印记功能 |
| 新增 | `rcad-pro/libs/rcad-coupling/src/coupling_interface.rs` | 耦合界面检测 |
| 新增 | `rcad-pro/libs/rcad-coupling/src/shared_topology.rs` | 共享拓扑分析 |

---

## 7. 验收标准

1. ✅ 开源版 `imprint_brep` 测试全部通过
2. ✅ 开源版 API 无破坏性变更
3. ⬜ 专业版 `imprint_assembly` 实现完成
4. ⬜ 专业版测试覆盖核心场景

---

*文档版本：1.0*
*最后更新：2026-04-16*
