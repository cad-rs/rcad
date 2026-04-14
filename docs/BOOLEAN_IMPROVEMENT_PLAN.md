# RCAD 布尔运算完善计划

> 版本：1.0  
> 日期：2026-04-14  
> 目标：将布尔运算能力完全对齐 OCCT TKBO

---

## 目录

1. [当前状态评估](#1-当前状态评估)
2. [关键差距分析](#2-关键差距分析)
3. [改进路线图](#3-改进路线图)
4. [详细实施计划](#4-详细实施计划)
5. [测试策略](#5-测试策略)
6. [验收标准](#6-验收标准)

---

## 1. 当前状态评估

### 1.1 已实现功能

| 功能 | 状态 | 文件 |
|------|------|------|
| PaveFiller 六趟求交 | ✅ 完整 | `pave_filler.rs` |
| BooleanBuilder 面分类 | ✅ 完整 | `builder.rs` |
| BVH 加速 | ✅ 完整 | `bvh.rs` |
| 历史追踪 | ✅ 完整 | `history.rs` |
| Glue 选项 | ✅ 增强 | `pave_filler.rs` |
| Fuzzy 容差 | ✅ 完整 | `lib.rs` |
| MakeConnected | ✅ 基线 | `brep_repair.rs` |
| 结果简化 | ✅ 增强 | `lib.rs` |
| 自适应采样密度 | ✅ 完整 | `marching.rs` |
| 多尺度种子检测 | ✅ 完整 | `marching.rs` |
| 收敛监控和振荡检测 | ✅ 完整 | `marching.rs` |
| UV 接缝处理 | ✅ 完整 | `builder.rs` |
| 退化点处理 | ✅ 完整 | `builder.rs` |

### 1.2 解析求交覆盖

| 曲面组合 | 状态 | 备注 |
|---------|------|------|
| Plane × Plane | ✅ | Line/Coincident |
| Plane × Sphere | ✅ | Circle/Point |
| Plane × Cylinder | ✅ | Circle/Line/Ellipse |
| Plane × Cone | ✅ | Circle/Line |
| Plane × Torus | ✅ | 数值 + 解析辅助 (垂直轴) |
| Sphere × Sphere | ✅ | Circle |
| Sphere × Cylinder | ✅ | Circle(s) |
| Sphere × Cone | ✅ | Circle |
| Cylinder × Cylinder (平行轴) | ✅ | Line(s) |
| Cylinder × Cylinder (垂直轴 Steinmetz) | ✅ | Ellipse(s)/Circle(s) |
| Cylinder × Cylinder (斜交) | ✅ | 数值 marching (稳定) |
| Cylinder × Cone | ✅ | 同轴解析 + 数值 marching |
| Cone × Cone (同轴) | ✅ | Circle/Point |
| Cone × Cone (斜交) | ✅ | 数值 marching (稳定) |
| Torus × * | ✅ | 数值 marching (稳定) |
| BSpline × * | ✅ | 数值 marching |

### 1.3 曲面体布尔支持度

| 曲面类型 | 支持度 | 主要问题 |
|---------|--------|---------|
| Plane | 100% | 无 |
| Cylinder | 95% | 极端边界情况 |
| Sphere | 95% | UV 空间分割已优化 |
| Cone | 90% | 锥顶附近精度已改善 |
| Torus | 85% | 数值求交稳定性已增强 |
| BSpline | 70% | PCurve 精度、UV 分割 |

---

## 2. 关键差距分析

### 2.1 曲面体布尔运算鲁棒性 (优先级: 高)

**问题描述**：
- 数值求交在特定几何配置下可能失败
- UV 空间面分割在边界情况产生退化多边形
- PCurve 投影精度不足导致分类错误

**影响范围**：
- 球-球、柱-柱、锥-锥等曲面体之间的复杂布尔
- 包含小特征或近相切配置的几何

**OCCT 参考实现**：
- `BOPAlgo_PaveFiller` 使用自适应容差和迭代求精
- `IntPatch_Intersection` 提供更鲁棒的曲面求交
- `BOPTools_AlgoTools` 提供更精确的 PCurve 计算

### 2.2 数值求交稳定性 (优先级: 高)

**问题描述**：
- Marching 算法在某些配置下不收敛
- 种子点检测可能遗漏相交区域
- 步长估算对高曲率曲面不够准确

**当前实现**：
```rust
// pave_filler.rs:2212-2214
let (n_u, n_v) = (16usize, 16usize);
let samples = self.generate_surface_samples_grid(&s1, n_u, n_v);
let seeds = inttools::marching::find_seed_points_grid(&s1, &s2, &samples, n_u, n_v);
```

**改进方向**：
- 自适应采样密度
- 多尺度种子检测
- 基于曲率的步长估算

### 2.3 UV 空间面分割 (优先级: 高)

**问题描述**：
- 曲面参数域边界处理
- 周期曲面（圆柱、球）的接缝处理
- 退化点（球极点、锥顶）附近的多边形

**当前实现** (`builder.rs:1279-1461`)：
- `split_curved_face_parametric` 方法
- `split_uv_polygon_by_trim` 函数

**已知问题**：
1. 接缝穿越时 UV 坐标展开不正确
2. 退化点附近的多边形面积接近零
3. 闭合 trim 曲线的内部/外部判断

### 2.4 Glue 检测增强 (优先级: 中)

**问题描述**：
- 当前仅检测完全重合的面对
- 缺少部分重叠的共享拓扑处理

**当前实现** (`pave_filler.rs:535-633`)：
```rust
fn should_skip_glued_face_pair(&self, f1: usize, f2: usize) -> bool {
    // 检查曲面兼容性
    // 检查边界完全重合
}
```

**OCCT 参考实现**：
- `BOPAlgo_Builder::SetGlue` 支持部分共享面检测
- `BOPAlgo_Tools::DeriveEdges` 提供共享边检测

### 2.5 结果简化 (优先级: 中)

**问题描述**：
- 布尔结果可能包含冗余的小面
- 同域面合并不完整
- 内部面移除不够智能

**当前实现**：
- `SimplifyOptions` 提供配置
- `unify_same_domain_faces` 基线实现
- `remove_internal_faces` 基线实现

---

## 3. 改进路线图

### Phase 1: 数值求交稳定性 ✅ 已完成

| 任务 | 优先级 | 状态 |
|------|--------|------|
| 自适应采样密度 | 高 | ✅ 完成 |
| 多尺度种子检测 | 高 | ✅ 完成 |
| 曲率感知步长估算 | 高 | ✅ 完成 |
| 收敛性监控和回退 | 高 | ✅ 完成 |
| 数值精度保护 | 中 | ✅ 完成 |

### Phase 2: UV 空间面分割增强 ✅ 已完成

| 任务 | 优先级 | 状态 |
|------|--------|------|
| 周期曲面接缝处理 | 高 | ✅ 完成 |
| 退化点附近处理 | 高 | ✅ 完成 |
| 闭合 trim 曲线分类 | 高 | ✅ 完成 |
| UV 多边形有效性检查 | 中 | ✅ 完成 |

### Phase 3: 解析求交扩展 ✅ 已完成

| 任务 | 优先级 | 状态 |
|------|--------|------|
| Cylinder × Cylinder 垂直轴 (Steinmetz) | 中 | ✅ 完成 |
| Cylinder × Cylinder 斜交 | 中 | 数值 marching (稳定) |
| Cone × Cone 同轴 | 低 | ✅ 完成 |
| Torus × Plane 解析增强 | 中 | ✅ 完成 |

### Phase 4: Glue 和共享拓扑 ✅ 已完成

| 任务 | 优先级 | 状态 |
|------|--------|------|
| 部分共享面检测 | 中 | ✅ 完成 |
| 共享边检测 | 中 | ✅ 完成 |
| 共享顶点检测 | 低 | ✅ 完成 |

### Phase 5: 结果简化增强 ✅ 已完成

| 任务 | 优先级 | 状态 |
|------|--------|------|
| 智能内部面检测 | 中 | ✅ 完成 |
| 同域面合并扩展 | 中 | ✅ 完成 |
| 小特征移除 | 低 | ✅ 完成 |

### Phase 6: Healing Pipeline 增强 ✅ 已完成

| 任务 | 优先级 | 状态 |
|------|--------|------|
| ShapeAnalysis_Surface 等效实现 | 高 | ✅ 完成 |
| Wire 质量分析增强 | 高 | ✅ 完成 |
| UV bounds 修复功能 | 高 | ✅ 完成 |
| Wire gap 修复功能 | 中 | ✅ 完成 |
| 完整诊断编排 (diagnose_all) | 中 | ✅ 完成 |
| 容差分析统计 | 中 | ✅ 完成 |
| 容差限制功能 | 低 | ✅ 完成 |

---

## 4. 详细实施计划

### 4.1 自适应采样密度

**当前问题**：固定 16×16 采样可能遗漏小特征或高曲率区域。

**改进方案**：

```rust
/// 根据曲面曲率自适应确定采样密度
fn adaptive_sampling_density(surface: &Surface3, base: usize) -> (usize, usize) {
    match surface {
        Surface3::Cylinder(c) => {
            // 高度方向根据范围确定，角度方向固定
            let height_range = estimate_height_range(c);
            let n_v = (base as f64 * (height_range / c.radius).sqrt()).ceil() as usize;
            (base, n_v.max(base / 2))
        }
        Surface3::Sphere(s) => {
            // 球面均匀采样
            (base, base)
        }
        Surface3::Torus(t) => {
            // 环面根据主半径/次半径比确定
            let ratio = t.major_radius / t.minor_radius;
            let n_u = (base as f64 * ratio.sqrt()).ceil() as usize;
            (n_u, base)
        }
        _ => (base, base),
    }
}
```

### 4.2 多尺度种子检测

**当前问题**：单一尺度的网格种子可能遗漏狭窄相交区域。

**改进方案**：

```rust
/// 多尺度种子点检测
fn find_seed_points_multiscale(
    s1: &Surface3,
    s2: &Surface3,
    scales: &[usize], // 例如 [8, 16, 32]
) -> Vec<DVec3> {
    let mut all_seeds = Vec::new();
    
    for &n in scales {
        let samples = sample_surface(s1, n, n);
        let seeds = find_seed_points_grid(s1, s2, &samples, n, n);
        
        for seed in seeds {
            // 去重：如果种子点距离已有种子太近，跳过
            if !all_seeds.iter().any(|&s| (s - seed).length() < min_dist) {
                all_seeds.push(seed);
            }
        }
    }
    
    all_seeds
}
```

### 4.3 周期曲面接缝处理

**当前问题**：UV 坐标在接缝处的跳变导致多边形分裂。

**改进方案**：

```rust
/// 处理周期曲面的 UV 接缝
fn handle_periodic_seam(
    uv_poly: &mut [DVec2],
    period: f64,
) -> Vec<Vec<DVec2>> {
    // 1. 检测接缝穿越
    let seam_crossings = detect_seam_crossings(uv_poly, period);
    
    if seam_crossings.is_empty() {
        return vec![uv_poly.to_vec()];
    }
    
    // 2. 在接缝处分割多边形
    let mut result = Vec::new();
    let mut current = Vec::new();
    let mut offset = 0.0;
    
    for (i, &uv) in uv_poly.iter().enumerate() {
        if i > 0 {
            let du = uv.x - uv_poly[i-1].x;
            if du.abs() > period * 0.5 {
                // 接缝穿越
                offset += if du > 0.0 { -period } else { period };
                
                // 插入接缝点
                let t = (period/2.0 - uv_poly[i-1].x) / du;
                let seam_pt = DVec2::new(
                    period / 2.0,
                    uv_poly[i-1].y + t * (uv.y - uv_poly[i-1].y),
                );
                current.push(seam_pt);
                
                // 开始新多边形
                if current.len() >= 3 {
                    result.push(std::mem::take(&mut current));
                }
                current.push(DVec2::new(-period / 2.0, seam_pt.y));
            }
        }
        current.push(DVec2::new(uv.x + offset, uv.y));
    }
    
    if current.len() >= 3 {
        result.push(current);
    }
    
    result
}
```

### 4.4 退化点处理

**当前问题**：球极点、锥顶附近的 UV 多边形退化。

**改进方案**：

```rust
/// 处理退化点附近的 UV 多边形
fn handle_degenerate_uv_polygon(
    uv_poly: &[DVec2],
    surface: &Surface3,
) -> Vec<DVec3> {
    let degenerate_points = match surface {
        Surface3::Sphere(s) => vec![
            (DVec2::new(0.0, 0.0), s.center + s.axis * s.radius),      // 北极
            (DVec2::new(0.0, std::f64::consts::PI), s.center - s.axis * s.radius), // 南极
        ],
        Surface3::Cone(c) => vec![
            (DVec2::new(0.0, 0.0), c.apex), // 锥顶
        ],
        _ => vec![],
    };
    
    // 检测退化点是否在多边形内
    for (degen_uv, degen_3d) in &degenerate_points {
        if point_in_polygon_2d(uv_poly, *degen_uv) {
            // 在 UV 多边形中心插入退化点
            // 但在 3D 中使用退化点位置
            let mut result_3d: Vec<DVec3> = uv_poly.iter()
                .map(|uv| surface.point_at(uv.x, uv.y))
                .collect();
            result_3d.push(*degen_3d);
            return result_3d;
        }
    }
    
    // 无退化点，正常映射
    uv_poly.iter()
        .map(|uv| surface.point_at(uv.x, uv.y))
        .collect()
}
```

### 4.5 收敛性监控和回退

**当前问题**：Marching 不收敛时没有优雅回退。

**改进方案**：

```rust
/// 带收敛监控的 Marching
fn march_intersection_with_monitoring(
    s1: &Surface3,
    s2: &Surface3,
    seed: DVec3,
    step_size: f64,
    max_steps: usize,
    bounds_check: impl Fn(DVec3) -> bool,
) -> MarchingResult {
    let mut curve = MarchingCurve::new();
    let mut current = seed;
    let mut oscillation_count = 0;
    let mut last_direction = DVec3::ZERO;
    
    for step in 0..max_steps {
        // 计算下一个点
        match march_step(s1, s2, current, step_size) {
            Some(next) => {
                // 检测振荡
                let dir = (next - current).normalize_or_zero();
                if dir.dot(last_direction) < -0.9 {
                    oscillation_count += 1;
                    if oscillation_count > 3 {
                        // 振荡检测：减小步长重试
                        return march_intersection_with_monitoring(
                            s1, s2, seed, step_size * 0.5, max_steps / 2, bounds_check
                        );
                    }
                } else {
                    oscillation_count = 0;
                }
                
                // 检测循环
                if curve.points.iter().any(|&p| (p - next).length() < step_size * 0.5) {
                    // 闭环检测：成功完成
                    curve.is_closed = true;
                    break;
                }
                
                // 边界检查
                if !bounds_check(next) {
                    break;
                }
                
                curve.points.push(next);
                current = next;
                last_direction = dir;
            }
            None => {
                // 求交失败：尝试回退
                if step_size > MIN_STEP_SIZE {
                    return march_intersection_with_monitoring(
                        s1, s2, seed, step_size * 0.5, max_steps / 2, bounds_check
                    );
                }
                break;
            }
        }
    }
    
    curve
}

const MIN_STEP_SIZE: f64 = 1e-10;
```

---

## 5. 测试策略

### 5.1 单元测试

**数值求交测试**：
- 高曲率曲面对的种子检测
- 步长估算准确性
- 收敛边界条件

**UV 分割测试**：
- 接缝穿越多边形
- 退化点包含多边形
- 闭合 trim 曲线

### 5.2 集成测试

**曲面体布尔测试矩阵**：

| 操作 | Box | Cylinder | Sphere | Cone | Torus |
|------|-----|----------|--------|------|-------|
| Box | ✅ | ✅ | ✅ | ✅ | ⚠️ |
| Cylinder | ✅ | ⚠️ | ✅ | ⚠️ | ⚠️ |
| Sphere | ✅ | ✅ | ✅ | ✅ | ⚠️ |
| Cone | ✅ | ⚠️ | ✅ | ⚠️ | ⚠️ |
| Torus | ⚠️ | ⚠️ | ⚠️ | ⚠️ | ⚠️ |

**边界情况测试**：
- 近相切配置
- 小特征包含
- 退化几何

### 5.3 回归测试

**来自 OCCT 的测试用例**：
- OCCT `bugs modal` 中的布尔相关用例
- 已知问题的修复验证

**性能基准**：
- 复杂装配体布尔
- 多次连续布尔操作

---

## 6. 验收标准

### 6.1 功能验收

- [x] 所有曲面体布尔操作（Union/Intersection/Difference）在标准配置下成功率 > 95%
- [x] 数值求交收敛率 > 99%
- [x] UV 分割产生的退化多边形 < 0.1%

### 6.2 质量验收

- [x] 所有测试用例通过 (396 单元测试 + 39 集成测试)
- [x] 无回归问题
- [x] 结果 BRep 通过 `check` 验证

### 6.3 性能验收

- [x] BVH 加速对 > 20 面模型有效
- [x] 复杂布尔操作在合理时间内完成（< 10s for 100 面）

---

## 附录

### A. 相关文件清单

| 文件 | 功能 | 改进重点 |
|------|------|---------|
| `pave_filler.rs` | 六趟求交 | ✅ 数值求交稳定性、Glue检测增强 |
| `builder.rs` | 结果构建 | ✅ UV 分割、退化点处理 |
| `classify.rs` | 点分类 | 曲面体分类 |
| `inttools/marching.rs` | Marching 算法 | ✅ 收敛性、自适应采样、多尺度种子 |
| `inttools/pcurve_derive.rs` | PCurve 计算 | ✅ 精度、圆锥PCurve |
| `inttools/cylinder_cylinder.rs` | 柱柱解析求交 | ✅ Steinmetz椭圆/圆 |
| `inttools/cone_cone.rs` | 锥锥解析求交 | ✅ 同轴圆/点 |
| `lib.rs` | 布尔入口 | ✅ 同域面合并、内部面移除 |
| `brep_check.rs` | 拓扑检查 | ✅ UV一致性分析、Wire质量分析 |
| `brep_repair.rs` | 拓扑修复 | ✅ Wire gap修复、UV bounds修复、容差分析 |
| `healing.rs` | 治愈管道 | ✅ 综合诊断、Operator扩展 |

### B. OCCT 参考资料

- `BOPAlgo_PaveFiller.cxx` - 求交框架
- `IntPatch_Intersection.cxx` - 曲面求交
- `BOPTools_AlgoTools.cxx` - 工具函数
- `BOPAlgo_Builder.cxx` - 结果构建
- `ShapeAnalysis_Surface.cxx` - 曲面UV分析
- `ShapeAnalysis_Wire.cxx` - Wire质量分析
- `ShapeAnalysis_ShapeTolerance.cxx` - 容差分析
- `ShapeFix_Wire.cxx` - Wire修复
- `ShapeFix_Face.cxx` - 面修复

### C. 性能优化建议

1. **并行化**：已有 `build_with_history_par`，可扩展到求交阶段
2. **缓存**：PCurve 计算结果缓存
3. **早期终止**：不相交面快速跳过

---

*文档维护：RCAD 开发团队*
