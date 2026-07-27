# DS OCCT Alignment — Handover 2026-07-27

## 当前状态

DS 字段数 **42 → 31**（已移除 11 个）

### 已完成的 Phase 2
所有 `let edge = &ds.edges[ei]` 绑定和 `ds.TYPE[IDX].FIELD` 内联访问 → `.get().unwrap()` / 访问器方法。所有外部代码不再直接索引 DS 并行数组。
（剩余的直接索引全部在 DS 内部实现文件 mod.rs/types.rs/topods_builder.rs 中，是预期的。）

### 已完成的 Phase 4
`my_images` / `my_origins` / `wire_images` / `shell_images` / `solid_images` — 死字段，直接从 DS 删除。

### 已完成的 Phase 3
| 字段 | 处理方式 |
|------|---------|
| `fuzzy_tol` | → PaveFiller + BOPDS_Iterator 参数 |
| `a_vertex_count` / `a_edge_count` / `a_face_count` | → 从 shape_info 计算的方法 |
| `ff_points` | 死字段，删除 |
| `section_edge_refs` | → PaveFiller |
| `same_domain_overlaps` | 死字段，删除 |

### 还剩 11 个冗余字段

**按优先级：**

1. `common_blocks`（9 文件引用）→ 移到 PaveFiller。需改 `bopds/tools::perform_common_blocks` 签名 + `brep_algo_api/section.rs` + PaveFiller 内部引用。尝试过一次但半途回退了。

2. `shared_topology` + `shape_sd`（10 文件引用）→ 移到 PaveFiller。`ds.detect_shared_topology()` 需改为独立函数。`BOPDS_Iterator` 也需要 shape_sd 来过滤同域形状对。

3. `intersection_curves`（21 文件引用）→ 移到 PaveFiller。Builder 也读取（用于 IC 边创建），需另传参。

4. `locations`（8 文件引用）→ 烤入 TShape。深层次 BRep 加载变更。

5. `vertices` / `edges` / `faces` + `*_shape_idx`（6 字段）→ 全代码库变更。最大任务，需专攻。

### 一键继续提示词

将以下文本粘贴到新 session：
