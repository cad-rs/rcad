# DS Parallel Array Removal — Handover 2026-07-27

## Progress

**Phase 0 ✅** — ShapeInfo新增4字段 + DS两个pool（face_info_pool + pave_blocks_pool），全部构造点已更新
**Phase 1 ✅** — 新增 vertex_count/edge_count/face_count 访问器
**Phase 2 ✅ ~90%** — 批量脚本 + 手动修复覆盖了大部分直接字段访问
**Phase 4 ✅** — 从 DS 删除 5 个 image 字段（my_images, my_origins, wire_images, shell_images, solid_images）

### Phase 4 详情

脱掉的字段：
| 字段 | 原因 |
|------|------|
| `my_images: Vec<Vec<usize>>` | 只被 `build_edge_images()` 写（死代码，从未调用） |
| `my_origins: Vec<usize>` | 只被 `build_edge_images()` 写（同上） |
| `wire_images` | 只被 `pipeline_dump.rs` 读（调试代码） |
| `shell_images` | 被设置但从未读取 |
| `solid_images` | 被设置但从未读取 |

同时删除：
- `build_edge_images()` 方法（死代码）
- `build_container_images()` → no-op
- `rebuild_wire_edges_ds()`、`wire_has_split_edges_ds()` 辅助函数
- `history::update_with_post_treat()` 中 `ds.my_images` 引用

**DS 字段数：42 → 37**

### 剩余工作

**Phase 2 剩余 ~111 处绑定模式，分布在 17 个文件：**

剩下 `let VAR = &ds.edges[EXPR]` / `let VAR = &ds.faces[EXPR]` 模式。
主要文件：
- `glue.rs` (14处)、`intersection.rs` (12处)、`classify.rs` (10处)
- `boptools/extra.rs` (12处)、`boptools/mod.rs` (18处)
- `interf.rs` (9处)、`filler_mod.rs` (8处)、`paves.rs` (7处)
- `shell_splitter.rs` (4处)、`builder.rs` (1处)、`make_blocks.rs` (3处)
- `posttreat.rs` (2处)、`checker_si.rs` (1处)

**Phase 3** — Phase 2 完成后，删除 `DSVertex` / `DSEdge` / `DSFace` 结构体，删除 `vertices / edges / faces / vertex_shape_idx / edge_shape_idx / face_shape_idx` 字段（再减 7 个字段 → DS 降到 ~30 字段）
