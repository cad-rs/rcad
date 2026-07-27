# DS Parallel Array Removal — Handover 2026-07-27

## Progress

**Phase 0 ✅** — ShapeInfo新增4字段 + DS两个pool
**Phase 1 ✅** — 新增 vertex_count/edge_count/face_count 访问器
**Phase 2 ✅ ~65%** — 大量直接字段访问已替换为访问器
**Phase 4 ✅** — 从 DS 删除 5 个 image 字段（42→37字段）

### 已完成的 Phase 2

通过批量脚本 + 手动修复覆盖了 ~200 处直接访问，包括：
- `ds.vertices[vi].point` → `ds.vertex_point(vi)`
- `ds.edges[ei].start_vertex/end_vertex/t_range/curve` → 访问器
- `ds.faces[fi].surface/boundary_edges/normal/origin/...` → 访问器
- `ds.faces[fi].face_info` → `ds.face_info(fi)/face_info_mut(fi)`
- `ds.vertices/edges/faces.len()` → `vertex_count/edge_count/face_count`
- 已修复文件：builder.rs, builder_face.rs, builder_utils.rs, edge_builders.rs,
  result_builder.rs, result_build_mod.rs, checker_si.rs, posttreat.rs, make_blocks.rs

### Phase 4 完成

从 DS 删除：`my_images`, `my_origins`, `wire_images`, `shell_images`, `solid_images`
DS 字段数：42 → 37

### 剩余工作

**Phase 2 剩余 105 处绑定模式，分布在 12 个文件：**

这些全是 `let VAR = &ds.TYPE[IDX]` → `VAR.FIELD` 模式。
难点：同一个变量名在不同作用域中绑定到不同的索引表达式（如 `df` = `fi` 在函数 A，
但 `df` = `dfi` 在函数 B）。需要逐个文件手动修复。

剩余文件：
- `filler_mod.rs` — df复用最复杂
- `glue.rs` — face1/face2多作用域
- `classify.rs` — face/edge多作用域
- `intersection.rs` — e/f单字符变量
- `interf.rs` — 同上
- `mod.rs` (pave_filler) — e/f单字符
- `paves.rs` — e/f单字符
- `shell_splitter.rs` — face多作用域
- `boptools/extra.rs` — edge多作用域
- `boptools/mod.rs` — e/f单字符
- `ff_intersect.rs` — edge/face
- `inttools/context.rs` — edge/face

**Phase 3** — Phase 2 完成后，删除 DSVertex/DSEdge/DSFace + 并行数组（再减 7 字段）
