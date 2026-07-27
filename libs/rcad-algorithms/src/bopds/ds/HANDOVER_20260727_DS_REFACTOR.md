# DS Parallel Array Removal — Handover 2026-07-27

## Progress

**Phase 0 ✅** — ShapeInfo新增4字段 + DS两个pool（face_info_pool + pave_blocks_pool），全部构造点已更新
**Phase 1 ✅** — 新增 vertex_count/edge_count/face_count 访问器
**Phase 2 ✅ 已完成 30+/35 个文件，~900/1076 处替换**

## 剩余工作

**Phase 2 剩余 ~176 处，分布在 5 个文件：**

1. **`result_builder.rs` (~80 处)** — `emit_wire_face` 函数内部有 `let face = &ds.faces[face_idx]` 绑定，大量使用 `face.normal` / `face.surface` / `face.location` / `face.geom_tol` / `face.source_face_idx` 等，需要用 `ds.face_*(face_idx)` 逐一替换

2. **`builder.rs` (~30 处)** — `build_draft_face`、`build_result_topo` 等函数接收 `face: &DSFace` 参数，需改为接收 `(ds: &DS, face_idx: usize)` 并用访问器

3. **`edge_builders.rs` (~20 处)** — 函数签名接收 `face: &DSFace` / `ds_edge: &DSEdge`，需改为 `(ds: &DS, fi: usize, ei: usize)` 

4. **`pipeline_dump.rs` (~10 处)** — 调试代码，迭代 `ds.vertices.iter()` / `ds.edges.iter()` / `ds.faces.iter()` 并读字段

5. **`integration_tests.rs` (~36 处)** — 测试代码，直接读字段

**Phase 3** — 上述全部完成后，删除 `DSVertex` / `DSEdge` / `DSFace` 结构体，删除 `vertices / edges / faces / vertex_shape_idx / edge_shape_idx / face_shape_idx` 字段

**Phase 4** — 将 `my_images / my_origins / wire_images / shell_images / solid_images` 从 DS 移到 `BooleanBuilder`
