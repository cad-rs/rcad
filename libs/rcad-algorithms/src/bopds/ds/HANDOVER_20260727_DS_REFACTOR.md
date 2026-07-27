# DS Parallel Array Removal — Handover 2026-07-27

## Progress

**Phase 0 ✅** — ShapeInfo新增4字段 + DS两个pool（face_info_pool + pave_blocks_pool），全部构造点已更新
**Phase 1 ✅** — 新增 vertex_count/edge_count/face_count 访问器
**Phase 2 ✅ ~90% 完成** — 批量脚本 + 手动修复覆盖了大部分直接字段访问

### 已完成的 Phase 2 工作

1. **batch 1 (fix_ds_access.py, 15 files)** — 批量替换简单模式：
   - `ds.vertices[vi].point` → `ds.vertex_point(vi)`
   - `ds.edges[ei].start_vertex/end_vertex` → `edge_start_vertex_ds/edge_end_vertex_ds`
   - `ds.edges[ei].t_range` → `ds.edge_range(ei)`
   - `ds.faces[fi].boundary_edges/boundary_verts` → `face_boundary_edges/face_boundary_verts`
   - `ds.faces[fi].origin/normal/location/natural_restriction` → accessors
   - `ds.vertices/edges/faces.len()` → `vertex_count/edge_count/face_count`

2. **batch 2 (fix_ds_access_2b.py, 12 files)** — surface/curve/face_info 模式：
   - `&ds.faces[fi].surface` → `ds.face_surface(fi).unwrap()`
   - `ds.faces[fi].surface.clone()` → `ds.face_surface(fi).cloned().unwrap()`
   - `&ds.edges[ei].curve` → `ds.edge_curve(ei).unwrap()`
   - `ds.faces[fi].face_info` → `ds.face_info(fi)/face_info_mut(fi)`
   - `ds.vertices[vi].geom_tol` → `ds.vertex_tolerance(vi)`

3. **手动修复：**
   - `edge_builders.rs`: 删除死函数 `build_sphere_seam_segments` 和 `build_cylinder_seam_segments`；`is_split_to_reverse` 改用访问器
   - `builder_face.rs`: `source_face_idx/natural_restriction/shapes_to_segments` 访问器迁移
   - `builder_utils.rs`: `is_tangent_face` + `build_edge_bounds` 改用访问器
   - `result_build_mod.rs`: edge binding 模式改用直接访问器
   - `builder.rs`: `build_draft_face` 改用 `face_boundary_edges/forwards/inner_boundary` 访问器
   - `result_builder.rs`: `emit_wire_face` 移除 face binding，改用 `face_normal/face_surface`
   - `mod.rs`: 新增 `face_boundary_edge_forwards()` 访问器

### 剩余工作

**Phase 2 剩余 ~111 处绑定模式，分布在 17 个文件：**

剩下的都是 `let VAR = &ds.edges[EXPR]` / `let VAR = &ds.faces[EXPR]` 这种先绑定再读写字段的模式。
每个文件需要单独处理。

主要文件列表：
- `glue.rs` (14处) — face/edge binding
- `intersection.rs` (12处) — edge/face binding
- `classify.rs` (10处) — face binding
- `boptools/extra.rs` (12处) — edge/face binding
- `boptools/mod.rs` (18处) — edge/face binding
- `interf.rs` (9处) — edge binding (pave_blocks)
- `filler_mod.rs` (8处) — edge/face binding
- `paves.rs` (7处) — edge/face binding
- `shell_splitter.rs` (4处) — face binding
- `builder.rs` (1处) — edge binding in closure
- `make_blocks.rs` (3处) — edge binding
- `posttreat.rs` (2处) — edge binding
- `checker_si.rs` (1处) — face binding
- 其他

**Phase 3** — 上述全部完成后，删除 `DSVertex` / `DSEdge` / `DSFace` 结构体，删除 `vertices / edges / faces / vertex_shape_idx / edge_shape_idx / face_shape_idx` 字段
**Phase 4** — 将 `my_images / my_origins / wire_images / shell_images / solid_images` 从 DS 移到 `BooleanBuilder`
