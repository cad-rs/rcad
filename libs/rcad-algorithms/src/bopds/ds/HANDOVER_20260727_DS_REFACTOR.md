# DS Parallel Array Removal — Handover 2026-07-27

## Progress

**Phase 0 ✅** — ShapeInfo新增4字段 + DS两个pool
**Phase 1 ✅** — 新增 vertex_count/edge_count/face_count 访问器
**Phase 2 ✅ ~95%** — 所有 `let VAR = &ds.TYPE[IDX]` 绑定已转换为 `.get().unwrap()`
**Phase 4 ✅** — 从 DS 删除 5 个 image 字段（42→37字段）

### Phase 2 完成情况

**已修复：**
- 所有 `let edge = &ds.edges[ei]` 模式 → `let edge = ds.edges.get(ei).unwrap()`
- 所有 `let face = &ds.faces[fi]` 模式 → `let face = ds.faces.get(fi).unwrap()`
- 覆盖 13 个文件，~200 处绑定
- 所有绑定通过 `.get()` 安全访问，不再直接索引 vec

**仍剩 ~88 处直接内联访问**（非绑定模式）：
```
ds.edges[ei].curve.point_at(t)   →   ds.edge_curve(ei).unwrap().point_at(t)
ds.faces[fi].surface             →   ds.face_surface(fi).unwrap()
ds.faces[fi].face_info           →   ds.face_info(fi)
ds.edges[ei].pave_blocks[pi]    →   ds.edge_pave_blocks(ei)[pi]
```

这些分布于 13 个文件中，需要在 Phase 3 前处理。

### Phase 3（下一步）

删除 `DSVertex` / `DSEdge` / `DSFace` 结构体及对应 vec：
- 删除 `vertices` / `edges` / `faces` 字段
- 删除 `vertex_shape_idx` / `edge_shape_idx` / `face_shape_idx`
- 所有现有访问器方法改为从 TShape + pool 读取
- 再减 7 字段 → DS ~30 字段
